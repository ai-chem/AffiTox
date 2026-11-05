#!/usr/bin/env python3
"""Cognate co-crystal ligand redocking RMSD anchor (CASF-style docking power).

Extracts drug-like co-crystal ligands from AffiTox input PDBs, redocks with
QVina2 and/or GNINA using the same receptor PDBQTs and LaBOX boxes as the
BindingDB panel, and reports heavy-atom symmetry-corrected RMSD vs the crystal
pose. DynamicBind / Boltz-2 cognate redocking is out of scope for this script
(GPU inference); PLAPT has no poses.

Outputs:
  analysis/tables/cognate_rmsd_redocking.csv
  analysis/tables/cognate_rmsd_summary.csv
  analysis/cognate_rmsd_work/  (refs, ligands, docked poses)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# Drug-like co-crystal ligand residue codes (exclude cofactors/solvents/glycans).
# None => no usable small-molecule cognate in the AffiTox surrogate PDB.
COGNATE_RESNAMES: Dict[str, Optional[str]] = {
    "1g5m": None,  # BCL2: Bak-peptide context; no drug-like HETATM in deposited PDB
    "2v5z": "SAG",  # safinamide (FAD is cofactor)
    "2z5x": "HRM",  # harmine (FAD is cofactor)
    "3eyg": "MI1",
    "3jy9": "JZH",
    "3lxk": "MI1",
    "3mjg": None,  # PDGFRB surrogate: glycans only (NAG/NDG)
    "4ase": "AV9",
    "4f65": "0S9",
    "4tz4": "LVY",
    "4zau": "YY3",
    "5jkv": "ASD",  # androstenedione substrate (HEM cofactor; 1PE solvent)
    "5mo4": "AY7",  # asciminib (allosteric ABL1 pose used in AffiTox)
    "6gqj": "F82",
    "6jok": "B49",
    "5lf3": "BO2",  # bortezomib in the chain-K site used for the panel box
    "7awe": "S5K",
    "7kk3": "2YQ",
    "11ue": None,  # PDGFRB kinase surrogate: no co-crystal HETATM in the deposited PDB
}

SKIP_NOTES: Dict[str, str] = {
    "1g5m": "No drug-like co-crystal HETATM in input/proteins/1g5m.pdb (Bak-peptide context).",
    "3mjg": "Surrogate PDB contains glycans (NAG/NDG) only; no drug-like cognate.",
    "11ue": "No HETATM records in input/proteins/11ue.pdb; cognate RMSD was not computed.",
}


def _repo_root() -> Path:
    return Path(os.environ.get("TOXAFFINITY_ROOT", Path(__file__).resolve().parents[2]))


def _het_centroid(lines: Sequence[str]) -> Tuple[float, float, float]:
    xs, ys, zs = [], [], []
    for line in lines:
        if not line.startswith("HETATM"):
            continue
        elem = (line[76:78].strip() if len(line) >= 78 else "") or line[12:16].strip()[:1]
        if elem.upper() == "H":
            continue
        xs.append(float(line[30:38]))
        ys.append(float(line[38:46]))
        zs.append(float(line[46:54]))
    if not xs:
        return (1e9, 1e9, 1e9)
    return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


def extract_hetatm(
    pdb_path: Path,
    resname: str,
    out_path: Path,
    prefer_center: Optional[Sequence[float]] = None,
) -> int:
    """Write HETATM lines for one residue copy of resname.

    Symmetrical assemblies repeat the ligand. A multi-residue PDBQT is not a
    valid Vina/GNINA ligand. When a docking-box center is given, the copy
    nearest that center is kept (the site the receptor PDBQT was prepared
    around). Otherwise chain A is preferred, then the lowest residue number.
    """
    by_chain_res: Dict[tuple, List[str]] = {}
    for line in pdb_path.read_text(errors="ignore").splitlines(True):
        if line.startswith("HETATM") and line[17:20].strip() == resname:
            chain = line[21].strip() or "_"
            resseq = line[22:26].strip() or "0"
            by_chain_res.setdefault((chain, resseq), []).append(line)
    if not by_chain_res:
        raise RuntimeError(f"No HETATM resname={resname} in {pdb_path}")
    if prefer_center is not None:
        cx, cy, cz = prefer_center

        def _dist(key: tuple) -> float:
            x, y, z = _het_centroid(by_chain_res[key])
            return (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2

        chain, resseq = min(by_chain_res, key=_dist)
    else:
        chains = {chain for chain, _resseq in by_chain_res}
        if "A" in chains:
            chain = "A"
        else:
            chain = max(
                chains,
                key=lambda c: sum(len(v) for (ch, _r), v in by_chain_res.items() if ch == c),
            )

        def _resseq_key(r: str) -> tuple:
            try:
                return (0, int(r))
            except ValueError:
                return (1, r)

        resseq = min((r for ch, r in by_chain_res if ch == chain), key=_resseq_key)
    lines = by_chain_res[(chain, resseq)]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(lines) + "END\n")
    return len(lines)


def cognate_box(
    ref_pdb: Path, padding: float = 10.0, min_size: float = 20.0
) -> Tuple[List[float], List[float]]:
    """Axis-aligned box around cognate heavy atoms (CASF-style redock box)."""
    xs, ys, zs = [], [], []
    for line in ref_pdb.read_text(errors="ignore").splitlines():
        if not line.startswith("HETATM"):
            continue
        elem = (line[76:78].strip() if len(line) >= 78 else "") or ""
        name = line[12:16].strip()
        if (elem or name[:1]).upper() == "H":
            continue
        xs.append(float(line[30:38]))
        ys.append(float(line[38:46]))
        zs.append(float(line[46:54]))
    if not xs:
        raise RuntimeError(f"No heavy atoms for box in {ref_pdb}")
    center = [
        0.5 * (min(xs) + max(xs)),
        0.5 * (min(ys) + max(ys)),
        0.5 * (min(zs) + max(zs)),
    ]
    size = [
        max(min_size, (max(xs) - min(xs)) + 2 * padding),
        max(min_size, (max(ys) - min(ys)) + 2 * padding),
        max(min_size, (max(zs) - min(zs)) + 2 * padding),
    ]
    return center, size


AD4_TO_ELEMENT = {
    "H": "H",
    "HD": "H",
    "HS": "H",
    "C": "C",
    "A": "C",
    "N": "N",
    "NA": "N",
    "NS": "N",
    "OA": "O",
    "OS": "O",
    "O": "O",
    "S": "S",
    "SA": "S",
    "P": "P",
    "F": "F",
    "CL": "CL",
    "BR": "BR",
    "I": "I",
    "MG": "MG",
    "MN": "MN",
    "ZN": "ZN",
    "CA": "CA",
    "FE": "FE",
}


def pdb_to_pdbqt(pdb_in: Path, pdbqt_out: Path, obabel: str = "obabel") -> None:
    pdbqt_out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [obabel, str(pdb_in), "-O", str(pdbqt_out), "-xh"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not pdbqt_out.exists():
        raise RuntimeError(
            f"obabel failed for {pdb_in}:\n{proc.stdout}\n{proc.stderr}"
        )
    # AutoDock Vina / QVina2 lack a boron atom type; type B as carbon for docking only.
    text = pdbqt_out.read_text()
    if re.search(r"\sB\s*$", text, flags=re.M) or " B " in text:
        fixed_lines = []
        for line in text.splitlines(True):
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 79:
                ad = line[77:79]
                if ad.strip().upper() == "B":
                    line = line[:77] + "C " + line[79:]
            fixed_lines.append(line)
        pdbqt_out.write_text("".join(fixed_lines))
        (pdbqt_out.parent / (pdbqt_out.stem + ".BORON_AS_CARBON")).write_text(
            "Boron AutoDock type B remapped to C for Vina-family docking.\n"
        )


def load_box(box_json: Path) -> Tuple[List[float], List[float]]:
    data = json.loads(box_json.read_text())
    return list(data["center"]), list(data["size"])


def run_qvina(
    binary: str,
    receptor: Path,
    ligand: Path,
    out_pdbqt: Path,
    center: Sequence[float],
    size: Sequence[float],
    exhaustiveness: int = 8,
    num_modes: int = 9,
    seed: int = 42,
) -> None:
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        binary,
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--center_x",
        str(center[0]),
        "--center_y",
        str(center[1]),
        "--center_z",
        str(center[2]),
        "--size_x",
        str(size[0]),
        "--size_y",
        str(size[1]),
        "--size_z",
        str(size[2]),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        str(num_modes),
        "--seed",
        str(seed),
        "--out",
        str(out_pdbqt),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not out_pdbqt.exists():
        raise RuntimeError(f"qvina failed:\n{proc.stdout}\n{proc.stderr}")


def run_gnina(
    binary: str,
    receptor: Path,
    ligand: Path,
    out_sdf: Path,
    center: Sequence[float],
    size: Sequence[float],
    exhaustiveness: int = 8,
    num_modes: int = 9,
    seed: int = 42,
    env: Optional[dict] = None,
) -> None:
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        binary,
        "-r",
        str(receptor),
        "-l",
        str(ligand),
        "--center_x",
        str(center[0]),
        "--center_y",
        str(center[1]),
        "--center_z",
        str(center[2]),
        "--size_x",
        str(size[0]),
        "--size_y",
        str(size[1]),
        "--size_z",
        str(size[2]),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        str(num_modes),
        "--seed",
        str(seed),
        "--no_gpu",
        "-o",
        str(out_sdf),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if proc.returncode != 0 or not out_sdf.exists():
        err = (proc.stderr or proc.stdout or "").strip().splitlines()
        brief = err[-1] if err else "unknown"
        raise RuntimeError(f"gnina failed: {brief}")


def _coords_from_pdb_het(pdb_path: Path) -> List[Tuple[str, Tuple[float, float, float]]]:
    atoms = []
    for line in pdb_path.read_text(errors="ignore").splitlines():
        if not line.startswith("HETATM"):
            continue
        elem = line[76:78].strip() if len(line) >= 78 else ""
        if not elem:
            name = line[12:16].strip()
            elem = re.sub(r"[0-9]", "", name)[:2] or "C"
        elem = elem.upper()
        if elem in {"H", "D"}:
            continue
        if elem == "CL":
            elem = "CL"
        elif elem == "BR":
            elem = "BR"
        atoms.append((elem, (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return atoms


def _coords_from_pdbqt_model1(pdbqt_path: Path) -> List[Tuple[str, Tuple[float, float, float]]]:
    atoms = []
    model_idx = 0
    take = False
    saw_model = False
    for line in pdbqt_path.read_text(errors="ignore").splitlines():
        if line.startswith("MODEL"):
            saw_model = True
            model_idx += 1
            take = model_idx == 1
            continue
        if line.startswith("ENDMDL"):
            if take:
                break
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if saw_model and not take:
            continue
        if not saw_model:
            take = True
        ad_type = line[77:79].strip().upper() if len(line) >= 79 else ""
        elem = AD4_TO_ELEMENT.get(ad_type, ad_type[:1] if ad_type else "C")
        if elem == "H":
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atoms.append((elem, (x, y, z)))
    if not atoms:
        raise RuntimeError(f"No heavy atoms in {pdbqt_path}")
    return atoms


def _coords_from_sdf_mol1(sdf_path: Path) -> List[Tuple[str, Tuple[float, float, float]]]:
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise RuntimeError("RDKit required for SDF parsing") from exc
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise RuntimeError(f"No molecule in {sdf_path}")
    atoms = []
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        atoms.append((atom.GetSymbol().upper(), (pos.x, pos.y, pos.z)))
    return atoms


def hungarian_rmsd(
    ref: Sequence[Tuple[str, Tuple[float, float, float]]],
    mob: Sequence[Tuple[str, Tuple[float, float, float]]],
) -> float:
    """Element-matched min-assignment RMSD (no Kabsch; poses share pocket frame)."""
    try:
        import numpy as np
        from scipy.optimize import linear_sum_assignment
    except ImportError as exc:
        raise RuntimeError("numpy/scipy required for RMSD") from exc

    def _norm_el(e: str) -> str:
        e = e.upper()
        # Vina-family boron→carbon remapping for docking
        if e == "B":
            return "C"
        return e

    ref_n = [(_norm_el(e), xyz) for e, xyz in ref]
    mob_n = [(_norm_el(e), xyz) for e, xyz in mob]
    if len(ref_n) != len(mob_n):
        raise RuntimeError(f"Atom-count mismatch: ref={len(ref_n)} mob={len(mob_n)}")

    elements = sorted({e for e, _ in ref_n})
    if set(e for e, _ in mob_n) != set(elements):
        raise RuntimeError(
            f"Element composition mismatch: ref={ {e:sum(1 for x,_ in ref_n if x==e) for e in elements} } "
            f"mob={ {e:sum(1 for x,_ in mob_n if x==e) for e in sorted({y for y,_ in mob_n})} }"
        )

    ref_coords_ordered: List[Tuple[float, float, float]] = []
    mob_coords_ordered: List[Tuple[float, float, float]] = []

    for el in elements:
        ref_i = [i for i, (e, _) in enumerate(ref_n) if e == el]
        mob_i = [i for i, (e, _) in enumerate(mob_n) if e == el]
        if len(ref_i) != len(mob_i):
            raise RuntimeError(f"Element count mismatch for {el}")
        R = np.array([ref_n[i][1] for i in ref_i], dtype=float)
        M = np.array([mob_n[i][1] for i in mob_i], dtype=float)
        diff = R[:, None, :] - M[None, :, :]
        cost = np.sum(diff * diff, axis=2)
        ri, mi = linear_sum_assignment(cost)
        for a, b in zip(ri, mi):
            ref_coords_ordered.append(tuple(R[a]))
            mob_coords_ordered.append(tuple(M[b]))

    R = np.asarray(ref_coords_ordered, dtype=float)
    M = np.asarray(mob_coords_ordered, dtype=float)
    return float(math.sqrt(np.mean(np.sum((R - M) ** 2, axis=1))))


def try_rdkit_best_rmsd(ref_pdb: Path, mobile_path: Path, mobile_fmt: str) -> Optional[float]:
    """Optional RDKit GetBestRMS when both molecules sanitize cleanly."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolAlign
    except ImportError:
        return None

    ref = Chem.MolFromPDBFile(str(ref_pdb), removeHs=True, sanitize=True)
    if ref is None:
        ref = Chem.MolFromPDBFile(str(ref_pdb), removeHs=True, sanitize=False)
        if ref is None:
            return None
        try:
            Chem.SanitizeMol(ref)
        except Exception:
            return None

    if mobile_fmt == "pdbqt":
        # Convert via obabel to SDF first
        tmp = mobile_path.with_suffix(".tmp.sdf")
        proc = subprocess.run(
            ["obabel", str(mobile_path), "-O", str(tmp), "-f", "1", "-l", "1"],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0 or not tmp.exists():
            return None
        mob = Chem.SDMolSupplier(str(tmp), removeHs=True)[0]
        tmp.unlink(missing_ok=True)
    else:
        mob = Chem.SDMolSupplier(str(mobile_path), removeHs=True)[0]

    if mob is None:
        return None
    try:
        return float(rdMolAlign.GetBestRMS(mob, ref))
    except Exception:
        try:
            return float(AllChem.CalcRMS(mob, ref))
        except Exception:
            return None


def summarize(rows: List[dict]) -> List[dict]:
    out = []
    for method in sorted({r["method"] for r in rows if r.get("rmsd_angstrom") not in (None, "")}):
        vals = [
            float(r["rmsd_angstrom"])
            for r in rows
            if r["method"] == method and r.get("rmsd_angstrom") not in (None, "")
        ]
        if not vals:
            continue
        n = len(vals)
        n_le2 = sum(1 for v in vals if v <= 2.0)
        vals_sorted = sorted(vals)
        mid = n // 2
        median = vals_sorted[mid] if n % 2 else 0.5 * (vals_sorted[mid - 1] + vals_sorted[mid])
        out.append(
            {
                "method": method,
                "n_targets": n,
                "n_rmsd_le_2A": n_le2,
                "pct_rmsd_le_2A": round(100.0 * n_le2 / n, 1),
                "median_rmsd_A": round(median, 3),
                "mean_rmsd_A": round(sum(vals) / n, 3),
                "min_rmsd_A": round(min(vals), 3),
                "max_rmsd_A": round(max(vals), 3),
            }
        )
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=None)
    ap.add_argument("--methods", nargs="+", default=["qvina", "gnina"])
    ap.add_argument("--targets", nargs="+", default=None)
    ap.add_argument("--qvina-bin", default=os.environ.get("QVINA_BIN", "qvina02"))
    ap.add_argument(
        "--gnina-bin",
        default=os.environ.get("GNINA_BIN", "gnina"),
    )
    ap.add_argument("--obabel", default="obabel")
    ap.add_argument("--skip-dock", action="store_true", help="Reuse existing docked poses")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument(
        "--merge-csv",
        action="store_true",
        help="Merge rows into existing cognate_rmsd_redocking.csv by (target,method)",
    )
    args = ap.parse_args(argv)

    root = args.root or _repo_root()
    work = root / "analysis" / "cognate_rmsd_work"
    refs = work / "refs"
    ligs = work / "ligands"
    docked = work / "docked"
    tables = root / "analysis" / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    # Ensure GNINA finds libcudart when launched outside its conda env.
    gnina_lib_raw = os.environ.get("GNINA_LIB", "")
    if gnina_lib_raw:
        gnina_lib = Path(gnina_lib_raw)
    else:
        conda_prefix = os.environ.get("CONDA_PREFIX", "")
        gnina_lib = Path(conda_prefix) / "lib" if conda_prefix else Path()
    run_env = os.environ.copy()
    if gnina_lib.is_dir():
        run_env["LD_LIBRARY_PATH"] = f"{gnina_lib}:{run_env.get('LD_LIBRARY_PATH', '')}"

    targets = args.targets or sorted(COGNATE_RESNAMES.keys())
    rows: List[dict] = []

    for target in targets:
        resname = COGNATE_RESNAMES.get(target)
        note = SKIP_NOTES.get(target, "")
        base = {
            "target": target,
            "cognate_resname": resname or "",
            "n_ref_heavy_atoms": "",
            "status": "",
            "note": note,
        }
        if resname is None:
            for method in args.methods:
                rows.append(
                    {
                        **base,
                        "method": method,
                        "rmsd_angstrom": "",
                        "rmsd_method": "",
                        "pose_path": "",
                        "status": "skipped_no_cognate",
                    }
                )
            continue

        pdb_src = root / "input" / "proteins" / f"{target}.pdb"
        if not pdb_src.exists():
            # try uppercase
            alt = root / "input" / "proteins" / f"{target.upper()}.pdb"
            pdb_src = alt if alt.exists() else pdb_src
        if not pdb_src.exists():
            for method in args.methods:
                rows.append(
                    {
                        **base,
                        "method": method,
                        "rmsd_angstrom": "",
                        "rmsd_method": "",
                        "pose_path": "",
                        "status": "missing_input_pdb",
                    }
                )
            continue

        ref_pdb = refs / f"{target}_{resname}.pdb"
        box_json = root / "processed" / "boxes" / f"{target}.json"
        prefer_center = None
        if box_json.exists():
            prefer_center = json.loads(box_json.read_text())["center"]
        try:
            n_lines = extract_hetatm(pdb_src, resname, ref_pdb, prefer_center=prefer_center)
        except Exception as exc:
            for method in args.methods:
                rows.append(
                    {
                        **base,
                        "method": method,
                        "rmsd_angstrom": "",
                        "rmsd_method": "",
                        "pose_path": "",
                        "status": f"extract_failed:{exc}",
                    }
                )
            continue

        ref_atoms = _coords_from_pdb_het(ref_pdb)
        base["n_ref_heavy_atoms"] = str(len(ref_atoms))

        lig_pdbqt = ligs / f"{target}_cognate.pdbqt"
        try:
            pdb_to_pdbqt(ref_pdb, lig_pdbqt, obabel=args.obabel)
        except Exception as exc:
            for method in args.methods:
                rows.append(
                    {
                        **base,
                        "method": method,
                        "rmsd_angstrom": "",
                        "rmsd_method": "",
                        "pose_path": "",
                        "status": f"pdbqt_failed:{exc}",
                    }
                )
            continue

        receptor = root / "processed" / "proteins" / f"{target}.pdbqt"
        if not receptor.exists():
            for method in args.methods:
                rows.append(
                    {
                        **base,
                        "method": method,
                        "rmsd_angstrom": "",
                        "rmsd_method": "",
                        "pose_path": "",
                        "status": "missing_receptor",
                    }
                )
            continue
        # Cognate-centered box (not the large AffiTox LaBOX volumes used for VS)
        center, size = cognate_box(ref_pdb, padding=10.0, min_size=20.0)

        for method in args.methods:
            row = {
                **base,
                "method": method,
                "rmsd_angstrom": "",
                "rmsd_method": "",
                "pose_path": "",
                "status": "",
            }
            try:
                if method == "qvina":
                    out = docked / "qvina" / f"{target}_out.pdbqt"
                    if not args.skip_dock or not out.exists():
                        qbin = shutil.which(args.qvina_bin) or args.qvina_bin
                        run_qvina(
                            qbin,
                            receptor,
                            lig_pdbqt,
                            out,
                            center,
                            size,
                            exhaustiveness=args.exhaustiveness,
                        )
                    mob = _coords_from_pdbqt_model1(out)
                    rmsd = hungarian_rmsd(ref_atoms, mob)
                    rd = try_rdkit_best_rmsd(ref_pdb, out, "pdbqt")
                    if rd is not None:
                        rmsd = rd
                        row["rmsd_method"] = "rdkit_GetBestRMS"
                    else:
                        row["rmsd_method"] = "hungarian_element_match"
                    row["pose_path"] = str(out.relative_to(root))
                    row["rmsd_angstrom"] = f"{rmsd:.3f}"
                    row["status"] = "ok"
                elif method == "gnina":
                    out = docked / "gnina" / f"{target}_out.sdf.gz"
                    # Prefer uncompressed sdf for RDKit
                    out_sdf = docked / "gnina" / f"{target}_out.sdf"
                    if not args.skip_dock or not out_sdf.exists():
                        run_gnina(
                            args.gnina_bin,
                            receptor,
                            lig_pdbqt,
                            out_sdf,
                            center,
                            size,
                            exhaustiveness=args.exhaustiveness,
                            env=run_env,
                        )
                    mob = _coords_from_sdf_mol1(out_sdf)
                    rmsd = hungarian_rmsd(ref_atoms, mob)
                    rd = try_rdkit_best_rmsd(ref_pdb, out_sdf, "sdf")
                    if rd is not None:
                        rmsd = rd
                        row["rmsd_method"] = "rdkit_GetBestRMS"
                    else:
                        row["rmsd_method"] = "hungarian_element_match"
                    row["pose_path"] = str(out_sdf.relative_to(root))
                    row["rmsd_angstrom"] = f"{rmsd:.3f}"
                    row["status"] = "ok"
                else:
                    row["status"] = f"unsupported_method:{method}"
            except Exception as exc:
                msg = str(exc).splitlines()[0][:200]
                row["status"] = f"error:{msg}"
            rows.append(row)
            print(
                f"[{target} {method}] status={row['status']} "
                f"rmsd={row['rmsd_angstrom']} atoms={base['n_ref_heavy_atoms']}",
                flush=True,
            )

    # Write detailed table (optionally merge with prior methods)
    detail_path = tables / "cognate_rmsd_redocking.csv"
    fieldnames = [
        "target",
        "cognate_resname",
        "method",
        "n_ref_heavy_atoms",
        "rmsd_angstrom",
        "rmsd_method",
        "status",
        "pose_path",
        "note",
    ]
    if args.merge_csv and detail_path.exists():
        import csv as _csv

        prev = list(_csv.DictReader(detail_path.open()))
        key = {(r.get("target"), r.get("method")) for r in rows}
        kept = [r for r in prev if (r.get("target"), r.get("method")) not in key]
        rows = kept + rows
        rows.sort(key=lambda r: (r.get("target", ""), r.get("method", "")))

    with detail_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})

    summary = summarize(rows)
    # annotate deferred methods
    for m in ("dynamicbind", "boltz2"):
        summary.append(
            {
                "method": m,
                "n_targets": "",
                "n_rmsd_le_2A": "",
                "pct_rmsd_le_2A": "",
                "median_rmsd_A": "",
                "mean_rmsd_A": "",
                "min_rmsd_A": "",
                "max_rmsd_A": "",
                "note": "Explicitly not computed in this session (GPU cognate redock deferred).",
            }
        )
    summary.append(
        {
            "method": "plapt",
            "n_targets": "",
            "n_rmsd_le_2A": "",
            "pct_rmsd_le_2A": "",
            "median_rmsd_A": "",
            "mean_rmsd_A": "",
            "min_rmsd_A": "",
            "max_rmsd_A": "",
            "note": "No poses; excluded from RMSD and PoseBusters.",
        }
    )
    sum_path = tables / "cognate_rmsd_summary.csv"
    sum_fields = [
        "method",
        "n_targets",
        "n_rmsd_le_2A",
        "pct_rmsd_le_2A",
        "median_rmsd_A",
        "mean_rmsd_A",
        "min_rmsd_A",
        "max_rmsd_A",
        "note",
    ]
    with sum_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sum_fields, extrasaction="ignore")
        w.writeheader()
        for r in summary:
            if "note" not in r:
                r["note"] = ""
            w.writerow(r)

    print(f"Wrote {detail_path}")
    print(f"Wrote {sum_path}")
    for s in summary:
        if s.get("n_targets"):
            print(
                f"  {s['method']}: n={s['n_targets']} "
                f"%≤2Å={s['pct_rmsd_le_2A']} median={s['median_rmsd_A']}"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
