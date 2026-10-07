#!/usr/bin/env python3
"""
fix_rotamers.py -- step 2 of the consistent FGF rebuild.

For any FGF residue still involved in a ring piercing after
prepare_fgf_core.py (as the ring owner or as the piercing atom), or with a
side-chain heavy atom closer than SEVERE A to a heavy atom of the scaffold or
of a non-adjacent FGF residue (e.g. a Lys NZ sitting on a heparin C1), scan its
side-chain chi angles on a 30-degree grid and keep the conformation with no
piercing and the fewest heavy-atom contacts < 2.5 A against the scaffold and
the rest of the FGF. Backbone is never moved. Rewrites fgf_PROC.pdb /
fgf_PROD.pdb in place and appends to core_report.txt.

Usage: python3 fix_rotamers.py            (all systems)
"""
import itertools
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
SEVERE = 1.8
CHIS = {  # (axis atom 1, axis atom 2, atoms moved by rotating about it)
    "LYS": [("CA", "CB", "CG CD CE NZ"), ("CB", "CG", "CD CE NZ"), ("CG", "CD", "CE NZ")],
    "ARG": [("CA", "CB", "CG CD NE CZ NH1 NH2"), ("CB", "CG", "CD NE CZ NH1 NH2"), ("CG", "CD", "NE CZ NH1 NH2")],
    "TYR": [("CA", "CB", "CG CD1 CD2 CE1 CE2 CZ OH"), ("CB", "CG", "CD1 CD2 CE1 CE2 CZ OH")],
    "PHE": [("CA", "CB", "CG CD1 CD2 CE1 CE2 CZ"), ("CB", "CG", "CD1 CD2 CE1 CE2 CZ")],
    "TRP": [("CA", "CB", "CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2"), ("CB", "CG", "CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2")],
    "HSD": [("CA", "CB", "CG ND1 CD2 CE1 NE2"), ("CB", "CG", "ND1 CD2 CE1 NE2")],
    "HSE": [("CA", "CB", "CG ND1 CD2 CE1 NE2"), ("CB", "CG", "ND1 CD2 CE1 NE2")],
    "GLN": [("CA", "CB", "CG CD OE1 NE2"), ("CB", "CG", "CD OE1 NE2")],
    "GLU": [("CA", "CB", "CG CD OE1 OE2"), ("CB", "CG", "CD OE1 OE2")],
    "ASN": [("CA", "CB", "CG OD1 ND2"), ("CB", "CG", "OD1 ND2")],
    "ASP": [("CA", "CB", "CG OD1 OD2"), ("CB", "CG", "OD1 OD2")],
    "MET": [("CA", "CB", "CG SD CE"), ("CB", "CG", "SD CE")],
    "LEU": [("CA", "CB", "CG CD1 CD2"), ("CB", "CG", "CD1 CD2")],
    "ILE": [("CA", "CB", "CG1 CG2 CD"), ("CB", "CG1", "CD")],
    "SER": [("CA", "CB", "OG")], "THR": [("CA", "CB", "OG1 CG2")], "VAL": [("CA", "CB", "CG1 CG2")],
    "CYS": [("CA", "CB", "SG")],
}
RINGS = {"PHE": "CG CD1 CD2 CE1 CE2 CZ", "TYR": "CG CD1 CD2 CE1 CE2 CZ", "TRP": "CD2 CE2 CE3 CZ2 CZ3 CH2",
         "HSD": "CG ND1 CD2 CE1 NE2", "HSE": "CG ND1 CD2 CE1 NE2", "PRO": "N CA CB CG CD"}
SUGAR = "C1 C2 C3 C4 C5 O5".split()


def rotate(xyz, a, b, ang):
    k = (b - a) / np.linalg.norm(b - a)
    v = xyz - b
    return b + v * math.cos(ang) + np.cross(k, v) * math.sin(ang) + np.outer(v @ k, k) * (1 - math.cos(ang))


def load(path):
    rows = []
    for l in open(path):
        if l.startswith("ATOM"):
            rows.append([l, l[72:76].strip(), int(l[22:26]), l[17:21].strip(), l[12:16].strip(),
                         np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])])
    return rows


def rings_of(rows):
    by = {}
    for r in rows:
        by.setdefault((r[1], r[2], r[3]), {})[r[4]] = r[5]
    out = []
    for (seg, rid, rn), at in by.items():
        names = RINGS.get(rn, "").split() or (SUGAR if seg.startswith("CAR") else [])
        if names and all(n in at for n in names):
            out.append(((seg, rid), np.mean([at[n] for n in names], axis=0)))
    return out


def score(res_atoms, key, env, env_rings):
    """(piercings, contacts) for one residue's side-chain atoms vs environment."""
    side = {n: x for n, x in res_atoms.items() if n not in ("N", "CA", "C", "O")}
    pierce = sum(1 for _, c in env_rings for x in side.values() if np.linalg.norm(c - x) < 1.6)
    sx_ = np.array(list(side.values()))
    if len(sx_):
        pierce += int(np.sum(np.linalg.norm(sx_[:, None, :] - env[None, :, :], axis=2) < SEVERE))
    names = RINGS.get(key[2], "").split()
    if names and all(n in res_atoms for n in names):
        c = np.mean([res_atoms[n] for n in names], axis=0)
        pierce += int(np.sum(np.linalg.norm(env - c, axis=1) < 1.6))
    sx = np.array(list(side.values()))
    contacts = int(np.sum(np.linalg.norm(sx[:, None, :] - env[None, :, :], axis=2) < 2.5)) if len(sx) else 0
    return pierce, contacts


def fix_system(sd):
    pre = sd / "precheck.pdb"
    out = subprocess.run([sys.executable, str(BASE.parent / "check_ring_piercing.py"), str(pre)],
                         capture_output=True, text=True).stdout.splitlines()[1:]
    targets = set()
    for line in out:  # "SEG:RESnnn ring pierced by SEG:RESnnn:ATOM (d A)"
        owner, piercer = line.split(" ring pierced by ")
        for tok in (owner.strip(), ":".join(piercer.split(":")[:2])):
            seg, rr = tok.split(":")
            if seg in ("PROC", "PROD"):
                targets.add((seg, int("".join(ch for ch in rr if ch.isdigit()))))
    allrows = load(pre)
    xyz = np.array([r[5] for r in allrows])
    for i, r in enumerate(allrows):
        if r[1] not in ("PROC", "PROD") or r[4] in ("N", "CA", "C", "O"):
            continue
        d = np.linalg.norm(xyz - r[5], axis=1)
        for j in np.nonzero(d < SEVERE)[0]:
            o = allrows[j]
            if o[1] == r[1] and abs(o[2] - r[2]) <= 1:
                continue
            targets.add((r[1], r[2]))
    if not targets:
        return []
    notes = []
    for seg in ("PROC", "PROD"):
        path = sd / f"fgf_{seg}.pdb"
        rows = load(path)
        for tseg, rid in sorted(targets):
            if tseg != seg:
                continue
            res = [r for r in rows if r[2] == rid]
            rn = res[0][3]
            if rn not in CHIS:
                notes.append(f"{seg}:{rn}{rid} has no chi definition; left for minimization")
                continue
            env_rows = [r for r in allrows if not (r[1] == seg and r[2] == rid)]
            env = np.array([r[5] for r in env_rows])
            env_rings = [c for k, c in rings_of(env_rows)]
            env_rings = [(None, c) for c in env_rings]
            orig = {r[4]: r[5].copy() for r in res}
            best = (score(orig, (seg, rid, rn), env, env_rings), orig)
            grid = [math.radians(a) for a in range(0, 360, 30)]
            for angs in itertools.product(grid, repeat=len(CHIS[rn])):
                at = {n: x.copy() for n, x in orig.items()}
                for (a, b, moving), ang in zip(CHIS[rn], angs):
                    mv = [n for n in moving.split() if n in at]
                    at.update(zip(mv, rotate(np.array([at[n] for n in mv]), at[a], at[b], ang)))
                sc = score(at, (seg, rid, rn), env, env_rings)
                if sc < best[0]:
                    best = (sc, at)
            for r in rows:
                if r[2] == rid:
                    r[5] = best[1][r[4]]
            notes.append(f"{seg}:{rn}{rid} rotamer scan -> piercings {best[0][0]}, contacts<2.5A {best[0][1]}")
        with open(path, "w") as fh:
            for l, *_rest, xyz in rows:
                fh.write(f"{l[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{l[54:]}")
            fh.write("END\n")
    with open(sd / "core_report.txt", "a") as fh:
        fh.write("".join(f"rotamer fix: {n}\n" for n in notes))
    return notes


def main():
    for tk in ("1fq9", "1e0o"):
        for sd in sorted((BASE / tk).iterdir()):
            for n in fix_system(sd):
                print(f"{tk}/{sd.name}: {n}")


if __name__ == "__main__":
    main()
