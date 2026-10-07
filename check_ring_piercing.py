#!/usr/bin/env python3
"""
check_ring_piercing.py PDB [PDB ...]

Geometric ring-piercing check on a built structure (e.g. CHARMM-GUI
step3_input.pdb), to run BEFORE MD. For every ring (PHE/TYR/TRP/HIS/PRO and
sugar pyranose rings), reports atoms of OTHER residues within 1.6 A of the
ring centroid -- the signature of a bond threaded through the ring, which MD
cannot undo (see qc_ring_piercing_v2.py for the MM/GBSA symptom).
"""
import math, sys
from collections import defaultdict

RINGS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"], "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"], "HSD": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HSE": ["CG", "ND1", "CD2", "CE1", "NE2"], "HSP": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "PRO": ["N", "CA", "CB", "CG", "CD"],
}
SUGAR_RING = ["C1", "C2", "C3", "C4", "C5", "O5"]
CUT = 1.6

for pdb in sys.argv[1:]:
    res = defaultdict(dict)
    order = []
    for l in open(pdb):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        seg = l[72:76].strip()
        if seg.startswith(("SOLV", "IONS", "WAT", "TIP", "POT", "CLA", "SOD")):
            continue
        key = (seg, int(l[22:26]), l[17:21].strip())
        if key not in res:
            order.append(key)
        res[key][l[12:16].strip()] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    grid = defaultdict(list)
    for k in order:
        for n, p in res[k].items():
            grid[tuple(int(c // 3) for c in p)].append((k, n, p))
    hits = []
    for k in order:
        names = RINGS.get(k[2]) or (SUGAR_RING if k[0].startswith("CAR") else None)
        if not names or not all(n in res[k] for n in names):
            continue
        c = [sum(res[k][n][i] for n in names) / len(names) for i in range(3)]
        g = tuple(int(x // 3) for x in c)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for k2, n2, p2 in grid.get((g[0] + dx, g[1] + dy, g[2] + dz), []):
                        if k2 != k and math.dist(c, p2) < CUT:
                            hits.append(f"{k[0]}:{k[2]}{k[1]} ring pierced by {k2[0]}:{k2[2]}{k2[1]}:{n2} ({math.dist(c, p2):.2f} A)")
    print(f"{pdb}: {'OK' if not hits else f'{len(hits)} piercing(s)'}")
    for h in hits:
        print("   ", h)
