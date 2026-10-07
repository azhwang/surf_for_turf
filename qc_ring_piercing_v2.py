#!/usr/bin/env python3
"""
qc_ring_piercing_v2.py -- flag v2 MM/GBSA runs with build artifacts.

Scans every namdrepN/mmpbsa_v2/*_run/FINAL_DECOMP_MMPBSA.dat and reports
residues whose DELTA van der Waals is > +THRESH kcal/mol. A real interface
residue contributes a few kcal/mol at most; hundreds means a side chain is
threaded through an aromatic ring (ring piercing from the CHARMM-GUI build),
which MD cannot undo and which dominates DELTA TOTAL.

Usage:
    python3 qc_ring_piercing_v2.py                 # original builds -> qc_ring_piercing_v2.csv
    python3 qc_ring_piercing_v2.py --set rebuild   # rebuild_v2/     -> qc_ring_piercing_rebuild.csv
"""
import argparse
import csv
from pathlib import Path

THRESH = 20.0
# (glob, function(path parts) -> (stoich, species, rep, split))
LAYOUTS = {
    "original": ("all_metrics_analysis*/1*_*/*/namdrep*/mmpbsa_v2/*_run/FINAL_DECOMP_MMPBSA.dat",
                 lambda p: ("2:2:1" if p[0].endswith("221") else "2:2:2", p[1].split("_", 1)[1], p[3], p[5])),
    "rebuild": ("rebuild_v2/1*/*/sim/namdrep[0-9]/mmpbsa_v2/*_run/FINAL_DECOMP_MMPBSA.dat",
                lambda p: ("2:2:1" if p[1] == "1e0o" else "2:2:2", p[2], p[4], p[6])),
}
OUT = {"original": "qc_ring_piercing_v2.csv", "rebuild": "qc_ring_piercing_rebuild.csv"}

ap = argparse.ArgumentParser()
ap.add_argument("--set", choices=LAYOUTS, default="original")
args = ap.parse_args()
pattern, parse = LAYOUTS[args.set]

out = []
for f in sorted(Path('.').glob(pattern)):
    rows = list(csv.reader(open(f)))
    try:
        start = next(i for i, r in enumerate(rows) if r and r[0] == 'DELTAS:')
    except StopIteration:
        continue
    stoich, species, rep, split = parse(f.parts)
    for r in rows[start + 4:]:
        if not r or not r[0].strip() or r[0].startswith(('Sidechain', 'Backbone', 'Residue')):
            break
        if float(r[5]) > THRESH:
            out.append(dict(stoich=stoich, species=species, rep=rep, split=split.replace('_run', ''),
                            residue=r[0].strip(), location=r[1].strip(),
                            dvdw=round(float(r[5]), 1), dtotal=round(float(r[17]), 1)))
with open(OUT[args.set], 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['stoich', 'species', 'rep', 'split', 'residue', 'location', 'dvdw', 'dtotal'])
    w.writeheader(); w.writerows(out)
print(f"{len(out)} flagged residues -> {OUT[args.set]}")
