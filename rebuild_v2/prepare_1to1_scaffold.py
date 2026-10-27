#!/usr/bin/env python3
"""
prepare_1to1_scaffold.py -- repaired FGFR scaffolds for the 1:1 FGF-FGFR systems.

The dedicated 1:1 crystal structures are the right templates for a 1:1 simulation, but their
receptor chains are incomplete:

    1EVT chain C (FGFR1)   147-359, 191 residues, missing 294-303 and 312-323
    1DJS chain A (FGFR2)   147-362, 206 residues, missing 297-306

Both gaps sit in the same D3 loop region, and 1FQ9 chain C resolves all of it (FGFR1 149-359,
211 residues, no gaps). FGFR1 and FGFR2 are collinear here and 1FQ9's numbering matches 1EVT
exactly, while 1DJS runs one residue ahead (FGFR1 i <-> FGFR2 i+1); across the FGFR2 gap the
donor loop is 8/10 identical (FGFR1 296-305 NGSKIGPDNL vs FGFR2 297-306 NGSKYGPDGL).

So each gap is closed with the structurally equivalent 1FQ9 loop, superposed on the host's own
flanking residues (ten on each side) so the two junctions close on the host backbone. Residues
identical between donor and host keep the donor side chain; the two substituted FGFR2 positions
and anything the donor lacks are left to psfgen's guesscoord. The same donor also repairs the
19 side chains 1EVT truncates at disordered surface Lys/Arg/Glu.

1DJS selenomethionines are converted back to methionine (MSE -> MET, SE -> SD).

Output: <out>/scaffold/receptor_<rec>.pdb (segid PROA) + scaffold_report.txt
"""
import sys
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
TMPL = BASE / "templates"
OUT = BASE / "onetoone" / "scaffold"

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
         "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "MSE": "M",
         "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
NHEAVY = {"ALA": 5, "ARG": 11, "ASN": 8, "ASP": 8, "CYS": 6, "GLN": 9, "GLU": 9, "GLY": 4,
          "HIS": 10, "ILE": 8, "LEU": 8, "LYS": 9, "MET": 8, "PHE": 11, "PRO": 7, "SER": 6,
          "THR": 7, "TRP": 14, "TYR": 12, "VAL": 7}
BB = ("N", "CA", "C", "O")
# receptor label -> (host pdb, host chain, seqres index of the first numbered residue,
#                   host resid - donor resid)
HOST = {"fgfr1": ("1EVT", "C", 6, 0), "fgfr2": ("1DJS", "A", 0, 1)}
DONOR = ("1FQ9", "C")
FLANK = 6


def read_chain(pdb, chain):
    """{resid: (resname, {atomname: xyz})} for the first conformer of one protein chain."""
    res = {}
    for l in (TMPL / f"{pdb}.pdb").read_text().splitlines():
        if l[:6] not in ("ATOM  ", "HETATM") or l[21] != chain or l[16] not in " A":
            continue
        rn, nm = l[17:20], l[12:16].strip()
        if rn not in THREE or nm.startswith("H") or nm in ("OXT",):
            continue
        if rn == "MSE":
            rn, nm = "MET", ("SD" if nm == "SE" else nm)
        r = int(l[22:26])
        if r not in res:
            res[r] = (rn, {})
        res[r][1][nm] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return res


def seqres(pdb, chain):
    s = []
    for l in (TMPL / f"{pdb}.pdb").read_text().splitlines():
        if l.startswith("SEQRES") and l[11] == chain:
            s += l[19:].split()
    return ["MET" if x == "MSE" else x for x in s]


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - pc @ R.T


def gaps(res):
    k = sorted(res)
    return [(a + 1, b - 1) for a, b in zip(k, k[1:]) if b > a + 1]


def build(rec, report):
    pdb, chain, sr_off, off = HOST[rec]
    host = read_chain(pdb, chain)
    donor = read_chain(*DONOR)
    sr = seqres(pdb, chain)
    lo, hi = min(host), max(host)
    report.append(f"== {rec}: host {pdb} chain {chain} {lo}-{hi} ({len(host)} res), "
                  f"donor {DONOR[0]} chain {DONOR[1]} ({len(donor)} res), "
                  f"host resid = donor resid + {off}")

    for g0, g1 in gaps(host):
        # Two separate superpositions, one per junction: a single rigid placement of the donor
        # loop cannot close both ends when the host and donor flanks differ in conformation
        # (1.7 A backbone RMSD for FGFR2). Taking the first half of the loop from the fit on
        # the N-side flank and the second half from the fit on the C-side flank leaves both
        # junctions on the host backbone and puts the one remaining break in the middle of the
        # loop, where the subsequent restrained minimisation can close it.
        mid = (g0 + g1) // 2
        fits, info = {}, []
        for side, rng in (("N", range(g0 - FLANK, g0)), ("C", range(g1 + 1, g1 + 1 + FLANK))):
            anc = [r for r in rng if r in host and r - off in donor]
            P = np.array([donor[r - off][1][a] for r in anc for a in BB])
            Q = np.array([host[r][1][a] for r in anc for a in BB])
            fits[side] = kabsch(P, Q)
            R, t = fits[side]
            info.append(f"{side}-side {len(anc)} res, backbone RMSD "
                        f"{np.sqrt(((P @ R.T + t - Q) ** 2).sum(1).mean()):.2f} A")
        copied = []
        for r in range(g0, g1 + 1):
            d = donor.get(r - off)
            want = sr[r - lo + sr_off]
            if d is None:
                report.append(f"   !! donor has no residue for {want}{r}")
                continue
            R, t = fits["N" if r <= mid else "C"]
            same = d[0] == want
            keep = set(d[1]) if same else (set(BB) | {"CB"}) & set(d[1])
            if want == "GLY":
                keep -= {"CB"}
            host[r] = (want, {a: d[1][a] @ R.T + t for a in keep})
            copied.append(f"{want}{r}" + ("" if same else f"/{d[0]}"))
        report.append(f"   gap {g0}-{g1}: {'; '.join(info)}; donor split at {mid}/{mid + 1}; "
                      f"grafted {' '.join(copied)}")

    # junction geometry
    k = sorted(host)
    bad = [(a, np.linalg.norm(host[a][1]["C"] - host[a + 1][1]["N"]))
           for a in k[:-1] if np.linalg.norm(host[a][1]["C"] - host[a + 1][1]["N"]) > 1.6]
    report.append("   peptide C-N bonds all <= 1.6 A" if not bad else
                  "   stretched C-N: " + ", ".join(f"{a}-{a+1} {d:.2f} A" for a, d in bad))
    assert not gaps(host), f"{rec} still has gaps: {gaps(host)}"

    # repair truncated side chains from the donor where the residue is identical and complete:
    # host backbone throughout, donor side chain superposed on the host residue's own N/CA/C
    fixed = []
    for r, (rn, at) in sorted(host.items()):
        if len(at) >= NHEAVY[rn] or r - off not in donor:
            continue
        d = donor[r - off]
        if d[0] != rn or len(d[1]) < NHEAVY[rn]:
            continue
        Rl, tl = kabsch(np.array([d[1][a] for a in ("N", "CA", "C")]),
                        np.array([at[a] for a in ("N", "CA", "C")]))
        host[r] = (rn, {a: (at[a] if a in BB else d[1][a] @ Rl.T + tl) for a in d[1]})
        fixed.append(f"{rn}{r}")
    left = [f"{rn}{r}" for r, (rn, at) in sorted(host.items()) if len(at) < NHEAVY[rn]]
    report.append(f"   side chains completed from donor ({len(fixed)}): {' '.join(fixed) or 'none'}")
    report.append(f"   still truncated, left to guesscoord ({len(left)}): {' '.join(left) or 'none'}")

    disu = [(a, b) for a in host for b in host
            if a < b and host[a][0] == host[b][0] == "CYS" and "SG" in host[a][1]
            and "SG" in host[b][1] and np.linalg.norm(host[a][1]["SG"] - host[b][1]["SG"]) < 2.5]
    report.append(f"   disulfides: {disu}")

    path = OUT / f"receptor_{rec}.pdb"
    with open(path, "w") as fh:
        n = 0
        for r, (rn, at) in sorted(host.items()):
            for a in sorted(at, key=lambda x: (BB.index(x) if x in BB else len(BB), x)):
                n += 1
                x, y, z = at[a]
                fh.write(f"ATOM  {n:5d} {a:<4s} {rn:<4s}X{r:4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00      PROA\n")
        fh.write("END\n")
    report.append(f"   wrote {path}: {len(host)} residues {min(host)}-{max(host)}, {n} heavy atoms")
    (OUT / f"disu_{rec}.txt").write_text("\n".join(f"{a} {b}" for a, b in disu) + "\n")
    return path


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = []
    for rec in HOST:
        build(rec, report)
    (OUT / "scaffold_report.txt").write_text("\n".join(report) + "\n")
    print("\n".join(report))


if __name__ == "__main__":
    sys.exit(main())
