#!/usr/bin/env python3
"""
prepare_1to1_fgf.py -- place each ortholog's FGF2 core on the 1:1 crystal FGF chain.

Same protocol as prepare_fgf_core.py (trim the AlphaFold model to the template-resolved core,
superpose by core CA, rebuild deviating N-terminal residues on the crystal backbone, take
crystal rotamers for conserved residues at the interface), with two changes for the 1:1 sets:

  * the placement targets are the 1:1 crystals' FGF chains, 1EVT chain A (the copy that packs
    against receptor chain C: 307 heavy-atom contacts within 4.5 A, against 9 for chain D) and
    1DJS chain B;
  * the scaffold used to decide which side chains sit at the interface is the repaired receptor
    alone, since these systems have no heparin and no second FGF2.

The FGF2 model for both receptor sets is taken from the 2:2:2 (1FQ9) CHARMM-GUI input of that
ortholog, so the two 1:1 sets differ only in the receptor and its crystal frame.

Output: <out>/onetoone/<rec>/<species>/fgf_PROC.pdb + core_report.txt
"""
import glob
import sys
import warnings
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices

from prepare_fgf_core import (ANCHOR_DEV, GRAFT_CUTOFF, SPECIES, THREE, anchor_residue, ca,
                              graft_sidechain, kabsch, read_residues, read_template_residues, seq)

warnings.filterwarnings("ignore")
BASE = Path(__file__).resolve().parent
ROOT = BASE.parent
OUT = BASE / "onetoone"
TEMPLATES = {"fgfr1": dict(pdb="1EVT", chain="A"), "fgfr2": dict(pdb="1DJS", chain="B")}
# Both crystals resolve more of the FGF1 N-terminus than 1E0O does (1EVT chain A 8-138, 1DJS
# chain B 6-140, against 1E0O chain A 10-138). Restricting both to the 1E0O span makes the
# placement target identical for the two receptor sets and gives the same FGF2 core as the
# 2:2:1 systems, so the 1:1 numbers are directly comparable with Table 1. The residues dropped
# are the start of the FGF1 N-terminal segment, where FGF1 and the FGF2 models diverge most:
# keeping them raised the core CA RMSD to 1.1-2.3 A and forced four N-terminal residues per
# system onto the crystal backbone.
MODEL_ROOT, MODEL_TMPL = "all_metrics_analysis", "1fq9"
SPAN = (10, 138)        # FGF1 residues, the span 1E0O resolves


def main():
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -10, -0.5
    aligner.end_gap_score = 0
    summary = []
    for rec, t in TEMPLATES.items():
        tres = [r for r in read_template_residues(BASE / "templates" / f"{t['pdb']}.pdb",
                                                  t["chain"]) if SPAN[0] <= r[0] <= SPAN[1]]
        scaf_xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                             for l in open(OUT / "scaffold" / f"receptor_{rec}.pdb")
                             if l.startswith("ATOM")])
        for sp in SPECIES:
            src = glob.glob(str(ROOT / MODEL_ROOT / f"{MODEL_TMPL}_{sp}" / "charmm-gui-*" /
                                "step1_pdbreader.pdb"))[0]
            model = read_residues(src, lambda l: l[72:76].strip() == "PROC"
                                  and not l[12:16].strip().startswith("H")
                                  and l[12:16].strip() not in ("OT2",))
            aln = aligner.align(seq(tres), seq(model))[0]
            pairs = [(tres[i0 + k], model[j0 + k])
                     for (i0, i1), (j0, _) in zip(*aln.aligned) for k in range(i1 - i0)]
            first = min(model.index(m) for _, m in pairs)
            last = max(model.index(m) for _, m in pairs)
            core = model[first:last + 1]
            core_ids = {r[0] for r in core}
            pairs = [(tr, m) for tr, m in pairs if m[0] in core_ids]
            P = np.array([ca(m) for _, m in pairs])
            Q = np.array([ca(tr) for tr, _ in pairs])
            R, tvec = kabsch(P, Q)
            rmsd = np.sqrt(((P @ R.T + tvec - Q) ** 2).sum(1).mean())
            ident = sum(THREE[tr[1]] == THREE[m[1]] for tr, m in pairs)
            placed = [(rid, rn, [(nm, xyz @ R.T + tvec) for nm, xyz in atoms])
                      for rid, rn, atoms in core]
            tmpl_of = {m[0]: tr for tr, m in pairs}
            anchored = []
            for i, res in enumerate(placed):
                tr = tmpl_of.get(res[0])
                if tr is None or np.linalg.norm(ca(res) - ca(tr)) <= ANCHOR_DEV:
                    break
                placed[i] = anchor_residue(res, tr)
                anchored.append(f"{res[1]}{res[0]}")
            grafted = []
            for i, res in enumerate(placed):
                tr = tmpl_of.get(res[0])
                if tr is None or THREE[tr[1]] != THREE[res[1]] or f"{res[1]}{res[0]}" in anchored:
                    continue
                xyz = np.array([a for _, a in res[2]])
                if np.min(np.linalg.norm(xyz[:, None, :] - scaf_xyz[None, :, :], axis=2)) > GRAFT_CUTOFF:
                    continue
                g = graft_sidechain(res, tr)
                if g:
                    placed[i] = g
                    grafted.append(f"{res[1]}{res[0]}")
            outdir = OUT / rec / sp
            outdir.mkdir(parents=True, exist_ok=True)
            with open(outdir / "fgf_PROC.pdb", "w") as fh:
                n = 0
                for rid, rn, atoms in placed:
                    for name, xyz in atoms:
                        n += 1
                        x, y, z = xyz
                        fh.write(f"ATOM  {n:5d} {name:<4s} {rn:<4s}X{rid:4d}    "
                                 f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00      PROC\n")
                fh.write("END\n")
            line = (f"PROC <- template {t['pdb']} chain {t['chain']}: core = model residues "
                    f"{core[0][0]}-{core[-1][0]} ({len(core)} res; trimmed N {first}, "
                    f"C {len(model) - last - 1}); {len(pairs)} aligned CA, identity "
                    f"{100 * ident / len(pairs):.0f}%, core CA RMSD to template {rmsd:.2f} A; "
                    f"N-terminal residues anchored: {' '.join(anchored) or 'none'}; "
                    f"{len(grafted)} conserved interface side chains from crystal: "
                    f"{' '.join(grafted) or 'none'}")
            (outdir / "core_report.txt").write_text(
                f"source model: {src} (segment PROC, heavy atoms)\n"
                f"full FGF: {len(model)} residues\n{line}\n")
            summary.append(f"{rec} {sp:18s} " + line.split(": core = ")[1])
    print("\n".join(summary))


if __name__ == "__main__":
    sys.exit(main())
