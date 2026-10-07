#!/usr/bin/env python3
"""
prepare_fgf_core.py -- step 1 of the consistent FGF rebuild.

For every species x template, take that species' FGF model (heavy atoms of
segment PROC in its existing CHARMM-GUI step1_pdbreader.pdb), trim it to the
residues that align to the crystal-resolved FGF core of the template, and
superpose the core onto each template FGF chain by core CA atoms only.

Why: the original builds placed the full-length AlphaFold FGF (including a
14-26 residue disordered N-terminal tail that no crystal structure resolves)
by whole-model alignment. The tail ended up threaded through FGFR rings in
several systems, and tail length differed by species. Trimming to the
template-resolved core and placing by the core makes every species
structurally equivalent to the crystal scaffold.

Scaffold: FGFR (PROA/PROB) + heparin (CAR*) from ONE clean build per template
(1FQ9: goldfish, heparin intact as 8+6 residues; 1E0O: chinook_salmon). The
scaffold coordinates are identical to the template frame (verified RMSD 0.000),
so the template FGF chains are used directly as placement targets:
    1FQ9: FGF chain A -> PROC (pairs with PROA), chain B -> PROD (pairs with PROB)
    1E0O: FGF chain A -> PROC,                  chain C -> PROD

Conserved interface side chains: for every core residue IDENTICAL to the
template and with a heavy atom within 6 A of the scaffold (FGFR/heparin), the
side chain is replaced by the crystal side chain, superposed on the model
residue's own N/CA/C (backbone and geometry stay the model's). AlphaFold
rotamers at conserved interface positions are not packed against the crystal
receptor/heparin and produced side-chain-into-ring overlaps (e.g. the conserved
heparin-binding Arg). Species-variable residues keep their AlphaFold side chains.

N-terminal anchor: the first 1-3 core residues of several models deviate
1.5-3 A (CA) from the crystal (they start the N-terminal strand, which contacts
FGFR; in the full-length models the now-removed tail pulled them out of
place). Contiguous N-terminal core residues with CA deviation > ANCHOR_DEV are
rebuilt on the crystal backbone (N/CA/C/O); side chain from the crystal if the
residue is conserved, otherwise the model side chain superposed on the crystal
N/CA/C. Internal loop deviations (FGF1 template vs FGF2-like models in 1E0O)
are real structural differences and are left alone.

Output per system: <out>/<tmpl>/<species>/fgf_PROC.pdb, fgf_PROD.pdb, core_report.txt
"""
import glob
import sys
import warnings
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices

warnings.filterwarnings("ignore")
BASE = Path(__file__).resolve().parent.parent          # ~/Desktop/Surf_for_Turf
OUT = Path(__file__).resolve().parent
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
         "HIS": "H", "HSD": "H", "HSE": "H", "HSP": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M",
         "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
TEMPLATES = {
    "1fq9": dict(pdb="1FQ9", root="all_metrics_analysis", fgf={"PROC": "A", "PROD": "B"}),
    "1e0o": dict(pdb="1E0O", root="all_metrics_analysis_221", fgf={"PROC": "A", "PROD": "C"}),
}
SCAFFOLD = {"1fq9": "all_metrics_analysis/1fq9_goldfish", "1e0o": "all_metrics_analysis_221/1e0o_chinook_salmon"}
GRAFT_CUTOFF = 6.0
ANCHOR_DEV = 1.5
CHARMM_NAME = {("ILE", "CD1"): "CD"}
SPECIES = ["human", "human_stab", "bovine_stab", "chinook_salmon", "goldfish", "japanese_ricefish", "zebrafish"]


def read_residues(path, select):
    """Ordered [(resid, resname, [(atomname, xyz), ...])] for atoms passing select(line)."""
    res, order = {}, []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")) or not select(l):
            continue
        rid = int(l[22:26])
        if rid not in res:
            res[rid] = (l[17:21].strip(), [])
            order.append(rid)
        res[rid][1].append((l[12:16].strip(), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return [(r, res[r][0], res[r][1]) for r in order]


def seq(residues):
    return "".join(THREE[rn] for _, rn, _ in residues)


def ca(residue):
    return next(xyz for n, xyz in residue[2] if n == "CA")


def kabsch(P, Q):
    """R, t minimising |P @ R.T + t - Q|."""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - pc @ R.T


def read_template_residues(path, chain):
    """Like read_residues but keeps only altloc ' ' or 'A' (first conformer)."""
    return read_residues(path, lambda l: l.startswith("ATOM") and l[21] == chain and l[16] in " A")


def graft_sidechain(model_res, tmpl_res):
    """Model residue with side chain (CB onward) taken from the template residue,
    superposed on N/CA/C. Returns None if the template residue is incomplete."""
    rid, rn, atoms = model_res
    m = dict(atoms)
    t = {CHARMM_NAME.get((rn if rn[:2] != "HS" else "HIS", n), n): xyz for n, xyz in tmpl_res[2]}
    side = [n for n, _ in atoms if n not in ("N", "CA", "C", "O", "OT1", "OT2")]
    if not side or not all(n in t for n in side) or not all(k in t for k in ("N", "CA", "C")):
        return None
    R, tv = kabsch(np.array([t[k] for k in ("N", "CA", "C")]), np.array([m[k] for k in ("N", "CA", "C")]))
    new = [(n, xyz if n not in side else t[n] @ R.T + tv) for n, xyz in atoms]
    return (rid, rn, new)


def anchor_residue(model_res, tmpl_res):
    """Model residue rebuilt on the template backbone: N/CA/C/O from the template;
    side chain from the template if conserved, else the model's, superposed on N/CA/C."""
    rid, rn, atoms = model_res
    t = {CHARMM_NAME.get((rn if rn[:2] != "HS" else "HIS", n), n): xyz for n, xyz in tmpl_res[2]}
    m = dict(atoms)
    conserved = THREE[tmpl_res[1]] == THREE[rn]
    R, tv = kabsch(np.array([m[k] for k in ("N", "CA", "C")]), np.array([t[k] for k in ("N", "CA", "C")]))
    new = []
    for n, xyz in atoms:
        if n in ("N", "CA", "C", "O"):
            new.append((n, t[n]))
        elif conserved and n in t:
            new.append((n, t[n]))
        else:
            new.append((n, xyz @ R.T + tv))
    return (rid, rn, new)


def main():
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -10, -0.5
    aligner.end_gap_score = 0
    summary = []
    for tk, t in TEMPLATES.items():
        tmpl_pdb = OUT / "templates" / f"{t['pdb']}.pdb"
        tchains = {seg: read_template_residues(tmpl_pdb, ch) for seg, ch in t["fgf"].items()}
        scaf_pdb = glob.glob(str(BASE / SCAFFOLD[tk] / "charmm-gui-*" / "step1_pdbreader.pdb"))[0]
        scaf_xyz = np.array([xyz for _, _, atoms in read_residues(
            scaf_pdb, lambda l: l[72:76].strip().startswith(("PROA", "PROB", "CAR"))
            and not l[12:16].strip().startswith("H")) for _, xyz in atoms])
        for sp in SPECIES:
            src = glob.glob(str(BASE / t["root"] / f"{tk}_{sp}" / "charmm-gui-*" / "step1_pdbreader.pdb"))[0]
            model = read_residues(src, lambda l: l[72:76].strip() == "PROC" and not l[12:16].strip().startswith("H")
                                  and l[12:16].strip() not in ("OT2",))
            mseq = seq(model)
            outdir = OUT / tk / sp
            outdir.mkdir(parents=True, exist_ok=True)
            report = [f"source model: {src} (segment PROC, heavy atoms)", f"full FGF: {len(model)} residues"]
            # both FGF copies get the same core: the model range aligned to BOTH template
            # FGF chains (1E0O chain C resolves one more N-terminal residue than chain A)
            all_pairs = {}
            for seg, tres in tchains.items():
                aln = aligner.align(seq(tres), mseq)[0]
                all_pairs[seg] = [(tres[i0 + k], model[j0 + k])
                                  for (i0, i1), (j0, _) in zip(*aln.aligned) for k in range(i1 - i0)]
            spans = [(min(model.index(m) for _, m in pr), max(model.index(m) for _, m in pr))
                     for pr in all_pairs.values()]
            first, last = max(sp[0] for sp in spans), min(sp[1] for sp in spans)
            core = model[first:last + 1]
            core_ids = {r[0] for r in core}
            for seg, tres in tchains.items():
                pairs = [(tr, m) for tr, m in all_pairs[seg] if m[0] in core_ids]
                P = np.array([ca(m) for _, m in pairs])
                Q = np.array([ca(tr) for tr, _ in pairs])
                R, tvec = kabsch(P, Q)
                rmsd = np.sqrt(((P @ R.T + tvec - Q) ** 2).sum(1).mean())
                ident = sum(THREE[tr[1]] == THREE[m[1]] for tr, m in pairs)
                placed = [(rid, rn, [(nm, xyz @ R.T + tvec) for nm, xyz in atoms]) for rid, rn, atoms in core]
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
                with open(outdir / f"fgf_{seg}.pdb", "w") as fh:
                    n = 0
                    for rid, rn, atoms in placed:
                        for name, xyz in atoms:
                            n += 1
                            x, y, z = xyz
                            fh.write(f"ATOM  {n:5d} {name:<4s} {rn:<4s}{'X':1s}{rid:4d}    "
                                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00      {seg:<4s}\n")
                    fh.write("END\n")
                report.append(
                    f"{seg} <- template {t['pdb']} chain {t['fgf'][seg]}: core = model residues "
                    f"{core[0][0]}-{core[-1][0]} ({len(core)} res; trimmed N {first}, C {len(model) - last - 1}); "
                    f"first {core[0][1]}, last {core[-1][1]}; {len(pairs)} aligned CA, identity "
                    f"{100 * ident / len(pairs):.0f}%, core CA RMSD to template {rmsd:.2f} A; "
                    f"N-terminal residues anchored to crystal backbone: {' '.join(anchored) or 'none'}; "
                    f"{len(grafted)} conserved interface side chains from crystal: {' '.join(grafted)}")
            (outdir / "core_report.txt").write_text("\n".join(report) + "\n")
            summary.append(f"{tk} {sp:18s} " + report[2].split(": core = ")[1])
    print("\n".join(summary))


if __name__ == "__main__":
    sys.exit(main())
