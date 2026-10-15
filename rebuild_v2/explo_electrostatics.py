#!/usr/bin/env python3
"""explo_electrostatics.py -- EXPLORATORY (2026-10-15). Static electrostatics of each FGF core
in the rebuilt starting structures (identical FGFR+heparin scaffold per template, so heparin
atoms sit at the same coordinates in every species and the numbers are directly comparable).
  - net charge of the FGF core (CHARMM charges from the PSF)
  - Debye-Hueckel screened Coulomb interaction of each FGF copy with heparin and with FGFR
    (eps = 80, 0.15 M -> Debye length 7.86 A; kcal/mol). Crude, comparative only.
  - goldfish vs human sequence differences, flagged by charge change and by contact
    (any heavy atom within 4.5 A of FGFR or heparin in the starting structure).
Writes ext_results/explo_electrostatics.txt."""
import math
from pathlib import Path
import numpy as np

B = Path(__file__).resolve().parent
SPECIES = ["human", "human_stab", "bovine_stab", "chinook_salmon", "goldfish", "japanese_ricefish", "zebrafish"]
AA = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', HSD='H', HSE='H', HSP='H',
          ILE='I', LEU='L', LYS='K', MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V')
CH = dict(K=1, R=1, D=-1, E=-1)
LD, EPS, KC = 7.86, 80.0, 332.0636


def load(d):
    psf = (d / "step3_input.psf").read_text().split("\n")
    i = next(k for k, l in enumerate(psf) if "!NATOM" in l)
    atoms = [l.split() for l in psf[i + 1:i + 1 + int(psf[i].split()[0])]]
    xyz = []
    for l in (d / "step3_input.pdb").read_text().split("\n"):
        if l.startswith(("ATOM", "HETATM")):
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    xyz = np.array(xyz)
    assert len(xyz) == len(atoms), d
    seg = np.array([a[1] for a in atoms]); resid = np.array([int(a[2]) for a in atoms])
    resn = np.array([a[3] for a in atoms]); name = np.array([a[4] for a in atoms]); q = np.array([float(a[6]) for a in atoms])
    return seg, resid, resn, name, q, xyz


def dh(qa, xa, qb, xb):
    r = np.sqrt(((xa[:, None, :] - xb[None, :, :]) ** 2).sum(-1))
    return float(KC / EPS * (qa[:, None] * qb[None, :] * np.exp(-r / LD) / r).sum())


def sysdir(tk, sp):
    for r in range(1, 6):
        d = B / tk / sp / "sim" / f"namdrep{r}"
        if (d / "step3_input.psf").exists() and (d / "step3_input.pdb").exists():
            return d
    for d in (B / tk / sp / "build",):
        if (d / "step3_input.psf").exists():
            return d
    return None


out = ["EXPLORATORY: static electrostatics of the FGF cores in the rebuilt starting structures",
       "E_hep / E_fgfr = Debye-Hueckel (eps 80, 0.15 M) interaction of each FGF copy with heparin / with both FGFRs, kcal/mol.",
       "Crude and comparative only (no desolvation, rigid start structure).", ""]
seqs, contact = {}, {}
for tk, label in (("1e0o", "2:2:1"), ("1fq9", "2:2:2")):
    out.append(f"== {label} ({tk})")
    out.append(f"  {'species':18s} {'core q':>7s} {'E_hep C':>8s} {'E_hep D':>8s} {'E_fgfr C':>9s} {'E_fgfr D':>9s}")
    for sp in SPECIES:
        d = sysdir(tk, sp)
        if d is None:
            out.append(f"  {sp:18s} (no structure locally)"); continue
        seg, resid, resn, name, q, xyz = load(d)
        hep = np.char.startswith(seg.astype(str), "CAR"); fgfr = np.isin(seg, ["PROA", "PROB"])
        row = []
        for c in ("PROC", "PROD"):
            m = seg == c
            row += [dh(q[m], xyz[m], q[hep], xyz[hep]), dh(q[m], xyz[m], q[fgfr], xyz[fgfr])]
        mC = seg == "PROC"
        out.append(f"  {sp:18s} {q[mC].sum():+7.1f} {row[0]:8.1f} {row[2]:8.1f} {row[1]:9.1f} {row[3]:9.1f}")
        # sequence + contacts of FGF C/D (positions in core numbering)
        rlist = list(dict.fromkeys(resid[mC]))
        seqs[(tk, sp)] = "".join(AA[resn[mC][resid[mC] == r][0]] for r in rlist)
        if sp in ("human", "goldfish"):
            heavy = ~np.char.startswith(name.astype(str), "H")
            env = {"FGFR": fgfr & heavy, "heparin": hep & heavy}
            cs = {}
            for c in ("PROC", "PROD"):
                rl = list(dict.fromkeys(resid[seg == c]))
                for k, r in enumerate(rl, 1):
                    a = (seg == c) & (resid == r) & heavy
                    for nm, e in env.items():
                        dmin = np.sqrt(((xyz[a][:, None] - xyz[e][None]) ** 2).sum(-1)).min()
                        if dmin < 4.5:
                            cs.setdefault(k, set()).add(f"{nm}({c[-1]})")
            contact[(tk, sp)] = cs
    out.append("")

for tk, label in (("1e0o", "2:2:1"), ("1fq9", "2:2:2")):
    h, g = seqs[(tk, "human")], seqs[(tk, "goldfish")]
    out.append(f"== human -> goldfish substitutions, {label} core numbering ({len(h)} residues; FGF2 numbering = core + {27 if tk=='1e0o' else 24})")
    dq = 0
    for k, (a, b) in enumerate(zip(h, g), 1):
        if a == b:
            continue
        dc = CH.get(b, 0) - CH.get(a, 0); dq += dc
        cts = sorted(contact[(tk, "human")].get(k, set()) | contact[(tk, "goldfish")].get(k, set()))
        out.append(f"  {a}{k}{b}  charge {dc:+d}  {'contacts ' + ', '.join(cts) if cts else ''}")
    out.append(f"  net charge change human -> goldfish (K/R/D/E only): {dq:+d}\n")
txt = "\n".join(out)
(B / "ext_results/explo_electrostatics.txt").write_text(txt + "\n")
print(txt)
