#!/usr/bin/env python3
"""explo_summary.py -- EXPLORATORY (2026-10-15), after analyze_exploratory.sh.
1. FGF-heparin (hep) and FGF-FGFR-without-heparin (fgfr) MM/GBSA energies, goldfish vs human
   per stoichiometry (Welch) and the species x stoichiometry interaction.
2. Per-residue decomposition (whole split, idecomp=1, 40 frames over 11-50 ns) of FGF residues:
   goldfish - human per core position and FGF copy; the total gap split into substituted vs
   conserved positions.
Writes ext_results/explo_summary.txt and explo_decomp_<tmpl>.csv."""
import math, re, statistics as st
from pathlib import Path
from scipy import stats

B = Path(__file__).resolve().parent
SEL = {"1e0o/human": [1, 2, 4], "1e0o/goldfish": [1, 2], "1fq9/human": [1, 2, 3], "1fq9/goldfish": [2, 3, 5]}
STO = {"1e0o": "2:2:1", "1fq9": "2:2:2"}
OFF = {"1e0o": 27, "1fq9": 24}   # core position -> FGF2 (155 aa) numbering
AA = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', HSD='H', HSE='H', HSP='H',
          ILE='I', LEU='L', LYS='K', MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V')


def dtot(f):
    m = re.search(r"Differences \(Complex - Receptor - Ligand\):.*?DELTA TOTAL\s+(-?\d+\.\d+)", Path(f).read_text(), re.S)
    return float(m.group(1)) if m else float("nan")


def welch(a, b):
    d = st.fmean(a) - st.fmean(b)
    va, vb = st.variance(a) / len(a), st.variance(b) / len(b)
    se = math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1)) if se else 1
    return d, stats.t.ppf(0.975, df) * se, stats.ttest_ind(a, b, equal_var=False).pvalue


def decomp(f, ncore):
    out = {}
    for l in Path(f).read_text().split("\n"):
        p = l.split(",")
        if len(p) > 17 and p[1].startswith("L "):
            k = int(p[1].split()[-1]); res = p[1].split()[1]
            c, pos = ("C", k) if k <= ncore else ("D", k - ncore)
            out[(c, pos)] = (AA.get(res, res), float(p[17]))
    return out


rows = []
L = ["EXPLORATORY analyses (not pre-registered). Replicate means over extension segments 11-50.", ""]
for split, lab in (("hep", "FGF-heparin (FGFR stripped)"), ("fgfr", "FGF-FGFR, heparin stripped"), ("whole", "whole assembly (reference)")):
    L.append(f"== {lab}")
    gaps = {}
    for tk in ("1e0o", "1fq9"):
        v = {}
        for sp in ("goldfish", "human"):
            v[sp] = []
            for r in SEL[f"{tk}/{sp}"]:
                f = B / tk / sp / "sim" / f"namdrep{r}" / f"ext_analysis/{split}_run/FINAL_RESULTS_{split}.dat"
                if f.exists():
                    v[sp].append(dtot(f)); rows.append((split, tk, sp, r, v[sp][-1]))
        if min(len(v["goldfish"]), len(v["human"])) < 2:
            L.append(f"  {STO[tk]}: goldfish {v['goldfish']} human {v['human']} (n<2)"); continue
        d, h, p = welch(v["goldfish"], v["human"])
        gaps[tk] = (v["goldfish"], v["human"])
        L.append(f"  {STO[tk]}  goldfish {st.fmean(v['goldfish']):8.1f}  human {st.fmean(v['human']):8.1f}  diff {d:+6.1f} [95% CI {d-h:+.1f}, {d+h:+.1f}]  p = {p:.2f}"
                 f"   reps g {[round(x,1) for x in v['goldfish']]} h {[round(x,1) for x in v['human']]}")
    if len(gaps) == 2:
        (g1, h1), (g2, h2) = gaps["1e0o"], gaps["1fq9"]
        c = (st.fmean(g1) - st.fmean(h1)) - (st.fmean(g2) - st.fmean(h2))
        se2 = [st.variance(x) / len(x) for x in (g1, h1, g2, h2)]; ns = [len(x) for x in (g1, h1, g2, h2)]
        df = sum(se2) ** 2 / sum(s * s / (n - 1) for s, n in zip(se2, ns))
        t = c / math.sqrt(sum(se2)); hw = stats.t.ppf(0.975, df) * math.sqrt(sum(se2))
        L.append(f"  interaction (g-h)2:2:1 - (g-h)2:2:2 = {c:+.1f} [95% CI {c-hw:+.1f}, {c+hw:+.1f}]  p = {2*stats.t.sf(abs(t), df):.2f}")
    L.append("")

with open(B / "ext_results/explo_energies.csv", "w") as fh:
    fh.write("split,tmpl,species,rep,dG\n")
    for r in rows:
        fh.write(",".join(map(str, r[:4])) + f",{r[4]:.2f}\n")

# --- decomposition ---
for tk in ("1e0o", "1fq9"):
    ncore = 127 if tk == "1e0o" else 129
    per = {"goldfish": [], "human": []}
    for sp in per:
        for r in SEL[f"{tk}/{sp}"]:
            f = B / tk / sp / "sim" / f"namdrep{r}" / "ext_analysis/decomp_run/FINAL_DECOMP_decomp.dat"
            if f.exists():
                per[sp].append(decomp(f, ncore))
    if min(len(per["goldfish"]), len(per["human"])) < 2:
        L.append(f"== decomposition {STO[tk]}: not enough replicates yet"); continue
    keys = sorted(per["human"][0])
    tab = []
    for k in keys:
        a = [x[k][1] for x in per["goldfish"]]; b = [x[k][1] for x in per["human"]]
        rg, rh = per["goldfish"][0][k][0], per["human"][0][k][0]
        d = st.fmean(a) - st.fmean(b)
        p = stats.ttest_ind(a, b, equal_var=False).pvalue if (st.variance(a) + st.variance(b)) > 0 else 1.0
        tab.append(dict(copy=k[0], pos=k[1], fgf2=k[1] + OFF[tk], human=rh, goldfish=rg, sub=rh != rg,
                        h=st.fmean(b), g=st.fmean(a), d=d, p=p))
    with open(B / f"ext_results/explo_decomp_{tk}.csv", "w") as fh:
        fh.write("copy,core_pos,fgf2_pos,human_res,goldfish_res,substituted,human_mean,goldfish_mean,g_minus_h,welch_p\n")
        for t in tab:
            fh.write(f"{t['copy']},{t['pos']},{t['fgf2']},{t['human']},{t['goldfish']},{int(t['sub'])},{t['h']:.2f},{t['g']:.2f},{t['d']:.2f},{t['p']:.3f}\n")
    L.append(f"== per-residue decomposition, {STO[tk]} ({tk}), FGF residues, whole-assembly split (n g {len(per['goldfish'])}, h {len(per['human'])})")
    for c in ("C", "D"):
        t = [x for x in tab if x["copy"] == c]
        s_sub = sum(x["d"] for x in t if x["sub"]); s_con = sum(x["d"] for x in t if not x["sub"])
        L.append(f"  FGF {c}: total g-h {s_sub+s_con:+.1f} = substituted positions {s_sub:+.1f} + conserved positions {s_con:+.1f}")
    L.append("  largest |goldfish - human| (kcal/mol; FGF2 numbering; * = substituted):")
    for x in sorted(tab, key=lambda x: -abs(x["d"]))[:15]:
        lab = f"{x['human']}{x['fgf2']}{x['goldfish'] if x['sub'] else ''}{'*' if x['sub'] else ''}"
        L.append(f"    FGF {x['copy']} {lab:8s} human {x['h']:7.2f}  goldfish {x['g']:7.2f}  diff {x['d']:+6.2f}  p = {x['p']:.2f}")
    L.append("")
txt = "\n".join(L)
(B / "ext_results/explo_summary.txt").write_text(txt + "\n")
print(txt)
