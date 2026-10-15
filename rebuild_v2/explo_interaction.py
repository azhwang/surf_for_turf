#!/usr/bin/env python3
"""explo_interaction.py -- EXPLORATORY (2026-10-15): species x stoichiometry interaction for
goldfish vs human, the form of the experimental claim (goldfish binds more strongly in 2:2:1,
no difference in 2:2:2). Two-way ANOVA (type III, sum contrasts) on replicate means, plus the
interaction contrast (g-h)_221 - (g-h)_222 with a Welch-Satterthwaite 95% CI.
Data: 1 ns rebuild set (mmpbsa_rebuild_results) and 50 ns extensions (ext_results)."""
import math, statistics as st
from pathlib import Path
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from scipy import stats

B = Path(__file__).resolve().parent
one = pd.read_excel(B.parent / "mmpbsa_rebuild_results/mmpbsa_v2_summary.xlsx", sheet_name="Replicates")
one = one[one.Species.isin(["goldfish", "human"])]
one = one.assign(stoich=one.Stoichiometry.str[:5], sp=one.Species, split=one.Split, dG=one["ΔG total"])
ext = pd.read_csv(B / "ext_results/replicate_table.csv")
ext = ext[ext.status == "complete"]
ext = pd.concat([ext.assign(split=s, dG=ext[f"dG_{s}"]) for s in ("iface1", "iface2", "whole")])
ext = ext.assign(stoich=ext.tmpl.map({"1e0o": "2:2:1", "1fq9": "2:2:2"}))

def contrast(d):
    cells = {(s, p): d[(d.stoich == s) & (d.sp == p)].dG.tolist() for s in ("2:2:1", "2:2:2") for p in ("goldfish", "human")}
    m = {k: st.fmean(v) for k, v in cells.items()}
    vn = {k: st.variance(v) / len(v) for k, v in cells.items()}
    c = (m["2:2:1", "goldfish"] - m["2:2:1", "human"]) - (m["2:2:2", "goldfish"] - m["2:2:2", "human"])
    se = math.sqrt(sum(vn.values()))
    df = sum(vn.values()) ** 2 / sum(vn[k] ** 2 / (len(cells[k]) - 1) for k in cells)
    h = stats.t.ppf(0.975, df) * se
    ns = {f"{s} {p[0]}": len(v) for (s, p), v in cells.items()}
    return m, c, h, ns

out = ["EXPLORATORY: species x stoichiometry interaction, goldfish vs human (not in the pre-registered plan)",
       "Interaction contrast = (goldfish - human) in 2:2:1  minus  (goldfish - human) in 2:2:2; the experimental",
       "claim predicts < 0. p_int from two-way ANOVA (type III).", ""]
for label, d in (("1 ns rebuild set", one), ("50 ns extensions, window 11-50 ns", ext)):
    out.append(f"== {label}")
    for s in ("iface1", "iface2", "whole"):
        x = d[d.split == s].dropna(subset=["dG"])
        fit = smf.ols("dG ~ C(sp, Sum) * C(stoich, Sum)", data=x).fit()
        a = anova_lm(fit, typ=3)
        m, c, h, ns = contrast(x)
        out.append(f"  {s:6s}  g-h 2:2:1 {m['2:2:1','goldfish']-m['2:2:1','human']:+6.1f}  g-h 2:2:2 {m['2:2:2','goldfish']-m['2:2:2','human']:+6.1f}"
                   f"  interaction {c:+6.1f} [95% CI {c-h:+.1f}, {c+h:+.1f}]  p_int = {a.loc['C(sp, Sum):C(stoich, Sum)','PR(>F)']:.2f}"
                   f"  p_species = {a.loc['C(sp, Sum)','PR(>F)']:.2f}   n {ns}")
    out.append("")
txt = "\n".join(out)
(B / "ext_results/explo_interaction.txt").write_text(txt + "\n")
print(txt)
