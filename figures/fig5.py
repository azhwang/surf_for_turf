#!/usr/bin/env python3
"""
fig5.py -- Fig. 5: the 50 ns extensions and the Arg90 finding.

  a  whole-assembly MM/GBSA over the 11-51 ns window, both assemblies (pre-specified test)
  b  heparin-site salt bridges over time, 2:2:1
  c  Arg90-heparin contacts per FGF2 copy
  d  Arg90 per-residue energy contribution per FGF2 copy
  e  per-residue goldfish - human difference along the sequence, 2:2:1

Writes figures/out/fig5.{svg,png}.
"""
import json
import math
import re

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

import figstyle as F

GROUPS = [("1e0o", "human"), ("1e0o", "goldfish"), ("1fq9", "human"), ("1fq9", "goldfish")]
GX = [0, 1, 2.5, 3.5]
R90POS = {"1e0o": 63, "1fq9": 66}       # core position of FGF2 Arg90 per template


def welch(a, b):
    a, b = list(a), list(b)
    d = float(np.mean(a) - np.mean(b))
    va, vb = np.var(a, ddof=1) / len(a), np.var(b, ddof=1) / len(b)
    se = math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    return d, float(stats.t.ppf(0.975, df) * se), float(stats.ttest_ind(a, b, equal_var=False).pvalue)


def r90_contacts(tk, sp, r):
    f = F.RB / tk / sp / "sim" / f"namdrep{r}" / "ext_analysis/r90.txt"
    if not f.exists():
        return None
    L = f.read_text().splitlines()
    return float(L[1].split()[3]), float(L[2].split()[3])


def r90_energy(tk, sp, r):
    f = F.RB / tk / sp / "sim" / f"namdrep{r}" / "ext_analysis/decomp_run/FINAL_DECOMP_decomp.dat"
    if not f.exists():
        return None
    pos, n = R90POS[tk], F.NCORE[tk]
    want = {pos: None, n + pos: None}
    for l in f.read_text().splitlines():
        p = l.split(",")
        if len(p) > 17 and p[1].startswith("L "):
            k = int(p[1].split()[-1])
            if k in want:
                want[k] = float(p[17])
    return (want[pos], want[n + pos]) if all(v is not None for v in want.values()) else None


def stoich_axis(ax):
    ax.set_xticks(GX, ["Human", "Goldfish", "Human", "Goldfish"], rotation=30, ha="right")
    ax.text(0.5, -0.46, "2:2:1", transform=ax.get_xaxis_transform(), ha="center",
            fontsize=7, color=F.INK)
    ax.text(3.0, -0.46, "2:2:2", transform=ax.get_xaxis_transform(), ha="center",
            fontsize=7, color=F.INK)
    ax.set_xlim(-0.55, 4.05)


def panel_energy(ax, ext, out):
    for (tk, sp), x in zip(GROUPS, GX):
        ys = ext[(ext.tmpl == tk) & (ext.sp == sp)].dG_whole.values
        F.dots(ax, x, ys, F.SP_COLOR[sp])
        F.meanbar(ax, x, ys)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.26 * (hi - lo))
    for tk, x0 in (("1e0o", 0.5), ("1fq9", 3.0)):
        g = ext[(ext.tmpl == tk) & (ext.sp == "goldfish")].dG_whole
        h = ext[(ext.tmpl == tk) & (ext.sp == "human")].dG_whole
        d, ci, p = welch(g, h)
        out[f"{tk}_whole"] = dict(diff=round(d, 1), ci=round(ci, 1), p=round(p, 2))
        ax.text(x0, 0.995, f"$\\Delta$ {d:+.0f} ± {ci:.0f}\n$p$ = {p:.2f}",
                transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6, color=F.INK)
    ax.set_ylabel("$\\Delta G_{bind}$, whole assembly\n(kcal mol$^{-1}$)")
    stoich_axis(ax)
    F.style(ax)


def panel_timeseries(ax, out):
    for sp in ("human", "goldfish"):
        for i, r in enumerate(F.EXT[f"1e0o/{sp}"]):
            d = F.timeseries("1e0o", sp, r)
            if d is None:
                continue
            ax.plot(d[9:, 1], F.running(d[:, 2]), color=F.SP_COLOR[sp], lw=0.8,
                    ls=["-", "--", ":"][i], alpha=0.95)
    ax.axvspan(1, 11, color=F.PALE, zorder=0)
    ax.axhline(14, color=F.MUTED, lw=0.6, ls="-.")
    ax.text(50.6, 14.4, "tight state $\\geq$ 14", ha="right", va="bottom", fontsize=6, color=F.MUTED)
    ax.text(6, ax.get_ylim()[1], "excluded", ha="center", va="top", fontsize=6, color=F.MUTED)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Heparin-site salt bridges\nFGF2 copy C")
    ax.set_xlim(1, 51)
    ax.legend(handles=[Line2D([], [], color=F.H, lw=1, label="Human (reps 1, 2, 4)"),
                       Line2D([], [], color=F.G, lw=1, label="Goldfish (reps 1–3)")],
              loc="lower right", ncol=2)
    ax.set_title("2:2:1, 1 ns running mean", color=F.INK)
    F.style(ax)


def panel_r90(ax, kind, out):
    vals = {}
    for (tk, sp), x in zip(GROUPS, GX):
        v = [(r90_contacts if kind == "contacts" else r90_energy)(tk, sp, r)
             for r in F.EXT[f"{tk}/{sp}"]]
        v = [q for q in v if q]
        vals[f"{tk}/{sp}"] = [[round(a, 1), round(b, 1)] for a, b in v]
        F.dots(ax, x - 0.14, [q[0] for q in v], F.SP_COLOR[sp], "o", jitter=0.05)
        F.dots(ax, x + 0.14, [q[1] for q in v], F.SP_COLOR[sp], "^", jitter=0.05)
    out[f"r90_{kind}"] = vals
    if kind == "contacts":
        ax.axhspan(5.5, 9.6, color=F.PALE, zorder=0)
        ax.axhline(7.5, color=F.MUTED, lw=0.6, ls="--")
        ax.text(4.0, 8.4, "engagement\nthreshold", ha="right", va="bottom", fontsize=6, color=F.MUTED)
        ax.set_ylabel("Arg90–heparin contacts\n(heavy atoms within 4.5 Å)")
        ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=F.INK, ms=4,
                                  label="FGF2 copy C"),
                           Line2D([], [], marker="^", ls="", mfc="white", mec=F.INK, ms=4,
                                  label="FGF2 copy D")], loc="upper left", handletextpad=0.2)
    else:
        ax.axhline(-5, color=F.MUTED, lw=0.6, ls="--")
        ax.set_ylabel("Arg90 contribution\n(kcal mol$^{-1}$)")
    stoich_axis(ax)
    F.style(ax)


def panel_decomp(ax, out):
    d = pd.read_csv(F.RB / "ext_results/explo_decomp_1e0o.csv")
    tot = d.groupby(["fgf2_pos", "substituted"]).g_minus_h.sum().reset_index()
    sub = tot[tot.substituted == 1]
    ax.axhline(0, color=F.INK, lw=0.5)
    ax.vlines(tot.fgf2_pos, 0, tot.g_minus_h, color=F.MUTED, lw=0.8)
    ax.vlines(sub.fgf2_pos, 0, sub.g_minus_h, color=F.G, lw=1.3)
    top = tot.reindex(tot.g_minus_h.abs().sort_values(ascending=False).index).head(4)
    lo, hi = tot.g_minus_h.min(), tot.g_minus_h.max()
    ax.set_ylim(lo - 0.22 * (hi - lo), hi + 0.30 * (hi - lo))
    seen = []
    for _, r in top.iterrows():
        lab = "Arg90" if int(r.fgf2_pos) == 90 else f"{int(r.fgf2_pos)}"
        dx = 9 if any(abs(r.fgf2_pos - s) < 6 for s in seen) else 0
        seen.append(r.fgf2_pos)
        ax.annotate(lab, (r.fgf2_pos, r.g_minus_h), textcoords="offset points",
                    xytext=(dx, -8 if r.g_minus_h < 0 else 4), ha="center", fontsize=6, color=F.INK)
    out["decomp_top"] = [[int(r.fgf2_pos), round(float(r.g_minus_h), 1)] for _, r in top.iterrows()]
    ax.set_xlabel("FGF2 residue")
    ax.set_ylabel("Goldfish − human\n(kcal mol$^{-1}$, both copies)")
    ax.set_title("2:2:1 per-residue difference", color=F.INK)
    ax.legend(handles=[Line2D([], [], color=F.G, lw=1.3, label="substituted"),
                       Line2D([], [], color=F.MUTED, lw=0.9, label="conserved")],
              loc="lower right", ncol=2)
    F.style(ax)


def main():
    out = {}
    ext = F.extensions()
    fig = F.plt.figure(figsize=(F.DOUBLE, 152 * F.MM))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.0, 0.95], hspace=1.15, wspace=0.60)

    ax = fig.add_subplot(gs[0, 0]); panel_energy(ax, ext, out); F.letter(ax, "a", x=-0.56)
    ax = fig.add_subplot(gs[0, 1:]); panel_timeseries(ax, out); F.letter(ax, "b", x=-0.095)
    ax = fig.add_subplot(gs[1, 0]); panel_r90(ax, "contacts", out); F.letter(ax, "c", x=-0.56)
    ax = fig.add_subplot(gs[1, 1]); panel_r90(ax, "energy", out); F.letter(ax, "d", x=-0.44)
    ax = fig.add_subplot(gs[2, :]); panel_decomp(ax, out); F.letter(ax, "e", x=-0.105)

    F.save(fig, "fig5")
    (F.OUT / "fig5_numbers.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
