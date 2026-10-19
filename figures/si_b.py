#!/usr/bin/env python3
"""
si_b.py -- Supplementary figures on the extensions and the exploratory analyses.

  S9   convergence: 1 ns vs 11-51 ns per replicate, and within- vs between-replicate spread
  S10  full extension time series, both assemblies
  S11  component energies: FGF2-heparin and FGF2-FGFR
  S12  per-residue goldfish - human difference, 2:2:2 (companion to Fig. 5e)
  S13  static electrostatics of the FGF2 cores
  S14  statistical power of the MM/GBSA comparison
"""
import json
import math
import re

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

import figstyle as F
from fig5 import welch

out = {}


def s9_convergence():
    one, ext = F.one_ns(), F.extensions()
    fig, axes = F.plt.subplots(1, 3, figsize=(F.DOUBLE, 62 * F.MM),
                               gridspec_kw=dict(wspace=0.42))
    for ax, tk in zip(axes[:2], ("1e0o", "1fq9")):
        for sp in ("human", "goldfish"):
            for r in F.EXT[f"{tk}/{sp}"]:
                a = one[(one.tk == tk) & (one.sp == sp) & (one.rep == r) & (one.split == "whole")].dG.values
                b = ext[(ext.tmpl == tk) & (ext.sp == sp) & (ext.rep == r)].dG_whole.values
                if len(a) and len(b):
                    ax.plot([0, 1], [a[0], b[0]], color=F.SP_COLOR[sp], lw=0.8, marker="o",
                            ms=3.2, mec=F.INK, mew=0.4)
        ax.set_xticks([0, 1], ["1 ns\n(10 frames)", "11–51 ns\n(80 frames)"])
        ax.set_xlim(-0.35, 1.35)
        ax.set_title(f"{F.TMPL[tk]} whole assembly", color=F.INK)
        ax.set_ylabel("$\\Delta G_{bind}$ (kcal mol$^{-1}$)")
        F.style(ax)
    axes[0].legend(handles=[Line2D([], [], color=F.H, lw=1, label="Human"),
                            Line2D([], [], color=F.G, lw=1, label="Goldfish")], loc="lower left")
    ax = axes[2]
    labels, within, between = [], [], []
    for tk in ("1e0o", "1fq9"):
        for sp in ("human", "goldfish"):
            d = ext[(ext.tmpl == tk) & (ext.sp == sp)]
            labels.append(f"{F.TMPL[tk]}\n{F.SP_LABEL[sp]}")
            within.append(float(d.blockSEM_whole.mean()))
            between.append(float(d.dG_whole.std(ddof=1)))
    out["within_mean"] = round(float(np.mean(within)), 1)
    out["between_mean"] = round(float(np.mean(between)), 1)
    x = np.arange(len(labels))
    ax.bar(x - 0.19, within, 0.36, color=F.GRID, edgecolor=F.INK, lw=0.4,
           label="within replicate (10 ns block SEM)")
    ax.bar(x + 0.19, between, 0.36, color=F.MUTED, edgecolor=F.INK, lw=0.4,
           label="between replicates (SD)")
    ax.set_xticks(x, labels, fontsize=6)
    ax.set_ylabel("Uncertainty (kcal mol$^{-1}$)")
    ax.legend(loc="upper left", fontsize=6)
    ax.set_title("Where the uncertainty sits", color=F.INK)
    F.style(ax)
    for a, l in zip(axes, "abc"):
        F.letter(a, l, x=-0.30)
    F.save(fig, "figS9")


def s10_timeseries():
    fig, axes = F.plt.subplots(2, 3, figsize=(F.DOUBLE, 104 * F.MM), sharex=True,
                               gridspec_kw=dict(hspace=0.42, wspace=0.40))
    cols = [(2, "Heparin-site salt bridges\nFGF2 copy C"),
            (3, "Heparin-site salt bridges\nFGF2 copy D"),
            (6, "FGFR–FGFR distance (Å)")]
    for r, tk in enumerate(("1e0o", "1fq9")):
        for c, (ci, lab) in enumerate(cols):
            ax = axes[r, c]
            for sp in ("human", "goldfish"):
                for i, rep in enumerate(F.EXT[f"{tk}/{sp}"]):
                    d = F.timeseries(tk, sp, rep)
                    if d is None:
                        continue
                    ax.plot(d[9:, 1], F.running(d[:, ci]), color=F.SP_COLOR[sp], lw=0.7,
                            ls=["-", "--", ":"][i])
            ax.axvspan(1, 11, color=F.PALE, zorder=0)
            ax.set_title(F.TMPL[tk], color=F.INK)
            ax.set_ylabel(lab, fontsize=6.5)
            ax.set_xlim(1, 51)
            if r == 1:
                ax.set_xlabel("Time (ns)")
            F.style(ax)
            F.letter(ax, "abcdef"[3 * r + c], x=-0.28)
    axes[0, 2].legend(handles=[Line2D([], [], color=F.H, lw=1, label="Human"),
                               Line2D([], [], color=F.G, lw=1, label="Goldfish")],
                      loc="upper right", ncol=2)
    F.save(fig, "figS10")


def s11_components():
    ext = F.extensions()
    fig, axes = F.plt.subplots(1, 2, figsize=(F.DOUBLE * 0.72, 62 * F.MM),
                               gridspec_kw=dict(wspace=0.46))
    groups = [("1e0o", "human"), ("1e0o", "goldfish"), ("1fq9", "human"), ("1fq9", "goldfish")]
    gx = [0, 1, 2.5, 3.5]
    for ax, col, lab in zip(axes, ("hep", "fgfr"),
                            ("FGF2–heparin\n(FGFR removed)", "FGF2–FGFR\n(heparin removed)")):
        for (tk, sp), x in zip(groups, gx):
            ys = ext[(ext.tmpl == tk) & (ext.sp == sp)][col].dropna().values
            F.dots(ax, x, ys, F.SP_COLOR[sp])
            F.meanbar(ax, x, ys)
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.30 * (hi - lo))
        for tk, x0 in (("1e0o", 0.5), ("1fq9", 3.0)):
            g = ext[(ext.tmpl == tk) & (ext.sp == "goldfish")][col].dropna()
            h = ext[(ext.tmpl == tk) & (ext.sp == "human")][col].dropna()
            d, ci, p = welch(g, h)
            out[f"{col}_{tk}"] = dict(diff=round(d, 1), p=round(p, 2))
            ax.text(x0, 0.995, f"$\\Delta$ {d:+.0f}\n$p$ = {p:.2f}", transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=6, color=F.INK)
        ax.set_xticks(gx, ["Human", "Goldfish", "Human", "Goldfish"], rotation=30, ha="right")
        ax.text(0.5, -0.46, "2:2:1", transform=ax.get_xaxis_transform(), ha="center", fontsize=7)
        ax.text(3.0, -0.46, "2:2:2", transform=ax.get_xaxis_transform(), ha="center", fontsize=7)
        ax.set_xlim(-0.55, 4.05)
        ax.set_title(lab, color=F.INK)
        ax.set_ylabel("$\\Delta G$ (kcal mol$^{-1}$)")
        F.style(ax)
    F.letter(axes[0], "a", x=-0.38); F.letter(axes[1], "b", x=-0.38)
    F.save(fig, "figS11")


def s12_decomp222():
    d = pd.read_csv(F.RB / "ext_results/explo_decomp_1fq9.csv")
    tot = d.groupby(["fgf2_pos", "substituted"]).g_minus_h.sum().reset_index()
    sub = tot[tot.substituted == 1]
    fig, ax = F.plt.subplots(figsize=(F.DOUBLE * 0.72, 54 * F.MM))
    ax.axhline(0, color=F.INK, lw=0.5)
    ax.vlines(tot.fgf2_pos, 0, tot.g_minus_h, color=F.MUTED, lw=0.8)
    ax.vlines(sub.fgf2_pos, 0, sub.g_minus_h, color=F.G, lw=1.3)
    top = tot.reindex(tot.g_minus_h.abs().sort_values(ascending=False).index).head(4)
    lo, hi = tot.g_minus_h.min(), tot.g_minus_h.max()
    ax.set_ylim(lo - 0.25 * (hi - lo), hi + 0.3 * (hi - lo))
    for _, r in top.iterrows():
        ax.annotate(f"{int(r.fgf2_pos)}", (r.fgf2_pos, r.g_minus_h), textcoords="offset points",
                    xytext=(0, -8 if r.g_minus_h < 0 else 4), ha="center", fontsize=6, color=F.INK)
    ax.set_xlabel("FGF2 residue")
    ax.set_ylabel("Goldfish − human\n(kcal mol$^{-1}$, both copies)")
    ax.set_title("2:2:2 per-residue difference (companion to Fig. 5e)", color=F.INK)
    ax.legend(handles=[Line2D([], [], color=F.G, lw=1.3, label="substituted"),
                       Line2D([], [], color=F.MUTED, lw=0.9, label="conserved")],
              loc="lower right", ncol=2)
    F.style(ax)
    F.save(fig, "figS12")


def s13_electrostatics():
    txt = (F.RB / "ext_results/explo_electrostatics.txt").read_text().split("\n")
    rec, tk = {}, None
    for l in txt:
        m = re.match(r"== (2:2:[12]) \((\w+)\)", l)
        if m:
            tk = m.group(2)
            continue
        p = l.split()
        if tk and len(p) == 6 and p[0] in F.SPECIES:
            rec[(tk, p[0])] = [float(v) for v in p[1:]]
    fig, axes = F.plt.subplots(1, 2, figsize=(F.DOUBLE * 0.82, 62 * F.MM),
                               gridspec_kw=dict(wspace=0.08, width_ratios=[1, 1.5]))
    ax = axes[0]
    for i, sp in enumerate(F.SPECIES):
        q = rec.get(("1fq9", sp), [np.nan])[0]
        ax.barh(i, q, color=F.SP_COLOR[sp], edgecolor=F.INK, lw=0.3, height=0.6)
    ax.set_yticks(range(len(F.SPECIES)), [F.SP_LABEL[s] for s in F.SPECIES])
    ax.invert_yaxis()
    ax.set_xlabel("Net charge of the FGF2 core ($e$)")
    F.style(ax, grid_y=False)
    ax.xaxis.grid(True, color=F.GRID, lw=0.4); ax.set_axisbelow(True)
    ax = axes[1]
    for j, tk in enumerate(("1e0o", "1fq9")):
        for i, sp in enumerate(F.SPECIES):
            v = rec.get((tk, sp))
            if v:
                ax.scatter([v[1] + v[2]], [i + (j - 0.5) * 0.28], s=17,
                           marker="o" if tk == "1e0o" else "s",
                           facecolor=F.SP_COLOR[sp], edgecolor=F.INK, lw=0.5, zorder=3)
    ax.set_yticks(range(len(F.SPECIES)), [""] * len(F.SPECIES))
    ax.invert_yaxis()
    ax.set_xlabel("Screened Coulomb FGF2–heparin energy,\nboth copies (kcal mol$^{-1}$)")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=F.INK, ms=4, label="2:2:1"),
                       Line2D([], [], marker="s", ls="", mfc="white", mec=F.INK, ms=4, label="2:2:2")],
              loc="lower right", ncol=2)
    F.style(ax, grid_y=False)
    ax.xaxis.grid(True, color=F.GRID, lw=0.4); ax.set_axisbelow(True)
    F.letter(axes[0], "a", x=-0.62); F.letter(axes[1], "b", x=-0.05)
    F.save(fig, "figS13")


def s14_power():
    one, ext = F.one_ns(), F.extensions()

    def pooled(vals):
        num = sum((len(v) - 1) * np.var(v, ddof=1) for v in vals)
        den = sum(len(v) - 1 for v in vals)
        return math.sqrt(num / den)

    sd1 = pooled([one[(one.tk == "1e0o") & (one.split == "whole") & (one.sp == s)].dG.values
                  for s in ("human", "goldfish")])
    sd50 = pooled([ext[(ext.tmpl == "1e0o") & (ext.sp == s)].dG_whole.values
                   for s in ("human", "goldfish")])
    out["sd_1ns"], out["sd_50ns"] = round(float(sd1), 1), round(float(sd50), 1)

    def mdd(n, sd):
        df, lo, hi = 2 * n - 2, 0.0, 12 * sd
        for _ in range(60):
            d = (lo + hi) / 2
            nc = d / (sd * math.sqrt(2 / n))
            tc = stats.t.ppf(0.975, df)
            power = 1 - stats.nct.cdf(tc, df, nc) + stats.nct.cdf(-tc, df, nc)
            lo, hi = (d, hi) if power < 0.8 else (lo, d)
        return d

    ns = np.arange(3, 41)
    fig, ax = F.plt.subplots(figsize=(F.SINGLE * 1.25, 60 * F.MM))
    # the two sampling lengths give the same between-replicate spread (SD 24.8 vs 24.9),
    # so the curves coincide; draw one and say so
    ax.plot(ns, [mdd(int(n), sd50) for n in ns], color=F.G, lw=1.3,
            label=f"replicate SD {sd50:.0f} kcal mol$^{{-1}}$\n(1 ns and 50 ns agree: {sd1:.0f} vs {sd50:.0f})")
    ax.axhline(1.0, color=F.INK, lw=0.8, ls="--")
    ax.text(40, 1.12, "experimental effect ≈ 1 kcal mol$^{-1}$", ha="right", fontsize=6, color=F.INK)
    for n, m in ((3, mdd(3, sd50)),):
        ax.scatter([n], [m], s=18, facecolor=F.G, edgecolor=F.INK, lw=0.5, zorder=4)
        ax.annotate(f"this study: n = 3,\nresolves {m:.0f} kcal mol$^{{-1}}$", (n, m), xytext=(12, -14),
                    textcoords="offset points", fontsize=6, color=F.INK)
    out["mdd_n3_50ns"] = round(float(mdd(3, sd50)), 1)
    ax.set_yscale("log")
    ax.set_xlabel("Replicates per species")
    ax.set_ylabel("Smallest difference detectable\nat 80% power (kcal mol$^{-1}$)")
    ax.set_title("2:2:1 whole-assembly MM/GBSA", color=F.INK)
    ax.legend(loc="upper right")
    F.style(ax)
    F.save(fig, "figS14")


if __name__ == "__main__":
    for fn in (s9_convergence, s10_timeseries, s11_components, s12_decomp222,
               s13_electrostatics, s14_power):
        fn()
    (F.OUT / "siB_numbers.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
