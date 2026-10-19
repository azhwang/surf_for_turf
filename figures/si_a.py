#!/usr/bin/env python3
"""
si_a.py -- Supplementary figures on the build and the 1 ns screen.

  S6  FGF2 trimming: residues kept vs trimmed, per ortholog and template
  S7  build quality control: largest single-residue dvdW per replicate, original vs rebuilt
  S8  per-interface binding energies for all orthologs

Writes figures/out/figS6.{svg,png}, figS7.*, figS8.*.
"""
import csv
import json
import re

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

import figstyle as F

out = {}


def s6_trimming():
    fig, axes = F.plt.subplots(1, 2, figsize=(F.DOUBLE, 62 * F.MM), sharey=True,
                               gridspec_kw=dict(wspace=0.08))
    for ax, tk in zip(axes, ("1fq9", "1e0o")):
        for i, sp in enumerate(F.SPECIES):
            rep = (F.RB / tk / sp / "core_report.txt").read_text()
            full = int(re.search(r"full FGF: (\d+)", rep).group(1))
            m = re.search(r"trimmed N (\d+), C (\d+)", rep)
            nt, ct = int(m.group(1)), int(m.group(2))
            core = full - nt - ct
            ax.barh(i, nt, color=F.GRID, height=0.6, edgecolor=F.INK, lw=0.3)
            ax.barh(i, core, left=nt, color=F.SP_COLOR[sp], height=0.6, edgecolor=F.INK, lw=0.3)
            if ct:
                ax.barh(i, ct, left=nt + core, color=F.GRID, height=0.6, edgecolor=F.INK, lw=0.3)
            ax.text(full + 2, i, f"{core} aa core", va="center", fontsize=6, color=F.INK)
        ax.set_yticks(range(len(F.SPECIES)), [F.SP_LABEL[s] for s in F.SPECIES])
        ax.invert_yaxis()
        ax.set_xlim(0, 185)
        ax.set_xlabel("Residues in the AlphaFold model")
        ax.set_title(f"{F.TMPL[tk]} ({tk.upper()})", color=F.INK)
        F.style(ax, grid_y=False)
    axes[0].legend(handles=[Line2D([], [], color=F.GRID, lw=5, label="trimmed (disordered tail)"),
                            Line2D([], [], color=F.MUTED, lw=5, label="core kept for simulation")],
                   loc="lower right", fontsize=6)
    F.letter(axes[0], "a", x=-0.40); F.letter(axes[1], "b", x=-0.06)
    F.save(fig, "figS6")


def max_dvdw_rows():
    pats = {"Original builds": ("all_metrics_analysis*/1*_*/*/namdrep*/mmpbsa_v2/*_run/FINAL_DECOMP_MMPBSA.dat",
                                lambda p: (p[-6].split("_", 1)[0], p[-6].split("_", 1)[1], p[-4])),
            "Rebuilt systems": ("rebuild_v2/1*/*/sim/namdrep[0-9]/mmpbsa_v2/*_run/FINAL_DECOMP_MMPBSA.dat",
                                lambda p: (p[-7], p[-6], p[-4]))}
    rows = []
    for setname, (pat, parse) in pats.items():
        for f in F.ROOT.glob(pat):
            r_ = list(csv.reader(open(f)))
            try:
                s0 = next(i for i, r in enumerate(r_) if r and r[0] == "DELTAS:")
            except StopIteration:
                continue
            mx = -1e9
            for r in r_[s0 + 4:]:
                if not r or not r[0].strip() or r[0].startswith(("Sidechain", "Backbone", "Residue")):
                    break
                mx = max(mx, float(r[5]))
            tk, sp, rep = parse(f.parts)
            if sp in F.SPECIES:
                rows.append((setname, tk, sp, rep, mx))
    return pd.DataFrame(rows, columns=["set", "tk", "sp", "rep", "mx"]).groupby(
        ["set", "tk", "sp", "rep"]).max().reset_index()


def s7_qc():
    df = max_dvdw_rows()
    out["qc_rebuild_max"] = round(float(df[df.set == "Rebuilt systems"].mx.max()), 1)
    out["qc_original_max"] = round(float(df[df.set == "Original builds"].mx.max()), 1)
    fig, axes = F.plt.subplots(1, 2, figsize=(F.DOUBLE, 66 * F.MM), sharey=True,
                               gridspec_kw=dict(wspace=0.06))
    for ax, tk in zip(axes, ("1fq9", "1e0o")):
        for i, sp in enumerate(F.SPECIES):
            for setname, dx, hollow in (("Original builds", -0.18, True), ("Rebuilt systems", 0.18, False)):
                v = df[(df.set == setname) & (df.tk == tk) & (df.sp == sp)].mx.clip(lower=0.12).values
                F.dots(ax, i + dx, v, F.SP_COLOR[sp], F.SP_MARK[sp], jitter=0.05, hollow=hollow)
        ax.axhline(4, color=F.MUTED, lw=0.6, ls="--")
        ax.set_yscale("log")
        ax.set_xticks(range(len(F.SPECIES)), [F.SP_LABEL[s] for s in F.SPECIES], rotation=38, ha="right")
        ax.set_title(f"{F.TMPL[tk]} ({tk.upper()})", color=F.INK)
        ax.set_xlim(-0.6, len(F.SPECIES) - 0.4)
        F.style(ax)
    axes[0].set_ylabel("Largest single-residue\n$\\Delta E_{vdW}$ (kcal mol$^{-1}$)")
    axes[0].text(-0.45, 4.6, "artefact-free threshold", fontsize=6, color=F.MUTED, va="bottom")
    axes[1].legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=F.INK, ms=4,
                                   label="original build"),
                            Line2D([], [], marker="o", ls="", mfc=F.INK, mec=F.INK, ms=4,
                                   label="rebuilt")], loc="upper right", ncol=2)
    F.letter(axes[0], "a", x=-0.165); F.letter(axes[1], "b", x=-0.05)
    F.save(fig, "figS7")


def s8_interfaces():
    one = F.one_ns()
    fig, axes = F.plt.subplots(2, 2, figsize=(F.DOUBLE, 118 * F.MM),
                               gridspec_kw=dict(hspace=0.55, wspace=0.22))
    for r, tk in enumerate(("1fq9", "1e0o")):
        for c, split in enumerate(("iface1", "iface2")):
            ax = axes[r, c]
            d = one[(one.split == split) & (one.tk == tk)]
            groups = []
            for j, sp in enumerate(F.SPECIES):
                ys = d[d.sp == sp].dG.values
                if len(ys):
                    F.dots(ax, j, ys, F.SP_COLOR[sp], F.SP_MARK[sp])
                    F.meanbar(ax, j, ys)
                    groups.append(ys)
                else:
                    ax.text(j, 0.5, "not yet run", transform=ax.get_xaxis_transform(), rotation=90,
                            ha="center", va="center", color=F.MUTED, fontsize=6)
            p = stats.f_oneway(*[g for g in groups if len(g) > 1]).pvalue
            out[f"anova_{tk}_{split}"] = round(float(p), 2)
            ax.set_title(f"{F.TMPL[tk]} interface {split[-1]}    ANOVA $p$ = {p:.2f}", color=F.INK)
            ax.set_xticks(range(len(F.SPECIES)),
                          [F.SP_LABEL[s] for s in F.SPECIES] if r == 1 else [""] * len(F.SPECIES),
                          rotation=38, ha="right")
            ax.set_xlim(-0.6, len(F.SPECIES) - 0.4)
            if c == 0:
                ax.set_ylabel("$\\Delta G_{bind}$ (kcal mol$^{-1}$)")
            F.style(ax)
            F.letter(ax, "abcd"[2 * r + c], x=-0.16)
    F.save(fig, "figS8")


if __name__ == "__main__":
    s6_trimming()
    s7_qc()
    s8_interfaces()
    (F.OUT / "siA_numbers.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
