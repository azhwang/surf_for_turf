#!/usr/bin/env python3
"""
fig4.py -- Fig. 4: the ortholog screen.

  a  the two simulated assemblies, Ca traces of the rebuilt goldfish systems
  b  2:2:2 whole-assembly MM/GBSA per ortholog (1 ns screen)
  c  2:2:1 whole-assembly MM/GBSA per ortholog
  d  sequence map: positions differing from human, with FGFR- and heparin-contacting
     positions marked underneath, Arg90 and the 86-97 segment annotated

Writes figures/out/fig4.{svg,png} and returns the numbers quoted in the legend.
"""
import json

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy import stats

import figstyle as F


def assembly_panel(ax, tk, title):
    """PyMOL cartoon render (render_structures.py); falls back to a Ca trace if absent."""
    png = F.OUT / f"struct_{'222' if tk == '1fq9' else '221'}.png"
    if png.exists():
        ax.imshow(F.plt.imread(png))
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    else:
        seg, rid, rn, nm, xyz = F.read_structure(F.RB / tk / "goldfish" / "precheck.pdb")
        ca = nm == "CA"
        centre = xyz[ca].mean(0)
        _, _, vt = np.linalg.svd(xyz[ca] - centre)
        P = (xyz - centre) @ vt[:2].T
        col = {"PROA": "#4c5a7a", "PROB": "#8c606d", "PROC": "#cb685f", "PROD": "#e5ad82"}
        for s, c in col.items():
            m = (seg == s) & ca
            if m.any():
                ax.plot(P[m, 0], P[m, 1], color=c, lw=0.9, zorder=2)
        for s in sorted({x for x in seg if x.startswith("CAR")}):
            m = seg == s
            cent = np.array([P[m & (rid == r)].mean(0) for r in dict.fromkeys(rid[m])])
            ax.plot(cent[:, 0], cent[:, 1], color=F.INK, lw=3.2, solid_capstyle="round", zorder=4)
        ax.set_aspect("equal")
        ax.axis("off")
    ax.set_title(title, loc="left", fontsize=7, color=F.INK, pad=1)


def energy_panel(ax, one, tk, label, out):
    d = one[(one.split == "whole") & (one.tk == tk)]
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
    out[f"anova_{tk}_whole"] = round(float(p), 2)
    out[f"n_{tk}"] = {sp: int((d.sp == sp).sum()) for sp in F.SPECIES}
    ax.set_xticks(range(len(F.SPECIES)), [F.SP_LABEL[s] for s in F.SPECIES], rotation=38, ha="right")
    ax.set_ylabel("$\\Delta G_{bind}$ (kcal mol$^{-1}$)")
    ax.set_title(f"{label}    one-way ANOVA $p$ = {p:.2f}", color=F.INK)
    ax.set_xlim(-0.6, len(F.SPECIES) - 0.4)
    F.style(ax)


def contacts_and_sequences(tk):
    """core sequences per ortholog, plus positions contacting FGFR / heparin in any ortholog."""
    seqs, fgfr_c, hep_c = {}, set(), set()
    for sp in F.SPECIES:
        seg, rid, rn, nm, xyz = F.read_structure(F.RB / tk / sp / "precheck.pdb")
        fgfr = np.isin(seg, ["PROA", "PROB"])
        hep = np.char.startswith(seg.astype(str), "CAR")
        for c in ("PROC",):          # one FGF copy, as in the manuscript text
            m = seg == c
            ids = list(dict.fromkeys(rid[m]))
            seqs[sp] = "".join(F.THREE[rn[m & (rid == r)][0]] for r in ids)
            for k, r in enumerate(ids, 1):
                a = xyz[m & (rid == r)]
                if np.sqrt(((a[:, None] - xyz[fgfr][None]) ** 2).sum(-1)).min() < 4.5:
                    fgfr_c.add(k)
                if np.sqrt(((a[:, None] - xyz[hep][None]) ** 2).sum(-1)).min() < 4.5:
                    hep_c.add(k)
    return seqs, sorted(fgfr_c), sorted(hep_c)


def sequence_panel(ax, out):
    tk = "1fq9"
    seqs, fc, hc = contacts_and_sequences(tk)
    ref, off = seqs["human"], F.OFF[tk]
    n = len(ref)
    ident = [k for k in range(1, n + 1) if len({seqs[s][k - 1] for s in F.SPECIES}) == 1]
    out.update(fgfr_contacts=len(fc), fgfr_identical=sum(k in ident for k in fc),
               hep_contacts=len(hc), hep_identical=sum(k in ident for k in hc),
               core_len=n, core_identity_pct=round(100 * len(ident) / n))
    for i, sp in enumerate(F.SPECIES):
        for k in range(1, n + 1):
            if seqs[sp][k - 1] != ref[k - 1]:
                ax.add_patch(Rectangle((k + off - 0.5, i - 0.40), 1, 0.80,
                                       color=F.SP_COLOR[sp], lw=0))
    for k in fc:
        ax.add_patch(Rectangle((k + off - 0.5, -1.05), 1, 0.42, color="#4c5a7a", lw=0))
    for k in hc:
        ax.add_patch(Rectangle((k + off - 0.5, -1.62), 1, 0.42, color=F.INK, lw=0))
    ax.text(off - 0.5, -0.84, "FGFR contact  ", ha="right", va="center", fontsize=6, color=F.INK)
    ax.text(off - 0.5, -1.41, "heparin contact  ", ha="right", va="center", fontsize=6, color=F.INK)
    ax.plot([86, 97], [-2.15, -2.15], color=F.SP_COLOR["goldfish"], lw=1.6, solid_capstyle="butt")
    ax.text(98.5, -2.15, "86–97 swapped segment", va="center", fontsize=6, color=F.INK)
    ax.annotate("Arg90", xy=(90, -1.70), xytext=(90, -2.75), ha="center", fontsize=6, color=F.INK,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=F.INK))
    ax.set_xlim(off + 0.5, off + n + 0.5)
    ax.set_ylim(len(F.SPECIES) - 0.4, -3.0)
    ax.set_yticks(range(len(F.SPECIES)), [F.SP_LABEL[s] for s in F.SPECIES])
    ax.set_xlabel("FGF2 residue (155-residue numbering)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_title(f"Differences from human FGF2 \u00b7 identical across all seven orthologs at "
                 f"{out['fgfr_identical']}/{out['fgfr_contacts']} FGFR-contacting and "
                 f"{out['hep_identical']}/{out['hep_contacts']} heparin-contacting positions",
                 color=F.INK, pad=4)


def main():
    out = {}
    fig = F.plt.figure(figsize=(F.DOUBLE, 158 * F.MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.25, 0.95, 1.05], hspace=0.40, wspace=0.30)

    a1, a2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    assembly_panel(a1, "1fq9", "2:2:2  (PDB 1FQ9)")
    assembly_panel(a2, "1e0o", "2:2:1  (PDB 1E0O)")
    F.letter(a1, "a", x=-0.10, y=1.0)
    a1.legend(handles=[Line2D([], [], color="#4c5a7a", lw=1.2, label="FGFR (chain A)"),
                       Line2D([], [], color="#8c606d", lw=1.2, label="FGFR (chain B)"),
                       Line2D([], [], color="#cb685f", lw=1.2, label="FGF2 (copy C)"),
                       Line2D([], [], color="#e5ad82", lw=1.2, label="FGF2 (copy D)"),
                       Line2D([], [], color=F.INK, lw=2.6, label="heparin")],
              loc="upper center", bbox_to_anchor=(1.08, 0.06), ncol=5, handlelength=1.3,
              columnspacing=1.4)

    ax = fig.add_subplot(gs[1, :])
    sequence_panel(ax, out)
    F.letter(ax, "d", x=-0.112, y=1.0)

    one = F.one_ns()
    for i, (tk, lab) in enumerate((("1fq9", "2:2:2 whole assembly"), ("1e0o", "2:2:1 whole assembly"))):
        ax = fig.add_subplot(gs[2, i])
        energy_panel(ax, one, tk, lab, out)
        F.letter(ax, "bc"[i], x=-0.145)

    F.save(fig, "fig4")
    (F.OUT / "fig4_numbers.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
