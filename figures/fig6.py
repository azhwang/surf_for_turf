#!/usr/bin/env python3
"""
fig6.py -- Fig. 6: the segment-swap chimeras.

  a  design: the FGF2 86-97 segment in human and goldfish, and the two chimeras
  b  Arg90 engagement in 2:2:1 for the two parents, with the pre-registered decision
     rule drawn and the chimera slots left open until those runs finish

Panel b fills in automatically once each chimera replicate has ext_analysis/r90.txt.
Writes figures/out/fig6.{svg,png}.
"""
import json

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import figstyle as F
from fig5 import r90_contacts

SEGMENT = range(86, 98)
CHIM = ["chim_human_gfseg", "chim_goldfish_hseg"]
CHIM_HOST = {"chim_human_gfseg": "human", "chim_goldfish_hseg": "goldfish"}


def design_panel(ax, out):
    tk = "1e0o"
    off = F.OFF[tk]
    h = F.core_sequence(tk, "human")
    g = F.core_sequence(tk, "goldfish")
    k0 = SEGMENT.start - off - 1
    hs, gs = h[k0:k0 + 12], g[k0:k0 + 12]
    out["segment_human"], out["segment_goldfish"] = hs, gs
    out["substitutions"] = [f"{a}{SEGMENT.start + i}{b}" for i, (a, b) in enumerate(zip(hs, gs)) if a != b]
    rows = [("Human", hs, F.H), ("Goldfish", gs, F.G),
            ("Human + goldfish 86–97", gs, F.H), ("Goldfish + human 86–97", hs, F.G)]
    for i, (lab, s, body) in enumerate(rows):
        y = -i
        ax.add_patch(Rectangle((-4.4, y - 0.30), 3.9, 0.60, color=body, lw=0))
        ax.text(-2.45, y, "core", ha="center", va="center", fontsize=6, color="white")
        for j, aa in enumerate(s):
            differs = hs[j] != gs[j]
            fc = (F.G if s[j] == gs[j] else F.H) if differs else "white"
            ax.add_patch(Rectangle((j, y - 0.30), 0.88, 0.60, facecolor=fc,
                                   edgecolor=F.INK, lw=0.4))
            ax.text(j + 0.44, y, aa, ha="center", va="center", fontsize=6.5,
                    color="white" if differs else F.INK)
        ax.text(-4.7, y, lab, ha="right", va="center", fontsize=6.5, color=F.INK)
    for j in range(12):
        ax.text(j + 0.44, 0.52, str(SEGMENT.start + j), ha="center", va="bottom",
                fontsize=5.5, color=F.MUTED, rotation=90)
    ax.annotate("Arg90 (conserved)", xy=(4.44, -3.34), xytext=(6.2, -3.95), fontsize=6,
                color=F.INK, ha="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=F.INK))
    ax.set_xlim(-12.5, 12.6)
    ax.set_ylim(-4.25, 1.35)
    ax.axis("off")
    ax.set_title("Segment swap, FGF2 86–97 (coloured = differs between species)", color=F.INK)


def result_panel(ax, out):
    xs = [0, 1, 2.4, 3.4]
    labels = ["Human", "Goldfish", "Human\n+ gf 86–97", "Goldfish\n+ h 86–97"]
    for x, sp in zip(xs[:2], ("human", "goldfish")):
        v = [r90_contacts("1e0o", sp, r) for r in F.EXT[f"1e0o/{sp}"]]
        v = [q for q in v if q]
        F.dots(ax, x - 0.13, [q[0] for q in v], F.SP_COLOR[sp], "o", jitter=0.04)
        F.dots(ax, x + 0.13, [q[1] for q in v], F.SP_COLOR[sp], "^", jitter=0.04)
        out[sp] = [[round(a, 1), round(b, 1)] for a, b in v]
    pending = []
    for x, name in zip(xs[2:], CHIM):
        vals = []
        for r in (1, 2, 3):
            q = r90_contacts("1e0o", name, r)
            if q:
                vals.append(q)
        col = F.SP_COLOR[CHIM_HOST[name]]
        if vals:
            F.dots(ax, x - 0.13, [q[0] for q in vals], col, "o", jitter=0.04)
            F.dots(ax, x + 0.13, [q[1] for q in vals], col, "^", jitter=0.04)
            out[name] = [[round(a, 1), round(b, 1)] for a, b in vals]
        else:
            pending.append(name)
            ax.add_patch(Rectangle((x - 0.33, 0.05), 0.66, 0.90, transform=ax.get_xaxis_transform(),
                                   fill=False, edgecolor=col, lw=0.7, ls=(0, (2, 2))))
            ax.text(x, 0.50, "simulations\nrunning", transform=ax.get_xaxis_transform(),
                    ha="center", va="center", fontsize=6, color=F.MUTED)
    out["pending"] = pending
    ax.axhspan(5.5, 9.6, color=F.PALE, zorder=0)
    ax.axhline(7.5, color=F.MUTED, lw=0.6, ls="--")
    ax.text(-0.5, 8.2, "engagement threshold", ha="left", va="bottom", fontsize=6, color=F.MUTED)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.16 * (hi - lo))
    ax.set_xticks(xs, labels, fontsize=6.5)
    ax.set_ylabel("Arg90–heparin contacts\n(heavy atoms within 4.5 Å)")
    ax.set_xlim(-0.6, 4.0)
    ax.set_title("2:2:1 · pre-registered test: the goldfish segment confers Arg90 engagement",
                 color=F.INK)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=F.INK, ms=4, label="copy C"),
                       Line2D([], [], marker="^", ls="", mfc="white", mec=F.INK, ms=4, label="copy D")],
              loc="upper center", handletextpad=0.2, ncol=2, bbox_to_anchor=(0.5, 1.0))
    F.style(ax)


def main():
    out = {}
    fig = F.plt.figure(figsize=(F.DOUBLE, 66 * F.MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1.0], wspace=0.28)
    ax = fig.add_subplot(gs[0, 0]); design_panel(ax, out); F.letter(ax, "a", x=0.0, y=0.97)
    ax = fig.add_subplot(gs[0, 1]); result_panel(ax, out); F.letter(ax, "b", x=-0.20)
    F.save(fig, "fig6")
    (F.OUT / "fig6_numbers.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
