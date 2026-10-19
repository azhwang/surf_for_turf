#!/usr/bin/env python3
"""
figstyle.py -- shared style and data loaders for the paper figures.

Colours: cmcrameri `lipari` throughout. The seven ortholog slots are taken from the ramp in an
order that alternates lightness, so adjacent orthologs stay distinguishable under deuteranopia
and protanopia (checked with the dataviz validator: worst adjacent pair dE 13.2 CVD, 22.6 normal).
Light slots carry a dark outline so they hold up against the page.

Type: Neue Haas Grotesk where installed, Helvetica Neue otherwise. SVG output keeps text as text
(`svg.fonttype = none`) so the font resolves on whatever machine opens it; the PNG/DOCX route
rasterises with whatever is installed here.

Sizes follow Nature: 89 mm single column, 183 mm double, 7 pt body text.
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cmcrameri import cm

ROOT = Path(__file__).resolve().parent.parent
RB = ROOT / "rebuild_v2"
OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)

MM = 1 / 25.4
SINGLE, DOUBLE = 89 * MM, 183 * MM

LIP = cm.lipari
INK = "#031326"        # lipari darkest: text, axes, heparin
MUTED = "#655e77"
GRID = "#ebcfa7"
PALE = "#fdf5da"       # lipari lightest: shaded bands

SPECIES = ["human", "human_stab", "bovine_stab", "chinook_salmon", "goldfish",
           "japanese_ricefish", "zebrafish"]
SP_LABEL = {"human": "Human", "human_stab": "Human STAB", "bovine_stab": "Bovine STAB",
            "chinook_salmon": "Chinook salmon", "goldfish": "Goldfish",
            "japanese_ricefish": "Japanese rice fish", "zebrafish": "Zebrafish",
            "chim_human_gfseg": "Human + gf 86–97", "chim_goldfish_hseg": "Goldfish + h 86–97"}
SP_COLOR = {"human": "#0f3354", "human_stab": "#8c606d", "bovine_stab": "#e5ad82",
            "chinook_salmon": "#4c5a7a", "goldfish": "#cb685f",
            "japanese_ricefish": "#e7c398", "zebrafish": "#655e77"}
SP_MARK = {"human": "o", "human_stab": "s", "bovine_stab": "D", "chinook_salmon": "^",
           "goldfish": "o", "japanese_ricefish": "v", "zebrafish": "P"}
H, G = SP_COLOR["human"], SP_COLOR["goldfish"]

TMPL = {"1e0o": "2:2:1", "1fq9": "2:2:2"}
OFF = {"1e0o": 27, "1fq9": 24}        # FGF2 (155-aa) number = core position + OFF
NCORE = {"1e0o": 127, "1fq9": 129}
EXT = {"1e0o/human": [1, 2, 4], "1e0o/goldfish": [1, 2, 3],
       "1fq9/human": [1, 2, 3], "1fq9/goldfish": [2, 3, 5]}
THREE = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H',
             HSD='H', HSE='H', HSP='H', ILE='I', LEU='L', LYS='K', MET='M', PHE='F', PRO='P',
             SER='S', THR='T', TRP='W', TYR='Y', VAL='V')

FONTS = ["Neue Haas Grotesk Text Pro", "Neue Haas Grotesk Display Pro", "Neue Haas Unica",
         "Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]
FONT_CSS = ("'Neue Haas Grotesk Text Pro','Neue Haas Grotesk Display Pro','Neue Haas Unica',"
            "'Helvetica Neue',Helvetica,Arial,sans-serif")

plt.rcParams.update({
    "font.family": FONTS, "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none", "figure.dpi": 150, "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03, "legend.frameon": False, "axes.titlelocation": "left",
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.edgecolor": INK,
})


def letter(ax, s, x=-0.16, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="left", color=INK)


def style(ax, grid_y=True):
    if grid_y:
        ax.yaxis.grid(True, color=GRID, lw=0.4)
        ax.set_axisbelow(True)


def save(fig, name, png_dpi=400):
    """SVG (text as text, for the journal) + PNG (for the Word document)."""
    svg = OUT / f"{name}.svg"
    fig.savefig(svg, facecolor="white")
    fig.savefig(OUT / f"{name}.png", facecolor="white", dpi=png_dpi)
    plt.close(fig)
    t = svg.read_text()
    t = re.sub(r"font-family:[^;\"']*", f"font-family:{FONT_CSS}", t)
    svg.write_text(t)
    return OUT / f"{name}.png"


def dots(ax, x, ys, color, marker="o", size=13, jitter=0.08, hollow=False, z=3):
    ys = list(ys)
    if not ys:
        return
    xs = x + np.linspace(-jitter, jitter, len(ys)) if len(ys) > 1 else [x]
    ax.scatter(xs, ys, s=size, marker=marker,
               facecolor="white" if hollow else color,
               edgecolor=color if hollow else INK,
               linewidth=0.9 if hollow else 0.5, zorder=z, clip_on=False)


def meanbar(ax, x, ys, color=INK, w=0.26):
    ys = list(ys)
    if not ys:
        return
    m = float(np.mean(ys))
    ax.plot([x - w, x + w], [m, m], color=color, lw=1.1, zorder=4, solid_capstyle="butt")
    if len(ys) > 1:
        s = float(np.std(ys, ddof=1))
        ax.plot([x, x], [m - s, m + s], color=color, lw=0.6, zorder=4)


# ---------------------------------------------------------------- data
def one_ns():
    """1 ns screen, all orthologs: tk, sp, rep, split, dG."""
    x = pd.read_excel(ROOT / "mmpbsa_rebuild_results/mmpbsa_v2_summary.xlsx", sheet_name="Replicates")
    x = x.assign(tk=np.where(x.Stoichiometry.str.startswith("2:2:1"), "1e0o", "1fq9"),
                 sp=x.Species, rep=x.Replicate.str.replace("namdrep", "").astype(int),
                 split=x.Split, dG=x["ΔG total"])
    return x[["tk", "sp", "rep", "split", "dG"]]


def extensions():
    """50 ns extensions, complete replicates, with the exploratory component energies joined."""
    t = pd.read_csv(RB / "ext_results/replicate_table.csv")
    t = t[t.status == "complete"].copy()
    e = (pd.read_csv(RB / "ext_results/explo_energies.csv")
         .pivot_table(index=["tmpl", "species", "rep"], columns="split", values="dG").reset_index())
    return t.merge(e, left_on=["tmpl", "sp", "rep"], right_on=["tmpl", "species", "rep"], how="left")


def timeseries(tk, sp, r):
    f = RB / tk / sp / "sim" / f"namdrep{r}" / "ext_analysis/timeseries.dat"
    return np.loadtxt(f, comments="#") if f.exists() else None


def running(y, k=10):
    return np.convolve(y, np.ones(k) / k, mode="valid")


def read_structure(path, heavy_only=True):
    """(segid, resid, resname, atomname, xyz) arrays from a PDB."""
    seg, rid, rn, nm, xyz = [], [], [], [], []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        name = l[12:16].strip()
        if heavy_only and name.startswith("H"):
            continue
        seg.append(l[72:76].strip()); rid.append(int(l[22:26])); rn.append(l[17:21].strip())
        nm.append(name); xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return (np.array(seg), np.array(rid), np.array(rn), np.array(nm), np.array(xyz))


def core_sequence(tk, sp):
    seg, rid, rn, _, _ = read_structure(RB / tk / sp / "fgf_PROC.pdb")
    return "".join(THREE[rn[rid == r][0]] for r in dict.fromkeys(rid))
