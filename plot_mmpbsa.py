#!/usr/bin/env python3
"""
MMPBSA cross-species comparison: replicate-averaged ΔG with error bars + ranges.
Usage: python plot_mmpbsa.py mmpbsa_summary_221.xlsx
"""
import re, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

XLSX = sys.argv[1] if len(sys.argv) > 1 else "mmpbsa_summary_221.xlsx"

# ---------- 1. LOAD & PARSE ----------
df = pd.read_excel(XLSX)
df.columns = [c.strip() for c in df.columns]

def parse_run(p):
    p = str(p)
    species = p.split("/")[0]
    rep = re.search(r"namdrep(\d+)", p)
    sect = re.search(r"(iface1|iface2|whole)_run", p)
    return pd.Series({
        "species": species,
        "rep": int(rep.group(1)) if rep else np.nan,
        "section": sect.group(1) if sect else "unknown",
    })

df = pd.concat([df, df["Run"].apply(parse_run)], axis=1)
df["species_label"] = (df["species"]
                       .str.replace(r"^1e0o_", "", regex=True)
                       .str.replace("_", " ").str.title())

DELTA = "Delta (Complex-Receptor-Ligand)"
SECT_ORDER = ["iface1", "iface2", "whole"]

# ---------- 2. AGGREGATE ----------
def agg(component, section_filter=DELTA):
    d = df[(df["Section"] == section_filter) & (df["Component"] == component)]
    g = (d.groupby(["species_label", "section"])["Average"]
           .agg(n="count", mean="mean", sd="std", min="min", max="max")
           .reset_index())
    g["sem"] = g["sd"] / np.sqrt(g["n"])
    # t-based 95% CI half-width; NaN-safe for n==1
    g["ci95"] = [r.sem * stats.t.ppf(0.975, r.n - 1) if r.n > 1 else np.nan
                 for r in g.itertuples()]
    g["section"] = pd.Categorical(g["section"], SECT_ORDER, ordered=True)
    return g.sort_values(["species_label", "section"])

tot = agg("DELTA TOTAL")
tot.to_csv("mmpbsa_deltatotal_summary.csv", index=False)
print(tot.to_string(index=False))

species = sorted(tot["species_label"].unique())
cmap = plt.get_cmap("tab10")
COL = {s: cmap(i) for i, s in enumerate(species)}

# ---------- FIG 1: grouped bars, SD error bars, per-rep points ----------
fig, ax = plt.subplots(figsize=(9, 5.5))
w = 0.8 / len(species)
x = np.arange(len(SECT_ORDER))
raw = df[(df["Section"] == DELTA) & (df["Component"] == "DELTA TOTAL")]

for i, sp in enumerate(species):
    sub = tot[tot["species_label"] == sp].set_index("section").reindex(SECT_ORDER)
    off = (i - (len(species) - 1) / 2) * w
    ax.bar(x + off, sub["mean"], w, yerr=sub["sd"], capsize=5,
           label=f"{sp} (n={int(sub['n'].max())})",
           color=COL[sp], edgecolor="k", lw=0.7,
           error_kw=dict(lw=1.4, ecolor="k"))
    # min–max range whiskers (thin, behind points)
    for xi, s in zip(x + off, SECT_ORDER):
        r = sub.loc[s]
        if pd.notna(r["min"]):
            ax.plot([xi, xi], [r["min"], r["max"]], color="0.25", lw=0.9, zorder=2)
    # individual replicates, jittered
    pts = raw[raw["species_label"] == sp]
    for s in SECT_ORDER:
        v = pts[pts["section"] == s]["Average"].values
        if len(v):
            xi = x[SECT_ORDER.index(s)] + off
            ax.scatter(np.full(len(v), xi) + np.random.uniform(-w*.22, w*.22, len(v)),
                       v, s=26, color="white", edgecolor="k", lw=.8, zorder=5)

ax.set_xticks(x); ax.set_xticklabels(SECT_ORDER)
ax.set_ylabel(r"$\Delta G_{\rm bind}$  (kcal mol$^{-1}$, MM/GBSA)")
ax.set_title("Interface binding energy by species\nbars = mean ± SD across reps; "
             "whiskers = min–max; points = individual reps")
ax.axhline(0, color="k", lw=.8)
ax.legend(frameon=False); ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig("fig1_deltaG_bars.png", dpi=300)

# ---------- FIG 2: box + swarm (shows full range honestly) ----------
fig, axes = plt.subplots(1, len(SECT_ORDER), figsize=(4*len(SECT_ORDER), 4.6),
                         sharey=True)
for axx, s in zip(np.atleast_1d(axes), SECT_ORDER):
    data = [raw[(raw.species_label == sp) & (raw.section == s)]["Average"].values
            for sp in species]
    bp = axx.boxplot(data, patch_artist=True, widths=.55, showmeans=True,
                     labels=[sp.split()[0] for sp in species])
    for patch, sp in zip(bp["boxes"], species):
        patch.set_facecolor(COL[sp]); patch.set_alpha(.45)
    for j, v in enumerate(data, start=1):
        axx.scatter(np.full(len(v), j) + np.random.uniform(-.09, .09, len(v)),
                    v, s=30, color="k", zorder=5)
    axx.set_title(s); axx.axhline(0, color="k", lw=.8)
    axx.spines[["top", "right"]].set_visible(False)
np.atleast_1d(axes)[0].set_ylabel(r"$\Delta G_{\rm bind}$ (kcal mol$^{-1}$)")
fig.suptitle("Replicate distributions per interface")
fig.tight_layout(); fig.savefig("fig2_box_swarm.png", dpi=300)

# ---------- FIG 3: component decomposition ----------
COMPS = ["VDWAALS", "EEL", "EGB", "ESURF", "DELTA TOTAL"]
fig, axes = plt.subplots(1, len(SECT_ORDER), figsize=(5*len(SECT_ORDER), 4.8),
                         sharey=False)
for axx, s in zip(np.atleast_1d(axes), SECT_ORDER):
    xc = np.arange(len(COMPS)); wc = .8/len(species)
    for i, sp in enumerate(species):
        m, e = [], []
        for c in COMPS:
            a = agg(c)
            r = a[(a.species_label == sp) & (a.section == s)]
            m.append(r["mean"].iat[0] if len(r) else np.nan)
            e.append(r["sd"].iat[0]   if len(r) else np.nan)
        axx.bar(xc + (i-(len(species)-1)/2)*wc, m, wc, yerr=e, capsize=3,
                color=COL[sp], edgecolor="k", lw=.6, label=sp)
    axx.set_xticks(xc); axx.set_xticklabels(COMPS, rotation=35, ha="right")
    axx.axhline(0, color="k", lw=.8); axx.set_title(s)
    axx.set_ylabel("kcal mol$^{-1}$"); axx.spines[["top","right"]].set_visible(False)
np.atleast_1d(axes)[0].legend(frameon=False, fontsize=8)
fig.suptitle("Energy-component decomposition (mean ± SD over reps)")
fig.tight_layout(); fig.savefig("fig3_components.png", dpi=300)

# ---------- FIG 4: forest plot with 95% CI ----------
fig, ax = plt.subplots(figsize=(7.5, 0.55*len(tot)+2))
lbl, yy = [], 0
for sp in species:
    for s in SECT_ORDER:
        r = tot[(tot.species_label == sp) & (tot.section == s)]
        if not len(r): continue
        r = r.iloc[0]
        err = r["ci95"] if pd.notna(r["ci95"]) else 0
        ax.errorbar(r["mean"], yy, xerr=err, fmt="o", ms=7, capsize=5,
                    color=COL[sp], mec="k")
        ax.plot([r["min"], r["max"]], [yy, yy], color=COL[sp], alpha=.35, lw=5,
                solid_capstyle="butt", zorder=1)
        lbl.append(f"{sp} · {s} (n={int(r['n'])})"); yy += 1
ax.set_yticks(range(len(lbl))); ax.set_yticklabels(lbl, fontsize=9)
ax.invert_yaxis(); ax.axvline(0, color="k", lw=.8)
ax.set_xlabel(r"$\Delta G_{\rm bind}$ (kcal mol$^{-1}$)")
ax.set_title("Point estimate ± 95% CI (t, n−1); shaded = observed range")
ax.spines[["top","right","left"]].set_visible(False)
fig.tight_layout(); fig.savefig("fig4_forest.png", dpi=300)

# ---------- 3. STATS ----------
print("\n--- Paired iface1 vs iface2 (within species, by rep) ---")
for sp in species:
    p = (raw[raw.species_label == sp]
         .pivot_table(index="rep", columns="section", values="Average"))
    if {"iface1","iface2"} <= set(p.columns):
        d = (p["iface1"] - p["iface2"]).dropna()
        if len(d) > 1:
            t, pv = stats.ttest_rel(p["iface1"].dropna(), p["iface2"].dropna())
            print(f"{sp}: Δ={d.mean():+.2f} ± {d.std(ddof=1):.2f}, t={t:.2f}, p={pv:.3f}")

print("\n--- Additivity: whole − (iface1+iface2) ---")
for sp in species:
    p = (raw[raw.species_label == sp]
         .pivot_table(index="rep", columns="section", values="Average"))
    if {"iface1","iface2","whole"} <= set(p.columns):
        d = (p["whole"] - p["iface1"] - p["iface2"]).dropna()
        if len(d) > 1:
            t, pv = stats.ttest_1samp(d, 0)
            print(f"{sp}: {d.mean():+.2f} ± {d.std(ddof=1):.2f}, t={t:.2f}, p={pv:.3f}")

print("\n--- Between-species (Welch, per section) ---")
for s in SECT_ORDER:
    for i in range(len(species)):
        for j in range(i+1, len(species)):
            a = raw[(raw.species_label==species[i]) & (raw.section==s)]["Average"]
            b = raw[(raw.species_label==species[j]) & (raw.section==s)]["Average"]
            if len(a) > 1 and len(b) > 1:
                t, pv = stats.ttest_ind(a, b, equal_var=False)
                print(f"{s}: {species[i]} vs {species[j]}  "
                      f"Δ={a.mean()-b.mean():+.2f}, t={t:.2f}, p={pv:.3f}")
plt.show()
