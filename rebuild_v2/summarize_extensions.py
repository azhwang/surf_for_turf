#!/usr/bin/env python3
"""
summarize_extensions.py -- apply ANALYSIS_PLAN_extensions.md to the analysed extensions.

Reads <rep>/ext_analysis/{timeseries.dat, <split>_run/FINAL_RESULTS_<split>.dat} for the
replicates selected in make_extensions.py and writes ext_results/:
  replicate_table.csv   per replicate: window means (extension segments 11-50) of every metric
  comparison.txt        goldfish vs human per stoichiometry: difference, 95% CI, Welch p,
                        and the pre-specified verdict on the v1 claim
  timeseries_<tmpl>.png heparin-site salt bridges and FGFR-FGFR distance over time
Primary analysis uses only replicates with all 50 segments (plan); others are listed.
"""
import math
import re
import statistics as st
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

BASE = Path(__file__).resolve().parent
OUT = BASE / "ext_results"
SEL = {"1e0o/human": [1, 2, 4], "1e0o/goldfish": [1, 2, 3], "1fq9/human": [1, 2, 3], "1fq9/goldfish": [2, 3, 5]}
STOICH = {"1e0o": "2:2:1", "1fq9": "2:2:2"}
WIN_START_NS = 11.0          # extension segments 11-50 -> time > 11 ns (time axis = 1 ns prod + ext)
TIGHT = 14                   # 2:2:1 interface-1 tight-state threshold (plan)


def ts(rep):
    rows = [l.split() for l in open(rep / "ext_analysis/timeseries.dat") if not l.startswith("#")]
    return [dict(t=float(r[1]), siteC=float(r[2]), siteD=float(r[3]), allC=float(r[4]),
                 allD=float(r[5]), dist=float(r[6])) for r in rows]


def gb(path):
    """(DELTA TOTAL avg, per-frame series) from a FINAL_RESULTS file."""
    txt = Path(path).read_text()
    m = re.search(r"Differences \(Complex - Receptor - Ligand\):.*?DELTA TOTAL\s+(-?\d+\.\d+)", txt, re.S)
    return float(m.group(1)) if m else None


def block_sem(csv_path, nblocks=4):
    """SEM of the 4 block means (10 ns each) of per-frame DELTA TOTAL, from MMPBSA -eo output."""
    p = Path(csv_path)
    if not p.exists():
        return float("nan")
    lines = p.read_text().splitlines()
    i = next((k for k, l in enumerate(lines) if l.startswith("DELTA Energy Terms")), None)
    if i is None:
        return float("nan")
    vals = []
    for l in lines[i + 2:]:
        if not l.strip() or not l[0].isdigit():
            break
        vals.append(float(l.split(",")[-1]))
    n = len(vals) // nblocks
    if n == 0:
        return float("nan")
    means = [st.fmean(vals[j * n:(j + 1) * n]) for j in range(nblocks)]
    return st.stdev(means) / math.sqrt(nblocks)


def main():
    OUT.mkdir(exist_ok=True)
    rows, complete = [], {}
    for sysname, reps in SEL.items():
        tk, sp = sysname.split("/")
        for r in reps:
            rep = BASE / sysname / "sim" / f"namdrep{r}"
            nseg = len(list(rep.glob("ext[0-9][0-9].dcd")))
            if not (rep / "ext_analysis/timeseries.dat").exists():
                rows.append(dict(tmpl=tk, sp=sp, rep=r, nseg=nseg, status="not analysed"))
                continue
            t = ts(rep)
            w = [x for x in t if x["t"] > WIN_START_NS]
            row = dict(tmpl=tk, sp=sp, rep=r, nseg=nseg, status="complete" if nseg == 50 else f"partial ({nseg} ns)")
            for k in ("siteC", "siteD", "allC", "allD", "dist"):
                row[k] = st.fmean(x[k] for x in w) if w else float("nan")
            if tk == "1e0o" and w:
                row["tight_frac"] = sum(x["siteC"] >= TIGHT for x in w) / len(w)
                states = [x["siteC"] >= TIGHT for x in w]
                row["switches"] = sum(a != b for a, b in zip(states, states[1:]))
            for s in ("iface1", "iface2", "whole"):
                f = rep / f"ext_analysis/{s}_run/FINAL_RESULTS_{s}.dat"
                row[f"dG_{s}"] = gb(f) if f.exists() else float("nan")
                row[f"blockSEM_{s}"] = block_sem(rep / f"ext_analysis/{s}_run/FRAMES_{s}.csv")
            rows.append(row)
            complete.setdefault((tk, sp), [])
            if nseg == 50:
                complete[(tk, sp)].append(row)

    keys = ["tmpl", "sp", "rep", "nseg", "status", "siteC", "siteD", "allC", "allD", "dist",
            "tight_frac", "switches", "dG_iface1", "dG_iface2", "dG_whole",
            "blockSEM_iface1", "blockSEM_iface2", "blockSEM_whole"]
    with open(OUT / "replicate_table.csv", "w") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(f"{r[k]:.2f}" if isinstance(r.get(k), float) else str(r.get(k, "")) for k in keys) + "\n")

    lines = ["Goldfish vs human, replicate means over extension segments 11-50 (complete replicates only)\n"]
    for tk in ("1e0o", "1fq9"):
        g, h = complete.get((tk, "goldfish"), []), complete.get((tk, "human"), [])
        lines.append(f"== {STOICH[tk]} ({tk}): n goldfish = {len(g)}, n human = {len(h)}")
        for k, label in (("dG_whole", "ΔG whole (kcal/mol)"), ("dG_iface1", "ΔG interface 1"), ("dG_iface2", "ΔG interface 2"),
                         ("siteC", "heparin-site salt bridges, FGF C"), ("siteD", "heparin-site salt bridges, FGF D"),
                         ("dist", "FGFR-FGFR distance (Å)"), ("tight_frac", "fraction tight (2:2:1)")):
            a = [x[k] for x in g if k in x and not math.isnan(x[k])]
            b = [x[k] for x in h if k in x and not math.isnan(x[k])]
            if len(a) < 2 or len(b) < 2:
                lines.append(f"  {label:38s} goldfish {[round(v,1) for v in a]}  human {[round(v,1) for v in b]}  (n<2, no test)")
                continue
            d = st.fmean(a) - st.fmean(b)
            se = math.sqrt(st.variance(a) / len(a) + st.variance(b) / len(b))
            res = stats.ttest_ind(a, b, equal_var=False)
            dfw = (st.variance(a)/len(a) + st.variance(b)/len(b))**2 / ((st.variance(a)/len(a))**2/(len(a)-1) + (st.variance(b)/len(b))**2/(len(b)-1)) if se > 0 else 1
            ci = stats.t.ppf(0.975, dfw) * se
            lines.append(f"  {label:38s} goldfish {st.fmean(a):8.1f}  human {st.fmean(b):8.1f}  diff {d:+7.1f} "
                         f"[95% CI {d-ci:+.1f}, {d+ci:+.1f}]  Welch p = {res.pvalue:.3f}")
            if tk == "1e0o" and k == "dG_whole":
                verdict = "SUPPORTED" if (d < 0 and res.pvalue < 0.05) else "NOT SUPPORTED (not resolved at this n)" if d < 0 else "NOT SUPPORTED (direction opposite or zero)"
                lines.append(f"  >> pre-specified test of the v1 claim (goldfish lower ΔG than human, 2:2:1 whole): {verdict}")
        lines.append("")
    pend = [f"{r['tmpl']}/{r['sp']} rep{r['rep']}: {r['status']}" for r in rows if r.get("status") != "complete"]
    if pend:
        lines.append("Not in the primary analysis yet: " + "; ".join(pend))
    (OUT / "comparison.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    for tk in ("1fq9", "1e0o"):
        reps = [(sysname, r) for sysname, rr in SEL.items() if sysname.startswith(tk) for r in rr]
        fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
        for sysname, r in reps:
            rep = BASE / sysname / "sim" / f"namdrep{r}"
            if not (rep / "ext_analysis/timeseries.dat").exists():
                continue
            t = ts(rep)
            sp = sysname.split("/")[1]
            col = "#D9822B" if sp == "goldfish" else "#2F4B7C"
            ls = {0: "-", 1: "--", 2: ":"}[[x for x in SEL[sysname]].index(r)]
            tt = [x["t"] for x in t]
            for ax, k in zip(axes, ("siteC", "siteD", "dist")):
                y = [x[k] for x in t]
                k2 = 10  # 1 ns running mean
                ysm = [st.fmean(y[max(0, i - k2 + 1):i + 1]) for i in range(len(y))]
                ax.plot(tt, ysm, color=col, ls=ls, lw=1.4, label=f"{sp} rep{r}")
        axes[0].set_ylabel("heparin-site salt bridges\nFGF C")
        axes[1].set_ylabel("heparin-site salt bridges\nFGF D")
        axes[2].set_ylabel("FGFR–FGFR distance (Å)")
        axes[2].set_xlabel("time (ns)")
        if tk == "1e0o":
            axes[0].axhline(TIGHT, color="grey", lw=0.8, ls="-.")
        for ax in axes:
            ax.axvspan(1, WIN_START_NS, color="#EEEEEE", zorder=0)
        axes[0].legend(fontsize=8, ncol=3)
        axes[0].set_title(f"{STOICH[tk]} ({tk.upper()}) extensions — 1 ns running mean; grey = excluded equilibration")
        fig.tight_layout()
        fig.savefig(OUT / f"timeseries_{tk}.png", dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    main()
