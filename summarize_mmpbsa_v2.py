#!/usr/bin/env python3
"""
summarize_mmpbsa_v2.py

Compile every v2 MM/GBSA result (namdrepN/mmpbsa_v2/*_run/FINAL_RESULTS_*.dat)
into one workbook plus per-split species comparison plots. Only v2 outputs
are read; v1 values appear solely in the "v1 vs v2" sheet for reference.

Usage (from ~/Desktop/Surf_for_Turf):
    python3 summarize_mmpbsa_v2.py -o mmpbsa_v2_results

Output:
    mmpbsa_v2_results/mmpbsa_v2_summary.xlsx
        Species Summary  -- mean / SEM / std / min / max / n of DELTA TOTAL
                            across replicates, per stoichiometry x split x species
        Replicates       -- one row per replicate x split: DELTA TOTAL and
                            VDW / EEL / EGB / ESURF (GB, Delta section)
        v1 vs v2         -- per replicate x split, where both exist
        Coverage + Notes -- what's missing / whole-only and why
    mmpbsa_v2_results/<stoich>_<split>_delta_total.png
"""

import argparse
import csv
import math
import re
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill

ROOTS = {"2:2:2 (1FQ9)": "all_metrics_analysis", "2:2:1 (1E0O)": "all_metrics_analysis_221"}
SPECIES_ORDER = ["human", "human_stab", "bovine_stab", "chinook_salmon",
                 "goldfish", "japanese_ricefish", "zebrafish"]
SPLITS = ["iface1", "iface2", "whole"]
COMPONENTS = ["VDWAALS", "EEL", "EGB", "ESURF", "DELTA G gas", "DELTA G solv", "DELTA TOTAL"]
FONT = "Arial"


def parse_delta(path):
    """{component: (avg, std, sem)} from the Differences section (first method = GB)."""
    out = {}
    in_diff = False
    for line in Path(path).read_text().splitlines():
        if line.startswith("Differences (Complex - Receptor - Ligand):"):
            in_diff = True
            continue
        if not in_diff:
            continue
        m = re.match(r"^(\S.*?)\s{2,}(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", line)
        if m and m.group(1).strip() in COMPONENTS:
            out[m.group(1).strip()] = tuple(float(m.group(i)) for i in (2, 3, 4))
            if m.group(1).strip() == "DELTA TOTAL":
                break
    return out


QC_FILE = {"original": "qc_ring_piercing_v2.csv", "rebuild": "qc_ring_piercing_rebuild.csv"}
REBUILD_ROOTS = {"2:2:2 (1FQ9)": "rebuild_v2/1fq9", "2:2:1 (1E0O)": "rebuild_v2/1e0o"}


def load_qc(base, which="original"):
    """{(stoich_label, species, split): 'residue list'} from the QC csv for this set."""
    f = base / QC_FILE[which]
    flags = {}
    if not f.exists():
        return flags
    label = {"2:2:2": "2:2:2 (1FQ9)", "2:2:1": "2:2:1 (1E0O)"}
    for r in csv.DictReader(open(f)):
        key = (label[r["stoich"]], r["species"], r["split"])
        res = f"{r['location'].split()[1]}{r['location'].split()[2]}({r['location'][0]})"
        flags.setdefault(key, set()).add(res)
    return {k: ", ".join(sorted(v)) for k, v in flags.items()}


def species_of(dirname):
    return dirname.split("_", 1)[1]           # "1fq9_human_stab" -> "human_stab"


def replicate_dirs(base, which):
    """Yield (stoich, species, replicate dir) for the original builds or rebuild_v2."""
    if which == "rebuild":
        for stoich, root in REBUILD_ROOTS.items():
            for sp_dir in sorted((base / root).iterdir()):
                for rep in sorted(sp_dir.glob("sim/namdrep[0-9]")):
                    yield stoich, sp_dir.name, rep
        return
    for stoich, root in ROOTS.items():
        for sp_dir in sorted((base / root).glob("1*_*")):
            if sp_dir.is_dir():
                for rep in sorted(sp_dir.glob("*/namdrep*")):
                    yield stoich, species_of(sp_dir.name), rep


def collect(base, which="original"):
    rows, coverage = [], []
    for stoich, sp, rep in replicate_dirs(base, which):
            if True:
                v2 = rep / "mmpbsa_v2"
                summ = v2 / "split_summary.txt"
                whole_only = summ.exists() and "SPLIT SKIPPED" in summ.read_text()
                have = [s for s in SPLITS if (v2 / f"{s}_run" / f"FINAL_RESULTS_{s}.dat").exists()]
                if not have:
                    reason = ("no production trajectory (step5_production.dcd)"
                              if not (rep / "step5_production.dcd").exists() else "MMPBSA not run")
                    coverage.append((stoich, sp, rep.name, "missing", reason))
                    continue
                coverage.append((stoich, sp, rep.name,
                                 "whole-only" if whole_only else "all splits",
                                 "ambiguous FGF-FGFR pairing" if whole_only else ""))
                for s in have:
                    d = parse_delta(v2 / f"{s}_run" / f"FINAL_RESULTS_{s}.dat")
                    v1f = rep / f"{s}_run" / f"FINAL_RESULTS_{s}.dat"
                    v1 = parse_delta(v1f).get("DELTA TOTAL", (None,))[0] if v1f.exists() else None
                    rows.append(dict(stoich=stoich, species=sp, rep=rep.name, split=s,
                                     d=d, v1=v1, path=str(v2 / f"{s}_run")))
    return rows, coverage


def stats(vals):
    n = len(vals)
    mean = statistics.fmean(vals)
    sd = statistics.stdev(vals) if n > 1 else float("nan")
    return mean, (sd / math.sqrt(n) if n > 1 else float("nan")), sd, min(vals), max(vals), n


def header(ws, cols):
    ws.append(cols)
    for c in ws[1]:
        c.font = Font(name=FONT, bold=True, color="FFFFFF")
        c.fill = PatternFill("solid", fgColor="4472C4")
        c.alignment = Alignment(horizontal="center", wrap_text=True)
    ws.freeze_panes = "A2"


def finish(ws, widths):
    for col, w in zip("ABCDEFGHIJKLMNOP", widths):
        ws.column_dimensions[col].width = w
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = Font(name=FONT)
            if isinstance(c.value, float):
                c.number_format = "0.00"


def ordered(keys):
    return sorted(keys, key=lambda k: (list(ROOTS).index(k[0]), SPLITS.index(k[1]),
                                       SPECIES_ORDER.index(k[2]) if k[2] in SPECIES_ORDER else 99))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=".")
    ap.add_argument("--set", choices=["original", "rebuild"], default="original",
                    help="original CHARMM-GUI builds, or the consistent rebuild in rebuild_v2/")
    ap.add_argument("-o", "--out", default=None,
                    help="default: mmpbsa_v2_results (original) / mmpbsa_rebuild_results (rebuild)")
    args = ap.parse_args()
    base = Path(args.base)
    out = Path(args.out or ("mmpbsa_v2_results" if args.set == "original" else "mmpbsa_rebuild_results"))
    out.mkdir(exist_ok=True)

    rows, coverage = collect(base, args.set)
    qc = load_qc(base, args.set)
    groups = {}
    for r in rows:
        if "DELTA TOTAL" in r["d"]:
            groups.setdefault((r["stoich"], r["split"], r["species"]), []).append(r)

    wb = Workbook()
    ws = wb.active
    ws.title = "Species Summary"
    header(ws, ["Stoichiometry", "Split", "Species", "Mean ΔG total (kcal/mol)",
                "SEM (across reps)", "Std dev", "Min", "Max", "n reps",
                "Mean VDW", "Mean EEL", "Mean EGB", "Mean ESURF", "QC"])
    for key in ordered(groups):
        g = groups[key]
        m, sem, sd, lo, hi, n = stats([r["d"]["DELTA TOTAL"][0] for r in g])
        comp = [statistics.fmean(r["d"][c][0] for r in g) for c in ("VDWAALS", "EEL", "EGB", "ESURF")]
        flag = qc.get((key[0], key[2], key[1]))
        ws.append([*key, m, sem, sd, lo, hi, n, *comp,
                   f"INVALID: ring piercing {flag}" if flag else "ok"])
        if flag:
            for c in ws[ws.max_row]:
                c.fill = PatternFill("solid", fgColor="F8CBAD")
    finish(ws, [16, 8, 18, 14, 11, 10, 10, 10, 7, 11, 11, 11, 11, 60])

    ws = wb.create_sheet("Replicates")
    header(ws, ["Stoichiometry", "Species", "Replicate", "Split", "ΔG total",
                "ΔG total SEM (frames)", "VDW", "EEL", "EGB", "ESURF", "ΔG gas", "ΔG solv", "Path"])
    for r in sorted(rows, key=lambda r: (list(ROOTS).index(r["stoich"]),
                                         SPECIES_ORDER.index(r["species"]), r["rep"], SPLITS.index(r["split"]))):
        d = r["d"]
        ws.append([r["stoich"], r["species"], r["rep"], r["split"],
                   d["DELTA TOTAL"][0], d["DELTA TOTAL"][2],
                   *(d[c][0] for c in ("VDWAALS", "EEL", "EGB", "ESURF", "DELTA G gas", "DELTA G solv")),
                   r["path"]])
    finish(ws, [16, 18, 10, 8, 10, 12, 10, 11, 11, 10, 11, 11, 70])

    ws = wb.create_sheet("v1 vs v2")
    header(ws, ["Stoichiometry", "Species", "Replicate", "Split", "v1 ΔG total (invalid split)",
                "v2 ΔG total", "v2 − v1"])
    for r in rows:
        if r["v1"] is not None:
            v2 = r["d"]["DELTA TOTAL"][0]
            ws.append([r["stoich"], r["species"], r["rep"], r["split"], r["v1"], v2, v2 - r["v1"]])
    finish(ws, [16, 18, 10, 8, 16, 12, 10])

    ws = wb.create_sheet("Coverage")
    header(ws, ["Stoichiometry", "Species", "Replicate", "Status", "Reason"])
    for c in coverage:
        ws.append(list(c))
    finish(ws, [16, 18, 10, 12, 50])

    ws = wb.create_sheet("Notes")
    for line in [
        "MM/GBSA v2 results (run_mmpbsa_pipeline_v2.sh). GB igb=5, saltcon 0.150 M, mbondi2, 10 frames/replicate.",
        "Ligand = FGF (PROC/PROD). Receptor = FGFR (PROA/PROB) + heparin. 2:2:1: shared heparin kept in both interfaces.",
        "ΔG total = MM + GB solvation enthalpy only; no entropy term. Use for relative ranking between species.",
        "iface1/iface2 and whole are not additive under GB. Compare species within one split, not across splits.",
        "'SEM (across reps)' is replicate-to-replicate SEM; the per-replicate 'SEM (frames)' understates true noise.",
        "2:2:2 bovine_stab, zebrafish, human_stab: each FGF contacts both FGFRs; per-interface split skipped (whole only).",
        "Those three 2:2:2 assemblies likely deviate from the 1FQ9 arrangement. Treat their whole-assembly values with caution.",
        "QC: rows marked INVALID have a side chain threaded through an FGF aromatic ring (build artifact; see qc_ring_piercing_v2.csv).",
        "Their ΔG is dominated by hundreds of kcal/mol of fake vdW repulsion and cannot be used. They are excluded from the plots.",
        "2:2:2 chinook_salmon: equilibration crashed in all 5 replicates (RATTLE failure), so there is no production data.",
        "v1 values in 'v1 vs v2' are from the old hardcoded-range pipeline and are invalid (see CLAUDE.md). Kept for reference only.",
    ]:
        ws.append([line])
    ws.column_dimensions["A"].width = 120
    for row in ws.iter_rows():
        for c in row:
            c.font = Font(name=FONT)

    xlsx = out / "mmpbsa_v2_summary.xlsx"
    wb.save(xlsx)
    print(f"wrote {xlsx}")

    # --- plots: one per stoichiometry x split, species on x ---
    for stoich in ROOTS:
        for split in SPLITS:
            keys = [k for k in ordered(groups) if k[0] == stoich and k[1] == split
                    and (k[0], k[2], k[1]) not in qc]
            dropped = [k[2] for k in ordered(groups) if k[0] == stoich and k[1] == split
                       and (k[0], k[2], k[1]) in qc]
            if not keys:
                continue
            fig, ax = plt.subplots(figsize=(8, 4.5))
            for x, k in enumerate(keys):
                vals = [r["d"]["DELTA TOTAL"][0] for r in groups[k]]
                m, sem, _, lo, hi, n = stats(vals)
                ax.bar(x, m, color="#4472C4", alpha=0.35, width=0.6)
                ax.vlines(x, lo, hi, color="#444444", lw=1)
                ax.errorbar(x, m, yerr=0 if math.isnan(sem) else sem, fmt="none",
                            ecolor="black", elinewidth=2, capsize=6)
                jit = [x + (i - (len(vals) - 1) / 2) * 0.06 for i in range(len(vals))]
                ax.scatter(jit, vals, s=14, color="#1F3864", zorder=3)
                ax.text(x, hi, f"n={n}", ha="center", va="bottom", fontsize=8)
            ax.set_xticks(range(len(keys)))
            ax.set_xticklabels([k[2] for k in keys], rotation=30, ha="right")
            ax.set_ylabel("ΔG total (kcal/mol, GB)")
            ax.axhline(0, color="grey", lw=0.5)
            ax.set_title(f"{stoich} — {split} — mean ± SEM across replicates (v2)")
            if dropped:
                ax.text(0.01, 0.01, "excluded (ring-piercing build artifact): " + ", ".join(dropped),
                        transform=ax.transAxes, fontsize=7, color="#C00000")
            fig.tight_layout()
            png = out / f"{stoich.split()[0].replace(':', '')}_{split}_delta_total.png"
            fig.savefig(png, dpi=150)
            plt.close(fig)
            print(f"wrote {png}")


if __name__ == "__main__":
    main()
