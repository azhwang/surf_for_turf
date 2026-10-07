#!/usr/bin/env python3
"""
analyze_mmpbsa.py

Parse and analyze AmberTools MMPBSA.py output files:
  - FINAL_RESULTS_MMPBSA.dat  (energy component summary: complex/receptor/ligand/delta)
  - FINAL_DECOMP_MMPBSA.dat   (per-residue energy decomposition)

Supports comparing multiple runs (e.g. iface1 vs iface2 vs whole-assembly,
or wild-type vs mutant) side by side.

Usage examples
--------------
# Summarize + plot one run's energy components
python analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat

# Compare DELTA TOTAL across several runs (bar chart with error bars)
python analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o compare_total.png

# Compare WT vs mutant DELTA TOTAL with a two-sample significance check
python analyze_mmpbsa.py ddg \
    --wt   wt_iface1_run/FINAL_RESULTS_iface1.dat \
    --mut  mut_iface1_run/FINAL_RESULTS_iface1.dat

# Per-residue decomposition plot (single run, residues near an interface)
python analyze_mmpbsa.py decomp iface1_run/FINAL_DECOMP_MMPBSA.dat \
    --resrange 1 60 -o decomp_iface1.png

# Per-residue decomposition WT vs mutant, highlighting one residue
python analyze_mmpbsa.py decomp-compare \
    --wt  wt_iface1_run/FINAL_DECOMP_MMPBSA.dat \
    --mut mut_iface1_run/FINAL_DECOMP_MMPBSA.dat \
    --highlight 45 \
    -o decomp_compare.png

Notes
-----
- All energies are kcal/mol, as reported by MMPBSA.py.
- "DELTA TOTAL" is the MM+solvation binding enthalpy estimate (no entropy
  term unless you ran &nmode / interaction-entropy -- see caveats printed
  by the script).
- Works for both GB- and PB-only result files, and for files containing
  both (it will parse whichever section(s) are present).
"""

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


# --------------------------------------------------------------------------
# Parsing: FINAL_RESULTS_MMPBSA.dat
# --------------------------------------------------------------------------

SECTION_HEADERS = ("Complex:", "Receptor:", "Ligand:", "Differences (Complex - Receptor - Ligand):")
METHOD_HEADERS = {"GENERALIZED BORN:": "GB", "POISSON BOLTZMANN:": "PB"}


@dataclass
class EnergyTable:
    """One energy component table (e.g. 'Complex' under 'GB')."""
    rows: dict = field(default_factory=dict)  # component -> (avg, std, sem)

    def get(self, component):
        return self.rows.get(component)


@dataclass
class MMPBSAResult:
    path: str
    receptor_mask: str = ""
    ligand_mask: str = ""
    nframes: float = None
    # method -> section -> EnergyTable
    tables: dict = field(default_factory=dict)

    def delta_total(self, method="GB"):
        """Return (avg, std, sem) for DELTA TOTAL under the given method."""
        sec = self.tables.get(method, {}).get("Differences (Complex - Receptor - Ligand):")
        if sec is None:
            return None
        return sec.get("DELTA TOTAL")


def _parse_energy_line(line):
    """Parse a line like 'EGB   -5095.1198   64.3752   20.3572' -> (name, avg, std, sem)."""
    parts = line.split()
    # component names can contain spaces (e.g. '1-4 VDW', 'G gas', 'DELTA TOTAL')
    # so parse from the right: last 3 tokens are numbers, rest is the name.
    if len(parts) < 4:
        return None
    try:
        sem = float(parts[-1])
        std = float(parts[-2])
        avg = float(parts[-3])
    except ValueError:
        return None
    name = " ".join(parts[:-3])
    return name, avg, std, sem


def parse_final_results(path):
    """Parse a MMPBSA.py FINAL_RESULTS_*.dat file into a MMPBSAResult."""
    path = Path(path)
    text = path.read_text()
    result = MMPBSAResult(path=str(path))

    m = re.search(r"Receptor mask:\s*\"([^\"]*)\"", text)
    if m:
        result.receptor_mask = m.group(1)
    m = re.search(r"Ligand mask:\s*\"([^\"]*)\"", text)
    if m:
        result.ligand_mask = m.group(1)
    m = re.search(r"Calculations performed using\s+([\d.]+)\s+complex frames", text)
    if m:
        result.nframes = float(m.group(1))

    lines = text.splitlines()
    current_method = None
    current_section = None

    for line in lines:
        stripped = line.strip()

        if stripped in METHOD_HEADERS:
            current_method = METHOD_HEADERS[stripped]
            result.tables.setdefault(current_method, {})
            current_section = None
            continue

        if stripped in SECTION_HEADERS:
            if current_method is None:
                # Some MMPBSA.py versions omit the GB/PB header if only one is run
                current_method = "GB"
                result.tables.setdefault(current_method, {})
            current_section = stripped
            result.tables[current_method][current_section] = EnergyTable()
            continue

        if current_method is not None and current_section is not None:
            if not stripped or stripped.startswith("-") or stripped.startswith("Energy Component"):
                continue
            parsed = _parse_energy_line(stripped)
            if parsed:
                name, avg, std, sem = parsed
                result.tables[current_method][current_section].rows[name] = (avg, std, sem)

    return result


# --------------------------------------------------------------------------
# Parsing: FINAL_DECOMP_MMPBSA.dat
# --------------------------------------------------------------------------

@dataclass
class DecompResult:
    path: str
    # residue_id (str, e.g. 'PROA:45') -> dict(component -> (avg, std, sem))
    complex_res: dict = field(default_factory=dict)


def parse_decomp(path):
    """
    Parse MMPBSA.py per-residue decomposition output (idecomp=1 or 2).

    Handles the 'Total Energy Decomposition' complex section, extracting
    each residue's Total avg/std/sem. Works with the CSV-like or
    whitespace/comma-delimited formats MMPBSA.py emits depending on version.
    """
    path = Path(path)
    text = path.read_text()
    result = DecompResult(path=str(path))

    # MMPBSA.py decomp output is typically comma-separated per-residue rows like:
    #   Resid 1,Resid 2,...,Internal,...,Total
    # or, for per-residue (not pairwise) idecomp=1/2 "Total Decomp" table:
    #   Residue, Location, Internal, van der Waals, Electrostatic, Polar Solv, Non-Polar Solv, TOTAL
    # We look for the "Total Energy Decomposition" / "Complex:" block and
    # grab lines that start with a residue-like token followed by numbers.

    in_complex_block = False
    header_cols = None

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        s = line.strip()

        if not s:
            continue

        if s.startswith("Complex"):
            in_complex_block = True
            header_cols = None
            continue
        if s.startswith("Receptor") or s.startswith("Ligand"):
            in_complex_block = False
            continue

        if not in_complex_block:
            continue

        # header line (contains commas and non-numeric first token)
        tokens = [t.strip() for t in s.split(",")] if "," in s else s.split()
        if header_cols is None:
            # detect header: first meaningful row containing "Total" or "Resid" or similar text
            if any(k in s for k in ("Total", "Resid", "Internal", "van der Waals")):
                header_cols = tokens
                continue

        if len(tokens) < 2:
            continue

        resid_token = tokens[0]
        # residue rows usually look like "PROA_45" / "45" / "1  ALA" etc. -- keep whatever it is
        nums = []
        for t in tokens[1:]:
            try:
                nums.append(float(t))
            except ValueError:
                nums.append(None)
        # heuristically: last 3 numeric columns of a decomp row are avg/std/sem of TOTAL
        numeric_vals = [n for n in nums if n is not None]
        if len(numeric_vals) >= 3:
            avg, std, sem = numeric_vals[-3], numeric_vals[-2], numeric_vals[-1]
            result.complex_res[resid_token] = {"TOTAL": (avg, std, sem)}

    return result


# --------------------------------------------------------------------------
# Reporting / plotting
# --------------------------------------------------------------------------

def print_summary(result: MMPBSAResult):
    print(f"\n=== {result.path} ===")
    if result.receptor_mask or result.ligand_mask:
        print(f"Receptor mask: {result.receptor_mask}   Ligand mask: {result.ligand_mask}")
    if result.nframes:
        print(f"Frames used: {result.nframes:.0f}")

    for method, sections in result.tables.items():
        diff = sections.get("Differences (Complex - Receptor - Ligand):")
        if diff is None:
            continue
        print(f"\n-- {method} --")
        for comp in ("DELTA G gas", "DELTA G solv", "DELTA TOTAL"):
            row = diff.get(comp)
            if row:
                avg, std, sem = row
                print(f"  {comp:15s}: {avg:9.3f}  (std {std:7.3f}, SEM {sem:6.3f})")

    print("\nCaveat: this is MM + solvation enthalpy only. No entropy term is")
    print("included unless you ran &nmode / interaction-entropy separately.")


def plot_energy_components(result: MMPBSAResult, method="GB", out=None):
    if plt is None:
        print("matplotlib not available; cannot plot.", file=sys.stderr)
        return
    diff = result.tables.get(method, {}).get("Differences (Complex - Receptor - Ligand):")
    if diff is None:
        print(f"No '{method}' Differences section found in {result.path}", file=sys.stderr)
        return

    components = ["VDWAALS", "EEL", "EGB" if method == "GB" else "EPB", "ESURF",
                  "DELTA G gas", "DELTA G solv", "DELTA TOTAL"]
    avgs, sems, labels = [], [], []
    for c in components:
        row = diff.get(c)
        if row is None:
            continue
        avgs.append(row[0])
        sems.append(row[2])
        labels.append(c)

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#4C72B0" if "DELTA" not in l else "#C44E52" for l in labels]
    ax.bar(labels, avgs, yerr=sems, capsize=4, color=colors)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_title(f"Energy components — {Path(result.path).stem} ({method})")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()

    if out:
        fig.savefig(out, dpi=150)
        print(f"Saved: {out}")
    else:
        out = str(Path(result.path).with_suffix(".png"))
        fig.savefig(out, dpi=150)
        print(f"Saved: {out}")
    plt.close(fig)


def compare_delta_total(labeled_paths, method="GB", out=None):
    """labeled_paths: list of (label, path) tuples."""
    labels, avgs, sems = [], [], []
    for label, path in labeled_paths:
        r = parse_final_results(path)
        row = r.delta_total(method=method)
        if row is None:
            print(f"Warning: no DELTA TOTAL ({method}) found in {path}", file=sys.stderr)
            continue
        avg, std, sem = row
        labels.append(label)
        avgs.append(avg)
        sems.append(sem)
        print(f"{label:15s}  DELTA TOTAL = {avg:9.3f} +/- {sem:6.3f} (SEM)  [{path}]")

    if plt is None:
        print("matplotlib not available; skipping plot.", file=sys.stderr)
        return

    fig, ax = plt.subplots(figsize=(1.2 * len(labels) + 2, 5))
    ax.bar(labels, avgs, yerr=sems, capsize=5, color="#4C72B0")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("ΔG (MM+solvation, kcal/mol)")
    ax.set_title("MM/GBSA binding enthalpy comparison")
    fig.tight_layout()

    out = out or "compare_delta_total.png"
    fig.savefig(out, dpi=150)
    print(f"Saved: {out}")
    plt.close(fig)


def wt_mut_ddg(wt_path, mut_path, method="GB"):
    """Compute WT vs mutant DELTA TOTAL, the difference (ddG), and a simple
    z-like significance estimate assuming independent SEMs."""
    wt = parse_final_results(wt_path)
    mut = parse_final_results(mut_path)

    wt_row = wt.delta_total(method=method)
    mut_row = mut.delta_total(method=method)
    if wt_row is None or mut_row is None:
        print("Could not find DELTA TOTAL in one or both files.", file=sys.stderr)
        return

    wt_avg, wt_std, wt_sem = wt_row
    mut_avg, mut_std, mut_sem = mut_row

    ddg = mut_avg - wt_avg
    combined_sem = np.sqrt(wt_sem**2 + mut_sem**2)
    # rough z-score; NOT a rigorous test given small, correlated MD samples --
    # treat as a heuristic for "is this difference plausibly above the noise floor"
    z = ddg / combined_sem if combined_sem > 0 else float("nan")

    print(f"\nWT  ({wt_path}):  DELTA TOTAL = {wt_avg:9.3f} +/- {wt_sem:6.3f} (SEM)")
    print(f"MUT ({mut_path}):  DELTA TOTAL = {mut_avg:9.3f} +/- {mut_sem:6.3f} (SEM)")
    print(f"\nddG (mut - wt) = {ddg:9.3f} kcal/mol")
    print(f"Combined SEM   = {combined_sem:9.3f} kcal/mol")
    print(f"ddG / SEM      = {z:9.2f}")

    print("\nInterpretation guide (heuristic, not a formal p-value):")
    print("  |ddG/SEM| < 1   -> difference is within noise; not resolved with this sampling")
    print("  |ddG/SEM| 1-2   -> suggestive, but weak given small/correlated MD sample sizes")
    print("  |ddG/SEM| > 2-3 -> more likely a real effect, but confirm with more frames")
    print("  Interpretation should ALWAYS be paired with per-residue decomposition")
    print("  at the mutated position to check the mechanism makes physical sense.")

    return ddg, combined_sem, z


def plot_decomp(decomp: DecompResult, resrange=None, highlight=None, out=None, title=None):
    if plt is None:
        print("matplotlib not available; cannot plot.", file=sys.stderr)
        return

    items = list(decomp.complex_res.items())
    if not items:
        print(f"No residue rows parsed from {decomp.path}. The decomposition file "
              f"format may differ from what this parser expects -- inspect the raw "
              f"file and adjust parse_decomp() if needed.", file=sys.stderr)
        return

    # try to sort numerically by residue id if possible
    def sort_key(kv):
        m = re.search(r"(\d+)", kv[0])
        return int(m.group(1)) if m else kv[0]
    items.sort(key=sort_key)

    if resrange:
        lo, hi = resrange
        def in_range(kv):
            m = re.search(r"(\d+)", kv[0])
            return m and lo <= int(m.group(1)) <= hi
        items = [kv for kv in items if in_range(kv)]

    labels = [k for k, _ in items]
    avgs = [v["TOTAL"][0] for _, v in items]
    sems = [v["TOTAL"][2] for _, v in items]

    colors = ["#4C72B0"] * len(labels)
    if highlight is not None:
        for i, lab in enumerate(labels):
            if str(highlight) in lab:
                colors[i] = "#C44E52"

    fig, ax = plt.subplots(figsize=(max(8, 0.25 * len(labels)), 5))
    ax.bar(labels, avgs, yerr=sems, capsize=2, color=colors)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Per-residue ΔG contribution (kcal/mol)")
    ax.set_title(title or f"Per-residue decomposition — {Path(decomp.path).stem}")
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    fig.tight_layout()

    out = out or str(Path(decomp.path).with_suffix(".png"))
    fig.savefig(out, dpi=150)
    print(f"Saved: {out}")
    plt.close(fig)


def plot_decomp_compare(wt: DecompResult, mut: DecompResult, highlight=None, out=None):
    if plt is None:
        print("matplotlib not available; cannot plot.", file=sys.stderr)
        return

    common = sorted(set(wt.complex_res) & set(mut.complex_res),
                     key=lambda k: int(re.search(r"(\d+)", k).group(1)) if re.search(r"(\d+)", k) else 0)
    if not common:
        print("No overlapping residue IDs found between WT and mutant decomposition "
              "files -- check that residue numbering matches between the two runs.",
              file=sys.stderr)
        return

    wt_avgs = [wt.complex_res[k]["TOTAL"][0] for k in common]
    mut_avgs = [mut.complex_res[k]["TOTAL"][0] for k in common]
    wt_sems = [wt.complex_res[k]["TOTAL"][2] for k in common]
    mut_sems = [mut.complex_res[k]["TOTAL"][2] for k in common]

    x = np.arange(len(common))
    width = 0.4

    fig, ax = plt.subplots(figsize=(max(10, 0.3 * len(common)), 5))
    ax.bar(x - width / 2, wt_avgs, width, yerr=wt_sems, capsize=2, label="WT", color="#4C72B0")
    ax.bar(x + width / 2, mut_avgs, width, yerr=mut_sems, capsize=2, label="Mutant", color="#DD8452")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(common, rotation=90, fontsize=7)
    ax.set_ylabel("Per-residue ΔG contribution (kcal/mol)")
    ax.set_title("Per-residue decomposition: WT vs Mutant")
    ax.legend()

    if highlight is not None:
        for i, lab in enumerate(common):
            if str(highlight) in lab:
                ax.axvspan(i - 0.5, i + 0.5, color="yellow", alpha=0.2)

    fig.tight_layout()
    out = out or "decomp_compare.png"
    fig.savefig(out, dpi=150)
    print(f"Saved: {out}")
    plt.close(fig)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_sum = sub.add_parser("summary", help="Print + plot one FINAL_RESULTS file")
    p_sum.add_argument("path")
    p_sum.add_argument("--method", default="GB", choices=["GB", "PB"])
    p_sum.add_argument("-o", "--out", default=None)

    p_cmp = sub.add_parser("compare", help="Compare DELTA TOTAL across several runs")
    p_cmp.add_argument("--label", nargs=2, action="append", required=True,
                        metavar=("LABEL", "PATH"),
                        help="Repeatable: --label NAME path/to/FINAL_RESULTS.dat")
    p_cmp.add_argument("--method", default="GB", choices=["GB", "PB"])
    p_cmp.add_argument("-o", "--out", default=None)

    p_ddg = sub.add_parser("ddg", help="WT vs mutant DELTA TOTAL comparison")
    p_ddg.add_argument("--wt", required=True)
    p_ddg.add_argument("--mut", required=True)
    p_ddg.add_argument("--method", default="GB", choices=["GB", "PB"])

    p_dec = sub.add_parser("decomp", help="Plot per-residue decomposition for one run")
    p_dec.add_argument("path")
    p_dec.add_argument("--resrange", nargs=2, type=int, default=None, metavar=("LO", "HI"))
    p_dec.add_argument("--highlight", default=None, help="Residue number to highlight")
    p_dec.add_argument("-o", "--out", default=None)

    p_decc = sub.add_parser("decomp-compare", help="Plot WT vs mutant per-residue decomposition")
    p_decc.add_argument("--wt", required=True)
    p_decc.add_argument("--mut", required=True)
    p_decc.add_argument("--highlight", default=None, help="Residue number to highlight")
    p_decc.add_argument("-o", "--out", default=None)

    args = parser.parse_args()

    if args.cmd == "summary":
        r = parse_final_results(args.path)
        print_summary(r)
        plot_energy_components(r, method=args.method, out=args.out)

    elif args.cmd == "compare":
        compare_delta_total([(lab, path) for lab, path in args.label],
                             method=args.method, out=args.out)

    elif args.cmd == "ddg":
        wt_mut_ddg(args.wt, args.mut, method=args.method)

    elif args.cmd == "decomp":
        d = parse_decomp(args.path)
        plot_decomp(d, resrange=args.resrange, highlight=args.highlight, out=args.out)

    elif args.cmd == "decomp-compare":
        wt = parse_decomp(args.wt)
        mut = parse_decomp(args.mut)
        plot_decomp_compare(wt, mut, highlight=args.highlight, out=args.out)


if __name__ == "__main__":
    main()
