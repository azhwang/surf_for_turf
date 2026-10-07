#!/usr/bin/env python3
"""
compare_species_replicates.py

Compare MM/GBSA binding energy (DELTA TOTAL) across multiple FGF species/
variants, where each species has several independent replicate runs.
Plots each species' across-replicate mean with SEM error bars, the full
min-max range across replicates, and the individual replicate values as
jittered points -- so you can see both the central estimate and how much
replicate-to-replicate variability there is, at a glance.

Each "replicate" here means one full FINAL_RESULTS_*.dat file (i.e. one
independent MD trajectory analyzed with MMPBSA.py) -- this script uses each
replicate's own frame-averaged DELTA TOTAL as a single data point, then
computes statistics ACROSS replicates. This is distinct from (and a better
measure of true variability than) the within-replicate SEM that MMPBSA.py
reports for frame-to-frame noise within one trajectory.

Usage
-----
python3 compare_species_replicates.py \\
    --species WT   wt/rep1/FINAL_RESULTS_iface1.dat wt/rep2/FINAL_RESULTS_iface1.dat \\
                   wt/rep3/FINAL_RESULTS_iface1.dat wt/rep4/FINAL_RESULTS_iface1.dat \\
                   wt/rep5/FINAL_RESULTS_iface1.dat \\
    --species R45A r45a/rep1/FINAL_RESULTS_iface1.dat r45a/rep2/FINAL_RESULTS_iface1.dat \\
                   r45a/rep3/FINAL_RESULTS_iface1.dat r45a/rep4/FINAL_RESULTS_iface1.dat \\
                   r45a/rep5/FINAL_RESULTS_iface1.dat \\
    -o species_replicate_compare.png

# Convenience: auto-discover replicates via a glob pattern instead of listing each path
python3 compare_species_replicates.py \\
    --species-glob WT   "wt/rep*/FINAL_RESULTS_iface1.dat" \\
    --species-glob R45A "r45a/rep*/FINAL_RESULTS_iface1.dat" \\
    -o species_replicate_compare.png

# Plot a different component or section instead of DELTA TOTAL
python3 compare_species_replicates.py \\
    --species WT wt/rep1/FINAL_RESULTS_iface1.dat wt/rep2/FINAL_RESULTS_iface1.dat \\
    --species MUT mut/rep1/FINAL_RESULTS_iface1.dat mut/rep2/FINAL_RESULTS_iface1.dat \\
    --component EGB --section Delta \\
    -o egb_compare.png

Notes
-----
- Bars/points show the MEAN of each replicate's frame-averaged value.
- Bold error bars = SEM across replicates (mean +/- SEM).
- Thin vertical line = full min-max range across replicates.
- Small dots = each individual replicate's value (jittered horizontally
  for visibility).
- A results table (mean, SEM, std, min, max, n) is printed to stdout, and
  a pairwise comparison (relative to the FIRST species given) is printed
  using a simple two-sample comparison of the replicate means.
"""

import argparse
import glob as globmod
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
    print("matplotlib is required: pip install matplotlib --break-system-packages",
          file=sys.stderr)
    raise


SECTION_HEADERS = {
    "Complex": "Complex:",
    "Receptor": "Receptor:",
    "Ligand": "Ligand:",
    "Delta": "Differences (Complex - Receptor - Ligand):",
}
METHOD_HEADERS = {"GENERALIZED BORN:": "GB", "POISSON BOLTZMANN:": "PB"}

SPECIES_COLORS = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
]


# --------------------------------------------------------------------------
# Parsing (self-contained)
# --------------------------------------------------------------------------

@dataclass
class EnergyTable:
    rows: dict = field(default_factory=dict)


@dataclass
class MMPBSAResult:
    path: str
    tables: dict = field(default_factory=dict)


def _parse_energy_line(line):
    parts = line.split()
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
    path = Path(path)
    text = path.read_text()
    result = MMPBSAResult(path=str(path))

    current_method = None
    current_section = None
    all_section_headers = set(SECTION_HEADERS.values())

    for raw_line in text.splitlines():
        stripped = raw_line.strip()

        if stripped in METHOD_HEADERS:
            current_method = METHOD_HEADERS[stripped]
            result.tables.setdefault(current_method, {})
            current_section = None
            continue

        if stripped in all_section_headers:
            if current_method is None:
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


def get_component_avg(result, method, section_key, component_name):
    section_header = SECTION_HEADERS[section_key]
    table = result.tables.get(method, {}).get(section_header)
    if table is None:
        return None
    if component_name in table.rows:
        return table.rows[component_name][0]
    normalized = component_name.replace("_", " ")
    row = table.rows.get(normalized)
    return row[0] if row else None


# --------------------------------------------------------------------------
# Stats + plotting
# --------------------------------------------------------------------------

@dataclass
class SpeciesStats:
    label: str
    values: list  # one value per replicate
    paths: list

    @property
    def n(self):
        return len(self.values)

    @property
    def mean(self):
        return float(np.mean(self.values)) if self.values else float("nan")

    @property
    def std(self):
        return float(np.std(self.values, ddof=1)) if self.n > 1 else 0.0

    @property
    def sem(self):
        return self.std / np.sqrt(self.n) if self.n > 1 else 0.0

    @property
    def vmin(self):
        return float(np.min(self.values)) if self.values else float("nan")

    @property
    def vmax(self):
        return float(np.max(self.values)) if self.values else float("nan")


def collect_species_stats(species_paths, method, section_key, component):
    """species_paths: list of (label, [path1, path2, ...])"""
    stats = []
    for label, paths in species_paths:
        values = []
        used_paths = []
        for p in paths:
            if not Path(p).is_file():
                print(f"Warning: file not found for {label}: {p}", file=sys.stderr)
                continue
            result = parse_final_results(p)
            val = get_component_avg(result, method, section_key, component)
            if val is None:
                print(f"Warning: component '{component}' not found in {p}", file=sys.stderr)
                continue
            values.append(val)
            used_paths.append(p)
        if not values:
            print(f"Warning: no usable replicates found for species '{label}' -- skipping.",
                  file=sys.stderr)
            continue
        stats.append(SpeciesStats(label=label, values=values, paths=used_paths))
    return stats


def print_stats_table(stats_list, component):
    print(f"\n{component.replace('_', ' ')} across replicates, kcal/mol:")
    header = f"{'Species':12s} {'n':>3s} {'mean':>10s} {'SEM':>8s} {'std':>8s} {'min':>10s} {'max':>10s}"
    print(header)
    print("-" * len(header))
    for s in stats_list:
        print(f"{s.label:12s} {s.n:3d} {s.mean:10.3f} {s.sem:8.3f} {s.std:8.3f} "
              f"{s.vmin:10.3f} {s.vmax:10.3f}")

    if len(stats_list) > 1:
        ref = stats_list[0]
        print(f"\nComparison relative to '{ref.label}' (first species given):")
        print(f"{'Species':12s} {'ddG (mean diff)':>18s} {'combined SEM':>14s} {'ddG/SEM':>10s}")
        for s in stats_list[1:]:
            ddg = s.mean - ref.mean
            combined_sem = np.sqrt(s.sem**2 + ref.sem**2)
            ratio = ddg / combined_sem if combined_sem > 0 else float("nan")
            print(f"{s.label:12s} {ddg:18.3f} {combined_sem:14.3f} {ratio:10.2f}")
        print("\n(ddG/SEM is a rough heuristic, not a formal statistical test --")
        print(" see prior discussion on interpreting this given limited/correlated sampling.)")


def plot_species_replicates(stats_list, component, method, section_key, out, title=None,
                             rng=None):
    n = len(stats_list)
    x = np.arange(n)

    fig, ax = plt.subplots(figsize=(max(6, 1.4 * n), 6))

    if rng is None:
        rng = np.random.default_rng(0)

    for i, s in enumerate(stats_list):
        color = SPECIES_COLORS[i % len(SPECIES_COLORS)]

        # thin min-max range line
        ax.plot([x[i], x[i]], [s.vmin, s.vmax], color=color, alpha=0.4, linewidth=2, zorder=1)

        # bold SEM error bar on the mean
        ax.errorbar(x[i], s.mean, yerr=s.sem, fmt="o", color=color,
                     markersize=9, capsize=6, elinewidth=2.2, markeredgecolor="black",
                     markeredgewidth=0.8, zorder=3)

        # individual replicate points, jittered
        jitter = rng.uniform(-0.10, 0.10, size=s.n)
        ax.scatter(x[i] + jitter, s.values, color=color, alpha=0.6, s=28,
                   edgecolor="black", linewidth=0.4, zorder=2)

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([s.label for s in stats_list])
    ax.set_ylabel(f"{component.replace('_', ' ')} (kcal/mol)")
    ax.set_title(title or f"{section_key} {component.replace('_', ' ')} by species "
                           f"({method}, mean +/- SEM across replicates, dots = individual reps)")
    fig.tight_layout()

    fig.savefig(out, dpi=150)
    print(f"\nSaved: {out}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--species", nargs="+", action="append", default=[],
                         metavar="LABEL PATH1 [PATH2 ...]",
                         help="Repeatable: --species LABEL path1 path2 path3 ... "
                              "(first token is the label, remaining tokens are "
                              "replicate FINAL_RESULTS_*.dat paths)")
    parser.add_argument("--species-glob", nargs=2, action="append", default=[],
                         metavar=("LABEL", "GLOB_PATTERN"),
                         help="Repeatable: --species-glob LABEL 'path/rep*/FINAL_RESULTS_*.dat' "
                              "(quote the glob so your shell doesn't expand it first)")
    parser.add_argument("--component", default="DELTA_TOTAL",
                         help="Energy component to compare (default: DELTA_TOTAL)")
    parser.add_argument("--section", choices=["Complex", "Receptor", "Ligand", "Delta"],
                         default="Delta")
    parser.add_argument("--method", choices=["GB", "PB"], default="GB")
    parser.add_argument("--title", default=None)
    parser.add_argument("-o", "--out", default="species_replicate_compare.png")
    args = parser.parse_args()

    species_paths = []
    for entry in args.species:
        label, paths = entry[0], entry[1:]
        if not paths:
            print(f"Warning: species '{label}' given with no replicate paths -- skipping.",
                  file=sys.stderr)
            continue
        species_paths.append((label, paths))

    for label, pattern in args.species_glob:
        paths = sorted(globmod.glob(pattern, recursive=True))
        if not paths:
            print(f"Warning: glob pattern for '{label}' matched no files: {pattern}",
                  file=sys.stderr)
            continue
        species_paths.append((label, paths))

    if not species_paths:
        print("No species specified. Use --species LABEL PATH... or "
              "--species-glob LABEL PATTERN.", file=sys.stderr)
        sys.exit(1)

    stats_list = collect_species_stats(species_paths, args.method, args.section, args.component)
    if not stats_list:
        print("No usable data found for any species.", file=sys.stderr)
        sys.exit(1)

    print_stats_table(stats_list, args.component)
    plot_species_replicates(stats_list, args.component, args.method, args.section,
                             args.out, title=args.title)


if __name__ == "__main__":
    main()
