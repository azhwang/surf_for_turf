#!/usr/bin/env python3
"""
export_mmpbsa_excel.py

Recursively find every MMPBSA.py FINAL_RESULTS_*.dat file under a root
directory (e.g. one containing per-replicate, per-variant, per-interface
run folders), parse all energy components, and compile everything into a
single Excel workbook for easy comparison/filtering/pivoting.

Usage
-----
# Auto-discover every FINAL_RESULTS_*.dat under a root folder
python3 export_mmpbsa_excel.py --root . -o mmpbsa_summary.xlsx

# Auto-discover, but only under specific subfolders (repeatable)
python3 export_mmpbsa_excel.py --root replicate1 --root replicate2 -o mmpbsa_summary.xlsx

# Explicit control over labels instead of auto-discovery
python3 export_mmpbsa_excel.py \
    --run WT_rep1_iface1  wt/rep1/iface1_run/FINAL_RESULTS_iface1.dat \
    --run WT_rep1_iface2  wt/rep1/iface2_run/FINAL_RESULTS_iface2.dat \
    --run MUT_rep1_iface1 mut/rep1/iface1_run/FINAL_RESULTS_iface1.dat \
    -o mmpbsa_summary.xlsx

Output workbook
---------------
- "All Components": long-format table, one row per (run, method, section,
  component) -- every number MMPBSA.py reported, easy to pivot/filter in
  Excel.
- "Delta Total Summary": one row per run with just DELTA TOTAL (avg/std/SEM)
  for each method (GB/PB) found, for quick side-by-side comparison.

Run labels are auto-derived from the file's path (parent-of-parent /
parent directory names, e.g. ".../WT/rep1/iface1_run/FINAL_RESULTS_..." ->
label "WT/rep1/iface1_run") when using --root, or taken directly from
--run when specified explicitly. Use --run for full control if your
directory naming doesn't cleanly convey variant/replicate/interface.
"""

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

FONT_NAME = "Arial"
SECTION_ORDER = ["Complex:", "Receptor:", "Ligand:", "Differences (Complex - Receptor - Ligand):"]
SECTION_LABEL = {
    "Complex:": "Complex",
    "Receptor:": "Receptor",
    "Ligand:": "Ligand",
    "Differences (Complex - Receptor - Ligand):": "Delta (Complex-Receptor-Ligand)",
}
METHOD_HEADERS = {"GENERALIZED BORN:": "GB", "POISSON BOLTZMANN:": "PB"}


# --------------------------------------------------------------------------
# Parsing (same logic as analyze_mmpbsa.py's parser, duplicated here so this
# script is self-contained and doesn't require another file to be present)
# --------------------------------------------------------------------------

@dataclass
class EnergyTable:
    rows: dict = field(default_factory=dict)  # component -> (avg, std, sem)


@dataclass
class MMPBSAResult:
    path: str
    receptor_mask: str = ""
    ligand_mask: str = ""
    nframes: float = None
    tables: dict = field(default_factory=dict)  # method -> section -> EnergyTable


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

    m = re.search(r'Receptor mask:\s*"([^"]*)"', text)
    if m:
        result.receptor_mask = m.group(1)
    m = re.search(r'Ligand mask:\s*"([^"]*)"', text)
    if m:
        result.ligand_mask = m.group(1)
    m = re.search(r"Calculations performed using\s+([\d.]+)\s+complex frames", text)
    if m:
        result.nframes = float(m.group(1))

    current_method = None
    current_section = None

    for raw_line in text.splitlines():
        stripped = raw_line.strip()

        if stripped in METHOD_HEADERS:
            current_method = METHOD_HEADERS[stripped]
            result.tables.setdefault(current_method, {})
            current_section = None
            continue

        if stripped in SECTION_ORDER:
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


# --------------------------------------------------------------------------
# Discovery + label derivation
# --------------------------------------------------------------------------

def discover_runs(roots):
    """Find every FINAL_RESULTS_*.dat under the given root(s), return
    list of (label, path) using the path relative to the root as label."""
    runs = []
    for root in roots:
        root = Path(root)
        for p in sorted(root.rglob("FINAL_RESULTS_*.dat")):
            try:
                rel = p.relative_to(root)
            except ValueError:
                rel = p
            label = str(rel.with_suffix(""))
            runs.append((label, str(p)))
    return runs


# --------------------------------------------------------------------------
# Excel writing
# --------------------------------------------------------------------------

def style_header_row(ws, row_idx, ncols):
    for c in range(1, ncols + 1):
        cell = ws.cell(row=row_idx, column=c)
        cell.font = Font(name=FONT_NAME, bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="4472C4")
        cell.alignment = Alignment(horizontal="center", vertical="center")


def autosize_columns(ws, min_width=10, max_width=45):
    widths = {}
    for row in ws.iter_rows():
        for cell in row:
            if cell.value is None:
                continue
            col = cell.column_letter
            length = len(str(cell.value))
            widths[col] = max(widths.get(col, 0), length)
    for col, width in widths.items():
        ws.column_dimensions[col].width = min(max(width + 2, min_width), max_width)


def write_all_components_sheet(wb, results):
    ws = wb.active
    ws.title = "All Components"

    headers = ["Run", "Path", "Frames", "Receptor Mask", "Ligand Mask",
               "Method", "Section", "Component", "Average", "Std Dev", "SEM"]
    ws.append(headers)
    style_header_row(ws, 1, len(headers))

    row_i = 2
    for label, result in results:
        for method, sections in result.tables.items():
            for section in SECTION_ORDER:
                table = sections.get(section)
                if table is None:
                    continue
                section_label = SECTION_LABEL.get(section, section)
                for component, (avg, std, sem) in table.rows.items():
                    ws.cell(row=row_i, column=1, value=label)
                    ws.cell(row=row_i, column=2, value=result.path)
                    ws.cell(row=row_i, column=3, value=result.nframes)
                    ws.cell(row=row_i, column=4, value=result.receptor_mask)
                    ws.cell(row=row_i, column=5, value=result.ligand_mask)
                    ws.cell(row=row_i, column=6, value=method)
                    ws.cell(row=row_i, column=7, value=section_label)
                    ws.cell(row=row_i, column=8, value=component)
                    ws.cell(row=row_i, column=9, value=round(avg, 4))
                    ws.cell(row=row_i, column=10, value=round(std, 4))
                    ws.cell(row=row_i, column=11, value=round(sem, 4))
                    for c in range(1, len(headers) + 1):
                        ws.cell(row=row_i, column=c).font = Font(name=FONT_NAME)
                    row_i += 1

    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(headers))}{row_i - 1}"
    autosize_columns(ws)
    return ws


def write_delta_summary_sheet(wb, results):
    ws = wb.create_sheet("Delta Total Summary")

    headers = ["Run", "Method", "DELTA G gas (avg)", "DELTA G gas (SEM)",
               "DELTA G solv (avg)", "DELTA G solv (SEM)",
               "DELTA TOTAL (avg)", "DELTA TOTAL (std)", "DELTA TOTAL (SEM)",
               "Frames"]
    ws.append(headers)
    style_header_row(ws, 1, len(headers))

    row_i = 2
    for label, result in results:
        for method, sections in result.tables.items():
            diff = sections.get("Differences (Complex - Receptor - Ligand):")
            if diff is None:
                continue
            gas = diff.rows.get("DELTA G gas")
            solv = diff.rows.get("DELTA G solv")
            total = diff.rows.get("DELTA TOTAL")
            if total is None:
                continue

            ws.cell(row=row_i, column=1, value=label)
            ws.cell(row=row_i, column=2, value=method)
            ws.cell(row=row_i, column=3, value=round(gas[0], 4) if gas else None)
            ws.cell(row=row_i, column=4, value=round(gas[2], 4) if gas else None)
            ws.cell(row=row_i, column=5, value=round(solv[0], 4) if solv else None)
            ws.cell(row=row_i, column=6, value=round(solv[2], 4) if solv else None)
            ws.cell(row=row_i, column=7, value=round(total[0], 4))
            ws.cell(row=row_i, column=8, value=round(total[1], 4))
            ws.cell(row=row_i, column=9, value=round(total[2], 4))
            ws.cell(row=row_i, column=10, value=result.nframes)
            for c in range(1, len(headers) + 1):
                ws.cell(row=row_i, column=c).font = Font(name=FONT_NAME)
            row_i += 1

    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(headers))}{row_i - 1}"
    autosize_columns(ws)
    return ws


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", action="append", default=[],
                         help="Root directory to search recursively for "
                              "FINAL_RESULTS_*.dat files. Repeatable.")
    parser.add_argument("--run", nargs=2, action="append", default=[],
                         metavar=("LABEL", "PATH"),
                         help="Explicit label + path to a FINAL_RESULTS_*.dat "
                              "file. Repeatable. Combine with --root if you "
                              "want both auto-discovered and explicit runs.")
    parser.add_argument("-o", "--out", default="mmpbsa_summary.xlsx")
    args = parser.parse_args()

    labeled_paths = []
    if args.root:
        labeled_paths.extend(discover_runs(args.root))
    if args.run:
        labeled_paths.extend([(label, path) for label, path in args.run])

    if not labeled_paths:
        print("No FINAL_RESULTS_*.dat files found. Use --root DIR to "
              "auto-discover, or --run LABEL PATH to specify explicitly.",
              file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(labeled_paths)} run(s):")
    results = []
    for label, path in labeled_paths:
        if not Path(path).is_file():
            print(f"  SKIP (not found): {label} -> {path}", file=sys.stderr)
            continue
        print(f"  {label:40s} -> {path}")
        results.append((label, parse_final_results(path)))

    if not results:
        print("No files could be parsed.", file=sys.stderr)
        sys.exit(1)

    wb = Workbook()
    write_all_components_sheet(wb, results)
    write_delta_summary_sheet(wb, results)

    wb.save(args.out)
    print(f"\nSaved: {args.out}")
    print(f"  'All Components'      -- {sum(1 for _ in wb['All Components'].iter_rows(min_row=2)) } data rows")
    print(f"  'Delta Total Summary' -- {sum(1 for _ in wb['Delta Total Summary'].iter_rows(min_row=2))} data rows")


if __name__ == "__main__":
    main()
