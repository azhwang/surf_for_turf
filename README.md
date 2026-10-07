# Surf for Turf: FGF–FGFR–heparin binding across species

MD + MM/GBSA pipeline comparing how FGF variants from different species (human,
human_stab, bovine_stab, chinook_salmon, goldfish, japanese_ricefish, zebrafish) bind
FGFR with heparin, in two stoichiometries: **2:2:2** (template 1FQ9, two heparins)
and **2:2:1** (template 1E0O, one shared heparin).

This repository holds the code, documentation, input files and summary results.
The simulation data itself (~134 GB of trajectories, restart files, solvated systems and
MM/GBSA intermediates) lives alongside it in this folder but is excluded by
`.gitignore`, which ignores everything by default and only lets in code, docs and
small results.

Project notes and history: **[CLAUDE.md](CLAUDE.md)**. Method notes: [docs/](docs/).

## Layout

| Path | What |
|---|---|
| `rebuild_v2/` | **Current systems.** Consistent rebuild of all 14 systems (FGF trimmed to the crystal-resolved core, placed on the crystal FGF, one FGFR+heparin scaffold per template). Build scripts, MD inputs, and the 50 ns extension setup. See [rebuild_v2/README.md](rebuild_v2/README.md). |
| `run_mmpbsa_pipeline_v2.sh`, `derive_mmpbsa_split.py` | MM/GBSA for one replicate. The receptor/ligand split is derived from each replicate's own PSF/PDB (FGF = ligand, heparin with FGFR). |
| `run_all_mmpbsa_v2.sh` | Batch runner (`JOBS=4` for parallel replicates; skips finished ones). |
| `qc_ring_piercing_v2.py`, `check_ring_piercing.py` | Build-artifact QC: energetic (per-residue ΔvdW) and geometric (pierced rings). |
| `summarize_mmpbsa_v2.py` | Spreadsheet + plots (`--set rebuild` → `mmpbsa_rebuild_results/`). |
| `mmpbsa_rebuild_results/` | **Current results** (rebuilt systems). |
| `mmpbsa_v2_results/` | Original CHARMM-GUI builds re-analysed with the v2 split. Several are marked INVALID because of ring piercing. |
| `all_*`, `all_metrics_analysis*`, `1to1/` | Run scripts from the original builds (provenance). |
| `run_mmpbsa_pipeline.sh` | v1 MM/GBSA pipeline. **Superseded**: hardcoded residue ranges were wrong for most systems. |

## Typical workflow

```bash
# MM/GBSA on all rebuilt replicates (skips finished ones), then QC and summary
JOBS=4 bash run_all_mmpbsa_v2.sh rebuild_v2/1fq9 rebuild_v2/1e0o
python3 qc_ring_piercing_v2.py --set rebuild
python3 summarize_mmpbsa_v2.py --set rebuild
```

Requires AmberTools23 (`~/anaconda3/envs/AmberTools23`). Rebuilding the systems also
needs VMD 2.0 (psfgen/solvate/autoionize) and Biopython. MD runs on the cluster with
NAMD 2.14 CUDA.
