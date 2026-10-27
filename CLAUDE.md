# FGF–FGFR Variant Binding (Surf for Turf)

## What this project is

Predicting how well FGF variants from different species, plus engineered
"stab" variants, bind FGF receptors. Heparin is included as a co-factor.
The workflow has five steps:

1. AlphaFold models of each FGF/FGFR variant.
2. Alignment onto crystal templates in PyMOL.
3. CHARMM-GUI system build using CHARMM36m and CHARMM36 carb/CGenFF.
4. NAMD MD, 5 replicates per system.
5. AmberTools23 `MMPBSA.py` (GB, igb=5), using `DELTA TOTAL` as the ranking
   metric.

Species/variants: human, human_stab, bovine_stab, chinook_salmon, goldfish,
japanese_ricefish, zebrafish.

## Stoichiometries / templates

| System | Template | Chains | Where |
|---|---|---|---|
| 1:1 FGF–FGFR | 1EVT (FGFR1), 1DJS (FGFR2) | – | `1to1FGFtoFGFR/`, `~/Desktop/Surf_for_Turf/1to1`, `all_11` |
| 2:2:2 (2 heparins) | 1FQ9 | PROA/PROB = FGFR, PROC/PROD = FGF, CAR* = heparin (2–4 fragments) | `~/Desktop/Surf_for_Turf/all_222`, `all_metrics_analysis/` |
| 2:2:1 (1 shared heparin) | 1E0O | PROA/PROB = FGFR, PROC/PROD = FGF, CARA = shared heparin | `~/Desktop/Surf_for_Turf/all_221`, `all_metrics_analysis_221/` |

PyMOL build recipes are in `Surf_for_Turf_7-8_Meeting.rtf` and
`master_notes_FGF_FGFR.rtf`.

## Where things live

Two folders share the name `Surf_for_Turf`; don't confuse them.
- `~/Desktop/LAMM/Surf_for_Turf/` (this folder): early 1:1 build work
  (psfgen/VMD scripts, first NAMD run in `first_success/`), AlphaFold models,
  docking results, results spreadsheets, and the methods draft
  (`Surf_for_Turf_Simulations_Methods.docx`).
- `~/Desktop/Surf_for_Turf/` (renamed from `~/Desktop/md_sims` on 2026-10-07;
  the git repository): the production CHARMM-GUI systems, NAMD replicate
  runs (`namdrep1`–`namdrep5`), and all MMPBSA scripts:
  - `run_mmpbsa_pipeline.sh`: the 2:2:2 pipeline. A copy is dropped into
    each replicate dir by `all_metrics_analysis*/run_md_analysis.sh`.
  - `~/Downloads/run_mmpbsa_pipeline_221.sh`: the 2:2:1 pipeline. It keeps
    the shared heparin in both interfaces. **Not currently used by anything;
    see Known issues.**
  - `analyze_mmpbsa.py`, `compare_species_replicates.py`,
    `export_mmpbsa_excel.py`, `plot_mmpbsa.py`: analysis and plotting.
- `docs/heparin_md_build_notes.md`: heparin parameterization, the
  CHARMM-GUI 2BNSF patch bug and its workarounds, and the minimization
  protocol.
- `docs/mmpbsa_pipeline_notes.md`: chamber conversion, strip masks,
  receptor/ligand splitting, and how to interpret results. **Read this
  before touching the MMPBSA scripts.**

## Environment

- AmberTools: `conda activate AmberTools23` (cpptraj, parmed,
  ante-MMPBSA.py, MMPBSA.py).
- MD runs on a cluster use NAMD 2.14 CUDA via slurm (`eq_slurm.sh`,
  `prod_slurm.sh`). Local runs use NAMD 3.0.2 (macOS ARM,
  `NAMD_3.0.2_MacOS-universal-multicore/`).
- macOS has BSD grep/sed. Don't use `grep -P`; use `awk`.

## Critical gotchas (details in docs/)

- Ions in the production systems are **POT/CLA**, not SOD, and chamber
  renames TIP3 to WAT. MMPBSA needs `strip_mask=":WAT,CLA,POT"` set
  explicitly in `&general`. The heparin notes' `autoionize -cation SOD`
  example does not match the production systems. Don't mix them.
- Include `toppar_all36_carb_imlab.str` in chamber whenever IdoA is present.
  Without it, the sulfated `CC3261` type has no bonded parameters.
- Chamber renumbers residues globally. **Residue ranges differ per species**
  because FGFR constructs and heparin segment counts differ. Derive ranges
  from each system's own PSF segids; never reuse one system's ranges for
  another:
  ```bash
  awk '/!NATOM/{f=1;next} f&&/!NBOND/{exit} f{k=$2" "$3; if(k!=last){n++; last=k; if($2!=seg){if(seg!="")printf "%s:%d-%d ",seg,start,prev; seg=$2; start=n} prev=n}} END{print ""}' step3_input.psf
  ```
- Group heparin with the receptor side. The ligand mask must be contiguous.
- Per-interface and whole-assembly ΔG are not additive under GB. Use
  per-interface values for relative comparisons.
- `DELTA TOTAL` has no entropy term.

## Known issues (found 2026-09-23, not yet fixed)

The MMPBSA residue ranges are hardcoded, and most systems were run with the
wrong ranges:

- **All 2:2:1 (1E0O) runs used the 2:2:2 script and its 1FQ9 ranges.** For
  example, in bovine_stab rep1 the true ranges are PROA 1-197, PROB 198-393,
  PROC 394-548, PROD 549-703, CARA 704-713. The generated `iface2_ligand`
  is 138 residues ending in heparin. **All 2:2:1 ΔG values
  (`mmpbsa_summary_221.xlsx`, `1E0O_221_Model_Results_Replicates.xlsx`) are
  invalid.** The standalone `run_mmpbsa_pipeline_221.sh` also has wrong
  ranges (its cofactor range is 713-720), so it can't simply be swapped in.
- **2:2:2 ranges match only goldfish exactly.** Chinook and ricefish are off
  by 1–2 residues at the FGFR/heparin ends. Human, human_stab and
  bovine_stab have FGFR chains 10+ residues longer, so the heparin and FGFR
  tails were misassigned. Zebrafish and bovine_stab have heparin split into
  3–4 CAR segments.
- v1 also wrote interface 2's `ante-MMPBSA.py -n` mask in full-system
  numbering, but `-n` is applied after `-s` strips and renumbers. So
  **iface2 was mis-split in every v1 run, goldfish included**. v1 also
  grouped heparin with FGF, because v1's "ligand" (PROA/PROB) is actually
  FGFR.
- **Segment identity (verified by sequence in all 70 replicates):
  PROA/PROB = FGFR, PROC/PROD = FGF, CAR\* = heparin.**

### Fix: v2 pipeline (2026-09-23), run alongside v1

- `~/Desktop/Surf_for_Turf/run_mmpbsa_pipeline_v2.sh` +
  `derive_mmpbsa_split.py`: derives the split per replicate from its own
  PSF/PDB. FGF is the ligand and heparin goes with FGFR. For 2:2:1 the
  shared heparin is kept in both interfaces. Component residue counts are
  verified before MMPBSA.py runs. Writes to `namdrepN/mmpbsa_v2/`; v1
  outputs are untouched.
- `run_all_mmpbsa_v2.sh <roots...>`: batch runner. It logs to
  `mmpbsa_v2_batch.log` and skips replicates that are already done.
- **1FQ9 bovine_stab and zebrafish (and possibly human_stab): each FGF's D3
  contacts both FGFRs' D3, at the start and at the end of MD.** They have no
  clean one-to-one interface, so v2 computes only the whole-assembly split
  for them (see `mmpbsa_v2/split_summary.txt`). Those assemblies probably
  deviate from the 1FQ9 arrangement; check them in VMD before comparing them
  with the other species.
- Results: `~/Desktop/Surf_for_Turf/mmpbsa_v2_results/` (`mmpbsa_v2_summary.xlsx`
  + one plot per stoichiometry × split), regenerated with
  `python3 summarize_mmpbsa_v2.py` after `python3 qc_ring_piercing_v2.py`.

### Ring piercing in the builds (found 2026-09-24)

- **The FGF N-terminal tail (roughly residues 4–26, from the AlphaFold
  models) is threaded through FGFR rings (roughly FGFR residues 250–350) in
  several builds.** MD can't undo it, and it adds +230 to +535 kcal/mol of
  fake vdW to one residue pair. Affected systems: 2:2:1 goldfish (Phe12 ↔
  Arg255), 2:2:1 zebrafish (Phe11 ↔ Ser252/Pro253), 2:2:2 human (Phe20 ↔
  Pro283), and 2:2:2 bovine_stab (Pro13 ↔ Leu257). Their ΔG is marked
  INVALID in the summary.
- **2:2:2 chinook_salmon has the same defect** (Phe12 pierced by FGFR Val316
  CB, 0.77 Å). That is why equilibration crashed with a RATTLE failure in
  all 5 replicates.
- QC tools: `check_ring_piercing.py step3_input.pdb` checks geometry. **Run
  it on every new build before MD.** Many input-structure hits are benign
  overlaps that minimization resolves, so confirm energetically with
  `qc_ring_piercing_v2.py` after MMPBSA. A healthy split has max residue
  ΔvdW ≤ ~4 kcal/mol; a pierced one is hundreds.
- **Fixed by the consistent rebuild below.** Decision (user, 2026-09-29): trim
  FGF to the residues the scaffold crystal structures resolve.

### Consistent rebuild: `~/Desktop/Surf_for_Turf/rebuild_v2/` (2026-09-29)

Rebuilds all 14 systems with one protocol. See `rebuild_v2/README.md`; the
single entry point is `run_rebuild.sh`.
- **FGF is trimmed to the template-resolved core** (1FQ9: 129 residues;
  1E0O: 127 residues, with the same span for both copies). This removes the
  14–26-residue disordered N-tail. Every FGF1/FGF2 crystal structure in the
  PDB leaves those residues out or unmodelled. The core carries all
  species/stab differences: fish cores are 75–79% identical to human, and all
  9–10 stab mutations sit in the core.
- **FGF placement:** the core is superposed onto the crystal FGF chains
  (1FQ9 A→PROC, B→PROD; 1E0O A→PROC, C→PROD). First core residues with CA
  deviation >1.5 Å are rebuilt on the crystal backbone. Conserved side chains
  within 6 Å of FGFR/heparin use crystal rotamers. A chi scan
  (`fix_rotamers.py`) fixes any remaining piercing or severe (<1.8 Å) clash.
  `precheck.py` gates the build on 0 piercings and 0 severe clashes.
- **Scaffold:** one per template for FGFR + heparin. 1FQ9 uses the goldfish
  build (heparin intact, 8+6 residues); 1E0O uses the chinook build. The
  scaffold coordinates equal the crystal frame exactly.
- **Result:** FGF–FGFR pairing is now unambiguous in all 14 systems (1FQ9
  about 70 vs 30 contacts, 1E0O about 24 vs 0). This includes 2:2:2
  bovine_stab, zebrafish and human_stab, which were whole-only before.
- Same box (124 Å / 129 Å), 0.15 M KCl, NAMD inputs, restraints and run
  lengths as the originals, plus an explicit per-replicate `seed`.
  `submit.sh` chains production on equilibration with `afterok`.
- **NAMD needs CHARMM-GUI's NAMD-format toppar.** That is `namdrepN/toppar`
  (or `namd/toppar`) in the original builds, NOT the top-level
  `charmm-gui-*/toppar`, which is the CHARMM version with `set nat ?NATC`
  lines. With the CHARMM version, every NAMD build (cluster 2.14 and local
  3.0.2) fails: "UNKNOWN PARAMETER IN CHARMM PARAMETER FILE". The rebuild
  first shipped the wrong one; fixed 2026-10-01. `sim/toppar` is now
  identical to the original runs', and the inputs read it as `../toppar/`
  (no symlinks, which the scp/SFTP transfer had broken).
- **MM/GBSA chamber needs the CHARMM-format toppar instead** (it has the
  `.rtf` files; the NAMD-format copy doesn't). Each rebuild `sim/` has
  `toppar_charmm/` for this. `run_mmpbsa_pipeline_v2.sh` picks the first of
  `$CHAMBER_TOPPAR`, `../toppar_charmm`, `../toppar`, `toppar` that contains
  `top_all36_prot.rtf`. `toppar_charmm/` does not need to go to the cluster.
- First rebuilt result (1FQ9 goldfish rep1, 2026-10-01): ΔG iface1 −171.8,
  iface2 −150.3, whole −382.2 kcal/mol. All are within the original goldfish
  replicate spread (that build wasn't pierced, so similar values are
  expected). Max residue ΔvdW is 1.8, so no artifacts.
- After MD: `run_all_mmpbsa_v2.sh rebuild_v2/1fq9 rebuild_v2/1e0o`, then
  `qc_ring_piercing_v2.py --set rebuild`, then
  `summarize_mmpbsa_v2.py --set rebuild` (writes `mmpbsa_rebuild_results/`,
  kept separate from the original-build results).
- `export_mmpbsa_excel.py --root` searches every subfolder, so run from a
  species root it picks up both v1 and v2. Point it at the `mmpbsa_v2`
  folders, or filter on that path.

### Rebuilt results (2026-10-02): `~/Desktop/Surf_for_Turf/mmpbsa_rebuild_results/`

- 61/70 replicates analysed (2:2:2 bovine_stab r5 was added 2026-10-04 and is
  now n = 5). Still missing: 2:2:1 zebrafish (all 5, not uploaded), 2:2:2
  human r5 / human_stab r1 / zebrafish r2 (not uploaded), and 2:2:1
  bovine_stab r4 (never ran). The summary and QC only read `namdrep[0-9]`,
  so leftovers like `namdrep5_old` are ignored. The
  cluster upload also dropped some `step3_input.psf/pdb`; those were restored
  from `rebuild_v2_old/*/build/ionized.*` after checksum checks.
- Clean: max single-residue ΔvdW is 2.4 kcal/mol across all 180 splits.
  Control: human_stab vs bovine_stab (one residue apart) differ by ≤2 kcal/mol.
- **No significant species differences anywhere** (2026-10-09, using the
  manuscript's statistics: one-way ANOVA + Tukey HSD on replicate means).
  2:2:2 whole p = 0.92; 2:2:1 whole p = 0.21; every split p > 0.18; no Tukey
  pair is significant. The earlier "ricefish binds more weakly" claim
  (uncorrected Welch p = 0.008) does NOT survive multiple-comparison
  correction (Tukey p = 0.22), so don't use it. The original builds' ricefish
  "strongest" (−269) was an artefact.
- Interface conservation (1FQ9/FGFR1 model, 4.5 Å contacts across all 7
  orthologs): 28/40 FGFR-contacting FGF positions identical (70%), 11/12
  heparin-contacting positions identical (92%); 65% across the whole core.
- Receptor species: modelled human FGFR2 (1E0O) is 97% identical to bovine.
  Human FGFR1 (1FQ9) is identical to bovine FGFR1 up to residue 312; the
  differences (313–352) sit in the IIIb/IIIc alternatively spliced exon (the
  bovine UniProt entry is likely the other isoform).
- Experimental claim (user): in 2:2:1, goldfish binds more strongly than
  human, and the effect disappears in 2:2:2. The simulations do NOT resolve
  this. Goldfish − human: 2:2:1 whole −4 (p = 0.8), 2:2:2 whole −5 (p = 0.7).
  Replicate SD is about 25 kcal/mol, against an experimental ΔΔG of about
  1 kcal/mol.
- Human 2:2:1 is bimodal at iface1: reps 1 and 4 at about −135, reps 2/3/5 at
  about −90. The ~44 kcal/mol gap is electrostatic: FGF's heparin-binding
  site (core L101/102/107/111/117 = FGF2 K128/R129/K134/K138/K144) engages
  the shared heparin much more in reps 1 and 4 (about 16 vs 10–12 Lys/Arg N to
  heparin O contacts within 4 Å; 233 vs ~175 heavy-atom contacts). The
  contacts still grow over the 1 ns, so they are not converged. Goldfish
  matches the low human state (7–10 contacts; same K/R conserved; R42→V).
  Scratch analysis: cpptraj nativecontacts on PROC residues
  411,445,504,505,510,514,520 vs CARA 394-403 sulfate/carboxylate O.
- **Extensions (set up 2026-10-05, `rebuild_v2/make_extensions.py`):** 50 ns
  for 2:2:1 human r1/r2/r4 and goldfish r1-3, and 2:2:2 human r1-3 and
  goldfish r2,r3,r5 (r1 lacked restart files on the cluster). They run in 1 ns self-resubmitting segments (`ext_slurm.sh`,
  `ext01..ext50`). The template uses `firsttimestep` plus a large `numsteps`,
  because `numsteps` is the ABSOLUTE last step in NAMD. Tested locally (2
  short segments chain correctly). Analysis after they finish still needs
  setting up: MM/GBSA over the concatenated `ext*.dcd`, plus heparin-contact
  time series.
- MM/GBSA runs locally with `JOBS=4` (8 swapped on this Mac):
  `JOBS=4 bash run_all_mmpbsa_v2.sh rebuild_v2/1fq9 rebuild_v2/1e0o`.

## Other open items

- `analyze_mmpbsa.py`'s decomposition parser has not been validated against
  a real `FINAL_DECOMP_MMPBSA.dat`. **This is the v1 parser only** — it is not
  on the path to any current result. The parser the extension and chimera
  analyses use, `decomp()` in `explo_summary.py`, WAS validated 2026-10-20
  (see below), so this open item does not gate the chimera analysis.
- The `nativecontacts` heparin-bridging check is still manual.
- The methods draft has placeholders (versions, temperature 310 K, frame
  count/interval).
- Only 10 frames per replicate. Treat single-replicate ΔG as preliminary.

### Extension results + exploratory analyses (2026-10-15)
- Pre-registered (ANALYSIS_PLAN_extensions.md): 2:2:1 whole ΔG goldfish − human −27
  (CI −96, +41; p = 0.28; n 2 vs 3) → **v1 claim NOT SUPPORTED (not resolved)**. 2:2:2 −19
  (p = 0.45). The 1 ns human two-state behaviour does not persist: all reps flicker
  (tight 1–26 % of frames), so it was a short-run artefact. Results in `rebuild_v2/ext_results/`.
- Exploratory (`analyze_exploratory.sh`, `explo_*.py`; label as exploratory in the paper; ~10
  tests, none survives correction): the species × stoichiometry interaction is null (whole p = 0.8).
  FGF–FGFR without heparin favours goldfish in BOTH stoichiometries (−15, p = 0.02; −18, p = 0.14),
  so it is not stoichiometry-specific. FGF–heparin interaction is −18 (right direction, p = 0.5).
  Net core charge is identical (+11), and static electrostatics are the same.
- Mechanistic lead: conserved **R90** (FGF2 155-aa numbering; known heparin-site residue)
  engages the shared 2:2:1 heparin in goldfish (both reps, both copies) but barely in human;
  per-residue −11.5 vs −0.5 to −2 kcal/mol. A falsifiable prediction for goldfish rep 3 was
  committed before its analysis: `ext_results/PREDICTION_goldfish_rep3.md`.
- **Goldfish rep 3 added (2026-10-16); final n = 3 vs 3.** Pre-registered 2:2:1 whole: goldfish − human
  −20 (CI −87, +46; p = 0.40) → v1 claim NOT SUPPORTED (not resolved). The species × stoichiometry interaction
  is null (whole p = 0.98). FGF–FGFR without heparin favours goldfish in both stoichiometries (2:2:1 −17,
  p = 0.01; 2:2:2 −18, p = 0.14).
  **The R90 prediction was confirmed in rep 3** (heparin contacts 18.9/15.8; −10.6/−9.4 kcal/mol). R90 engages the
  heparin in 3/3 goldfish and 0/3 human in 2:2:1 and contributes the same in both species in 2:2:2, so it is a
  stoichiometry-specific difference in heparin-binding mode. The net FGF–heparin energy is compensated (−6, p = 0.77).

### Segment-swap chimeras (set up 2026-10-16; not yet run)
- `rebuild_v2/make_chimeras.py` builds `chim_human_gfseg` (human core + goldfish FGF2 86–97:
  K86N E87S L92I A93G S94T K95R C96R V97T) and `chim_goldfish_hseg` (the reverse), in 1E0O and 1FQ9,
  3 replicates each. Host backbone everywhere, donor side chains superposed; precheck clean, builds
  match the parents (identical close contacts, identical local NAMD start energies).
- `build_all.py`: `chim_*` systems get 3 replicates and seeds +500. `precheck.py` now skips non-system
  entries (.DS_Store crashed it).
- Each replicate has `submit_chain.sh` (eq → prod → 50 × 1 ns ext chain); top level
  `submit_chimeras.sh [1e0o|1fq9]`. Upload to the cluster `extension_test/` folder (where the extensions ran).
- Pre-registered plan: `rebuild_v2/ANALYSIS_PLAN_chimeras.md`. Primary metric: Arg90–heparin contacts
  (`analyze_r90.sh`, validated against the parents; reference in `ext_results/r90_parent_reference.txt`).
  Engaged copy = contacts > 7.5 and Arg90 decomposition < −5 kcal/mol. H1: human_gfseg engaged in ≥ 2/3;
  H2: goldfish_hseg engaged in ≤ 1/3 (2:2:1).

#### Run provenance (2026-10-20)
- Duplicate `ext_slurm.sh` jobs were queued on three chains (1e0o/chim_goldfish_hseg r1 and r3,
  1fq9/chim_human_gfseg r1) — leftovers of the 2026-10-16 accidental resubmission and/or the
  script's `sbatch` retry loop. Cancelled 2026-10-20 before they ran. Two jobs on one chain both
  pick the first missing `extNN.coor` at start time, so they race on the same segment files;
  each chain must hold exactly ONE queued job. Check with
  `squeue -u azhwang -h -o "%.10i %.40j %Z" | grep ext_ | sort -k3` — expect 12 lines, 12 paths.
- Earlier concurrent duplicates did overlap on nine segments (doubled "End of program" in
  ext06/08/09/10/11 across five replicates). Trajectory integrity test: within a replicate every
  `ext*.dcd` must be the same size (10 frames at `dcdfreq 50000`). Only one file failed it —
  **1fq9/chim_goldfish_hseg r1 `ext10.dcd` has 8 frames of 10** (17,472,916 vs 21,841,076 bytes;
  182,000 atoms -> 276-byte header + 2,184,080/frame). Segment 10 is outside the pre-specified
  11–50 window, the DCD is never read by the dynamics, and `ext10.coor` is complete and
  continuous (ext09 ends 5,000,000; ext10 runs 5,000,000->5,500,000; ext11 starts 5,500,000), so
  the analysis is unaffected and the file was left as is. Mention in the Methods.
- A truncated log does NOT imply bad data: the parent `1e0o/goldfish` r3 `ext33.log` has no
  "End of program" because the transfer to the Mac cut it, while its `.dcd` matches all 49 others
  byte for byte. Judge segments by the `.dcd`/`.coor`, not the log.
- Measured throughput before the 2026-10-20 maintenance window: 6.3-9.4 ns/day per chain
  (mean 8.0), i.e. roughly one 1 ns segment per 3 h with almost no queue wait. `sched_mit_buehler_gpu`
  has only 3 GPU nodes (node1228-1230), and each chain holds one job at a time, so it has to win a
  slot 40 more times; competing submissions cost it multiplicatively, not additively.

#### Decomposition parser validated (2026-10-20)
`decomp()` in `explo_summary.py` reads `float(p[17])` from
`ext_analysis/decomp_run/FINAL_DECOMP_decomp.dat`. Column 17 IS `TOTAL / Avg.` Header line 7 is
`Residue,Location,Internal,,,van der Waals,,,Electrostatic,,,Polar Solvation,,,Non-Polar Solv.,,,TOTAL,,`
and each named block is Avg./Std.Dev./SEM, so 2 = Internal, 5 = vdW, 8 = Electrostatic,
11 = Polar Solvation, 14 = Non-Polar, **17 = TOTAL**, 18 = SD, 19 = SEM. Checked independently:
Internal + vdW + Elec + PolarSolv + NonPolar equals column 17 on 254/254 ligand residues of
1e0o/goldfish r3, zero mismatches — an identity a shifted index could not reproduce. Re-run the
check on any new file with:
```bash
python3 -c 'import sys;from pathlib import Path;b=0;n=0
for l in Path(sys.argv[1]).read_text().split("\n"):
 p=l.split(",")
 if len(p)>17 and p[1].startswith("L "):
  n+=1; b+= abs(sum(float(p[k]) for k in (2,5,8,11,14))-float(p[17]))>1e-6
print(f"{n-b}/{n} rows consistent")' <FINAL_DECOMP_decomp.dat>
```
- **Arg90 is FGF2 155-aa numbering = core position 63 (1E0O) / 66 (1FQ9), NOT core 90** (core 90 is
  a Ser; the 1E0O core starts at model residue 19, offset +27 — the same offset as
  core L101/102/107/111/117 = FGF2 K128/R129/K134/K138/K144). `analyze_r90.sh` line 6 records it.
  Looking up core 90 by mistake returns +0.04 and looks like a null result.
- The published goldfish r3 values reproduce exactly from column 17: copy C `L 63 = ARG` **-10.60**,
  copy D `L 190 = ARG` **-9.40** kcal/mol. The ranking is chemically coherent — the most negative
  ligand residues are core 102 (-17.6), 101 (-13.8), 117 (-13.7) and 111 copy D (-12.6), i.e. the
  canonical heparin-binding set, with Arg90 fifth. So the pre-registered criterion's second half
  (decomposition < -5 kcal/mol) rests on a verified column.

### 1:1 FGF2–FGFR systems (rebuilt from 1EVT/1DJS 2026-10-19; not yet run)

`rebuild_v2/run_1to1.sh` -> `rebuild_v2/onetoone/{fgfr1,fgfr2}/<species>/`. 7 species x 2
receptors x 5 replicates = 70 runs. No heparin (user decision 2026-10-19).

**Built on the dedicated 1:1 crystals** (user decision 2026-10-19), not on the 2:2:2 / 2:2:1
solutes. A 2:2 receptor conformation is the wrong starting point for a 1:1 run: 1DJS and 1E0O
FGFR2 differ by a large D2–D3 hinge rotation (whole-chain CA RMSD 12.15 A; D2 1.04, D3 1.74)
and engage FGF very differently (575 vs 244 heavy-atom contacts within 4.5 A). 1EVT and the
1FQ9 FGFR1 are equivalent (0.87 A over 189 CA), so for FGFR1 the two routes agree; using both
crystals keeps one protocol for both receptors. The first version of `build_1to1.py`, which
deleted PROB/PROD/CAR* from the assembly builds, is superseded.

Three steps, all run by `run_1to1.sh`:

1. `prepare_1to1_scaffold.py` -> `onetoone/scaffold/receptor_{fgfr1,fgfr2}.pdb`.
   Both crystal receptors have a gap in the same D3 loop: 1EVT chain C 294–303 and 312–323
   (chain D is no better), 1DJS chain A 297–306. **1E0O cannot donate the 1DJS loop — its own
   gap, 294–308, spans it.** 1FQ9 chain C resolves all of it (FGFR1 149–359, no gaps) and is
   the donor for all three: its numbering matches 1EVT exactly and runs one residue behind
   1DJS (FGFR1 i <-> FGFR2 i+1), and across the FGFR2 gap the loop is 8/10 identical
   (FGFR1 296–305 NGSKIGPDNL vs FGFR2 297–306 NGSKYGPDGL; only I300->Y301 and N304->G305
   differ, and those two keep only backbone + CB from the donor). One rigid placement cannot
   close both junctions, so each loop is fitted twice, once on each 6-residue flank, and split
   in the middle; both junctions then sit on the host backbone and one mild mid-loop break is
   left (1.96 A FGFR1, 2.89 A FGFR2) for the minimisation in step 3 to close.
   Also: MSE -> MET for 1DJS, and 16 side chains that 1EVT truncates at disordered surface
   Lys/Arg/Glu completed from the donor (host backbone, donor side chain on the host N/CA/C).
   Four side chains (R148 K195 E294 R329 of FGFR1, Y301 of FGFR2) are left to `guesscoord`.
2. `prepare_1to1_fgf.py` -> `onetoone/<rec>/<sp>/fgf_PROC.pdb`. Same protocol as
   `prepare_fgf_core.py`, with the 1:1 crystals' FGF chains as the placement targets (1EVT
   chain A, the copy that packs against receptor chain C: 307 contacts vs 9 for chain D; and
   1DJS chain B) and the repaired receptor alone as the scaffold for the interface side-chain
   graft. **Both targets are restricted to FGF1 residues 10–138, the span 1E0O resolves**, so
   all 14 systems get the same 127-residue FGF2 core as the 2:2:1 set and the numbers are
   directly comparable with Table 1. Core CA RMSD to template 0.79–0.95 A and no N-terminal
   anchoring anywhere. Using the full 1DJS span (6–140) instead gave 131–132-residue cores,
   1.1–2.3 A RMSD and four anchored residues per system, because FGF1 and the FGF2 models
   diverge at the start of the N-terminal segment. The FGF2 model is taken from the 1FQ9
   CHARMM-GUI input for both receptor sets, so the two sets differ only in the receptor.
3. `build_1to1.py`. psfgen builds PROA (FGFR) + PROC (FGF2) with both D3 disulfides patched
   (FGFR1 178–230, 277–341; FGFR2 179–231, 278–342) and `guesscoord` for the rest. Then NAMD
   3.0.2 relaxes, in implicit solvent (GBIS, 0.15 M, 2000 steps), the grafted loops plus
   everything within 5 A of them plus any residue in an FGFR–FGF heavy-atom contact < 2.8 A
   (54–86 residues per system); **everything else is held fixed at its crystal position.** That
   closes the mid-loop break and the few tight interface contacts the rigid placement leaves
   (worst 1.11 A in fgfr2/goldfish -> 2.75 A). Finally re-solvated at 12 A padding and 0.15 M
   KCl.

- **Box sizing changed:** the solute is first rotated onto its principal axes
  (`orient()`), then the cube is sized from the longest principal extent + 2 x PAD. In the
  1EVT/1DJS frames the long axis runs diagonally, so an axis-aligned measurement understates
  it. The complex is ~101 A end to end, so the first version's 103 A boxes were undersized and
  would have let it approach its own periodic image once it rotated. Now 125 A / ~187.5k atoms
  (FGFR1) and 126 A / ~192.1k atoms (FGFR2) — the same cost per ns as the 2:2:2 runs
  (124 A / 204.8k), not cheaper.
- **Gotcha:** `orient()` aligns the solute to its principal axes before solvating. The first
  version used `np.linalg.eigh`'s eigenvector matrix as given; it is orthonormal but can be
  left-handed, and for all seven FGFR1 systems `det(V) = -1`, so the solute was **reflected
  into D-amino acids**. NAMD minimised and ran those systems with entirely normal energies and
  flagged nothing. `orient()` now negates a column when the determinant is negative, and
  `chirality_violations()` asserts zero D-residues on both the oriented solute and the final
  ionised system. Check the sign of (N-CA) x (C-CA) . (CB-CA) per residue after any rigid-body
  transform built from an eigendecomposition.
- QC on all 14: no peptide gaps, 0 ring piercings, closest FGFR–FGF contact 2.65–2.75 A,
  closest contact to a grafted loop 2.68–2.71 A. Interface contacts within 4.5 A are 388
  (fgfr1/goldfish) and 352 (fgfr2/goldfish) against 307 and 443 in the crystals. All 140
  `seed` values distinct (30000 + 10000*receptor + 100*species + 10*rep, +1 for production).
  Local NAMD 3.0.2 smoke test on fgfr2/goldfish minimises and runs cleanly.
- Upload `onetoone/` to the cluster and `bash submit_1to1.sh [fgfr1|fgfr2]`. 4.9 GB.

### A pre-existing defect in the 2:2:1 FGFR2 (found 2026-10-19)

CHARMM-GUI left the FGFR2 294–308 loop unmodelled in the 1E0O builds, and psfgen bonded
residue 293 to 309 across the gap. In `1e0o/*/build/solute.pdb` that peptide bond is
**7.75 A**; minimisation closes it to 1.37 A, so the two gap-flanking residues end up
covalently fused. It stays closed through 50 ns.

The distortion is local: superposing PROA D3 (250–360) before and after minimisation gives
0.40 A CA RMSD (D2 0.29, whole chain 0.35), and only residues 309 (1.7 A) and 293 (1.3 A) move
more than 1 A. Both sit 19–20 A from FGF. So the 2:2:1 ΔG values are not invalidated, but the
receptor carries a non-physical link and this is part of why 1E0O and 1DJS FGFR2 differ so
much. Worth a sentence in the Methods. The 1FQ9 builds are unaffected (FGFR1 chain C has no
gaps).
