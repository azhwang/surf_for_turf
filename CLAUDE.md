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
  a real `FINAL_DECOMP_MMPBSA.dat`.
- The `nativecontacts` heparin-bridging check is still manual.
- The methods draft has placeholders (versions, temperature 310 K, frame
  count/interval).
- Only 10 frames per replicate. Treat single-replicate ΔG as preliminary.
