# MM/GBSA Pipeline for FGF–FGFR–Heparin Assemblies — Session Notes

Distilled from an extended debugging/build session. System: CHARMM-GUI-built
FGF–FGFR–heparin complexes (2:2:2 and 2:2:1 stoichiometries), MD run in NAMD
2.14, binding free energy computed with AmberTools23 MMPBSA.py on macOS
(conda env `AmberTools23`).

Companion scripts referenced below (already delivered earlier in this
conversation — recreate/fetch them in the new session):
- `run_mmpbsa_pipeline.sh` — full 2:2:2 pipeline
- `run_mmpbsa_pipeline_221.sh` — full 2:2:1 pipeline (shared single heparin)
- `analyze_mmpbsa.py` — parses FINAL_RESULTS/FINAL_DECOMP, plots summary/compare/ddg/decomp
- `compare_species_energies.py` — grouped bar chart of energy components across species/variants
- `compare_species_replicates.py` — cross-replicate mean/SEM/range comparison across species
- `export_mmpbsa_excel.py` — compiles all FINAL_RESULTS_*.dat under a root into one Excel workbook

---

## 1. Environment

```bash
conda activate AmberTools23
which cpptraj   # verify env is active before running the pipeline
```

Tool paths were under `~/anaconda3/envs/AmberTools23/bin/` (cpptraj, sander,
MMPBSA.py, ante-MMPBSA.py, parmed). macOS ships BSD grep/sed — **GNU-only
flags like `grep -P` fail** ("invalid option -- P"); use `awk` instead
wherever the pipeline needs to extract numbers from text output.

---

## 2. Inputs

- `step3_input.psf` / `step3_input.pdb` — CHARMM-GUI-built system (confirm
  these are the actual files referenced by `structure`/`coordinates` in the
  NAMD production config — not necessarily the highest step-number files)
- `step5_production.dcd` (+ any additional segments) and its NAMD log
  (needed for the atom-count sanity check)
- `toppar/` — full CHARMM-GUI toppar directory (large; only a subset is
  actually needed for chamber conversion, see below)

System composition: protein (FGF + FGFR chains) + carbohydrate (heparin,
containing iduronic acid, `AIDOA` residue) + CHARMM36m/CHARMM36
carb/CGenFF force fields. No lipid/membrane component.

---

## 3. Chamber conversion (PSF/toppar → Amber prmtop)

CHARMM-GUI's `toppar.str` streams the ENTIRE toppar library regardless of
what's actually in the system — don't feed chamber the whole thing (causes
duplicate-parameter conflicts). Use an explicit, minimal file list matching
the system's actual composition:

```bash
parmed
```
```
chamber -top toppar/top_all36_prot.rtf \
        -top toppar/top_all36_carb.rtf \
        -top toppar/top_all36_cgenff.rtf \
        -param toppar/par_all36m_prot.prm \
        -param toppar/par_all36_carb.prm \
        -param toppar/par_all36_cgenff.prm \
        -str toppar/toppar_all36_moreions.str \
        -str toppar/toppar_water_ions.str \
        -str toppar/toppar_ions_won.str \
        -str toppar/toppar_all36_carb_glycopeptide.str \
        -str toppar/toppar_all36_carb_imlab.str \
        -psf step3_input.psf \
        -crd step3_input.pdb

outparm complex_chamber.parm7
```
(omit `outparm`'s second/`.rst7` argument unless coordinates were also
loaded in the same ParmEd session — otherwise it errors with
`Got 1 coordinates for N atoms`)

**Key fix discovered:** iduronic acid (`AIDOA`, a GAG/heparin sugar) uses
ring carbon types `CC3161`/`CC3163` normally, but a **sulfated/patched
copy** of the same residue can carry a different type (`CC3261` for C4)
that has no bonded parameters in the base `par_all36_carb.prm`. The
`CC3261`–`CC3163` bond (and related angle/dihedral terms) exist only in
`toppar_all36_carb_imlab.str` — **this stream must be included whenever
IdoA/GAG-type sugars are in the system**, or chamber fails with:
```
ChamberError: Problem assigning parameters to PSF: Missing bond type for ...
```
If a similar "missing bond/angle type" error recurs for other residues,
diagnose by: (1) get the atom's real type from the PSF directly
(`grep "SEGID " file.psf | awk '$3==N && $4=="RESNAME"'`), NOT from the
`.rtf`'s default listing — a patch can retype an atom; (2) grep every
`.str`/`.rtf`/`.prm` for that atom type to find which stream actually
defines its bonded parameters.

**Atom count sanity check** (must match NAMD's own count exactly):
```bash
cpptraj complex_chamber.parm7
```
```
parminfo
```
Compare to the NAMD log's `Info: N ATOMS` line (portable extraction:
`awk '/ATOMS/ {print $2; exit}' namd.log` — not `grep -P`).

---

## 4. Solvent/ion residue naming — a recurring gotcha

CHARMM-GUI names: water = `TIP3`, ions = `CLA` (Cl⁻), `POT` (K⁺).

**Chamber conversion renames water `TIP3` → `WAT`**, but leaves `CLA`/`POT`
as-is (asymmetric — easy to miss). Verify post-conversion names with:
```bash
cpptraj complex_chamber.parm7
```
```
resinfo :1-<N> out full_resinfo.txt
```
```bash
awk '{print $2}' full_resinfo.txt | sort -u
```

**Correct strip mask for this system: `:WAT,CLA,POT`** (not `:TIP3,CLA,POT`).

**MMPBSA.py also strips solvent internally**, using its OWN hardcoded
default mask (`:WAT,Cl*,CIO,Cs+,IB,K*,Li+,MG*,Na+,Rb+,CS,RB,NA,F,CL`) —
this default does **not** include `CLA`/`POT` as literal strings, so it
silently fails to strip your ions even though your `-s` mask to
`ante-MMPBSA.py` was correct. Symptom: sander fails with `FATAL: NATOM
mismatch in coord and topology files`, and the internal trajectory
(`_MMPBSA_complex.pdb`/`.mdcrd.0`) has thousands more atoms than the
`-cp` topology. **Fix: explicitly override MMPBSA.py's internal strip
mask** in the `&general` block:
```
strip_mask=":WAT,CLA,POT",
```

---

## 5. Trajectory processing

ASCII `.mdcrd` trajectory output can corrupt box-line parsing
(`expect only 3 or 6 box coords, got 7`) — **always use NetCDF** instead:

```bash
cpptraj
```
```
parm complex_chamber.parm7
trajin step5_production.dcd
autoimage
trajout complex_wet.nc netcdf
run
```

**Do NOT pre-strip solvent from this trajectory.** MMPBSA.py's `-y` must
be the SAME atom count/composition as `-sp` (solvated topology) — it does
its own internal stripping to match `-cp`/`-rp`/`-lp`. Feeding it an
already-stripped trajectory throws `Number of atoms in NetCDF file (...)
does not match number in associated topology (...)`.

---

## 6. Component topology construction (receptor/ligand split)

`ante-MMPBSA.py` **requires ligand residues to be contiguous** in the
complex topology (`PrmtopError: Couldn't predict mask... ligand residues
must be sequential`). A non-covalent cofactor (heparin) sitting physically
"between" ligand and receptor chains must be grouped entirely with ONE
side (convention used here: **receptor**), never split across both, or
the mask will end up non-contiguous.

Residue numbering: **chamber renumbers residues sequentially across the
WHOLE system** — PSF-local residue numbers are NOT preserved. Determine
true global ranges per chain/molecule with:
```bash
cpptraj complex_chamber.parm7
```
```
resinfo ^1 out mol1.dat
resinfo ^2 out mol2.dat
... (one per molecule number)
```
then match each molecule number to a chain by inspecting residue names at
block start/end (e.g. sequence identity for protein chains, `AIDOA`/`AGLC`
for heparin).

### 2:2:2 system (each interface has its own heparin copy)
- Interface 1 = ligand(FGF copy1) vs. receptor(FGFR copy1) + cofactor(heparin1)
- Interface 2 = ligand(FGF copy2) vs. receptor(FGFR copy2) + cofactor(heparin2)
- When building each interface, strip the ENTIRE other interface (its FGF,
  FGFR, AND its heparin) from `-s`.
- Confirmed via `nativecontacts` that heparin can still bridge BOTH FGF
  copies even when clearly interfacing mainly with one FGFR copy — treat
  per-interface numbers as approximations, not fully independent.

### 2:2:1 system (single heparin shared between both interfaces)
- The one heparin chain must be **retained in BOTH** per-interface
  complexes (never stripped from either) — it is not split or exclusively
  assigned.
- This is a bigger non-additivity approximation than 2:2:2; run a
  `nativecontacts` scan (heparin vs. both ligand chains, both receptor
  chains) before trusting the per-interface split's magnitude.

### Critical: `-sp` (solvated topology) must exactly match `-cp` minus only solvent
Reusing the FULL system's solvated topology (`complex_chamber.parm7`) as
`-sp` for a PER-INTERFACE run is wrong — it has extra chains `-cp` doesn't.
Must build an interface-specific solvated topology AND a matching
interface-specific solvated trajectory:

```bash
parmed
```
```
parm complex_chamber.parm7
strip :<other-interface-chains>
outparm iface1_solvated.parm7
```
```bash
cpptraj
```
```
parm complex_chamber.parm7
trajin complex_wet.nc
strip :<other-interface-chains>
trajout iface1_wet.nc netcdf
run
```
(Only the WHOLE-assembly run can reuse `complex_chamber.parm7`/
`complex_wet.nc` directly as `-sp`/`-y`, since it strips only solvent.)

### Full ante-MMPBSA.py command template
```bash
ante-MMPBSA.py -p complex_chamber.parm7 \
  -c iface1_complex.parm7 -r iface1_receptor.parm7 -l iface1_ligand.parm7 \
  -s ':WAT,CLA,POT,<other-interface-chains>' \
  -n ':<ligand-chain-range>' \
  --radii mbondi2
```

---

## 7. MMPBSA.py execution

```
&general
  startframe=1, endframe=<N>, interval=1, verbose=2,
  strip_mask=":WAT,CLA,POT",
/
&gb
  igb=5, saltcon=0.150,
/
&decomp
  idecomp=2, dec_verbose=1,
  print_res="all"
/
```
**`print_res="within 4"` syntax is NOT supported in AmberTools23** — fails
with `SelectionError: Invalid selection! Integers expected.` Use `"all"`
or an explicit residue range instead.

```bash
mkdir -p iface1_run && cd iface1_run
MMPBSA.py -O -i ../mmpbsa.in -o FINAL_RESULTS_iface1.dat \
  -sp ../iface1_solvated.parm7 \
  -cp ../iface1_complex.parm7 \
  -rp ../iface1_receptor.parm7 \
  -lp ../iface1_ligand.parm7 \
  -y ../iface1_wet.nc
```
Repeat for `iface2_run` (its own solvated topology/trajectory) and
`whole_run` (reuses `complex_chamber.parm7`/`complex_wet.nc` as `-sp`/`-y`).

If debugging a failure: check `<run>/_MMPBSA_complex_gb.mdout.0` (sander's
actual error, e.g. `FATAL: NATOM mismatch...`) and
`<run>/_MMPBSA_normal_traj_cpptraj.out` (shows the literal strip mask
cpptraj actually applied — the fastest way to see whether `strip_mask`
took effect and whether atom counts line up with the target `-cp`).

---

## 8. Interpreting results

- **Use `DELTA TOTAL`** (Differences: Complex − Receptor − Ligand section)
  as the binding-affinity comparison metric across FGF species/variants.
  It's MM + GB solvation enthalpy — **no entropy term** (no `&nmode`/
  interaction-entropy run) — standard simplification for relative
  comparisons between close variants (entropy assumed to roughly cancel).
- Individual components (VDWAALS, EEL, EGB, ESURF) explain *why* DELTA
  TOTAL moved, not what "binding efficiency" means by themselves. EEL
  (favorable) and EGB (desolvation penalty, unfavorable) commonly
  partially cancel for polar/charged interfaces like this one.
- **Per-interface DELTA TOTAL can be negative while the whole-assembly
  DELTA TOTAL is positive.** Not a bug — known non-additivity of GB
  decomposition across different topology splits:
  1. GB Born radii are computed per-topology (not pairwise-additive) —
     stripping the other interface's chains changes the retained atoms'
     effective desolvation, making per-interface numbers look more
     favorable than their "true" contribution within the intact assembly.
  2. Whole-assembly calc captures FGF–FGF and FGFR–FGFR contacts (direct
     or heparin-bridged) that per-interface calcs are structurally blind to.
  3. (2:2:1 specific) the shared heparin's favorable contacts get
     effectively "double-counted" across both per-interface numbers but
     appear only once in the whole-assembly split.
  → Treat whole-assembly DELTA TOTAL as the more complete absolute
  estimate; treat per-interface numbers as better for *relative*
  comparisons (interface 1 vs 2, WT vs mutant at a given interface) than
  as standalone absolute binding energies.
- With only ~10 frames per replicate, treat single-replicate ΔG as
  preliminary. For WT-vs-mutant comparisons at a direct interface
  substitution site, per-residue decomposition (`idecomp=2`) at the
  mutated residue is the strongest corroborating evidence — check it
  moves in a chemically sensible direction alongside any ΔΔG shift.
  Compare replicate-to-replicate spread (5 reps recommended) via
  `compare_species_replicates.py`'s ddG/SEM heuristic — not a formal
  p-value, but a reasonable noise-floor gut check given correlated MD
  sampling.

---

## 9. Known outstanding items / good next steps in a code session

- Formalize the `nativecontacts` bridging check into the pipeline scripts
  (currently documented as manual steps / printed instructions in
  `run_mmpbsa_pipeline_221.sh`, not auto-executed).
- Consider whether to add `&nmode` or interaction-entropy calculation if
  reviewers/publication requires an entropy-corrected ΔG rather than
  enthalpy-only.
- `analyze_mmpbsa.py`'s `FINAL_DECOMP_MMPBSA.dat` parser has not yet been
  validated against a real decomposition output file (only the summary
  `FINAL_RESULTS_*.dat` parser has been tested against real data) — verify
  and patch once a real decomp file is available.
- Methods-section draft (`mmpbsa_methods.md`, delivered earlier) still has
  bracketed placeholders (NAMD/AmberTools versions, temperature, exact
  frame count/interval used) to fill in once finalized.
