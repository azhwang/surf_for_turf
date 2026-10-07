#!/usr/bin/env bash
#
# run_mmpbsa_pipeline.sh
#
# Full MM/GBSA pipeline: CHARMM-GUI PSF/PDB + NAMD DCD  ->  MMPBSA.py results
#
# Covers, in order:
#   1. chamber conversion (PSF/toppar -> Amber-format prmtop)
#   2. atom-count sanity check against the NAMD log
#   3. PBC-wrapped, solvated NetCDF trajectory from the DCD
#   4. per-interface (and whole-assembly) component topologies + trajectories
#   5. MMPBSA.py runs (GB, with per-residue decomposition) for each split
#
# This encodes the *specific* fixes discovered while debugging this pipeline:
#   - chamber needs an explicit, minimal list of toppar files (not the whole
#     toppar.str, which streams the entire CHARMM36 library and causes
#     duplicate-parameter conflicts)
#   - CHARMM-GUI's TIP3/CLA/POT naming: TIP3 -> WAT after chamber conversion,
#     but CLA/POT are left as-is -- MMPBSA.py's *own* internal default
#     strip_mask does NOT know CLA/POT, so it must be overridden explicitly
#   - MMPBSA.py's -sp (solvated) topology/trajectory must describe the EXACT
#     same atom set as -cp/-rp/-lp minus only solvent -- for a multi-chain
#     assembly split into interfaces, that means building an
#     interface-specific *solvated* topology and trajectory too, not just
#     reusing the full assembly's solvated files
#   - ante-MMPBSA.py requires ligand residues to be contiguous in the
#     complex topology -- non-covalent cofactors (e.g. heparin) that
#     "bridge" both sides must be grouped entirely with one side (here:
#     grouped with the receptor)
#
# Adjust the CONFIG section below for your system, then run:
#   bash run_mmpbsa_pipeline.sh
#
# Requires: AmberTools (parmed, ante-MMPBSA.py, cpptraj, MMPBSA.py, sander)
# on your PATH.

set -euo pipefail

# ============================================================================
# CONFIG -- edit these for your system
# ============================================================================

# --- Input files ---
PSF="step3_input.psf"
PDB="step3_input.pdb"
DCD="step5_production.dcd"          # add more DCDs below in the trajin list if you have several segments
NAMD_LOG="step5_production.log"     # used only for the atom-count sanity check; set to "" to skip
TOPPAR_DIR="../toppar"

# --- chamber conversion: explicit file list (edit for YOUR system's composition) ---
# This example is for a protein + carbohydrate(glycan/heparin) + CGenFF-ligand system.
# Trim/extend based on what your toppar.str actually streams AND what residues
# your PSF actually contains (see the "grep residue name in *.str/*.rtf" steps
# from the interactive debugging session if you hit a missing-parameter error).
CHAMBER_TOP_FILES=(
    "$TOPPAR_DIR/top_all36_prot.rtf"
    "$TOPPAR_DIR/top_all36_carb.rtf"
    "$TOPPAR_DIR/top_all36_cgenff.rtf"
)
CHAMBER_PARAM_FILES=(
    "$TOPPAR_DIR/par_all36m_prot.prm"
    "$TOPPAR_DIR/par_all36_carb.prm"
    "$TOPPAR_DIR/par_all36_cgenff.prm"
)
CHAMBER_STR_FILES=(
    "$TOPPAR_DIR/toppar_all36_moreions.str"
    "$TOPPAR_DIR/toppar_water_ions.str"
    "$TOPPAR_DIR/toppar_ions_won.str"
    "$TOPPAR_DIR/toppar_all36_carb_glycopeptide.str"
    "$TOPPAR_DIR/toppar_all36_carb_imlab.str"       # needed for e.g. iduronic acid / GAG sugars
    # "$TOPPAR_DIR/cam.str"                          # uncomment + set if you have a CGenFF ligand/custom residue
)

# --- solvent/ion residue names AFTER chamber conversion ---
# NOTE: chamber typically renames CHARMM's TIP3 water residue to WAT.
# Ion names (CLA, POT, SOD, etc.) are usually left as CHARMM-style names.
# Verify with: cpptraj complex_chamber.parm7 -> resinfo, or awk '{print $4}' step3_input.psf | sort -u
STRIP_MASK=":WAT,CLA,POT"

# --- system decomposition: global residue ranges in the CHAMBER topology ---
# Find these with cpptraj's `resinfo ^N` per molecule number (see debugging
# session) -- chamber renumbers residues sequentially across the WHOLE
# system, so these are NOT the same as PSF-local residue numbers.
# Example below: 2:2:2 FGF-FGFR-heparin assembly.
LIGAND1_RANGE="1-211"          # e.g. PROA (FGF copy 1)
RECEPTOR1_RANGE="408-552"      # e.g. PROC (FGFR copy 1)
COFACTOR1_RANGE="698-705"      # e.g. CARA (heparin bound at interface 1) -- grouped with receptor

LIGAND2_RANGE="212-407"        # e.g. PROB (FGF copy 2)
RECEPTOR2_RANGE="553-697"      # e.g. PROD (FGFR copy 2)
COFACTOR2_RANGE="706-711"      # e.g. CARB (heparin bound at interface 2) -- grouped with receptor

# whole-assembly split: both ligands vs. both receptors+cofactors
LIGAND_WHOLE_RANGE="1-407"
# (receptor side = everything else non-solvent; ante-MMPBSA.py infers this from -n and -s)

# --- MMPBSA settings ---
START_FRAME=1
END_FRAME=10          # set to your actual frame count; see README notes on convergence
INTERVAL=1
SALTCON=0.150         # match your system's ionic strength from CHARMM-GUI solvation step
GB_MODEL=5            # igb=5 (or 8) are common choices
DECOMP_PRINT_RES="all"  # AmberTools23 does not support "within N" syntax here

# --- radii for GB/PB ---
RADII="mbondi2"

# ============================================================================
# END CONFIG -- shouldn't need to edit below this line
# ============================================================================

echo "=============================================="
echo " Step 1: chamber conversion (PSF/toppar -> Amber prmtop)"
echo "=============================================="

CHAMBER_ARGS=()
for f in "${CHAMBER_TOP_FILES[@]}";   do CHAMBER_ARGS+=(-top "$f");   done
for f in "${CHAMBER_PARAM_FILES[@]}"; do CHAMBER_ARGS+=(-param "$f"); done
for f in "${CHAMBER_STR_FILES[@]}";   do CHAMBER_ARGS+=(-str "$f");   done

if [[ -f complex_chamber.parm7 ]]; then
    echo "complex_chamber.parm7 already exists, skipping chamber conversion."
else
    parmed <<EOF
chamber ${CHAMBER_ARGS[@]} -psf $PSF -crd $PDB
outparm complex_chamber.parm7
EOF
fi

echo
echo "=============================================="
echo " Step 2: atom count sanity check"
echo "=============================================="

CHAMBER_NATOM=$(cpptraj complex_chamber.parm7 <<'EOF' | awk '/contains/ {print $(NF-1)}'
parminfo
EOF
)
echo "complex_chamber.parm7 atom count: $CHAMBER_NATOM"

if [[ -n "$NAMD_LOG" && -f "$NAMD_LOG" ]]; then
    NAMD_NATOM=$(awk '/ATOMS/ {print $2; exit}' "$NAMD_LOG" || true)
    if [[ -n "$NAMD_NATOM" ]]; then
        echo "NAMD log atom count:            $NAMD_NATOM"
        if [[ "$CHAMBER_NATOM" != "$NAMD_NATOM" ]]; then
            echo "WARNING: atom count mismatch between chamber topology and NAMD log!"
            echo "         Double check you used the correct PSF/PDB (matching the"
            echo "         'structure'/'coordinates' lines in your NAMD config), and"
            echo "         that the chamber toppar file list didn't drop/add atoms."
            echo "         Continuing anyway, but results may be invalid."
        else
            echo "Atom counts match -- good."
        fi
    else
        echo "Could not find an 'N ATOMS' line in $NAMD_LOG; skipping automated check."
    fi
else
    echo "NAMD_LOG not set or not found; skipping automated atom-count check."
    echo "Manually verify complex_chamber.parm7's atom count matches your NAMD log's"
    echo "'Info: N ATOMS' line before trusting downstream results."
fi

echo
echo "=============================================="
echo " Step 3: PBC-wrapped, solvated NetCDF trajectory"
echo "=============================================="

if [[ -f complex_wet.nc ]]; then
    echo "complex_wet.nc already exists, skipping."
else
    cpptraj <<EOF
parm complex_chamber.parm7
trajin $DCD
autoimage
trajout complex_wet.nc netcdf
run
EOF
fi

echo
echo "=============================================="
echo " Step 4: per-interface + whole-assembly component topologies"
echo "=============================================="

build_interface () {
    local name="$1" lig_range="$2" rec_range="$3" cofactor_range="$4" strip_other="$5"

    echo "--- building interface: $name ---"

    # dry component topologies (complex / receptor / ligand)
    if [[ ! -f "${name}_complex.parm7" ]]; then
        ante-MMPBSA.py -p complex_chamber.parm7 \
            -c "${name}_complex.parm7" \
            -r "${name}_receptor.parm7" \
            -l "${name}_ligand.parm7" \
            -s "$STRIP_MASK,$strip_other" \
            -n ":$lig_range" \
            --radii "$RADII"
    fi

    # solvated topology matching -cp minus only solvent (needed because -sp
    # must be an exact atom-superset of -cp; the full assembly's solvated
    # topology is NOT a match once other chains are removed from -cp)
    if [[ ! -f "${name}_solvated.parm7" ]]; then
        parmed <<EOF
parm complex_chamber.parm7
strip :$strip_other
outparm ${name}_solvated.parm7
EOF
    fi

    # matching solvated trajectory
    if [[ ! -f "${name}_wet.nc" ]]; then
        cpptraj <<EOF
parm complex_chamber.parm7
trajin complex_wet.nc
strip :$strip_other
trajout ${name}_wet.nc netcdf
run
EOF
    fi
}

# Interface 1: ligand1 vs (receptor1 + cofactor1); strip everything from interface 2
build_interface "iface1" "$LIGAND1_RANGE" "$RECEPTOR1_RANGE" "$COFACTOR1_RANGE" \
    "$LIGAND2_RANGE,$RECEPTOR2_RANGE,$COFACTOR2_RANGE"

# Interface 2: ligand2 vs (receptor2 + cofactor2); strip everything from interface 1
build_interface "iface2" "$LIGAND2_RANGE" "$RECEPTOR2_RANGE" "$COFACTOR2_RANGE" \
    "$LIGAND1_RANGE,$RECEPTOR1_RANGE,$COFACTOR1_RANGE"

# Whole assembly: both ligands vs everything else (receptors + cofactors); nothing extra to strip
echo "--- building whole-assembly split ---"
if [[ ! -f whole_complex.parm7 ]]; then
    ante-MMPBSA.py -p complex_chamber.parm7 \
        -c whole_complex.parm7 \
        -r whole_receptor.parm7 \
        -l whole_ligand.parm7 \
        -s "$STRIP_MASK" \
        -n ":$LIGAND_WHOLE_RANGE" \
        --radii "$RADII"
fi
# whole assembly's solvated topology/trajectory = the full system minus only solvent,
# i.e. complex_chamber.parm7 / complex_wet.nc themselves -- no extra stripping needed.

echo
echo "=============================================="
echo " Step 5: write mmpbsa.in"
echo "=============================================="

cat > mmpbsa.in << EOF
&general
  startframe=$START_FRAME, endframe=$END_FRAME, interval=$INTERVAL, verbose=2,
  strip_mask="$STRIP_MASK",
/
&gb
  igb=$GB_MODEL, saltcon=$SALTCON,
/
&decomp
  idecomp=2, dec_verbose=1,
  print_res="$DECOMP_PRINT_RES"
/
EOF
cat mmpbsa.in

echo
echo "=============================================="
echo " Step 6: run MMPBSA.py for each split"
echo "=============================================="

run_mmpbsa () {
    local name="$1" sp="$2" y="$3"
    echo "--- running MMPBSA.py: $name ---"
    mkdir -p "${name}_run"
    (
        cd "${name}_run"
        MMPBSA.py -O -i ../mmpbsa.in -o "FINAL_RESULTS_${name}.dat" \
            -sp "../$sp" \
            -cp "../${name}_complex.parm7" \
            -rp "../${name}_receptor.parm7" \
            -lp "../${name}_ligand.parm7" \
            -y "../$y"
    )
    echo "Done: ${name}_run/FINAL_RESULTS_${name}.dat"
}

run_mmpbsa "iface1" "iface1_solvated.parm7" "iface1_wet.nc"
run_mmpbsa "iface2" "iface2_solvated.parm7" "iface2_wet.nc"
run_mmpbsa "whole"  "complex_chamber.parm7" "complex_wet.nc"

echo
echo "=============================================="
echo " Pipeline complete."
echo "=============================================="
echo "Results:"
echo "  iface1_run/FINAL_RESULTS_iface1.dat"
echo "  iface2_run/FINAL_RESULTS_iface2.dat"
echo "  whole_run/FINAL_RESULTS_whole.dat"
echo
echo "Decomposition files (per-residue):"
echo "  iface1_run/FINAL_DECOMP_MMPBSA.dat  (etc.)"
echo
echo "Analyze with analyze_mmpbsa.py, e.g.:"
echo "  python3 analyze_mmpbsa.py compare \\"
echo "      --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \\"
echo "      --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \\"
echo "      --label whole  whole_run/FINAL_RESULTS_whole.dat"
echo
echo "REMINDER: with END_FRAME=$END_FRAME frames, treat results as preliminary."
echo "Check SEM in FINAL_RESULTS_*.dat and per-residue decomposition before"
echo "drawing conclusions, especially for small differences between variants."
