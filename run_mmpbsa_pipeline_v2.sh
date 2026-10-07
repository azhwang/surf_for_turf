#!/usr/bin/env bash
#
# run_mmpbsa_pipeline_v2.sh
#
# Corrected MM/GBSA pipeline for FGF-FGFR-heparin assemblies (2:2:2 and
# 2:2:1). Runs inside one replicate directory (namdrepN/) and writes
# everything to ./mmpbsa_v2/, leaving v1 outputs untouched for
# reproducibility.
#
# Differences from run_mmpbsa_pipeline.sh (v1):
#   - No hardcoded residue ranges. derive_mmpbsa_split.py reads THIS
#     replicate's PSF/PDB, identifies FGF/FGFR by sequence and heparin by
#     CAR* segid, pairs chains by contacts, and aborts if anything is
#     ambiguous. (v1 reused 1FQ9-goldfish ranges for every system.)
#   - ante-MMPBSA.py's -n ligand mask is written in the numbering it
#     actually sees (after the other interface and solvent are stripped).
#     v1 used full-system numbering, which mis-split interface 2 everywhere.
#   - FGF is the ligand. Heparin is grouped with the FGFR (receptor) side.
#     For 2:2:1 the single shared heparin is kept in both interfaces.
#   - Each interface's dry topologies are built from that interface's own
#     solvated topology, so -sp/-cp/-rp/-lp are consistent by construction.
#   - After ante-MMPBSA.py, component residue counts are verified against
#     the derived split before MMPBSA.py is run.
#   - The solvent/ion strip mask is derived from the PSF (WAT + actual ions).
#   - If FGF-FGFR pairing is ambiguous (each FGF contacts both receptors, as
#     in the 1FQ9 bovine_stab/zebrafish builds), the per-interface split is
#     skipped and only the whole-assembly split is computed.
#
# Usage (from a namdrepN/ directory, with AmberTools on PATH):
#   bash ~/Desktop/Surf_for_Turf/run_mmpbsa_pipeline_v2.sh
#
# Reuses ../complex_chamber.parm7 and complex_wet.nc from the v1 run in the
# same directory if present (they don't depend on the split); otherwise
# builds them.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPLIT_HELPER="$HERE/derive_mmpbsa_split.py"
OUT_DIR="mmpbsa_v2"

# ============================================================================
# CONFIG
# ============================================================================
PSF="step3_input.psf"
PDB="step3_input.pdb"
DCD="step5_production.dcd"
NAMD_LOG="step5_production.log"     # atom-count sanity check; "" to skip
# chamber needs CHARMM-format toppar (with .rtf files). In CHARMM-GUI builds that is
# ../toppar (namdrepN/toppar is the NAMD-format copy, without .rtf); in rebuild_v2 it is
# ../toppar_charmm. Override with CHAMBER_TOPPAR=<dir>.
TOPPAR_DIR=""
for cand in "${CHAMBER_TOPPAR:-}" ../toppar_charmm toppar_charmm ../toppar toppar; do
    if [[ -n "$cand" && -f "$cand/top_all36_prot.rtf" ]]; then TOPPAR_DIR="$cand"; break; fi
done
if [[ -z "$TOPPAR_DIR" ]]; then
    echo "ERROR: no CHARMM-format toppar (with top_all36_prot.rtf) found; set CHAMBER_TOPPAR" >&2
    exit 1
fi

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
    "$TOPPAR_DIR/toppar_all36_carb_imlab.str"       # required for IdoA (CC3261 type)
)

# MMPBSA settings -- identical to v1 so results are comparable
START_FRAME=1
END_FRAME=10
INTERVAL=1
SALTCON=0.150
GB_MODEL=5
DECOMP_PRINT_RES="all"
RADII="mbondi2"
# ============================================================================

echo "== Step 1: chamber topology + wet trajectory (shared with v1) =="

if [[ ! -f complex_chamber.parm7 ]]; then
    CHAMBER_ARGS=()
    for f in "${CHAMBER_TOP_FILES[@]}";   do CHAMBER_ARGS+=(-top "$f");   done
    for f in "${CHAMBER_PARAM_FILES[@]}"; do CHAMBER_ARGS+=(-param "$f"); done
    for f in "${CHAMBER_STR_FILES[@]}";   do CHAMBER_ARGS+=(-str "$f");   done
    parmed <<EOF
chamber ${CHAMBER_ARGS[@]} -psf $PSF -crd $PDB
outparm complex_chamber.parm7
EOF
else
    echo "reusing complex_chamber.parm7"
fi

CHAMBER_NATOM=$(cpptraj complex_chamber.parm7 <<'EOF' | awk '/contains/ {print $(NF-1)}'
parminfo
EOF
)
echo "complex_chamber.parm7 atoms: $CHAMBER_NATOM"
if [[ -n "$NAMD_LOG" && -f "$NAMD_LOG" ]]; then
    NAMD_NATOM=$(awk '/ATOMS/ {print $2; exit}' "$NAMD_LOG" || true)
    if [[ -n "$NAMD_NATOM" && "$NAMD_NATOM" != "$CHAMBER_NATOM" ]]; then
        echo "ERROR: atom count mismatch (NAMD log $NAMD_NATOM vs chamber $CHAMBER_NATOM)" >&2
        exit 1
    fi
fi

if [[ ! -f complex_wet.nc ]]; then
    cpptraj <<EOF
parm complex_chamber.parm7
trajin $DCD
autoimage
trajout complex_wet.nc netcdf
run
EOF
else
    echo "reusing complex_wet.nc"
fi

echo
echo "== Step 2: derive split from this replicate's PSF/PDB =="

mkdir -p "$OUT_DIR"
SPLIT_VARS=$(python3 "$SPLIT_HELPER" --psf "$PSF" --pdb "$PDB" \
    --prmtop complex_chamber.parm7 --summary "$OUT_DIR/split_summary.txt")
eval "$SPLIT_VARS"
printf '%s\n' "$SPLIT_VARS" > "$OUT_DIR/split_vars.sh"

cd "$OUT_DIR"

echo
echo "== Step 3: component topologies + trajectories =="

build_interface () {
    local name="$1" other_mask="$2" lig_mask="$3"
    echo "--- $name ---"
    if [[ ! -f "${name}_solvated.parm7" ]]; then
        parmed <<EOF
parm ../complex_chamber.parm7
strip $other_mask
outparm ${name}_solvated.parm7
EOF
    fi
    if [[ ! -f "${name}_wet.nc" ]]; then
        cpptraj <<EOF
parm ../complex_chamber.parm7
trajin ../complex_wet.nc
strip $other_mask
trajout ${name}_wet.nc netcdf
run
EOF
    fi
    if [[ ! -f "${name}_complex.parm7" ]]; then
        ante-MMPBSA.py -p "${name}_solvated.parm7" \
            -c "${name}_complex.parm7" \
            -r "${name}_receptor.parm7" \
            -l "${name}_ligand.parm7" \
            -s "$STRIP_MASK" \
            -n "$lig_mask" \
            --radii "$RADII"
    fi
}

if [[ "$IFACE_SPLIT" == "1" ]]; then
    build_interface iface1 "$IFACE1_OTHER_MASK" "$IFACE1_LIGAND_MASK"
    build_interface iface2 "$IFACE2_OTHER_MASK" "$IFACE2_LIGAND_MASK"
    SPLITS=(iface1 iface2 whole)
else
    echo "NOTE: FGF-FGFR pairing is ambiguous in this structure -- per-interface"
    echo "      split skipped, whole-assembly only (see split_summary.txt)."
    SPLITS=(whole)
fi

echo "--- whole ---"
if [[ ! -f whole_complex.parm7 ]]; then
    ante-MMPBSA.py -p ../complex_chamber.parm7 \
        -c whole_complex.parm7 -r whole_receptor.parm7 -l whole_ligand.parm7 \
        -s "$STRIP_MASK" -n "$WHOLE_LIGAND_MASK" --radii "$RADII"
fi

echo
echo "== Step 4: verify component residue counts =="
python3 "$SPLIT_HELPER" --psf "../$PSF" --pdb "../$PDB" --check .

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

echo
echo "== Step 5: MMPBSA.py =="

run_mmpbsa () {
    local name="$1" sp="$2" y="$3"
    if [[ -f "${name}_run/FINAL_RESULTS_${name}.dat" ]]; then
        echo "--- $name: already done, skipping ---"
        return
    fi
    echo "--- $name ---"
    mkdir -p "${name}_run"
    (
        cd "${name}_run"
        MMPBSA.py -O -i ../mmpbsa.in -o "FINAL_RESULTS_${name}.dat" \
            -sp "$sp" \
            -cp "../${name}_complex.parm7" \
            -rp "../${name}_receptor.parm7" \
            -lp "../${name}_ligand.parm7" \
            -y "$y"
    )
}

if [[ "$IFACE_SPLIT" == "1" ]]; then
    run_mmpbsa iface1 ../iface1_solvated.parm7 ../iface1_wet.nc
    run_mmpbsa iface2 ../iface2_solvated.parm7 ../iface2_wet.nc
fi
run_mmpbsa whole  ../../complex_chamber.parm7 ../../complex_wet.nc

echo
echo "== Done: $(pwd) =="
cat split_summary.txt
for n in "${SPLITS[@]}"; do
    printf '%-7s ' "$n"; grep -m1 "DELTA TOTAL" "${n}_run/FINAL_RESULTS_${n}.dat"
done
