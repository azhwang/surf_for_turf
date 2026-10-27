#!/usr/bin/env bash
#
# fetch_and_analyze_chimeras.sh [replicate ...] -- run the pre-registered chimera analysis
# (ANALYSIS_PLAN_chimeras.md) one replicate at a time, staging trajectories through local disk.
#
# Why staged: each replicate is ~1.1 GB of ext*.dcd and all 50 segments are needed (the analysis
# scripts slice win=segs[10:], i.e. they drop the first TEN FILES PRESENT, so a partial download
# would silently shift the window off the pre-specified 11-50). Twelve replicates at once is
# ~13 GB of trajectory plus ~3 GB of MM/GBSA intermediates; this keeps the peak near 1.5 GB by
# deleting each replicate's trajectories once its analysis outputs are written.
#
# The 2:2:1 (1e0o) replicates come first: H1 and H2 are defined on the 2:2:1 assembly, so those
# six decide the pre-registered result. The 1fq9 six are the companion 2:2:2 set.
#
# Cluster access costs a two-step auth per connection, so fetching is separable from analysing:
# --fetch pulls a whole template's six replicates in ONE rsync (one auth), and the analysis then
# runs with no cluster access at all. A replicate whose 50 segments are already on disk is never
# re-fetched, so the two modes compose.
#
# Usage:  conda activate AmberTools23
#         bash fetch_and_analyze_chimeras.sh --fetch 1e0o     # one auth, ~7 GB, no analysis
#         bash fetch_and_analyze_chimeras.sh 1e0o             # analyse those six, no cluster access
#         bash fetch_and_analyze_chimeras.sh --fetch 1fq9     # second auth
#         bash fetch_and_analyze_chimeras.sh 1fq9
#         bash fetch_and_analyze_chimeras.sh                 # all 12, fetching per replicate
#         bash fetch_and_analyze_chimeras.sh 1e0o/chim_human_gfseg/sim/namdrep1
set -uo pipefail
cd "$(dirname "$0")"

HOST=azhwang@orcd-login.mit.edu
REMOTE=Surf_for_Turf/extension_test
JOBS=${JOBS:-4}
SURF="$(cd .. && pwd)"                       # ~/Desktop/Surf_for_Turf, where the pipelines live
export PATH="$HOME/anaconda3/envs/AmberTools23/bin:$PATH"
export AMBERHOME="$HOME/anaconda3/envs/AmberTools23"
MIN_FREE_GB=${MIN_FREE_GB:-6}

FETCH_ONLY=0
if [ "${1:-}" = "--fetch" ]; then FETCH_ONLY=1; shift; fi

fetch() {   # fetch <replicate dir>; retries, since a ~1.2 GB transfer can stall and drop
  local REP="$1" try
  mkdir -p "$REP"
  for try in 1 2 3; do
    # all 50 segments, the inputs, and the 1 ns production trajectory the pipeline needs.
    # .coor/.vel/.restart are NOT fetched: nothing here continues the run.
    # --partial keeps a half-sent file so a retry resumes it; the ssh keepalives stop the
    # connection being dropped while rsync is busy building its file list on the far end.
    rsync -az --progress --partial --timeout=180 \
      -e 'ssh -o ServerAliveInterval=20 -o ServerAliveCountMax=10 -o TCPKeepAlive=yes' \
      --include='ext[0-9][0-9].dcd' --include='step3_input.psf' --include='step3_input.pdb' \
      --include='step5_production.dcd' --include='step3_pbcsetup.str' --exclude='*' \
      "$HOST:$REMOTE/$REP/" "$REP/" && return 0
    echo "   rsync attempt $try failed for $REP" >&2
    [ "$try" -lt 3 ] && { echo "   retrying in 30 s"; sleep 30; }
  done
  return 1
}

fetch_batch() {   # fetch_batch <template>  -- every still-incomplete replicate in ONE rsync,
  # so the cluster is contacted once and the two-step auth is answered once. The include list is
  # built from exactly the replicates that still need data, anchored to the transfer root so the
  # parent extension directories (goldfish/, human/) and the .coor/.vel/.log files stay behind.
  local tk="$1" rep inc=() need=()
  for rep in "${ALL[@]}"; do
    [[ $rep == $tk/* ]] || continue
    have_all_segments "$rep" && continue
    need+=("$rep")
    local sub="${rep#$tk/}"                       # chim_*/sim/namdrepN
    inc+=(--include="/${sub%/*}/" --include="/$sub/"
          --include="/$sub/ext[0-9][0-9].dcd"
          --include="/$sub/step3_input.psf" --include="/$sub/step3_input.pdb"
          --include="/$sub/step5_production.dcd")
  done
  if [ ${#need[@]} -eq 0 ]; then echo "   all $tk replicates already complete"; return 0; fi
  echo "   one rsync for ${#need[@]} replicate(s): ${need[*]#$tk/}"
  local try
  for try in 1 2 3; do
    rsync -az --progress --partial --prune-empty-dirs --timeout=180 \
      -e 'ssh -o ServerAliveInterval=20 -o ServerAliveCountMax=10 -o TCPKeepAlive=yes' \
      --include="/chim_*/" --include="/chim_*/sim/" --include="/chim_*/sim/step3_pbcsetup.str" \
      "${inc[@]}" --exclude='*' \
      "$HOST:$REMOTE/$tk/" "$tk/" && return 0
    echo "   attempt $try failed" >&2
    [ "$try" -lt 3 ] && { echo "   retrying in 30 s (this needs the auth again)"; sleep 30; }
  done
  return 1
}

have_all_segments() { [ "$(ls "$1"/ext[0-9][0-9].dcd 2>/dev/null | wc -l | tr -d ' ')" -eq 50 ]; }

ALL=()
for sys in 1e0o/chim_human_gfseg 1e0o/chim_goldfish_hseg 1fq9/chim_human_gfseg 1fq9/chim_goldfish_hseg; do
  for r in 1 2 3; do ALL+=("$sys/sim/namdrep$r"); done
done
REPS=()
if [ $# -eq 0 ]; then
  REPS=("${ALL[@]}")
else
  for a in "$@"; do
    case "$a" in
      1e0o|1fq9) for r in "${ALL[@]}"; do [[ $r == $a/* ]] && REPS+=("$r"); done ;;
      *)         REPS+=("$a") ;;
    esac
  done
fi

FAILED=()
if [ "$FETCH_ONLY" = 1 ]; then
  # one rsync per template when whole templates were asked for, else one per replicate
  tks=(); only_tmpl=1
  for a in "$@"; do case "$a" in 1e0o|1fq9) tks+=("$a") ;; *) only_tmpl=0 ;; esac; done
  [ $# -eq 0 ] && { tks=(1e0o 1fq9); only_tmpl=1; }
  if [ "$only_tmpl" = 1 ]; then
    for tk in "${tks[@]}"; do
      echo "== $tk"; fetch_batch "$tk" || FAILED+=("$tk")
    done
  else
    for REP in "${REPS[@]}"; do
      if have_all_segments "$REP"; then echo "   $REP already complete, skipping"; continue; fi
      echo "== $REP"; fetch "$REP" || { echo "!! gave up on $REP, continuing" >&2; FAILED+=("$REP"); }
    done
  fi
  [ ${#FAILED[@]} -eq 0 ] || { echo "not fetched: ${FAILED[*]} -- re-run --fetch to retry just these"; }
  echo "fetched; now run without --fetch to analyse (no cluster access needed)"
  exit 0
fi

for REP in "${REPS[@]}"; do
  echo "================ $REP  ($(date '+%F %T'))"
  free=$(df -g . | awk 'NR==2{print $4}')
  if [ "$free" -lt "$MIN_FREE_GB" ]; then
    echo "!! only ${free} GB free, need ${MIN_FREE_GB} -- stopping" >&2; exit 1
  fi

  if [ -f "$REP/ext_analysis/r90.txt" ] && [ -f "$REP/ext_analysis/decomp_run/FINAL_DECOMP_decomp.dat" ]; then
    echo "   already analysed, skipping"; continue
  fi

  if have_all_segments "$REP"; then
    echo "   segments already on disk, no cluster access needed"
  else
    fetch "$REP" || { echo "!! rsync failed for $REP" >&2; exit 1; }
  fi

  n=$(ls "$REP"/ext[0-9][0-9].dcd 2>/dev/null | wc -l | tr -d ' ')
  if [ "$n" -ne 50 ]; then echo "!! $REP has $n segments, expected 50 -- skipping" >&2; continue; fi
  u=$(ls -l "$REP"/ext[0-9][0-9].dcd | awk '{print $5}' | sort -u | wc -l | tr -d ' ')
  [ "$u" -eq 1 ] || echo "   note: $u distinct dcd sizes in $REP (expected 1) -- check before trusting the window"

  # 1. the 1 ns v2 MM/GBSA, which produces complex_chamber.parm7 and the splits.
  # Call the pipeline directly rather than run_all_mmpbsa_v2.sh: that wrapper takes a SPECIES
  # root and locates replicates with `find -mindepth 2 -maxdepth 3 -name namdrep*`, so handing
  # it a replicate directory matches nothing and it exits 0 without doing anything.
  if [ -f "$REP/mmpbsa_v2/whole_run/FINAL_RESULTS_whole.dat" ]; then
    echo "   step 1 already done"
  else
    ( cd "$REP" && bash "$SURF/run_mmpbsa_pipeline_v2.sh" > mmpbsa_v2.log 2>&1 ) || {
      echo "!! step 1 failed for $REP; last errors:" >&2
      grep -m5 -E "SPLIT ERROR|MISMATCH|Error|error" "$REP/mmpbsa_v2.log" >&2 || tail -5 "$REP/mmpbsa_v2.log" >&2
      exit 1; }
  fi
  [ -f "$REP/complex_chamber.parm7" ] || {
    echo "!! step 1 left no complex_chamber.parm7 in $REP -- see $REP/mmpbsa_v2.log" >&2; exit 1; }
  # 2. primary metric: Arg90-heparin contacts (cpptraj only)
  bash analyze_r90.sh "$REP" || { echo "!! analyze_r90 failed for $REP" >&2; exit 1; }
  # 3. planned MM/GBSA over the 11-50 window, needed by step 4
  bash analyze_extension.sh "$REP" || { echo "!! analyze_extension failed for $REP" >&2; exit 1; }
  # 4. per-residue decomposition: the second half of the engagement criterion
  bash analyze_exploratory.sh "$REP" || { echo "!! analyze_exploratory failed for $REP" >&2; exit 1; }

  # keep the analysis outputs and the topology; drop the trajectories
  rm -f "$REP"/ext[0-9][0-9].dcd "$REP"/step5_production.dcd "$REP"/complex_wet.nc
  echo "   done, $(du -sh "$REP" | cut -f1) retained"
done

echo "================ all requested replicates analysed"
