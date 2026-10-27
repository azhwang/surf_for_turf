#!/usr/bin/env bash
#
# queue_chimeras.sh [replicate|template ...] -- run the chimera analysis over many replicates,
# PAR at a time, each through fetch_and_analyze_chimeras.sh.
#
# Concurrency is limited by MEMORY, not CPU: each replicate's ParmEd/ante-MMPBSA steps load the
# full ~205k-atom solvated system and take several GB, and the project notes record that 8 at
# once swapped on this 64 GB Mac. PAR=3 leaves room for one replicate already in flight.
#
# Already-analysed replicates (r90.txt + decomp_run present) are skipped, so this is safe to
# re-run after a failure or a reboot. Per-replicate log: chim_<tmpl>_<system>_<rep>.log
#
# Usage:  bash queue_chimeras.sh 1e0o              # the six 2:2:1 replicates
#         bash queue_chimeras.sh                   # all 12
#         PAR=2 bash queue_chimeras.sh 1e0o        # gentler on memory
set -uo pipefail
cd "$(dirname "$0")"
PAR=${PAR:-3}

ALL=()
for sys in 1e0o/chim_human_gfseg 1e0o/chim_goldfish_hseg 1fq9/chim_human_gfseg 1fq9/chim_goldfish_hseg; do
  for r in 1 2 3; do ALL+=("$sys/sim/namdrep$r"); done
done

REPS=()
if [ $# -eq 0 ]; then REPS=("${ALL[@]}"); else
  for a in "$@"; do
    case "$a" in
      1e0o|1fq9) for r in "${ALL[@]}"; do [[ $r == $a/* ]] && REPS+=("$r"); done ;;
      *) REPS+=("$a") ;;
    esac
  done
fi

QUEUE=()
for r in "${REPS[@]}"; do
  if [ -f "$r/ext_analysis/r90.txt" ] && [ -f "$r/ext_analysis/decomp_run/FINAL_DECOMP_decomp.dat" ]; then
    echo "skip (done)      $r"; continue
  fi
  if [ "$(ls "$r"/ext[0-9][0-9].dcd 2>/dev/null | wc -l | tr -d ' ')" -ne 50 ]; then
    echo "skip (no data)   $r  -- fetch it first"; continue
  fi
  # a replicate currently being processed by another invocation
  if pgrep -f "fetch_and_analyze_chimeras.sh $r\$" >/dev/null 2>&1; then
    echo "skip (running)   $r"; continue
  fi
  QUEUE+=("$r"); echo "queued           $r"
done

[ ${#QUEUE[@]} -eq 0 ] && { echo "nothing to do"; exit 0; }
echo
echo "running ${#QUEUE[@]} replicate(s), $PAR at a time -- $(date '+%F %T')"
printf '%s\n' "${QUEUE[@]}" | xargs -P "$PAR" -I{} bash -c '
  rep="$1"; log="chim_$(echo "$rep" | sed "s|/sim/|_|; s|/|_|g").log"
  echo "START $rep -> $log"
  if bash fetch_and_analyze_chimeras.sh "$rep" > "$log" 2>&1; then
    echo "OK    $rep"
  else
    echo "FAIL  $rep  (see $log)"; tail -3 "$log" | sed "s/^/        /"
  fi' _ {}

echo
echo "=== summary $(date '+%F %T')"
for r in "${REPS[@]}"; do
  printf "%-44s r90 %-3s decomp %-3s\n" "$r" \
    "$([ -f "$r/ext_analysis/r90.txt" ] && echo yes || echo NO)" \
    "$([ -f "$r/ext_analysis/decomp_run/FINAL_DECOMP_decomp.dat" ] && echo yes || echo NO)"
done
