#!/bin/bash
# Run ON THE CLUSTER from the copied rebuild_v2/ folder. Submits equilibration for every
# replicate of every system, with production chained on --dependency=afterok (production
# never starts unless its equilibration finished cleanly).
#   bash submit_all.sh                 # all 14 systems x 5 replicates
#   bash submit_all.sh 1e0o/goldfish   # one system
cd "$(dirname "$0")"
targets=("$@"); [ ${#targets[@]} -eq 0 ] && targets=(1fq9/* 1e0o/*)
for t in "${targets[@]}"; do
  for rep in "$t"/sim/namdrep*; do
    [ -f "$rep/submit.sh" ] || continue
    echo "== $rep"; bash "$rep/submit.sh"
  done
done
