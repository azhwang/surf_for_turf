#!/bin/bash
# Run ON THE CLUSTER from rebuild_v2/. For every chimera replicate: equilibration -> 1 ns production
# -> self-resubmitting 50 x 1 ns extension chain (each step only if the previous one finished cleanly).
#   bash submit_chimeras.sh            # all 4 chimera systems x 3 replicates (12 chains)
#   bash submit_chimeras.sh 1e0o       # 2:2:1 only (the primary test; 6 chains)
cd "$(dirname "$0")"
tmpls=("$@"); [ ${#tmpls[@]} -eq 0 ] && tmpls=(1e0o 1fq9)
for t in "${tmpls[@]}"; do
  for rep in "$t"/chim_*/sim/namdrep*; do
    [ -f "$rep/submit_chain.sh" ] || continue
    echo "== $rep"; bash "$rep/submit_chain.sh"
  done
done
