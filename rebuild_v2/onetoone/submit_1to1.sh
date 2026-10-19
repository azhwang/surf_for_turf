#!/bin/bash
# Run ON THE CLUSTER from onetoone/. Submits equilibration for every replicate, with
# production chained on --dependency=afterok.
#   bash submit_1to1.sh            # both receptors, 14 systems x 5 replicates
#   bash submit_1to1.sh fgfr1      # one receptor only
cd "$(dirname "$0")"
targets=("$@"); [ ${#targets[@]} -eq 0 ] && targets=(fgfr1 fgfr2)
for t in "${targets[@]}"; do
  for rep in "$t"/*/sim/namdrep*; do
    [ -f "$rep/submit.sh" ] || continue
    echo "== $rep"; bash "$rep/submit.sh"
  done
done
