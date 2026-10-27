#!/bin/bash
# equilibration, then production only if equilibration finished cleanly
cd "$(dirname "$0")"
jid=$(sbatch --parsable eq_slurm.sh)
sbatch --dependency=afterok:$jid prod_slurm.sh
