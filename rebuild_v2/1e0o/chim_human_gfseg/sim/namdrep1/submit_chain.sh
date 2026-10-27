#!/bin/bash
# equilibration -> 1 ns production -> 50 x 1 ns extension chain, each only if the previous step finished cleanly
cd "$(dirname "$0")"
e=$(sbatch --parsable eq_slurm.sh)
p=$(sbatch --parsable --dependency=afterok:$e prod_slurm.sh)
sbatch --dependency=afterok:$p ext_slurm.sh
