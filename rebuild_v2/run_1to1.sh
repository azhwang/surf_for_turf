#!/bin/bash
# Build the 1:1 FGF2-FGFR systems from the dedicated 1:1 crystals, in order.
#   1. prepare_1to1_scaffold.py  close the D3 loop gaps in 1EVT chain C and 1DJS chain A
#                                with the equivalent loop of the complete 1FQ9 chain C
#   2. prepare_1to1_fgf.py       trim each ortholog's FGF2 model to the 1E0O-span core and
#                                superpose it on the crystal FGF chain
#   3. build_1to1.py             psfgen, local relaxation of the grafted loops and of the few
#                                tight interface contacts, solvation, ionisation, replicates
# Then upload onetoone/ to the cluster and run `bash submit_1to1.sh [fgfr1|fgfr2]` there.
set -e
cd "$(dirname "$0")"
python3 prepare_1to1_scaffold.py
python3 prepare_1to1_fgf.py
python3 build_1to1.py "$@"
