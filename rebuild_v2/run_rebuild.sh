#!/usr/bin/env bash
# Consistent FGF rebuild, all steps in order. Re-running is safe (everything is regenerated).
#   1. prepare_fgf_core.py  trim FGF to the template-resolved core, place on crystal FGF,
#                           anchor deviating N-terminal residues, crystal side chains at
#                           conserved interface positions
#   2. precheck.py          write combined structures for the checks
#   3. fix_rotamers.py      chi scan for any FGF side chain still in a pierced ring
#   4. precheck.py          must report 0 piercings everywhere, or stop
#   5. build_all.py         psfgen + solvate + ionize + NAMD/slurm replicate folders
set -euo pipefail
cd "$(dirname "$0")"
python3 prepare_fgf_core.py
python3 precheck.py > /dev/null || true
python3 fix_rotamers.py
python3 precheck.py || { echo "piercings remain -- not building" >&2; exit 1; }
python3 build_all.py "$@"
