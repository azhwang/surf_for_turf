#!/usr/bin/env bash
#
# run_all_mmpbsa_v2.sh -- run run_mmpbsa_pipeline_v2.sh in every replicate
# under the given roots. A replicate whose split is ambiguous or fails is
# logged and skipped; the batch continues.
#
# Usage:
#   bash run_all_mmpbsa_v2.sh all_metrics_analysis all_metrics_analysis_221
#   JOBS=8 bash run_all_mmpbsa_v2.sh rebuild_v2/1fq9 rebuild_v2/1e0o    # 8 replicates at once
#
# Each replicate is single-threaded, but the ParmEd/ante-MMPBSA steps load the
# full solvated system (~180-205k atoms) and take several GB each. On the 64 GB
# Mac, JOBS=8 ran out of memory and swapped (slower than serial); use JOBS=4.
# Default JOBS=1.
# A replicate being processed holds <namdrepN>/mmpbsa_v2.lock (an atomic mkdir),
# so overlapping batches never work on the same replicate; replicates that
# already have whole_run/FINAL_RESULTS_whole.dat are skipped.
#
# Per-replicate log: <namdrepN>/mmpbsa_v2.log
# Batch summary:     mmpbsa_v2_batch.log (in the current directory)

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$HOME/anaconda3/envs/AmberTools23/bin:$PATH"
export AMBERHOME="$HOME/anaconda3/envs/AmberTools23"
export BATCH_LOG="$PWD/mmpbsa_v2_batch.log"
export PIPELINE="$HERE/run_mmpbsa_pipeline_v2.sh"
JOBS="${JOBS:-1}"

run_one () {
    local rep="$1"
    if [[ ! -f "$rep/step3_input.psf" || ! -f "$rep/step3_input.pdb" || ! -f "$rep/step5_production.dcd" ]]; then
        echo "$(date '+%F %T') SKIP (missing step3_input.psf/pdb or step5_production.dcd) $rep" | tee -a "$BATCH_LOG"
        return 0
    fi
    if [[ -f "$rep/mmpbsa_v2/whole_run/FINAL_RESULTS_whole.dat" ]]; then
        echo "$(date '+%F %T') SKIP (done)  $rep" | tee -a "$BATCH_LOG"
        return 0
    fi
    if ! mkdir "$rep/mmpbsa_v2.lock" 2>/dev/null; then
        echo "$(date '+%F %T') SKIP (locked) $rep" | tee -a "$BATCH_LOG"
        return 0
    fi
    echo "$(date '+%F %T') START        $rep" | tee -a "$BATCH_LOG"
    if (cd "$rep" && bash "$PIPELINE" > mmpbsa_v2.log 2>&1); then
        echo "$(date '+%F %T') OK           $rep" | tee -a "$BATCH_LOG"
    else
        reason=$(grep -m1 -E "SPLIT ERROR|ERROR|MISMATCH|Error" "$rep/mmpbsa_v2.log" || echo "see log")
        echo "$(date '+%F %T') FAILED       $rep :: $reason" | tee -a "$BATCH_LOG"
    fi
    rmdir "$rep/mmpbsa_v2.lock"
}
export -f run_one

for root in "$@"; do
    find "$root" -mindepth 2 -maxdepth 3 -type d -name 'namdrep*'
done | sort | xargs -P "$JOBS" -I{} bash -c 'run_one "$1"' _ {}
