#!/bin/bash
# Run ON THE CLUSTER from rebuild_v2/. Starts the self-resubmitting extension
# chain in every prepared replicate (or only the given ones, e.g. 1e0o/human/sim/namdrep1).
cd "$(dirname "$0")"
reps=("$@"); [ ${#reps[@]} -eq 0 ] && reps=($(ls -d */*/sim/namdrep*/ext_slurm.sh | xargs -n1 dirname))
for r in "${reps[@]}"; do
    if squeue -u "$USER" -h -o "%j" | grep -qx "$(awk -F= '/job-name/{print $2}' "$r/ext_slurm.sh")"; then
        echo "== $r: already queued/running, skipped"; continue
    fi
    echo "== $r"; (cd "$r" && sbatch ext_slurm.sh)
done
