#!/bin/bash
#SBATCH -p sched_mit_buehler_gpu
#SBATCH --job-name=ext_221_goldfish_r3
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=slurmx-%j.out
#
# Self-resubmitting production extension (50 x 1 ns segments).
# Runs the next missing segment, then submits itself again until extNN.coor
# exists for all segments. Safe to resubmit by hand at any time (`sbatch ext_slurm.sh`):
# it always resumes from the last completed segment.

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
TOTAL=50
SEG_STEPS=500000
BASE_SEED=203100
NAMD=/home/azhwang/Surf_for_Turf/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2

n=1
while [ -f "$(printf 'ext%02d.coor' $n)" ]; do n=$((n+1)); done
if [ $n -gt $TOTAL ]; then echo "all $TOTAL segments done"; exit 0; fi
out=$(printf 'ext%02d' $n)
if [ $n -eq 1 ]; then prev=step5_production; else prev=$(printf 'ext%02d' $((n-1))); fi
first=$((500000 + (n-1)*SEG_STEPS))

sed -e "s/__INPUT__/$prev/" -e "s/__OUTPUT__/$out/" -e "s/__STEPS__/$SEG_STEPS/g" \
    -e "s/__FIRST__/$first/" -e "s/__SEED__/$((BASE_SEED + n))/" ext_template.inp > $out.inp

echo "$(date '+%F %T') segment $out (from $prev) on $(hostname)"
$NAMD +p${SLURM_CPUS_PER_TASK:-4} $out.inp > $out.log
if grep -q "End of program" $out.log && [ -f $out.coor ]; then
    echo "$(date '+%F %T') $out done"
    [ $n -lt $TOTAL ] && sbatch ext_slurm.sh
else
    echo "$(date '+%F %T') $out FAILED -- see $out.log; chain stopped" >&2
    exit 1
fi
