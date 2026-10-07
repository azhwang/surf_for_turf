#!/usr/bin/env python3
"""
make_extensions.py -- set up long production extensions for selected replicates.

Purpose: the 1 ns x 5 replicate runs cannot resolve the experimental
goldfish-vs-human difference (2:2:1 vs 2:2:2), and human 2:2:1 splits into a
tight and a loose heparin-binding state that had not converged after 1 ns.
These extensions continue chosen replicates to TOTAL_NS of production so the
heparin-binding state and the MM/GBSA energy can converge.

Each replicate continues from its step5_production.coor/vel/xsc in 1 ns
segments (ext01 ... extNN), with the same NPT settings as step5_production.inp,
a frame every 100 ps, and a per-segment Langevin seed. ext_slurm.sh is
self-resubmitting: each job finds the next missing segment (no extNN.coor yet),
runs it, and submits the following one, so a segment killed by the walltime
just reruns and the chain stops by itself at TOTAL_NS.

Files written into each chosen namdrepN/: ext_template.inp, ext_slurm.sh
Top level: submit_extensions.sh, ext_status.sh

Usage: python3 make_extensions.py
"""
import re
from pathlib import Path

BASE = Path(__file__).resolve().parent
TOTAL_NS = 50
SEG_STEPS = 500000                       # 1 ns at 2 fs
NAMD = "/home/azhwang/Surf_for_Turf/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2"
# human 2:2:1: reps 1 and 4 were in the tight heparin-binding state, rep 2 in the loose one
SELECTION = {
    "1e0o/human": [1, 2, 4],
    "1e0o/goldfish": [1, 2, 3],
    "1fq9/human": [1, 2, 3],
    "1fq9/goldfish": [2, 3, 5],
}
LABEL = {"1e0o": "221", "1fq9": "222"}

SLURM = """#!/bin/bash
#SBATCH -p sched_mit_buehler_gpu
#SBATCH --job-name={job}
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=slurmx-%j.out
#
# Self-resubmitting production extension ({total} x 1 ns segments).
# Runs the next missing segment, then submits itself again until extNN.coor
# exists for all segments. Safe to resubmit by hand at any time (`sbatch ext_slurm.sh`):
# it always resumes from the last completed segment.

cd "${{SLURM_SUBMIT_DIR:-$(dirname "$0")}}"
TOTAL={total}
SEG_STEPS={seg}
BASE_SEED={seed}
NAMD={namd}

n=1
while [ -f "$(printf 'ext%02d.coor' $n)" ]; do n=$((n+1)); done
if [ $n -gt $TOTAL ]; then echo "all $TOTAL segments done"; exit 0; fi
out=$(printf 'ext%02d' $n)
if [ $n -eq 1 ]; then prev=step5_production; else prev=$(printf 'ext%02d' $((n-1))); fi
first=$((500000 + (n-1)*SEG_STEPS))

sed -e "s/__INPUT__/$prev/" -e "s/__OUTPUT__/$out/" -e "s/__STEPS__/$SEG_STEPS/g" \\
    -e "s/__FIRST__/$first/" -e "s/__SEED__/$((BASE_SEED + n))/" ext_template.inp > $out.inp

echo "$(date '+%F %T') segment $out (from $prev) on $(hostname)"
$NAMD +p${{SLURM_CPUS_PER_TASK:-4}} $out.inp > $out.log
if grep -q "End of program" $out.log && [ -f $out.coor ]; then
    echo "$(date '+%F %T') $out done"
    [ $n -lt $TOTAL ] && sbatch ext_slurm.sh
else
    echo "$(date '+%F %T') $out FAILED -- see $out.log; chain stopped" >&2
    exit 1
fi
"""

SUBMIT_ALL = """#!/bin/bash
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
"""

STATUS = """#!/bin/bash
# Progress of every extension chain: completed 1 ns segments and last segment timing.
cd "$(dirname "$0")"
for r in $(ls -d */*/sim/namdrep*/ext_slurm.sh 2>/dev/null | xargs -n1 dirname); do
    done_n=$(ls $r/ext[0-9][0-9].coor 2>/dev/null | wc -l | tr -d ' ')
    # timing of the last FINISHED segment (the newest log is usually the one still running)
    last=$( [ "$done_n" -gt 0 ] && printf '%s/ext%02d.log' "$r" "$done_n" )
    wall=$( [ -n "$last" ] && grep -h WallClock "$last" | awk '{printf "%.1f h", $2/3600}')
    fail=$(grep -l "FAILED" $r/slurmx-*.out 2>/dev/null | wc -l | tr -d ' ')
    printf "%-36s %3s ns done   last segment wall: %-8s %s\\n" "$r" "$done_n" "${wall:--}" "$( [ "$fail" != 0 ] && echo "FAILURES: $fail" )"
done
"""


def make_template(prod_inp: str) -> str:
    t = prod_inp
    subs = [
        (r"(?m)^outputName\s+step5_production;", "outputName              __OUTPUT__;"),
        (r"(?m)^set inputname\s+step4_equilibration;", "set inputname           __INPUT__;\nfirsttimestep           __FIRST__;"),
        (r"(?m)^seed\s+\d+;", "seed                    __SEED__;"),
        # numsteps is the ABSOLUTE last step in NAMD, so with firsttimestep > 0 it must be large;
        # the segment length is set by `run` (same pattern as CHARMM-GUI's step4 input)
        (r"(?m)^numsteps\s+\d+;", "numsteps                90000000;"),
        (r"(?m)^run\s+\d+;", "run                     __STEPS__;"),
    ]
    for pat, rep in subs:
        t, n = re.subn(pat, rep, t)
        if n != 1:
            raise SystemExit(f"template substitution failed for {pat!r} ({n} matches)")
    return "# Production extension segment, generated from step5_production.inp by make_extensions.py\n" + t


def main():
    made = []
    for system, reps in SELECTION.items():
        tk, sp = system.split("/")
        for r in reps:
            d = BASE / system / "sim" / f"namdrep{r}"
            if not (d / "step5_production.inp").exists():
                raise SystemExit(f"{d}/step5_production.inp missing")
            # the restart files are read on the cluster; a local copy can be missing after
            # an incomplete download, so only warn here
            for f in ("step5_production.coor", "step5_production.vel", "step5_production.xsc"):
                if not (d / f).exists():
                    print(f"WARNING: {d.relative_to(BASE)}/{f} not present locally -- make sure it exists on the cluster")
            prod = (d / "step5_production.inp").read_text()
            seed = int(re.search(r"(?m)^seed\s+(\d+);", prod).group(1)) * 100
            (d / "ext_template.inp").write_text(make_template(prod))
            job = f"ext_{LABEL[tk]}_{sp}_r{r}"
            (d / "ext_slurm.sh").write_text(SLURM.format(job=job, total=TOTAL_NS, seg=SEG_STEPS, seed=seed, namd=NAMD))
            (d / "ext_slurm.sh").chmod(0o755)
            made.append(f"{system}/sim/namdrep{r}  (job {job}, base seed {seed})")
    for name, txt in (("submit_extensions.sh", SUBMIT_ALL), ("ext_status.sh", STATUS)):
        (BASE / name).write_text(txt)
        (BASE / name).chmod(0o755)
    print("\n".join(made))


if __name__ == "__main__":
    main()
