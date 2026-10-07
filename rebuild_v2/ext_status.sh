#!/bin/bash
# Progress of every extension chain: completed 1 ns segments and last segment timing.
cd "$(dirname "$0")"
for r in $(ls -d */*/sim/namdrep*/ext_slurm.sh 2>/dev/null | xargs -n1 dirname); do
    done_n=$(ls $r/ext[0-9][0-9].coor 2>/dev/null | wc -l | tr -d ' ')
    # timing of the last FINISHED segment (the newest log is usually the one still running)
    last=$( [ "$done_n" -gt 0 ] && printf '%s/ext%02d.log' "$r" "$done_n" )
    wall=$( [ -n "$last" ] && grep -h WallClock "$last" | awk '{printf "%.1f h", $2/3600}')
    fail=$(grep -l "FAILED" $r/slurmx-*.out 2>/dev/null | wc -l | tr -d ' ')
    printf "%-36s %3s ns done   last segment wall: %-8s %s\n" "$r" "$done_n" "${wall:--}" "$( [ "$fail" != 0 ] && echo "FAILURES: $fail" )"
done
