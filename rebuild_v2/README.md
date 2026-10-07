# Consistent FGF–FGFR–heparin rebuild (rebuild_v2)

Rebuilds all 14 systems (7 species × 2:2:2/1FQ9 and 2:2:1/1E0O) with one protocol,
fixing the ring-piercing and inconsistency problems in the original CHARMM-GUI builds.

## What changed vs the original builds
- **FGF trimmed to the template-resolved core** (the residues 1FQ9/1E0O actually model):
  129 residues for 1FQ9, 127 for 1E0O, identical span for both FGF copies. Removes the
  14–26-residue disordered N-terminal tail (unresolved in every FGF crystal structure),
  which was threaded through FGFR rings in several builds and differed in length by species.
- **FGF placed by its core onto the crystal FGF chains**, not by whole-model alignment.
  First 1–3 core residues that deviated >1.5 Å are rebuilt on the crystal backbone;
  conserved interface side chains (within 6 Å of FGFR/heparin) use crystal rotamers.
  Result: FGF–FGFR pairing matches the crystal in every system.
- **One scaffold per template** for FGFR + heparin (1FQ9: goldfish build, heparin intact
  8+6 residues; 1E0O: chinook build). The original bovine_stab/zebrafish 2:2:2 builds had
  heparin broken into 3–4 fragments.
- Same box (1FQ9 124 Å, 1E0O 129 Å cubic), 0.15 M KCl (POT/CLA), same NAMD inputs,
  restraints and run lengths as the originals; explicit per-replicate `seed`.
- `submit.sh` chains production on equilibration (`--dependency=afterok`).

## Layout
    <tmpl>/<species>/fgf_PROC.pdb, fgf_PROD.pdb, core_report.txt   placed FGF cores
    <tmpl>/<species>/build/     psfgen/solvate/autoionize inputs, logs, ionized.psf/pdb
    <tmpl>/<species>/sim/       toppar/, step3_pbcsetup.str, namdrep1..5/  <- copy to cluster

## Rebuild from scratch
    ./run_rebuild.sh            # prepare -> precheck -> fix_rotamers -> precheck -> build

## Run MD (on the cluster)
    rsync -a --include='*/' --include='sim/***' --exclude='*' rebuild_v2/ cluster:~/Surf_for_Turf/rebuild_v2/
    ssh cluster; cd ~/Surf_for_Turf/rebuild_v2; bash submit_all.sh

## After MD
Copy the namdrep*/ folders back, then run MM/GBSA and QC exactly as for v2:
    bash run_all_mmpbsa_v2.sh rebuild_v2/1fq9 rebuild_v2/1e0o
    python3 qc_ring_piercing_v2.py; python3 summarize_mmpbsa_v2.py

## Long production extensions (2026-10-05)
Why: 1 ns x 5 replicates cannot resolve the experimental goldfish-vs-human difference, and human
2:2:1 shows a tight vs loose heparin-binding state (reps 1,4 vs 2,3,5) that had not converged.

Selected replicates (`make_extensions.py` -> `SELECTION`):
2:2:1 human r1, r2, r4 (both states) | 2:2:1 goldfish r1-3 | 2:2:2 human r1-3 | 2:2:2 goldfish r2,r3,r5 (r1 restart files incomplete on the cluster)

Each continues from `step5_production.coor/vel/xsc` in 1 ns segments `ext01..ext50` (50 ns),
same NPT settings, frame every 100 ps, `+p4`. `ext_slurm.sh` is self-resubmitting and always
resumes from the last completed segment (a walltime kill just reruns that segment).

On the cluster, from rebuild_v2/:
    bash submit_extensions.sh        # start all 12 chains (skips ones already queued)
    bash ext_status.sh               # ns completed per run + wall time of the last segment
Restart one chain by hand: `cd <rep> && sbatch ext_slurm.sh`.
