#!/bin/bash
#SBATCH -p sched_mit_buehler_gpu
#SBATCH --job-name=11r2_bovine_stab_r1_eq
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=slurme-%j.out

/home/azhwang/Surf_for_Turf/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2 +p${SLURM_CPUS_PER_TASK:-4} step4_equilibration.inp > eq.log
