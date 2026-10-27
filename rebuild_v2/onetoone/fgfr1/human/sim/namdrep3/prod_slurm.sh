#!/bin/bash
#SBATCH -p sched_mit_buehler_gpu
#SBATCH --job-name=11r1_human_r3_pr
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=slurmp-%j.out

/home/azhwang/Surf_for_Turf/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2 +p${SLURM_CPUS_PER_TASK:-4} step5_production.inp > prod.log
