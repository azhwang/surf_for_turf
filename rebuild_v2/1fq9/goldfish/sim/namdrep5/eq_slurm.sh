#!/bin/bash
#SBATCH -p sched_mit_buehler_gpu     # partition name
#SBATCH --job-name=222_goldfish_r5_eq
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=slurme-%j.out

/home/azhwang/Surf_for_Turf/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2 step4_equilibration.inp > eq.log
