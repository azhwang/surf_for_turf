#!/bin/bash
#SBATCH -p sched_mit_buehler_gpu     # partition name
#SBATCH --job-name=222_bovine_testeq # name for your job
#SBATCH --gres=gpu:1                 # if you need GPUs
#SBATCH --ntasks=1                   # number of tasks (often 1 for serial jobs)
#SBATCH --cpus-per-task=4            # CPU cores per task
#SBATCH --mem=64G                    # memory per node
#SBATCH --time=06:00:00              # max walltime (HH:MM:SS)
#SBATCH --output=slurme-test-%j.out  # output file (%j = job ID) to capture logs for debugging

/home/azhwang/Surf_for_Turf/NAMD_3.0.2_Linux-x86_64-multicore-CUDA/namd3 step4_equilibration.inp > eq.log
