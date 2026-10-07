for d in */; do
    echo "Processing directory: $d"
    # Your commands here, e.g., cd "$d" && ... && cd ..
    cd $d
    sbatch md-slurm.sh
    sbatch md-slurm_rep2.sh
    sbatch md-slurm_rep3.sh
    sbatch md-slurm_rep4.sh
    sbatch md-slurm_rep5.sh
    cd ..
done
