#!/bin/bash
#SBATCH -J gpu_run_link_peaks      # Job name
#SBATCH -o gpu_run_link_peaks.log       # Name of stdout output file
#SBATCH -e gpu_run_link_peaks.log        # Name of stderr error file
#SBATCH -p gpu_p   # Specify the partition name
#SBATCH -n 2                # Number of tasks (usually cores) to run
#SBATCH --qos=gpu_long
#SBATCH -t 48:00:00
#SBATCH --mem 128G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/icb/dominik.klein/mambaforge/envs/r_env


# Run R script
Rscript run_link_peaks.R
