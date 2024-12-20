#!/bin/bash
#SBATCH -J run_chromvar      # Job name
#SBATCH -o run_chromvar.log       # Name of stdout output file
#SBATCH -e run_chromvar.log        # Name of stderr error file
#SBATCH -p cpu_p   # Specify the partition name
#SBATCH -n 4                # Number of tasks (usually cores) to run
#SBATCH --qos cpu_normal
#SBATCH -t 16:00:00
#SBATCH --mem 128G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/icb/dominik.klein/mambaforge/envs/r_env


# Run R script
Rscript run_chromVar.R
