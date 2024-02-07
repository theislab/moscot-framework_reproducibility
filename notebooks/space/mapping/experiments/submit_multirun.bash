#!/bin/bash

#SBATCH -o "slurm_%j.out"
#SBATCH -e "slurm_%j.err"
#SBATCH -J cpu_sub
#SBATCH --qos cpu_normal
#SBATCH -p cpu_p
#SBATCH -t 1:00:00
#SBATCH --mem 20G

source ~/.bashrc
conda activate moscot

unset SLURM_CPU_BIND
python experiment.py --multirun launcher=icb