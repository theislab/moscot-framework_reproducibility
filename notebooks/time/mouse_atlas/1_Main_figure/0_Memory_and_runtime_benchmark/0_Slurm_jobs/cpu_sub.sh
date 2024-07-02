#!/bin/bash
#SBATCH --job-name=Moscot_job
#SBATCH --output=/home/icb/manuel.gander/Moscot_job.txt
#SBATCH --error=/home/icb/manuel.gander/error.txt
#SBATCH --nice=10000
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
# SBATCH --nodelist=cpusrv05


#
k="$1"
method="$2"
tau="$3"

source jax2/bin/activate

python benchmarker.py --param1 "$k" --param2 "$method" --param3 "$tau"
