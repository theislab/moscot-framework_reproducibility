#!/bin/bash
#SBATCH --job-name=gp
#SBATCH --output=/home/icb/manuel.gander/gb.txt
#SBATCH --error=/home/icb/manuel.gander/gb_error.txt
#SBATCH --nice=10000
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
# SBATCH --mem=50G
# SBATCH --cpus-per-task=12
#SBATCH --nodelist=supergpu08


k="$1"
method="$2"
tau="$3"

echo "$tau"
source jax2/bin/activate

python benchmarker.py --param1 "$k" --param2 "$method" --param3 "$tau"
