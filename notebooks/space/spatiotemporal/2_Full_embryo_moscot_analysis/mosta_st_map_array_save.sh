#!/bin/bash
#SBATCH --job-name=save_rank_moscot_arrayjob
#SBATCH --array=4-7
#SBATCH -c1
#SBATCH --gres=gpu:a5000:1
#SBATCH --time=4-0
#SBATCH --mem=50gb
#SBATCH --output=/logs/save_rank_moscot_%a.log


PROJECT_DIR=""
DATA_DIR=""

# Specify the path to the config file
config=${DATA_DIR}/save_config.txt

# Extract the alpha name for the current $SLURM_ARRAY_TASK_ID
alpha=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Extract the epsilon for the current $SLURM_ARRAY_TASK_ID
epsilon=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

# Extract the rank for the current $SLURM_ARRAY_TASK_ID
rank=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)

# Extract the tps_couple for the current $SLURM_ARRAY_TASK_ID
tps_couple=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)

# Extract the gamma for the current $SLURM_ARRAY_TASK_ID
gamma=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)

# Extract the cost for the current $SLURM_ARRAY_TASK_ID
cost=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)

# Extract tau_a for the current $SLURM_ARRAY_TASK_ID
tau_a=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $config)

# Extract tau_a for the current $SLURM_ARRAY_TASK_ID
tau_b=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $9}' $config)


source ""

echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample number is ${tps_couple}."

python3 ${PROJECT_DIR}/mosta_st_map_accuracies_save.py --alpha ${alpha} --epsilon ${epsilon} --rank ${rank} --tps_couple ${tps_couple} --gamma ${gamma} --cost ${cost} --tau_a ${tau_b} --tau_a ${tau_b}