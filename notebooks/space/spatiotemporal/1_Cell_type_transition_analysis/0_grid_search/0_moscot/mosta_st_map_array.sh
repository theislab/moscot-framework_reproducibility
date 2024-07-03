#!/bin/bash
#SBATCH --job-name=unbalanced_moscot_arrayjob
#SBATCH --array=0-41%8
#SBATCH -c1
#SBATCH --gres=gpu:a5000:1
#SBATCH --time=4-0
#SBATCH --mem=50gb
#SBATCH --output=/logs/unbalanced_moscot_%a.log


PROJECT_DIR=""
DATA_DIR=""

# Specify the path to the config file
config=${DATA_DIR}/grid_config.txt

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


source ""

echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample is ${tps_couple} and the alpha is ${alpha}."

python3 ${PROJECT_DIR}/mosta_st_map_accuracies.py --alpha ${alpha} --epsilon ${epsilon} --rank ${rank} --tps_couple ${tps_couple} --gamma ${gamma} --cost ${cost}