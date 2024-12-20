#! /bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/icb/dominik.klein/mambaforge/envs/moscot_rev_new

python $(dirname "$0")/run_stability_analysis.py \
    "$1" \
    "$2" \
    "$3" \
    "$4" \
