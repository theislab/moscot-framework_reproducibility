#!/bin/bash
methods=("WOT" "moscot_CPU" "moscot_GPU" "moscot_LR_CPU" "moscot_LR_GPU")
methods=("WOT")
methods=("moscot_CPU")
ks=(25 50 75 100 125 150 175 200 225 250 275)

taus=(0.9)

for tau in "${taus[@]}"; do
    for k in "${ks[@]}"; do
        for method in "${methods[@]}"; do
            echo "$k"
            echo "$rank"
            echo "$tau"
            if [[ $method == *"GPU"* ]]; then
                echo 'GPU'
                sbatch gpu_sub.sh "$k" "$method" "$tau"
            else
                echo 'CPU'
                sbatch cpu_sub.sh "$k" "$method" "$tau"
            fi
        done
    done
done
