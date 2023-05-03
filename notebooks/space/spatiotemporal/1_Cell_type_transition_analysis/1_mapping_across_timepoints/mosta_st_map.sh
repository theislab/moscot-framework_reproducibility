#!/bin/bash

PROJECT_DIR='./'

alpha=$1
epsilon=$2
rank=$3
tps_couple=$4
save_res=$5

source

python3 ${PROJECT_DIR}/mosta_st_map_transitions.py --alpha ${alpha} --epsilon ${epsilon} --rank ${rank} --tps_couple ${tps_couple} --save_res ${save_res}