#!/bin/bash

PROJECT_DIR=""

tps_couple=$1
random_state=$2

source ""

python3 ${PROJECT_DIR}/run_paste2.py --tps_couple ${tps_couple} --random_state ${random_state}