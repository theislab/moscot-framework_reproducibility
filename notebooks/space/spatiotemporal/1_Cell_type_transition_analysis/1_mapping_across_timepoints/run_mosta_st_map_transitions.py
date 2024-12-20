from subprocess import Popen
import os
import sys


DATA_DIR = '/cs/labs/mornitzan/zoe.piran/research/projects/moscot_framework_analysis/data/mosta/'


alphas = [0.4, 0.99, 0.4, 0.8, 0.8, 0.6, 0.99]
epsilon = 1e-3
rank = 500 
save_res = 0

NUM_TPS = 7

for i in range(NUM_TPS):
    cmdline0 = ['sbatch', '--gres=gpu:a5000:1', '-c1', '--time=1-0', '--mem=300gb',
                f'--output={DATA_DIR}/output_res/mosta-full-{i}-{alphas[i]}-{epsilon}_{rank}.log',
                f'--job-name=mosta-st-{i}-{alphas[i]}-{epsilon}-{rank}',
                'mosta_st_map.sh', str(alphas[i]), str(epsilon), str(rank), str(i), str(save_res)
               ]
    print(' '.join(cmdline0))
    Popen(cmdline0)
