from subprocess import Popen
import os
import sys


DATA_DIR = './'


alphas = [0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.99]
epsilon = 1e-3
ranks = [500]

NUM_TPS = 7

for i in range(NUM_TPS):
    for rank in ranks:
        for alpha in alphas:
            cmdline0 = ['sbatch', '--gres=gpu:a5000:1', '-c1', '--time=1-0', '--mem=300gb',
                        f'--output={DATA_DIR}/output/mosta-st-{i}-{alpha}-{epsilon}_{rank}_{alpha}.log',
                        f'--job-name=mosta-st-{i}-{alpha}-{epsilon}-{rank}-{alpha}',
                        'mosta_st_map.sh', str(alpha), str(epsilon), str(rank), str(i)
                        ]

            print(' '.join(cmdline0))
            Popen(cmdline0)
