from subprocess import Popen


DATA_DIR = './'


alpha = 0.5
epsilon = 1e-3
rank_high = 10000
rank_mid = 5000

NUM_TPS = 7
tps_alphas = ["0.4", "0.99", "0.4", "0.8", "0.8", "0.6", "0.5"]

for i in range(NUM_TPS):
    alpha = tps_alphas[i]
    rank = rank_high if i < 6 else rank_mid
    cmdline0 = ['sbatch', '--gres=gpu:a5000:1', '-c1', '--time=1-0', '--mem=100gb',
                f'--output={DATA_DIR}/output_brain/mosta-st-{i}-{alpha}-{epsilon}_{rank}_{alpha}_100.log',
                f'--job-name=mosta-st-{i}-{alpha}-{epsilon}-{rank}-{alpha}',
                'mosta_st_map_brain.sh', str(alpha), str(epsilon), str(rank), str(i)
               ]

    print(' '.join(cmdline0))
    Popen(cmdline0)
