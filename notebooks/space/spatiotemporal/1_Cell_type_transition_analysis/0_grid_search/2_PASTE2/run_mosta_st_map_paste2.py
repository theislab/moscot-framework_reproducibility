from subprocess import Popen



DATA_DIR = ''
random_states = [23860, 50349, 36489, 59707, 38128, 25295, 49142, 12102, 30139, 4698]
NUM_TPS = 7

for i in range(NUM_TPS):
    for random_state in random_states:
        cmdline0 = ['sbatch','--gres=gpu:a5000:1', '-c1', '--time=7-0', '--mem=20gb',
                    f'--output={DATA_DIR}/output_paste2_seeds/mosta-st-{i}-{random_state}.log',
                    f'--job-name=mosta-paste2-{i}-{random_state}',
                    'mosta_st_map.sh', str(i), str(random_state)
                    ]
        
        print(' '.join(cmdline0))
        Popen(cmdline0)
