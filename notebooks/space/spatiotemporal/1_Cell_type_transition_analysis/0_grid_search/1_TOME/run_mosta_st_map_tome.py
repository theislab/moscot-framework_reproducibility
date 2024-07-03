from subprocess import Popen


DATA_DIR = ""

NUM_TPS = 7

for i in range(NUM_TPS):
    cmdline0 = ['sbatch', '-c1', '--time=7-0', '--mem=20gb',
                f'--output={DATA_DIR}/output_tome/mosta-st-{i}.log',
                f'--job-name=mosta-tome-{i}',
                'mosta_st_map.sh', str(i)
                ]
    
    print(' '.join(cmdline0))
    Popen(cmdline0)
