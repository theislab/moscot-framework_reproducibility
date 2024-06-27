import argparse
# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Description of your script')
# Define the command-line arguments
parser.add_argument('--param1', type=str, help='Description of param1')
parser.add_argument('--param2', type=str, help='Description of param2')
parser.add_argument('--param3', type=str, help='Decr')
# Parse the arguments from the command line
args = parser.parse_args()
# Access the parameter values
k = int(args.param1)
method = args.param2
tau=float(args.param3)

import os
os.environ["WANDB__SERVICE_WAIT"]="10000"

import numpy as np
import wandb
import scheduler as sch

print(method)
print(tau)
batch_size=2500
rs="0"
benchmark_result, ot_result=sch.benchmark(k=k, method=method, tau=tau, batch_size=batch_size, rs=rs)
print(len(ot_result))

max_memory=max(benchmark_result)

if method=='WOT':
    print('WOT')
    acc0, acc1, apoptosis_rate, t1, time1=ot_result
    iterations=np.NaN
    rank=-1
    
else:
    acc0, acc1, apoptosis_rate, t1, time1, iterations, batch_size, rs=ot_result

if 'LR' in method:
    gamma=500
    rank=2000
    tau=tau/10
else:
    gamma=np.NaN
    rank=-1

print(method)
print(k)
print(max_memory)

wandb.init(project=f"memory_runtime_benchmark_run_8")
wandb.log({"Method":method, "Accuracy_Curated": acc0, 'Accuracy_Germ':acc1, 'Iteration':iterations, 'Gamma':gamma, 'Rank':rank,'Max_memory':max_memory, 'Apoptosis_rate':apoptosis_rate, 'Sinkhorn_time':t1, 'Evaluation_time':time1, 'tau1':tau, 'k':k, "batch_size":batch_size, "seed":rs})
wandb.finish()


