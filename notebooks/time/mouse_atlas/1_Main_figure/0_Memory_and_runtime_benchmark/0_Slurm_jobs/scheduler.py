import sys
sys.path.append('/home/icb/manuel.gander')
import c3
import c4

from typing import Callable
from functools import wraps, partial
from memory_profiler import memory_usage

def benchmark_memory(func: Callable) -> Callable:
    @wraps(func)
    def wrapper(*args, **kwargs):
        return usage((func, args, kwargs))
    usage = partial(memory_usage, interval=0.01, include_children=False, multiprocess=False, retval=True)
    return wrapper


def benchmark(k=25, method='moscot', tau=0.8, batch_size=1000, rs="0"):
    if method=='WOT':
        benchmarked = benchmark_memory(c3.solve_WOT)
        benchmark_result, ot_result=benchmarked(k=k, tau1=tau, rs=rs)
    else:
        benchmarked = benchmark_memory(c4.solve_moscot)
        if 'LR' in method:
             benchmark_result, ot_result=benchmarked(k=k, rank=2000, tau1=tau, batch_size=batch_size, rs=rs)
        else:
             benchmark_result, ot_result=benchmarked(k=k, rank=-1, tau1=tau, batch_size=batch_size, rs=rs)
    return(benchmark_result, ot_result)


