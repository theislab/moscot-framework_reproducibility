from logging import getLogger
from pathlib import Path
from typing import Any, Dict, Literal, Optional
from pandas.api.types import is_numeric_dtype

import numpy as np
import seml
from sacred import Experiment
from sacred.run import Run

import sys

sys.path.insert(0, "../../../../")  
sys.path.insert(0, "../")  
#try:
#    import cr_utils as utils
#except RuntimeError as e:
#    if "jaxlib" not in str(e):
#        raise
        
logger = getLogger()

ex = Experiment()
seml.setup_logger(ex)


@ex.post_run_hook
def collect_stats(_run: Run) -> None:
    seml.collect_exp_stats(_run)


@ex.config
def config():
    overwrite = None
    db_collection = None
    if db_collection is not None:
        ex.observers.append(
            seml.create_mongodb_observer(db_collection, overwrite=overwrite)
        )


@ex.automain
def benchmark(
    epsilon: Literal[0.05,0.005],
    lambda_1: Literal[1,0.1],
    lambda_2: Literal[50],
    method: Literal['WOT', 'offline'],
    size: Literal[25,50,75,100],
) -> Dict[str, Any]:
    
    import scanpy as sc
    import time
    import numpy as np
    
    from typing import Callable
    from functools import wraps, partial

    from memory_profiler import memory_usage

    def benchmark_memory(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            return usage((func, args, kwargs))

        mp = False
        usage = partial(memory_usage, interval=0.01, include_children=mp, multiprocess=mp, retval=True)

        return wrapper
    
    
    if method=='WOT':   
        import wot
        import sklearn
        
        Path="/home/icb/manuel.gander/moscotTime_Reproducibility/Data"
        # get the anndata object
        adata=sc.read(f'{Path}/Subsampled/adata_{size}k.h5ad')
        
        A0=adata[adata.obs['day']==11.5]
        A1=adata[adata.obs['day']==12.5]

        cost_matrix = sklearn.metrics.pairwise.pairwise_distances(A0.obsm['X_pcaS'], A1.obsm['X_pcaS'],
                                                                      metric='sqeuclidean', n_jobs=-1)
        cost_matrix=cost_matrix.astype(np.float64)
        # matrix is reascaled in the ot-model with median, so rescaling is not needed at this step

        ot_model = wot.ot.OTModel(adata, day_field="day", epsilon=epsilon, lambda1=lambda_1, lambda2=lambda_2)

        starting_time=time.time()
        benchmarked = benchmark_memory(ot_model.compute_transport_map)
        benchmark_result, ot_result = benchmarked(11.5, 12.5, cost_matrix=cost_matrix)
        D={}
        D[f'{size}_time']=time.time()-starting_time
        D[size]=benchmark_result

        np.save(f'{Path}/Benchmark/{method}_{size}_CPU.npy', D, allow_pickle=True)
        
    elif method=='offline':
    
        # Setting jax to use float64, makes the algorithm more stable
        from jax.config import config
        config.update("jax_enable_x64", True)

        #import moscot
        from moscot.problems.time._lineage import TemporalProblem
        from moscot.backends.ott._solver import SinkhornSolver

        import anndata

        tau1=lambda_1/(lambda_1+epsilon)
        tau2=lambda_2/(lambda_2+epsilon)

        Path="/home/icb/manuel.gander/moscotTime_Reproducibility/Data"

        # get the anndata object
        adata=sc.read(f'{Path}/Subsampled/adata_{size}k.h5ad')

        tp=TemporalProblem(adata)
        tp.score_genes_for_marginals(gene_set_proliferation='mouse',  gene_set_apoptosis='mouse')
        tp = tp.prepare('day', joint_attr='X_pcaS')

        starting_time=time.time()
        benchmarked = benchmark_memory(tp.solve)
        benchmark_result, ot_result = benchmarked(batch_size=3*10**6, epsilon=epsilon, tau_a=tau1, tau_b=tau2, scale_cost="mean",  max_iterations=10**5)
        
        D={}
        D[f'{size}_time']=time.time()-starting_time

        D[size]=benchmark_result

        np.save(f'{Path}/Benchmark/{method}_{size}_CPU.npy', D, allow_pickle=True)

    #return D