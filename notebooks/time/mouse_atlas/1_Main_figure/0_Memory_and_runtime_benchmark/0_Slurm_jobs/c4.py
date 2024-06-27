import jax
jax.config.update('jax_platform_name', 'cpu')
#from jax.config import config
#config.update("jax_enable_x64", True)
import numpy as np
import anndata
import pandas as pd
import scanpy as sc
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
import os
import time
#os.environ["JAX_PLATFORM_NAME"] = "cpu"
from moscot.problems.time._lineage import TemporalProblem
import warnings
import jax
import jax.numpy as jnp
from ott.geometry import pointcloud
from ott.problems.linear import linear_problem
from ott.solvers.linear import sinkhorn, sinkhorn_lr
import sys,os
sys.path.append('/home/icb/manuel.gander/mouse_atlas/notebook')
import scripts as scr
import c2
warnings.simplefilter(action='ignore', category=FutureWarning)
sc.settings.verbosity = 0
from memory_profiler import profile



def solve_moscot(k=25, rank=-1, tau1=0.9, batch_size=1000, rs="0"):

    if rank==-1:
        epsilon=0.005
        #tau1=0.8
    else:
        epsilon=0.0001
        gamma=500
        iterations=1000
        tau1=tau1/10
    tau2=0.99995

    Path="/home/icb/manuel.gander/mouse_atlas/data/Sub2"
    print(k)
    adata=sc.read(f"{Path}/adata_{k}k_{rs}.h5ad")
    del adata.raw
    adata.obs['day']=adata.obs['day'].astype('category')

    day0,day1=sorted(set(adata.obs['day']))
    inds0=list(adata[adata.obs['day']==day0].obs.index)
    inds1=list(adata[adata.obs['day']==day1].obs.index)

    tp=TemporalProblem(adata)
    tp.score_genes_for_marginals(gene_set_proliferation='mouse',  gene_set_apoptosis='mouse')
    tp = tp.prepare('day', joint_attr=f'X_pcaS')
 
    t0=time.time()
    if rank==-1:
        result=tp.solve(batch_size=batch_size, epsilon=epsilon, tau_a=tau1, tau_b=tau2, rank=rank)
        iterations=float(tp[(day0, day1)].solution._output.n_iters)
    else:
        inners=max(int(iterations/20),1)
        result=tp.solve(batch_size=batch_size, epsilon=epsilon, tau_a=tau1, tau_b=tau2, max_iterations=iterations, rank=rank, threshold=0, inner_iterations=inners, gamma=gamma)
    t1=time.time()-t0
    print('Sinkhorn done')
    
    time0=time.time()
    gr=tp[(day0, day1)].solution.a
    gr=gr/gr.mean()*2600000/1100000
    cell_dying=np.sum((1-gr[gr<1]))
    apoptosis_rate=float(cell_dying/len(gr))

    cell_states0={'cell_state': list(set(adata[adata.obs['day']==day0].obs['cell_state']))}
    cell_states1={'cell_state': list(set(adata[adata.obs['day']==day1].obs['cell_state']))}
    CT=tp.cell_transition(day0, day1, cell_states0, cell_states1, batch_size=batch_size)

    ev=c2.evaluate_using_curated_transitions(CT)
    acc0=list(ev['Accuracy'])[0]

    ev=c2.evaluate_using_germ_layers(CT)
    acc1=list(ev['Accuracy'])[0]
    time1=time.time()-time0
    #acc0, acc1, apoptosis_rate, time1, iteration=0,0,0,0,0
    return(acc0, acc1, apoptosis_rate, t1/3600, time1/3600, iterations, batch_size, rs)
