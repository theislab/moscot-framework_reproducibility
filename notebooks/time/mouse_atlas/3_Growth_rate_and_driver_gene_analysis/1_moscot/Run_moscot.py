import numpy as np
import anndata
import pandas as pd
import scanpy as sc
import scipy
from moscot.problems.time._lineage import TemporalProblem

Path="/home/mgander/moscot-framework_reproducibility/data/mouse_temporal"
ts=['E3.5', 'E4.5', 'E5.25', 'E5.5', 'E6.25', 'E6.5', 'E6.75', 'E7.0', 'E7.25', 'E7.5', 'E7.75', 'E8.0', 'E8.25', 'E8.5a', 'E8.5b', 'E9.5', 'E10.5', 'E11.5', 'E12.5', 'E13.5']


i=17
ts0=ts[i]
ts1=ts[i+1]
print(f'{ts0}_{ts1}')
print('------------------------')

adata=sc.read(f"{Path}/anndatas/Comb_anndatas/adata_{ts0}_{ts1}.h5ad")
# Raw count matrix not needed here, but causes problems in "score_genes_for_marginals"
del adata.raw
    
#sc.pp.subsample(adata, n_obs=150000, random_state=0)

tp=TemporalProblem(adata)
tp.score_genes_for_marginals(gene_set_proliferation='mouse',  gene_set_apoptosis='mouse')
tp = tp.prepare('day', joint_attr=f'X_pcaS')

batch_size=25*10**2
eps=0.005
tau1=0.87
tau2=1

result=tp.solve(batch_size=batch_size, epsilon=eps, tau_a=tau1, tau_b=tau2, scale_cost="mean", max_iterations=10**5)

tp.save(f'{Path}/moscot_maps/', f'{ts0}_{ts1}', overwrite=True)