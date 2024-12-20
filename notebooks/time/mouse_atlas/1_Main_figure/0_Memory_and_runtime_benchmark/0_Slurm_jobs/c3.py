import wot
import numpy as np
import anndata
import pandas as pd
import scanpy as sc
import time
from moscot.problems.time._lineage import TemporalProblem
import sys
sys.path.append('/home/icb/manuel.gander/mouse_atlas/notebook')
import c2





def solve_WOT(k, tau1=0.8, rs="0"):
    Path="/home/icb/manuel.gander/mouse_atlas/data/Sub2" 
    adata=sc.read(f"{Path}/adata_{k}k_{rs}.h5ad")
    del adata.raw
    # Compute intiial growth rates the same as in moscot
    adata.obs['day']=adata.obs['day'].astype('category') 
    adata.obs['cell_growth_rate']=compute_initial_growth_rates(adata)
    adata=anndata.AnnData(X=adata.obsm['X_pcaS'], obs=adata.obs)
    day0,day1=sorted(set(adata.obs['day']))
    eps=0.005
    #eps=1
    #tau1=0.8
    lam1=tau1*eps/(1-tau1)
    tau2=0.99995
    lam2=tau2*eps/(1-tau2)
    
    # compute WOT mapping
    ot_model = wot.ot.OTModel(adata, epsilon = eps, lambda1 = lam1, lambda2 = lam2, growth_iters=1)
    print('Running WOT')
    t0=time.time()
    tmap = ot_model.compute_transport_map(day0, day1)
    t1=time.time()-t0

    # evaluate
    time0=time.time()
    tmap=tmap.X
    tmap=tmap/tmap.sum()
    gr=tmap.sum(1)
    gr=gr/gr.mean()*2600000/1100000
    cell_dying=np.sum((1-gr[gr<1]))
    apoptosis_rate=float(cell_dying/len(gr))
    CT=compute_cell_type_transitions(tmap, adata, day0, day1)
    ev=c2.evaluate_using_curated_transitions(CT)
    acc0=list(ev['Accuracy'])[0]
    ev=c2.evaluate_using_germ_layers(CT)
    acc1=list(ev['Accuracy'])[0]
    time1=time.time()-time0
   
    return(acc0, acc1, apoptosis_rate, t1/3600, time1/3600)


def compute_initial_growth_rates(adata):
    # Get prior growth rates
    tp=TemporalProblem(adata)
    tp=tp.score_genes_for_marginals(gene_set_proliferation='mouse', gene_set_apoptosis='mouse')
    tp = tp.prepare('day', joint_attr=f'X_pcaS') 
    gr=tp.prior_growth_rates.loc[adata.obs.index]['prior_growth_rates'].values
    return(gr)


def compute_cell_type_transitions(tmap, adata, day0, day1):
    a0=adata[adata.obs['day']==day0].copy()
    a1=adata[adata.obs['day']==day1].copy()
    cs0=np.array(a0.obs['cell_state'])
    cs1=np.array(a1.obs['cell_state'])
    cellstates0=sorted(set(cs0))
    cellstates1=sorted(set(cs1))

    CT=np.zeros((len(cellstates0), len(cellstates1)))
    for j,c1 in enumerate(cellstates1):
        wh1=np.where(cs1==c1)[0]
        tmap0=tmap[:,wh1]
        for i,c0 in enumerate(cellstates0):
            wh0=np.where(cs0==c0)[0]
            CT[i,j]=tmap0[wh0,:].sum()
    CT=pd.DataFrame(data=CT, index=cellstates0, columns=cellstates1)
    return(CT)
