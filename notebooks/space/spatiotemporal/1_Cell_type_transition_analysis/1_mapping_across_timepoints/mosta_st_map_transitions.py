from os.path import exists
import sys
import numpy as np
import argparse
import pickle
from typing import List


import networkx as nx

import scanpy as sc
import pandas as pd

from jax.config import config
config.update("jax_enable_x64", True)


from moscot.problems.time import TemporalProblem
from moscot.problems.spatiotemporal import SpatioTemporalProblem

from scipy.sparse import hstack as sparse_hstack, csr_matrix

DATA_DIR = './'


def sparsify(stp, threshold = None, batch_size: int = 1024) -> csr_matrix:
    tmaps = []
    for batch in range(0, stp.shape[1], batch_size):
        x = np.eye(stp.shape[1], min(batch_size, stp.shape[1] - batch), -(min(batch, stp.shape[1])))
        res = stp.pull(x, scale_by_marginals=False)  # tmap @ indicator_vectors
        tmaps.append(res)
    tmat = np.hstack(tmaps)
    if threshold is not None:
        threshold = min(np.max(tmat_) for tmat_ in tmat) if threshold == "Auto" else threshold
        tmat[tmat < threshold] = 0
        return csr_matrix(tmat)


def solve(args):
    ## load data
    gamma = 10
    adata = sc.read(DATA_DIR + 'mouse_embryo_all_stage_renormalized.h5ad')
    
    adata.obs["time"] = adata.obs["timepoint"].copy()
    adata.obs["time"] = adata.obs["time"].cat.rename_categories(lambda x: x.split("E")[1]).astype(float)
    
    tps_numeric = adata.obs["time"].unique()
    tps_couple = [[i, i+1] for i in tps_numeric[:-1]]
    tps_c  = tps_couple[args.tps_couple]
    start, end = tps_c
    
    adata_c = adata[adata.obs["time"].isin(tps_c)].copy()
    sc.pp.pca(adata_c, use_highly_variable=True)
    adata_c.obsm["X_pca_30"] = adata_c.obsm["X_pca"][:, :30]
    
    adata_c.obs["transition"] = adata_c.obs_names
    adata_c.obs.loc[adata_c.obs["time"].isin([start]), "transition"] = \
    adata_c[adata_c.obs["time"].isin([start])].obs["annotation"]
    adata_c.obs["transition"] = adata_c.obs["transition"].astype("category")

    scale_cost= "max_cost"
    print(f"eval {args.alpha}")

    if args.alpha == 0:
        stp = TemporalProblem(adata=adata_c)
        stp = stp.score_genes_for_marginals(gene_set_proliferation="mouse", gene_set_apoptosis="mouse")
        
        stp = stp.prepare(
            time_key="time",
            epsilon=args.epsilon,
            joint_attr="X_pca_30"
        )
        
        stp = stp.solve(
                epsilon=args.epsilon, 
                rank=args.rank,
                gamma=gamma,
                scale_cost=scale_cost,
                batch_size=1024
        )
        
    else:
        stp = SpatioTemporalProblem(adata=adata_c)
        stp = stp.score_genes_for_marginals(gene_set_proliferation="mouse", gene_set_apoptosis="mouse")
            
        stp = stp.prepare(
                time_key="time",
                spatial_key="spatial",
                joint_attr="X_pca_30",
            )
            
        stp = stp.solve(
                alpha=args.alpha, 
                epsilon=args.epsilon, 
                rank=args.rank, 
                gamma=gamma,
                scale_cost=scale_cost
            )
    if stp.solutions[list(stp.solutions.keys())[0]].converged:
        
        fname_stp = DATA_DIR + f"output_res/mouse_embryo_eps_{args.epsilon}_rank_{args.rank}_gamma_{gamma}_alpha_{args.alpha}_tp_{args.tps_couple}_stp.pkl"
        with open(fname_stp, "wb") as handle:
            pickle.dump(stp.solutions[list(stp.solutions.keys())[0]], handle)
            
        if args.save_res:

            df_push = stp.push(
                    source=start,
                    target=end,
                    data="annotation",
                    subset="Heart",
                    return_data=True
                )
            fname_push = DATA_DIR + f"output_res/mouse_embryo_eps_{args.epsilon}_rank_{args.rank}_gamma_{gamma}_alpha_{args.alpha}_tp_{args.tps_couple}_heart_push.pkl"
            pd.DataFrame(df_push[end], index=stp[start, end].adata_tgt.obs_names).to_pickle(fname_push)

            tmat = sparsify(stp=stp.solutions[list(stp.solutions.keys())[0]], threshold="Auto")
            fname_tmat = DATA_DIR + f"output_res/mouse_embryo_eps_{args.epsilon}_rank_{args.rank}_gamma_{gamma}_alpha_{args.alpha}_tp_{args.tps_couple}_stp.pkl"
            with open(fname_tmat, "wb") as handle:
                pickle.dump(tmat, handle)

    return
   
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--alpha',
                        type=float,
                        required=False,
                        default=0.5)
    parser.add_argument('--epsilon',
                        type=float,
                        required=False,
                        default=0.001)
    parser.add_argument('--rank',
                        type=int,
                        required=False,
                        default=None)
    parser.add_argument('--tps_couple',
                        type=int,
                        required=False,
                        default=0)
    parser.add_argument('--save_res',
                        type=int,
                        required=False,
                        default=0)
    args = parser.parse_args()

    solve(args)
        
    
    
    
    
    