import sys
import argparse

import numpy as np
import scanpy as sc

from jax.config import config
config.update("jax_enable_x64", True)

from moscot.problems.time import TemporalProblem
from moscot.problems.spatiotemporal import SpatioTemporalProblem

from scipy.sparse import csr_matrix
import pickle

DATA_DIR = "./"


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
    adata = sc.read(DATA_DIR + "mouse_embryo_brain_normalized.h5ad")
    tps_numeric = adata.obs["time"].unique()
    tps_couple = [[i, i+1] for i in tps_numeric[:-1]]
    tps_c = tps_couple[args.tps_couple]
    start, end = tps_c
    
    adata_c = adata[adata.obs["time"].isin(tps_c)].copy()
    
    adata_c.obs["transition"] = adata_c.obs_names
    if end == 16.5:
        adata_c.obs.loc[adata_c.obs["time"].isin([end]), "transition"] = adata_c[adata_c.obs["time"].isin([end])].obs["annotation"]
    
    adata_c.obs["transition"] = adata_c.obs["transition"].astype("category") 
    
    gamma = 100
    scale_cost = "max_cost"
    
    if args.alpha == 0:
        stp = TemporalProblem(adata=adata_c).prepare(
            time_key="time",
            epsilon=args.epsilon)
        stp = stp.score_genes_for_marginals(gene_set_proliferation="mouse", gene_set_apoptosis="mouse")
        stp = stp.solve(
            epsilon=args.epsilon, 
            rank=args.rank,
            gamma=gamma,
            scale_cost=scale_cost,
        )
    else:
        stp = SpatioTemporalProblem(adata=adata_c).prepare(
            time_key="time",
            spatial_key="spatial",
        )
        
        stp = stp.score_genes_for_marginals(gene_set_proliferation="mouse", gene_set_apoptosis="mouse")
        stp = stp.solve(
            alpha=args.alpha, 
            epsilon=args.epsilon,
            rank=args.rank, 
            gamma=gamma,
            scale_cost=scale_cost
        )
    
    if stp.solutions[list(stp.solutions.keys())[0]].converged:
        df_tmat_forward = stp.cell_transition(
            source=start,
            target=end,
            source_groups="transition",
            target_groups="transition",
            forward=True
        )
        
        df_tmat_back = stp.cell_transition(
            source=start,
            target=end,
            source_groups="transition",
            target_groups="transition",
            forward=False
        )
        
        df_tmat_sparsify = sparsify(stp=stp[(start, end)].solution, threshold="Auto")
        
        df_tmat_forward.to_csv(DATA_DIR + f"output_brain/brain_eps_{args.epsilon}_rank_{args.rank}_gamma_{gamma}_alpha_{args.alpha}_tp_{args.tps_couple}_forward.csv")
        df_tmat_back.to_csv(DATA_DIR + f"output_brain/brain_eps_{args.epsilon}_rank_{args.rank}_gamma_{gamma}_alpha_{args.alpha}_tp_{args.tps_couple}_back.csv")
        fname = DATA_DIR +f"output_brain/brain_eps_{args.epsilon}_rank_{args.rank}_gamma_{gamma}_alpha_{args.alpha}_tp_{args.tps_couple}_sparse.pkl"
        with open(fname, "wb") as handle:
            pickle.dump(df_tmat_sparsify, handle)
    else:
        print(f"eps_{args.epsilon}_rank_{args.rank}_gamma_{gamma}_alpha_{args.alpha}_tp_{args.tps_couple} didn't converge")

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
    args = parser.parse_args()

    solve(args)
        
    
    
    
    
    