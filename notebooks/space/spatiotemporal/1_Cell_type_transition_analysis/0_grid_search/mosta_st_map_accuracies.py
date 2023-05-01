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


def set_transition_graph():
    transitions = pd.read_csv(DATA_DIR  + 'MOSTA_curated_transitions.csv', sep='\t')
    del transitions['Unnamed: 3']
    del transitions['/']

    G = nx.DiGraph()

    for n in transitions['Cell_type']:
        G.add_node(n)

    for i in range(len(transitions)):
        cell_type=transitions['Cell_type'][i]
        desc=transitions['Known_descendants'][i].split(', ')
        for des in desc:
            G.add_edge(cell_type, des)

        # Add transition to self
        G.add_edge(cell_type, cell_type)
    return G


def accuracy(df, G):
    correct = 0
    wrong = 0
    for cell_i in df.index:
        for cell_j in df.columns:
            if (cell_i, cell_j) in list(G.edges):
                correct += df[cell_j][cell_i]
            else:
                wrong +=  df[cell_j][cell_i]
                            
    return correct / (correct + wrong)


def solve(args):
    ## load data
    gamma = 10
    G = set_transition_graph()
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
    
    accuracies = {}
    accuracies_marg = {}
    probs = {}
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
        df_ = stp.cell_transition(
                source=start,
                target=end,
                source_groups="annotation",
                target_groups="annotation",
                forward=True,
            )
        accuracies_marg[key] = accuracy(df_, G)
    else:
        accuracies_marg[key] = -1
    res = accuracies_marg[key]
    print(f"for {key} accuracy is: {res}")
    fname = DATA_DIR + f"output/mouse_embryo_eps_{args.epsilon}_rank_{args.rank}_gamma_{gamma}_alpha_{args.alpha}_tp_{args.tps_couple}.csv"

    accuracies[args.alpha] = accuracies_marg
    pd.DataFrame(accuracies).to_csv(fname)

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
        
    
    
    
    
    