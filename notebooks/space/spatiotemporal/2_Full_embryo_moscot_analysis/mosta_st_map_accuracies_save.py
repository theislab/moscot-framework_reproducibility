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


DATA_DIR = ""


import ott
print(ott.__file__)

def set_transition_graph():
    transitions = pd.read_csv(DATA_DIR + 'MOSTA_curated_transitions.csv', sep='\t')
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
    gamma = args.gamma
    tau_a = args.tau_a
    tau_b = args.tau_b
    G = set_transition_graph()
    
    tps_couple = [
        [9.5, 10.5],
        [10.5, 11.5],
        [11.5, 12.5],
        [12.5, 13.5],
        [13.5, 14.5],
        [14.5, 15.5],
        [15.5, 16.5],
        [9.5, 11.5]
    ]
    tps_c = tps_couple[args.tps_couple]
    start, end = tps_c
    
    start_str = str(start).replace(".", "_")
    end_str = str(end).replace(".", "_")
    print(f"Loading `mouse_embryo_{start_str}_{end_str}_renormalized.h5ad`")
    adata_c = sc.read_h5ad(DATA_DIR + f"mouse_embryo_{start_str}_{end_str}_renormalized.h5ad")

    scale_cost = "max_cost"
    print(f"eval {args.alpha}")

    accuracies_taus = []
    min_iterations = 500
    max_iterations = 10_000
    inner_iterations = 2_000
    if args.alpha == 0:
        stp = TemporalProblem(adata=adata_c)
        stp = stp.score_genes_for_marginals(gene_set_proliferation="mouse", gene_set_apoptosis="mouse")
        stp = stp.prepare(
                        time_key="time",
                        epsilon=args.epsilon,
                        joint_attr="X_pca_30",
                        cost=args.cost,
                    )
        stp = stp.solve(
                        epsilon=args.epsilon, 
                        rank=args.rank,
                        gamma=gamma,
                        scale_cost=scale_cost,
                        tau_a=tau_a,
                        tau_b=tau_b,
                        min_iterations=min_iterations,
                        max_iterations=max_iterations,
                        inner_iterations=inner_iterations,
                        batch_size=1024
                    )
            
    else:
        stp = SpatioTemporalProblem(adata=adata_c)
        stp = stp.score_genes_for_marginals(gene_set_proliferation="mouse", gene_set_apoptosis="mouse")
                        
        stp = stp.prepare(
                        time_key="time",
                        spatial_key="spatial",
                        joint_attr="X_pca_30",
                        cost=args.cost,
                    )
                        
        stp = stp.solve(
                        alpha=args.alpha, 
                        epsilon=args.epsilon, 
                        rank=args.rank, 
                        gamma=gamma,
                        scale_cost=scale_cost,
                        tau_a=tau_a,
                        tau_b=tau_b,
                        min_iterations=min_iterations,
                        max_iterations=max_iterations,
                    )
    df_ = stp.cell_transition(
                            source=start,
                            target=end,
                            source_groups="annotation",
                            target_groups="annotation",
                            forward=True,
                        )
    df_.divide(df_.sum(axis=0), axis="columns")
    df_ = df_.fillna(0.0)
    accuracy_val = accuracy(df_, G)
                
    converged = stp.solutions[list(stp.solutions.keys())[0]].converged
    costs = stp.solutions[list(stp.solutions.keys())[0]]._costs
    n_iters = len(np.where(costs != -1)[0])
    costs_last = str(costs[np.where(costs != -1)[0]][-2]) + "_" +str(costs[np.where(costs != -1)[0]][-1])
    accuracies_taus.append(pd.DataFrame({
                    "tau_a": [tau_a], 
                    "tau_b": [tau_b], 
                    "accuracy": [accuracy_val], 
                    "scale_cost": [scale_cost], 
                    "converged":[converged], 
                    "costs_last": [costs_last],
                    "n_iters":[n_iters]
                })) 
                
    print(f"for {tau_a}-{tau_b} accuracy={accuracy_val}. n_iters={n_iters}, converged={converged}, costs_last={costs_last}")
    fname_acc = DATA_DIR + f"output_liver/mouse_embryo_eps_{args.epsilon}_rank_{args.rank}_gamma_{gamma}_alpha_{args.alpha}_tp_{start_str}_{end_str}_cost_{args.cost}.csv"
    fname_liver = DATA_DIR + f"output_liver/liver_pull_eps_{args.epsilon}_rank_{args.rank}_gamma_{gamma}_alpha_{args.alpha}_tp_{start_str}_{end_str}_cost_{args.cost}.csv"
    fname_tfs = DATA_DIR + f"output_liver/tfs_push_eps_{args.epsilon}_rank_{args.rank}_gamma_{gamma}_alpha_{args.alpha}_tp_{start_str}_{end_str}_cost_{args.cost}.csv"
    accuracies_df = pd.concat(accuracies_taus).reset_index().drop(columns=["index"])
    accuracies_df.to_csv(fname_acc)

    stp.save(DATA_DIR + f"output_liver/stp_eps_{args.epsilon}_rank_{args.rank}_gamma_{gamma}_alpha_{args.alpha}_tp_{start_str}_{end_str}_cost_{args.cost}.pkl", overwrite=True)
    ct = "Liver"

    stp.pull(source=start, target=end, data="annotation", subset=ct, key_added=f"{ct}_pull", normalize=False)
    df_liver = stp.compute_feature_correlation(
        obs_key=f"{ct}_pull",
        annotation={"time": [start, end]}
    )
    df_liver.to_csv(fname_liver)

    tfs = pd.read_csv(DATA_DIR + "Mus_musculus_TF.txt", delimiter="\t")
    stp.adata.var["TF"] = stp.adata.var_names.isin(tfs["Symbol"])

    ## push & correlate tfs
    tfs_dfs = []
    for tf in stp.adata.var_names[stp.adata.var["TF"]]:
        stp.adata.obs[tf] = stp.adata[:, tf].X.A.copy()
        stp.push(source=start, target=end, data=tf, key_added=f"{tf}_push", normalize=False)
        del stp.adata.obs[tf]
        tfs_dfs.append(
            stp.compute_feature_correlation(obs_key=f"{tf}_push", annotation={"time": [start, end]}))
        del stp.adata.obs[f"{tf}_push"]

    tfs_df = pd.concat(tfs_dfs, axis=1)
    tfs_df.to_csv(fname_tfs)

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
                        default=-1)
    parser.add_argument('--tps_couple',
                        type=int,
                        required=False,
                        default=0)
    parser.add_argument('--gamma',
                        type=float,
                        required=False,
                        default=None)
    parser.add_argument('--cost',
                        type=str,
                        required=False,
                        default="sq_euclidean")
    parser.add_argument('--tau_a',
                        type=float,
                        required=False,
                        default=1.0)
    parser.add_argument('--tau_b',
                        type=float,
                        required=False,
                        default=1.0)
    args = parser.parse_args()

    solve(args)
        
    
    
    
    
    