
import sys

import anndata
import argparse
import numpy as np
from pandas import DataFrame
from sklearn.neighbors import NearestNeighbors
import scanpy as sc
import pandas as pd
from typing import Tuple, Literal
from numpy.typing import ArrayLike

import networkx as nx

from paste2 import PASTE2, model_selection, projection

DATA_DIR =  ""

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
                wrong += df[cell_j][cell_i]

    return correct / (correct + wrong)

def paste2_agg(
        adata: anndata.AnnData,
        use_rep: str = None,
        s: float = 0.8,
        max_cells: int = 2000,
        nb_dummies: int = 1,
        dissimilarity: Literal["glmpca", "kl", "euclidean"] = "glmpca",
        random_state: int = 0,
        return_adata: bool = False
) -> Tuple:
    """
    finding ancestor node for each node
    :param adata: the input adata (n_cells X n_features)
    :param use_rep: If none, uses slice.X to calculate dissimilarity between spots, otherwise uses the representation given by slice.obsm[use_rep]
    :param s: overlap fraction (float)
    :param max_cells: number of max cells (int)
    :param nb_dummies: int, optional, default:1
        number of reservoir points to be added (to avoid numerical
        instabilities, increase its value if an error is raised)
    :param dissimilarity: Expression dissimilarity measure: 'kl' or 'euclidean' or 'glmpca'. Default is glmpca.
    """

    n_cells, _ = adata.shape
    source, target = np.unique(adata.obs["time"])
    adata_1 = adata[adata.obs["time"].isin([source])].copy()
    adata_2 = adata[adata.obs["time"].isin([target])].copy()

    if adata_1.shape[0] > max_cells:
        sc.pp.subsample(adata_1, n_obs=max_cells, random_state=random_state)
    if adata_2.shape[0] > max_cells:
        sc.pp.subsample(adata_2, n_obs=max_cells, random_state=random_state)

    # keep stated consistent across reps
    pi = PASTE2.partial_pairwise_align(adata_1, adata_2, s=s, use_rep=use_rep, nb_dummies=nb_dummies, dissimilarity=dissimilarity)
    df_res = pd.DataFrame(pi, index=adata_1.obs["annotation"].values, columns=adata_2.obs["annotation"].values)
    df_res = df_res.reset_index().groupby("index").sum()
    df_res = df_res.T.reset_index().groupby("index").sum()

    if return_adata:
        return df_res, pi, adata_1, adata_2
    
    return df_res, pi


def solve(args):
    ## load data
    G = set_transition_graph()
    tps_couple = [
        [9.5, 10.5],
        [10.5, 11.5],
        [11.5, 12.5],
        [12.5, 13.5],
        [13.5, 14.5],
        [14.5, 15.5],
        [15.5, 16.5]
    ]
    tps_c = tps_couple[args.tps_couple]
    start, end = tps_c
    
    start_str = str(start).replace(".", "_")
    end_str = str(end).replace(".", "_")
    print(f"Loading `mouse_embryo_{start_str}_{end_str}_renormalized.h5ad`")
    adata_c = sc.read_h5ad(DATA_DIR + f"mouse_embryo_{start_str}_{end_str}_renormalized.h5ad")

    res, _ = paste2_agg(
        adata=adata_c,
        s=args.s,
        max_cells=args.max_cells,
        random_state=args.random_state
    )

    res.divide(res.sum(axis=0), axis="columns")
    res = res.fillna(0.0)
    accuracy_val = accuracy(res, G)
    print(f"for {tps_c} accuracy is: {accuracy_val}")
    fname = DATA_DIR + f"output_paste2/mouse_embryo_tp_{args.tps_couple}_{args.random_state}.csv"
    pd.DataFrame({"tp": [args.tps_couple], "start": [start], "end": [end], "accuracy": [accuracy_val]}).to_csv(fname)
    return 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tps_couple',
                        type=int,
                        required=False,
                        default=0)
    parser.add_argument('--max_cells',
                        type=int,
                        required=False,
                        default=2000)
    parser.add_argument('--random_state',
                        type=int,
                        required=False,
                        default=0)
    parser.add_argument('--s',
                        type=float,
                        required=False,
                        default=0.8)
    args = parser.parse_args()

    solve(args)