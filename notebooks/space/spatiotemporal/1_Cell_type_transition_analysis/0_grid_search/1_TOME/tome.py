import os
import sys

import argparse
import pickle

import numpy as np
from sklearn.neighbors import NearestNeighbors
import scanpy as sc
import pandas as pd
from typing import Tuple
from numpy.typing import ArrayLike

import networkx as nx


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

def createlineage_knn(
        emb: ArrayLike,
        time: ArrayLike,
        annotation: ArrayLike,
        n_reps: int=500,
        ratio_remove_cells: float = 0.2,
        n_neighbors: int = 5
) -> Tuple:
    """
    finding ancestor node for each node
    :param emb: the embedding for k-nn computation (n_cells X n_emb)
    :param time: an array of time point per-cell (n_cells)
    :param annotation: an array of cell state annotations per-cell (n_cells)
    :param n_reps: number of reps (int)
    :param ratio_remove_cells: fraction of cells to remove in each rep (float)
    :param n_neighbors: number of neighbors (int)
    :return: List of n=reps mappings from target to source states.
    """
    res = []

    states = dict()

    n_cells, _ = emb.shape
    source, target = np.unique(time)
    ratio_keep_cells = 1.0-ratio_remove_cells

    # keep stated consistent across reps
    for tp in [source, target]:
        states[tp] = np.unique(annotation[time == tp])

    res_full = np.zeros((n_reps, len(states[target]), len(states[source])))

    for i in range(n_reps):

        sampling_idx = np.random.choice(n_cells, int(n_cells*ratio_keep_cells), replace=False)
        time_rep = time[sampling_idx]

        embs = dict()
        annts = dict()
        for tp in [source, target]:
            embs[tp] = emb[sampling_idx, :][time_rep==tp, :]
            annts[tp] = annotation[sampling_idx][time_rep == tp]

        neigh = NearestNeighbors(n_neighbors=n_neighbors, algorithm="kd_tree").fit(embs[source])
        k_neighbors = neigh.kneighbors(embs[target], n_neighbors=n_neighbors, return_distance=False)

        k_neighbors_ann = []
        for k in range(n_neighbors):
            k_neighbors_ann.append(annts[source][k_neighbors[:, k]][:, np.newaxis])

        k_neighbors_ann = np.concatenate(k_neighbors_ann, axis=1)
        df_res = pd.DataFrame(index=states[target], columns=states[source])

        for state_tgt in states[target]:
            k_neighbors_state = k_neighbors_ann[annts[target] == state_tgt, ]
            for state_src in states[source]:
                df_res.loc[state_tgt, state_src] = np.sum(k_neighbors_state == state_src)

        df_norm = df_res.divide(df_res.sum(axis=1), axis="index")
        res.append(df_norm)
        res_full[i, :] = df_norm.values

    # output summarized mappings from source (tp0) to target (tp1)
    res_median = pd.DataFrame(np.median(res_full, axis=0), index=states[target], columns=states[source]).T
    res_std = pd.DataFrame(np.std(res_full, axis=0), index=states[target], columns=states[source]).T # sqrt(var)

    return res_median, res_std


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

    emb = adata_c.obsm["X_pca_30"]
    time = adata_c.obs["time"].values.astype(float)
    annotation = adata_c.obs["annotation"].values.astype(str)

    res_median, _ = createlineage_knn(
        emb = emb,
        time=time,
        annotation=annotation,
        n_reps=args.n_reps,
        ratio_remove_cells=args.ratio_remove_cells,
        n_neighbors=args.n_neighbors
    )

    res_median.divide(res_median.sum(axis=0), axis="columns")
    res_median = res_median.fillna(0.0)
    accuracy_val = accuracy(res_median, G)
    print(f"for {tps_c} accuracy is: {accuracy_val}")
    fname = DATA_DIR + f"output_tome/mouse_embryo_tp_{args.tps_couple}.csv"
    pd.DataFrame({"tp":[args.tps_couple], "start": [start], "end": [end], "accuracy": [accuracy_val]}).to_csv(fname)
    return 


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tps_couple',
                        type=int,
                        required=False,
                        default=0)
    parser.add_argument('--n_reps',
                        type=int,
                        required=False,
                        default=500)
    parser.add_argument('--ratio_remove_cells',
                        type=float,
                        required=False,
                        default=0.2)
    parser.add_argument('--n_neighbors',
                        type=int,
                        required=False,
                        default=5)
    args = parser.parse_args()

    solve(args)