import jax
from jax import config
config.update("jax_enable_x64", True)

import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import moscot
from moscot.problems.time import TemporalProblem
import moscot.plotting as mpl
import pandas as pd
import os
import muon
from ott.geometry import pointcloud, geometry
from sklearn.preprocessing import StandardScaler
import networkx as nx
import itertools
import anndata
from mudata import MuData
import jax.numpy as jnp
from typing import Dict, Tuple
from ott import tools
from tqdm import tqdm
import jax
import sys

output_dir = "/lustre/groups/ml01/workspace/moscot_paper/pancreas_revision/stability_analysis"
mudata = muon.read("/lustre/groups/ml01/workspace/moscot_paper/pancreas_revision/mudata_with_annotation_all.h5mu")

endocrine_celltypes = [
    "Ngn3 low",
    "Ngn3 high",
    "Ngn3 high cycling",
    "Fev+",
    "Fev+ Alpha",
    "Fev+ Beta",
    "Fev+ Delta",
    "Eps. progenitors",
    "Alpha",
    "Beta",
    "Delta",
    "Epsilon"
]

mudata = mudata[mudata.obs["cell_type"].isin(endocrine_celltypes)].copy()


arguments = sys.argv
emb_1 = arguments[1] if arguments[1] != "None" else None
emb_2 = arguments[2] if arguments[2] != "None" else None
cost_1 = arguments[3] if arguments[3] != "None" else None
cost_2 = int(arguments[4]) if arguments[4] != "None" else None

EMB = "embedding"

def adapt_time(x):
        if x["stage"]=="E14.5":
            return 14.5
        if x["stage"]=="E15.5":
            return 15.5
        if x["stage"]=="E16.5":
            return 16.5
        raise ValueError
    
def create_adata(mudata: MuData, rna_embedding: str, atac_embedding: str) -> anndata.AnnData:

    adata = mudata["rna"]
    adata.obs["cell_type_refined"] = mudata.obs["cell_type_refined"]
    adata.obs['time'] = adata.obs.apply(adapt_time, axis=1).astype("category")
    if rna_embedding == "X_MultiVI":
        adata.obsm[EMB] = mudata.obsm[rna_embedding].copy()
        return adata
    else:
        if rna_embedding is not None:
            rna_emb = mudata["rna"].obsm[rna_embedding]
            rna_emb_scaled = rna_emb #StandardScaler().fit_transform(rna_emb)
        if atac_embedding is not None:
            atac_emb = mudata["atac"].obsm[atac_embedding]
            atac_emb_scaled = atac_emb # StandardScaler().fit_transform(atac_emb)
        if rna_embedding is not None and atac_embedding is not None:
            emb = np.concatenate((rna_emb_scaled, atac_emb_scaled), axis=1)
        elif rna_embedding is None and atac_embedding is not None:
            emb = atac_emb
        elif rna_embedding is not None and atac_embedding is None:
            emb = rna_emb
        else:
            raise NotImplementedError
    
        adata.obsm[EMB] = emb
        return adata


def create_graphs(adata: anndata.AnnData, n_neighbors: int) -> Dict[Tuple, pd.DataFrame]:
    dfs = {}
    batch_column = "time"
    unique_batches = [14.5, 15.5, 16.5]
    for i in range(len(unique_batches) - 1):
        batch1 = unique_batches[i]
        batch2 = unique_batches[i + 1]
    
        indices = np.where((adata.obs[batch_column] == batch1) | (adata.obs[batch_column] == batch2))[0]
        adata_subset = adata[indices]
        sc.pp.neighbors(adata_subset, use_rep=EMB, n_neighbors=n_neighbors)
        G = nx.from_numpy_array(adata_subset.obsp["connectivities"].A)
        assert nx.is_connected(G)
    
        dfs[(batch1, batch2)] = (
            pd.DataFrame(
                index=adata_subset.obs_names,
                columns=adata_subset.obs_names,
                data=adata_subset.obsp["connectivities"].A.astype("float64"),
            )
        )
    return dfs

cm = jnp.ones((144, 144)) - jnp.eye(144)

def compute_metrics(df_reference: jax.Array, df: pd.DataFrame, emb_0: str, emb_1: str, cost_0: str, cost_1: str) -> pd.DataFrame:
    
    sink_div = tools.sinkhorn_divergence.sinkhorn_divergence(geometry.Geometry, cost_matrix=(cm,cm,cm), a=df_reference.values.flatten(), b=df.values.flatten(), epsilon=1e-2).divergence
    eps_from_eps_prog = df.loc["Eps. progenitors", "Epsilon"]
    delta_from_fev_delta = df.loc["Fev+ Delta", "Delta"]
    fev_delta_from_eps_prog = df.loc["Eps. progenitors", "Fev+ Delta"]
    eps_from_fev_delta = df.loc["Fev+ Delta", "Epsilon"]
    beta_from_fev_beta = df.loc["Fev+ Beta", "Beta"]
    delta_from_ngn3_low = df.loc["Ngn3 low", "Delta"]
    print(emb_0 , emb_1, cost_0, cost_1, sink_div)
    print(eps_from_eps_prog, delta_from_fev_delta, fev_delta_from_eps_prog, eps_from_fev_delta)
    data = [[str(emb_0),str(emb_1), str(cost_0), str(cost_1), sink_div, eps_from_eps_prog, delta_from_fev_delta, fev_delta_from_eps_prog, eps_from_fev_delta, beta_from_fev_beta, delta_from_ngn3_low]]
    
    return pd.DataFrame(data=data, columns=["emb_0", "emb_1", "cost_0", "cost_1", "sink_div", "eps_from_eps_prog", "delta_from_fev_delta", "fev_delta_from_eps_prog", "eps_from_fev_delta", "beta_from_fev_beta", "delta_from_ngn3_low"])
                            
order_cell_types = endocrine_celltypes


tp_reference = TemporalProblem.load("/lustre/groups/ml01/workspace/moscot_paper/pancreas_revision/plots/OT_encodrine_analysis/TemporalProblem.pkl")
metrics_early = pd.DataFrame(columns=["emb_0", "emb_1", "cost_0", "cost_1", "sink_div", "eps_from_eps_prog", "delta_from_fev_delta", "fev_delta_from_eps_prog", "eps_from_fev_delta", "beta_from_fev_beta", "delta_from_ngn3_low"])
metrics_late= pd.DataFrame(columns=["emb_0", "emb_1", "cost_0", "cost_1", "sink_div", "eps_from_eps_prog", "delta_from_fev_delta", "fev_delta_from_eps_prog", "eps_from_fev_delta", "beta_from_fev_beta", "delta_from_ngn3_low"])
            

reference_tmap_early = tp_reference.cell_transition(14.5, 15.5, {"cell_type": order_cell_types}, {"cell_type": order_cell_types}, forward=False)
reference_tmap_late = tp_reference.cell_transition(15.5, 16.5, {"cell_type": order_cell_types}, {"cell_type": order_cell_types}, forward=False)

adata = create_adata(mudata, emb_1, emb_2)

tp = TemporalProblem(adata)
tp = tp.prepare("time", joint_attr=EMB, cost=cost_1 if cost_1!="geodesic" else "sq_euclidean")

if cost_1=="geodesic":
    dfs = create_graphs(adata, cost_2)
    tp[14.5, 15.5].set_graph_xy((dfs[14.5, 15.5]).astype("float"), t=100.0)
    tp[15.5, 16.5].set_graph_xy((dfs[15.5, 16.5]).astype("float"), t=100.0)

tp = tp.solve(max_iterations=1e7)

df_early = tp.cell_transition(14.5, 15.5, {"cell_type": order_cell_types}, {"cell_type": order_cell_types}, forward=False)
df_late = tp.cell_transition(15.5, 16.5, {"cell_type": order_cell_types}, {"cell_type": order_cell_types}, forward=False)
metrics_early = pd.concat((metrics_early, compute_metrics(reference_tmap_early, df_early,emb_1, emb_2,  str(cost_1), str(cost_2))))
metrics_late = pd.concat((metrics_late, compute_metrics(reference_tmap_late, df_late, emb_1, emb_2, str(cost_1), str(cost_2))))

metrics_early.to_csv(os.path.join(output_dir, f"stability_metrics_early_{emb_1}_{emb_2}_{cost_1}_{cost_2}_new2.csv"))
metrics_late.to_csv(os.path.join(output_dir, f"stability_metrics_late_{emb_1}_{emb_2}_{cost_1}_{cost_2}_new2.csv"))
