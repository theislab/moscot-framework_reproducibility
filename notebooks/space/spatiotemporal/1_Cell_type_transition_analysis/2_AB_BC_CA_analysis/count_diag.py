# %%
import moscot as mt
import pandas as pd
import numpy as np
import gc
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
from ott.geometry import pointcloud
from ott.problems.linear import linear_problem
from ott.solvers.linear import sinkhorn
from ott.tools import sinkhorn_divergence
import jax.numpy as jnp
import pandas as pd

# %%
adata = mt.datasets.mosta()
adata.obs["time"] = pd.Categorical(adata.obs["time"])
# %%
stp = mt.problems.SpatioTemporalProblem(adata)
# %%
stp = stp.score_genes_for_marginals(
    gene_set_proliferation="mouse", gene_set_apoptosis="mouse"
)
# %%
stp = stp.prepare(
    time_key="time",
    spatial_key="spatial",
    joint_attr=None,
    callback="local-pca",
    policy="explicit",
    subset=[(9.5, 10.5), (10.5, 11.5), (11.5, 9.5)],
)
# %%
result = {}
for epsilon in [1e-2, 5e-2, 1e-1, 5e-1]:
    stp.solve(epsilon=epsilon)
    print(f"done epsilon {epsilon}")
    t0 = stp[(9.5, 10.5)].solution.transport_matrix
    t1 = stp[(10.5, 11.5)].solution.transport_matrix
    t2 = stp[(11.5, 9.5)].solution.transport_matrix
    tmat = t0 @ t1 @ t2
    tmat = np.fliplr(tmat)
    del t0, t1, t2
    gc.collect()
    diagonal_entries = np.diag(tmat)
    comparison = tmat >= diagonal_entries[:, np.newaxis]
    count_larger_than_diagonal = np.sum(comparison, axis=1)
    print(count_larger_than_diagonal.mean())
    result[epsilon] = count_larger_than_diagonal.mean()

# tmat = np.diag(np.eye(len(diagonal_entries)) / (len(diagonal_entries) ** 2))
# comparison = tmat >= diagonal_entries[:, np.newaxis]
# count_larger_than_diagonal = np.sum(comparison, axis=1)
# print(count_larger_than_diagonal.mean())
# result[np.nan] = count_larger_than_diagonal.mean()
result["outer"] = ((len(tmat) - 1) * len(tmat)) / (len(tmat))
print(result)
df = pd.DataFrame(result, index=[0]).melt()
df.columns = ["epsilon", "mean_count_diag"]
df.to_csv("mean_count_diag.csv")
