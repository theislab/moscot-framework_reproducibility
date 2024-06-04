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
stp.solve(epsilon=0, rank=300, gamma=100)
# %%
t0 = stp[(9.5, 10.5)].solution.transport_matrix
t1 = stp[(10.5, 11.5)].solution.transport_matrix
t2 = stp[(11.5, 9.5)].solution.transport_matrix
tmat = t0 @ t1 @ t2
# del t0, t1, t2
# gc.collect()
# plt.imshow(total)
# %%
tmat = np.fliplr(tmat)
sc.pp.pca(adata)
feat = adata[adata.obs.time == 9.5].obsm["X_pca"]
idx = np.random.choice(len(feat), size=(50,))
feat = feat[idx].copy()
N = len(feat)
B = np.repeat(feat, N, axis=0)
A = np.tile(feat, (N, 1))
# pc = np.concatenate((A, B), axis=1)
# %%

# A_sub, B_sub = A[idx].copy(), B[idx].copy()
total_sub = tmat[idx, ...][..., idx].flatten().copy()
eye_sub = (np.eye(len(idx)) / len(idx)).flatten()


# %%
def sink_divergence(x, y, a, b, epsilon=0.01) -> jnp.ndarray:
    """Sink divergence."""
    ot = sinkhorn_divergence.sinkhorn_divergence(
        pointcloud.PointCloud,
        x=x,
        y=y,
        epsilon=epsilon,
        static_b=True,
        a=a,
        b=b,
    )
    return ot.divergence


print(sink_divergence(A, B, total_sub, eye_sub, epsilon=0.01))
# %%

# %%
idx = np.random.choice(len(feat), size=(50,))
shuffled_total = tmat[idx, ...][..., idx].flatten().copy()
# shuffled_total = np.repeat(1 / (len(total_sub) ** 2), len(total_sub))
print(sink_divergence(A, B, shuffled_total, eye_sub, epsilon=0.01))
