# %%
import moscot as mt
import pandas as pd
import numpy as np
import anndata as ad
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
import mplscience

mplscience.set_style()
# %%
adata = mt.datasets.mosta()
adata.obs["time"] = pd.Categorical(adata.obs["time"])
adata_12 = sc.read("/lustre/groups/ml01/workspace/moscot_paper/mosta/mosta_E12.5.h5ad")
adata_12.obs["time"] = pd.Categorical(adata_12.obs["time"])
adata = ad.concat([adata, adata_12])
print(adata.obs["time"].unique())
# %%

# %%
# stp = stp.score_genes_for_marginals(
#     gene_set_proliferation="mouse", gene_set_apoptosis="mouse"
# )
# %%
final_df = []
for timepoint, subset in zip(
    [
        # "10.5-11.5_11.5-12.5_12.5-10.5",
        "9.5-10.5_10.5-11.5_11.5-9.5",
    ],
    [
        # [(10.5, 11.5), (11.5, 12.5), (12.5, 10.5)],
        [(9.5, 10.5), (10.5, 11.5), (11.5, 9.5)],
    ],
):
    stp = mt.problems.SpatioTemporalProblem(adata)
    stp = stp.prepare(
        time_key="time",
        spatial_key="spatial",
        joint_attr=None,
        callback="local-pca",
        policy="explicit",
        subset=subset,
    )

    result = {}
    for epsilon in [0.05, 0.1, 1, 10]:
        stp.solve(epsilon=epsilon)
        print(f"done epsilon {epsilon}")
        print(subset)
        t0 = stp[subset[0]].solution.transport_matrix
        t1 = stp[subset[1]].solution.transport_matrix
        t2 = stp[subset[2]].solution.transport_matrix
        tmat = t0 @ t1 @ t2
        tmat = np.fliplr(tmat)
        del t0, t1, t2
        gc.collect()
        diagonal_entries = np.diag(tmat)
        comparison = tmat >= diagonal_entries[:, np.newaxis]
        count_larger_than_diagonal = np.sum(comparison, axis=1)
        print(count_larger_than_diagonal.mean())
        result[epsilon] = count_larger_than_diagonal.mean()

        if epsilon == 0.01:
            plt.imshow(tmat, vmax=1e-18)
            plt.savefig(f"tmat_{timepoint}.png")

    result["outer"] = ((len(tmat) - 1) * len(tmat)) / (len(tmat))
    print(result)
    df = pd.DataFrame(result, index=[0]).melt()
    df.columns = ["epsilon", "mean_count_diag"]
    df["timepoints"] = timepoint
    final_df.append(df)

df = pd.concat(final_df)
df.to_csv("mean_count_diag_2.csv")
