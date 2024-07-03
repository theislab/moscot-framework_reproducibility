# %%
import moscot as mt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
# %%
df = pd.read_csv("./mean_count_diag.csv", index_col=0)
df["mean_count_diag"] = df["mean_count_diag"]/df["mean_count_diag"].max()
df = pd.read_csv("./mean_count_diag_2.csv", index_col=0)
df["mean_count_diag"] = df["mean_count_diag"]/df["mean_count_diag"].max()
# df = pd.concat([df, df2])
df["timepoints"] = df["timepoints"].replace(
    {
        # "10.5-11.5_11.5-12.5_12.5-10.5":"A[10.5-11.5] B[11.5-12.5] C[12.5-10.5]",
        "9.5-10.5_10.5-11.5_11.5-9.5":"A[9.5-10.5] B[10.5-11.5] C[11.5-9.5]"
        }
    )
# %%
import seaborn as sns
import mplscience
mplscience.set_style()
g = sns.scatterplot(data=df, x="epsilon", y="mean_count_diag", hue="timepoints", s=100)
g.set(xlabel="epsilon", ylabel="Fraction of entries equal or larger \n than diagonal in transport matrix", title="A[E9.5-E10.5] B[E10.5-E11.5] C[E11.5-E9.5 \n")
g.legend(loc='center left', bbox_to_anchor=(1, 0.5)).remove()
# %%
import squidpy as sq
adata = mt.datasets.mosta()
adata.obs["time"] = pd.Categorical(adata.obs["time"])
adata.obs["x"] = adata.obsm["spatial"][:,0]
adata.obs["y"] = adata.obsm["spatial"][:,1]
adata.obs['color_index'] = np.arange(len(adata))
# %%
# sq.pl.spatial_scatter(adata[adata.obs.time==9.5], shape=None, size=10, color="x")
# sq.pl.spatial_scatter(adata[adata.obs.time==9.5], shape=None, size=10, color="y")
sq.pl.spatial_scatter(adata[adata.obs.time==9.5], shape=None, size=10, color="color_index", title="Color index")

# %%
