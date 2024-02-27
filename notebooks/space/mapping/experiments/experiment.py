from pathlib import Path

from omegaconf import DictConfig
from hydra.core.hydra_config import HydraConfig
import hydra
import pandas as pd
from anndata import AnnData
import os


def process_adata(cfg: DictConfig) -> tuple[AnnData, AnnData, AnnData]:
    import scanpy as sc
    from sklearn.preprocessing import StandardScaler
    import numpy as np

    scaler = StandardScaler()
    path_data = Path(cfg.paths.path_data)
    adata_sc = sc.read(path_data / cfg.paths.adata_sc)
    adata_p = sc.read(path_data / cfg.paths.adata_p)
    adata_spatial = sc.read(path_data / cfg.paths.adata_sp)

    adata_sc.obsm["cite"] = adata_p.X.copy()
    adata_sc.obsm["X_pca_cite"] = adata_p.obsm["X_pca"].copy()

    adata_spatial.X = adata_spatial.layers["normalized"].copy()
    sc.pp.pca(adata_spatial, n_comps=50)
    spatial = adata_spatial.obsm["spatial"]
    adata_spatial.obsm["spatial_scaled"] = (spatial - spatial.mean()) / spatial.std()
    adata_spatial.obsm["X_pca_scaled"] = scaler.fit_transform(
        adata_spatial.obsm["X_pca"]
    )
    adata_spatial.obsm["X_pca_spatial_scaled"] = np.concatenate(
        [adata_spatial.obsm["X_pca_scaled"], adata_spatial.obsm["spatial_scaled"]],
        axis=1,
    )

    return adata_sc, adata_spatial, adata_p


def _corr_results(true_df: pd.DataFrame, pred_df: pd.DataFrame) -> pd.DataFrame:
    pred_df = pred_df[true_df.columns].copy()
    corr_pearson = pred_df.corrwith(true_df, method="pearson")
    corr_spearman = pred_df.corrwith(true_df, method="spearman")
    out = pd.concat([corr_pearson, corr_spearman], axis=1)
    out.columns = ["pearson", "spearman"]
    return out


def benchmark(cfg):
    from datetime import datetime
    import numpy as np
    import scanpy as sc
    from ott.geometry import costs, pointcloud
    from ott.problems.linear import linear_problem
    from ott.solvers.linear import univariate
    import jax.numpy as jnp
    from jax.config import config

    config.update("jax_enable_x64", True)

    unique_id = HydraConfig.get().job.override_dirname
    path_results = Path(cfg.paths.path_results) / unique_id

    os.makedirs(path_results, exist_ok=True)

    adata_sc, adata_spatial, adata_p = process_adata(cfg)
    var_names = list(
        set(adata_sc.var_names.values).intersection(adata_spatial.var_names.values)
    )

    from moscot.problems.space import MappingProblem

    # sc.pp.subsample(adata_sc, fraction=0.1)
    # sc.pp.subsample(adata_spatial, fraction=0.1)
    mp = MappingProblem(adata_sc, adata_spatial).prepare(
        sc_attr={"attr": "obsm", "key": cfg.moscot.sc_attr},
        joint_attr={"attr": "X"},
        spatial_key="X_pca_scaled",
        normalize_spatial=False,
        var_names=var_names,
        cost=cfg.moscot.cost,
    )

    start_time = datetime.now()

    mp = mp.solve(
        alpha=cfg.moscot.alpha,
        rank=cfg.moscot.rank,
        epsilon=cfg.moscot.epsilon,
        gamma=cfg.moscot.gamma,
        initializer=cfg.moscot.initializer,
        tau_a=cfg.moscot.tau_a,
        tau_b=1,
        min_iterations=10_000,
        max_iterations=100_000,
        threshold=1e-20,
    )
    end_time = datetime.now()

    mp.problems[("src", "tgt")].solution.plot_costs(save=path_results / "costs.png")

    # impute protein
    prot = adata_p[adata_sc.obs_names].X.A
    spatial_prot = mp.problems[("src", "tgt")].solution.pull(
        prot,
        scale_by_marginals=True,
    )
    adata_spatial_prot = AnnData(
        np.array(spatial_prot),
        obs=adata_spatial.obs.copy(),
        obsm=adata_spatial.obsm.copy(),
    )
    adata_spatial_prot.var_names = adata_p.var_names.copy()

    # impute gexp
    genes = genes = ["Gja5", "Adgrg6", "Mmrn1", "Epcam", "Spp1"]
    adata_genes_pred = mp.impute(var_names=genes, device="cpu")
    adata_spatial_prot.obsm["pred_genes"] = sc.get.obs_df(adata_genes_pred, keys=genes)

    # gexp = adata_sc[adata_sc.obs_names, genes].X
    # spatial_gexp = mp.problems[("src", "tgt")].solution.pull(
    #     gexp,
    #     scale_by_marginals=True,
    # )
    # adata_spatial_gexp = AnnData(
    #     np.array(spatial_gexp),
    #     obs=adata_spatial.obs.copy(),
    #     obsm=adata_spatial.obsm.copy(),
    # )
    # adata_spatial_gexp.var_names = genes

    # evaluate imputed gexp
    adata_pred = mp.impute(var_names=var_names, device="cpu")

    pred_df = sc.get.obs_df(adata_pred, keys=var_names)
    true_df = sc.get.obs_df(adata_spatial, keys=var_names)
    np.testing.assert_array_equal(pred_df.columns, true_df.columns)
    corr_results = _corr_results(true_df, pred_df)
    corr_results.to_csv(path_results / "corr_results.csv")


    # impute annotation
    dummy = pd.get_dummies(adata_sc.obs["annot"])
    temp = mp[("src", "tgt")].pull(dummy, scale_by_marginals=True)
    clusters = pd.Categorical([dummy.columns[i] for i in np.array(temp.argmax(1))])
    adata_spatial_prot.obs["celltype_mapped"] = clusters
    
    results = {"time": end_time - start_time}

    # get results
    for corr in ["pearson", "spearman"]:
        results[f"corr_results_mean_{corr}"] = corr_results[corr].mean()
        results[f"corr_results_median_{corr}"] = corr_results[corr].median()
        results[f"corr_results_var_{corr}"] = corr_results[corr].var()
        results[f"corr_results_max_{corr}"] = corr_results[corr].max()
        results[f"corr_results_min_{corr}"] = corr_results[corr].min()

    # get protein results
    celltype_protein_target = {
        "Endothelial cells": "ESAM",
        "Endothelial cells": "CD300LG",
        "Endothelial cells": "CD38",
        "Kupffer cells": "CD86",
        "Kupffer cells": "Tim4",
        "Kupffer cells": "Folate",
    }

    solver = univariate.UnivariateSolver()
    for celltype, protein in celltype_protein_target.items():
        ct = adata_spatial_prot.obs["celltype_mapped"].isin([celltype]).values
        nct = ~ct
        x = adata_spatial_prot[ct, protein].X
        y = adata_spatial_prot[nct, protein].X
        if len(x) > 100 and len(y) > 100:
            a = jnp.ones(x.shape) / len(x)
            b = jnp.ones(y.shape) / len(y)
            geom = pointcloud.PointCloud(x, y, cost_fn=costs.PNormP(1.0))
            prob = linear_problem.LinearProblem(geom=geom, a=a, b=b)
            ott_d = solver(prob).ot_costs[0].item()
            results[f"wd_{celltype}_{protein}"] = ott_d
        else:
            results[f"wd_{celltype}_{protein}"] = np.nan

    pd.DataFrame(results, index=[0]).to_csv(path_results / "results.csv")
    adata_spatial_prot.write(path_results / "adata_spatial_prot.h5ad")


@hydra.main(config_path=".", config_name="config.yaml", version_base="1.2")
def main(cfg: DictConfig) -> None:
    sweep_params = HydraConfig.get().job.override_dirname
    print("SWEEP PARAMS", sweep_params)
    benchmark(cfg)


if __name__ == "__main__":
    main()

# python experiment.py +method.name=moscot +method.epsilon=1,10 +method.lambda=10 --multirun
