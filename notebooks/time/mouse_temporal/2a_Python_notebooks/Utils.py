import numpy as np
import anndata
import pandas as pd
import scanpy as sc
import scipy

Path="/home/icb/manuel.gander/moscotTime_Reproducibility/Data"

def load_adata(Path, ts0, ts1):
    
    # Ignore warning here
    import warnings
    warnings.filterwarnings('ignore')
    sc.settings.verbosity = 0
    A0=sc.read(f"{Path}/anndatas/adata_{ts0}.h5ad")
    A1=sc.read(f"{Path}/anndatas/adata_{ts1}.h5ad")

    if ts0=='E8.5a':
        A0.obs['day']=0
    adata=A0.concatenate(A1, join='outer', index_unique=None, batch_key=None).copy()
    
    
    # For E8.5b to E9.5, use the recomputed (better) integration (uses more hvgs=features)
    if ts0=='E8.5b':
        k='_new'
    else:
        k=''

    # Load the representation from Seurat integration
    PCA=pd.read_csv(f"{Path}/Seurat_Representations/{ts0}_{ts1}_pca{k}.csv", sep= ",", index_col='Unnamed: 0')
    adata.obsm['X_pcaS']=PCA.loc[list(adata.obs['cellID'])].values

    UMAP=pd.read_csv(f"{Path}/Seurat_Representations/{ts0}_{ts1}_umap3{k}.csv", sep= ",", index_col='Unnamed: 0')
    adata.obsm['X_umap3']=UMAP.loc[list(adata.obs['cellID'])].values
    
    # Put the gene names into adata.var.index for scoring genes later
    adata.var['index']=[str(a) for a in adata.var['gene_names']]
    adata.var=adata.var.set_index('index')
    adata.var_names_make_unique()
    
    return(adata)


def growth_rates_to_apoptosis_ratio(growth_rates, ts0, ts1):
    ts=['E3.5', 'E4.5', 'E5.25', 'E5.5', 'E6.25', 'E6.5', 'E6.75', 'E7.0', 'E7.25', 'E7.5', 'E7.75', 'E8.0', 'E8.25', 'E8.5a', 'E8.5b', 'E9.5', 'E10.5', 'E11.5', 'E12.5', 'E13.5']
    
    # I got these cell numbers from http://tome.gs.washington.edu/, and they got them from the experiments or 
    # for E8.5b they estimated it themselves (as it was their experiment)
    
    cells=[32, 80, 100, 120, 400, 660, 1720, 4500, 8200, 15000, 30000, 60000, 73000, 90000, 90000, 200000, 1100000, 2600000, 6000000, 13000000]
    Cell_number_dict={}
    for i in range(20):
        Cell_number_dict[ts[i]]=cells[i]
    
    growth_rates=growth_rates/growth_rates.mean()
    cellular_growth_rates=growth_rates*Cell_number_dict[ts1]/Cell_number_dict[ts0]
    apoptotic_cells=cellular_growth_rates[cellular_growth_rates<1]
    sum_apoptotic_cells=(1-apoptotic_cells).sum()
    perc_apoptotic_cells=sum_apoptotic_cells/len(growth_rates)
    
    return(perc_apoptotic_cells, cellular_growth_rates)
