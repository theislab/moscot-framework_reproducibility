## Moscot enables atlas-scale temporal mapping of mouse embryogenesis
We applied moscot to a 1.7 M. cells mouse developental atlas. The data was preprocessed, annotated and partially created/ sequenced deeper by Qiu et al. (https://www.nature.com/articles/s41588-022-01018-x)
For reproducibility, the analysis was divided into 5 folders:


## Notebook folders





<details>
    <summary>0_Data_preparation  </summary>
Contains Seurat integration and conversion of the downloaded .RDS files (http://tome.gs.washington.edu/) into anndata objects.

#### 0_Integration_notebooks:
  
  Runs TOME on the integrated data as done by Qiu et al. It contains the following notebooks:
  
   * ```MG_05-01-2023_Seurat_Integartion.ipynb```: Performs Seurat's anchor based batch correction analogous to Qiu et al, using using code obtained from https://github.com/ChengxiangQiu/tome_code
   * ```MG_05-01-2023_Seurat_Integartion_E8.5b-E9.5_Redone.ipynb```: Performs the same integration, but using 3000 hvgs instead of 2000hvgs since integration with 2000 hvgs was not able to separate nerual crest and allantois sufficiently

#### 1_Seurat_object_to_anndata_notebooks:

  Transforms the downloaded .RDS objects into anndata objects, which are then concatenated and the intefration result is added.
  
   * ```MG_05-01-2023_Seurat_object_to_anndata.ipynb```: Runs SeuratDisk/Data to transform .RDS into anndata objects
   * ```MG_05-01-2023_Ensemble_to_gene_symbol.ipynb```: Uses Biomart to construct a dictionary translating ENSEMBL IDs to gene symbols
   * ```MG_05-01-2023_Fix_anndata_annotations.ipynb```: Metadata is not transformed correctly by SeuratDisk/Data. This is fixed here, and further annotations are added.
   * ```MG_05-01-2023_Concatenate_time_pair_anndatas.ipynb```: Anndatas of adjacent time points are concatenated and the latent representation obtained from the integration is added
   * ```MG_05-01-2023_adata_to_obs.ipynb```: Saves the Anndata annotation, which is needed when defining growth rates in the case when TOME has been run on the data where extraembryonic tissues have been removed.
</details>







<details>
    <summary>1_Cell_type_transition_analysis  </summary>
Both moscot and TOME are run on the same latent representation to obtain cell type transition rates, which are then evaluated.

#### 0_TOME:
  
  Performs integration as done by Qiu et al. It contains the following notebooks:
  
   * ```MG_05-01-2023_TOME_Maps_for_cell_type_transitions```: Runs TOME as in  https://github.com/ChengxiangQiu/tome_code.

#### 1_moscot:

  Runs moscot on the same representation as used in TOME  
   * ```Run_moscot.py```: Python script running moscot saving the resulting solution.
   * ```MG_05-01-2023_Check_growth_rates.ipynb```: Loads the calculated solution to inspect growth/apoptisis rates.
   * ```MG_05-01-2023_moscot_transport_matrix_to_cell_type_transitions.ipynb```: Used the moscot solutions to compute cell type transition rates.
   
   
#### 2_Validation:

  Evaluating the transitions obtained from TOME and moscot  
   * ```MG_05-01-2023_Evaluation_of_cell_type_transitions.ipynb```: Uses curated transitions and germ layer annotation (Supplementary Table 1) to calculate validation scores.
</details>




<details>
    <summary>2_Growth_rate_and_driver_gene_analysis  </summary>
To get a more detailed view of transitions on the cell level we extend the kNN-approach intruduced by to to cell-level TOME (clTOME), which is then compared to moscot. For this analyis, extraembryonic tissues (inlcuding Blood progenitors and Primitive erythroid cells until E8.5) have been excluded for gastulation and organogenesis.


#### 0_clTOME:
  
  Performs integration as done by Qiu et al. It contains the following notebooks:
  
   * ```MG_05-01-2023_TOME_Maps_for_growth_rate_and_driver_genes_analysis.ipynb```: Saves the identified neirest neighors obtained while running TOME.
   * ```MG_05-01-2023_Transforming_Identified_Neigbors_to_Transport_Matrix.ipynb```: Takes the neirest neighors and shapes them into a sparse matrix.
   * ```MG_05-01-2023_TOME_transport_matrix_to_growth_rates.ipynb```: Uses the neirest neighbor matrix to calculate growth rates.
   * ```MG_05-01-2023_TOME_transport_matrix_to_pulls.ipynb```: Uses the neirest neighbor matrix to calculate pulls of selected cell types.
   

#### 1_moscot:

  Runs moscot on the same representation as used in TOME  
   * ```Run_moscot.py```: Python script running moscot saving the resulting solution.
   * ```MG_05-01-2023_Check_growth_rates.ipynb```: Loads the calculated solution to inspect growth/apoptisis rates.
   * ```MG_05-01-2023_moscot_transport_matrix_to_growth_rates.ipynb```: Used the moscot solutions to compute growth rates.
   * ```MG_05-01-2023_moscot_transport_matrix_to_pulls.ipynb```: Used the moscot solutions to compute pulls of selected cell types.
   
   
#### 2_Validation:

  Evaluates obtained growth rates and cell type pulls
  
   * 0_scVI_computations: Contains 1 notebook running scVI to obtain scVI normalized gene expression.
   * 1_Driver_gene_correlations: Contains 4 notebooks calculating correlation of scVI normalized gene expression to cell type pulls.
   * 2_Growth_rate_correlations: Contains 2 notebooks used to calculate correlation of obtained growth rates to cell cycle scores.
   
</details>







<details>
    <summary>3_Memory_and_runtime_benchmark  </summary>
Benchmarking memory consumption and running time of WOT, moscot and moscot low rank.


#### 0_Subsampling:
  
  Subsamples cells from the biggest time pair into anndata objects.
  
   * ```MG_05-01-2023_E11.5_subsampling```: Subsamples such that earlier and later time point both contain the same amount of cells, which increases in steps of 25,000, starting form 0, up to 275,000 cells.

#### 1_Scripts:
  Contains python scripts and yaml_files with which the benchmark has been performed. For each yaml file the exists the corresponding python scirpt (e.g. bm_CPU_offline.yml and run_cpu_offline.py).

</details>


<details>
    <summary>4_Figures  </summary>
Notebooks to create plots and figures

   * ```MG_05-01-2023_Memory_and_runtime_benchmark.ipynb```: Plots result of memory and runtime benchmark.
   * ```MG_05-01-2023_Cell_type_transition_accuracy.ipynb```: Plots result of cell type transition analysis.
   * ```MG_05-01-2023_E8_UMAPs.ipynb```: Plots UMAPS of growth rates, pulls and gene expression of E8.0 to E8.25 data.
   * ```MG_05-01-2023_Plotting_marker_gene_correlation.ipynb```: Plots result of driver gene correlations.
   * ```MG_05-01-2023_Illustrate_growth_rates.ipynb```: Plots computed growth rates for all time pairs.

</details>
