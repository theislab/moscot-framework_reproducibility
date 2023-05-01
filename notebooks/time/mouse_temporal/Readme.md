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
Contains Seurat integration and conversion of the downloaded .RDS files (http://tome.gs.washington.edu/) into anndata objects.

#### 0_TOME:
  
  Performs integration as done by Qiu et al. It contains the following notebooks:
  
   * ```MG_05-01-2023_TOME_Maps_for_cell_type_transitions```: Runs TOME as in  https://github.com/ChengxiangQiu/tome_code.

#### 1_moscot:

  Runs moscot on the same representation as used in TOME  
   * ```Run_moscot.py```: Python script running moscot saving the resulting solution.
   * ```MG_05-01-2023_Check_growth_rates.ipynb```: Loads the calculated solution to inspect growth/apoptisis rates.
   * ```MG_05-01-2023_moscot_transport_matrix_to_cell_type_transitions.ipynb```: Used the moscot solutions to compute cell type transition rates.
   
   
#### 2_Validation:

  Runs moscot on the same representation as used in TOME  
   * ```MG_05-01-2023_Evaluation_of_cell_type_transitions.ipynb```: Uses curated transitions and germ layer annotation (Supplementary Table 1) to calculate validation scores.
</details>








<details>
    <summary>1_TOME/clTOME_computations</summary>
    
We ran TOME on the integrated data. Default TOME output are cell type transitions. To see if TOME's strategy also results in reasonable coupling on the single cell level we save the neirest neighobrs TOME identifies and transform it into a coupling/transport matrix.
    
   * ```MG_XXX_Running_TOME```: Runs TOME using code obtianed from https://github.com/ChengxiangQiu/tome_code. In addition to saving the cell type transitions, we also saved TOME's identified neirest neighbors, which will be used in cell-level (cl)TOME
   * ```MG_XXX_Transforming_Identified_Neigbors_to_Transport_Matrix.ipynb```: Notebooks where moscot is applied to the data
   * ```MG_XXX_TOME_transport_matrix_to_growth_rates.ipynb```: Notebooks where moscot is applied to the data
   * ```MG_XXX_TOME_transport_matrix_to_pulls.ipynb```: Notebooks where moscot is applied to the data
    
</details>

<details>
    <summary>2_moscot_computations</summary>
    
We ran moscot on the integrated data, afterwards we extracted cell type trasitions, grwoth rates and pulls of specific cell types for later evaluation.
    
   * ```MG_XXX_Running_moscot```: Notebooks where moscot is applied to the data
   * ```MG_XXX_moscot_map_to_cell_type_transitions```: Notebooks where moscot is applied to the data
   * ```MG_XXX_moscot_map_to_growth_rates```: Notebooks where moscot is applied to the data
   * ```MG_XXX_moscot_map_to_pull```: Notebooks where moscot is applied to the data
</details>


<details>
    <summary>3_Evaluation</summary>

We evaluated TOME/clTOME and moscot output using 3 different metrics: 
    
   * ```MG_XXX_Evalution_of_germ_layer_and_cell_type_transitions```: Notebooks where moscot is applied to the data
   * ```MG_XXX_Evaluation_of_growth_rates```: Notebooks where moscot is applied to the data
   * ```MG_XXX_Running_scVI```: scVI was used to infer gene expression using "get_normalized_genes"
   * ```MG_XXX_Evaluation_of_driver_gene_correlations```: Notebooks where moscot is applied to the data
    
</details>


<details>
    <summary>4_Memory_and_runtime_benchmark</summary>
    
We used the time pair with the most cells (E11.5 --> E12.5 with 455,124 cells --> 292,726 cells) and subsampled it such that E11.5 and E12.5 contain the same amount of cells x, where x was chosen in steps of 25,000. Memory and runtime of moscot and low rank version of moscot where compared to WaddingtonOT (PMID: 30712874). Supplementary Table XXX containes the results of this benchmark.
    
   * ```yaml_files_and_stuff```: Notebooks where moscot is applied to the data
   * ```MG_XXX_Evaluation_of_growth_rates```: Notebooks where moscot is applied to the data
    
</details>

<details>
    <summary>5_Figure_creation</summary>
    
For the main figure:
    
   * ```MG_XXX_Memory_and_runtime_plot```: Notebooks where moscot is applied to the data
   * ```MG_XXX_Germ_layer_and_cell_type_transition_plot```: Notebooks where moscot is applied to the data
   * ```MG_XXX_E8.0_to_E8.25_UMAPS:and_growth_rate_histogram```: Notebooks where moscot is applied to the data
   * ```MG_XXX_Driver_gene_correlation_plots```: Notebooks where moscot is applied to the data
    
For the supplementary figure:
    
   * ```MG_XXX_All_growth_rate_histograms```: Notebooks where moscot is applied to the data
   * ```MG_XXX_Growth_rate_to_cell_cycle_score_correlations```: Notebooks where moscot is applied to the data
    
</details>




## Rules
1. Have an **empty line** after the `</summary>` tag or markdown/code blocks will not render.
1. Have an **empty line** after each `</details>` tag if you have multiple collapsible sections.
