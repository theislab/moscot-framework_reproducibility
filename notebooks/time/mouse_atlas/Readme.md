## Recovering trajectories of mouse embryogenesis on a developmental atlas using moscot
We applied moscot to a 1.7M cells mouse developental atlas [[1]](https://www.nature.com/articles/s41588-022-01018-x). First we conducted a memory and runtime analyisis comparing moscot to Waddington OT [[2]](https://www.sciencedirect.com/science/article/pii/S009286741930039X?via%3Dihub) and moscot low-rank approaches. Then we benchmarked the performance of moscot and TOME [[1]](https://www.nature.com/articles/s41588-022-01018-x) using 3 validation metrics: cell type transition accuracy, cellular growth rates and driver gene correlations. 

## Notebook folders

### 0_Data_preparation
Contains Seurat integration and conversion of the downloaded .RDS files (http://tome.gs.washington.edu/) into anndata objects. We share the created anndatas on figshare (https://figshare.com/articles/dataset/Mouse_embryogenesis_atlas_-_pariwise_anndatas/25040537)

<details>
    <summary>0_Integration_notebooks </summary>
    &nbsp; 
    
   * ```MG_05-01-2023_Seurat_Integartion.ipynb```: Performs Seurat's anchor based batch correction analogous to Qiu et al, using using code obtained from https://github.com/ChengxiangQiu/tome_code.
   * ```MG_05-01-2023_Seurat_Integartion_E8.5b-E9.5_Redone.ipynb```: Performs the same integration, but using 3000 hvgs instead of 2000 hvgs since integration with 2000 hvgs was not able to separate neural crest and allantois sufficiently.
</details>

<details>
    <summary>1_Seurat_object_to_anndata_notebooks </summary>
    &nbsp;  

  Transforms the downloaded .RDS objects into anndata objects, which are then concatenated and the integration emebdding is added.
  
   * ```MG_05-01-2023_Seurat_object_to_anndata.ipynb```: Runs SeuratDisk/Data to transform .RDS into anndata objects.
   * ```MG_05-01-2023_Ensemble_to_gene_symbol.ipynb```: Uses Biomart to construct a dictionary translating ENSEMBL IDs to gene symbols.
   * ```MG_05-01-2023_Fix_anndata_annotations.ipynb```: Metadata is not transformed correctly by SeuratDisk/Data. This is fixed here, and additional annotations are added.
   * ```MG_05-01-2023_Concatenate_time_pair_anndatas.ipynb```: Anndatas of adjacent time points are concatenated and the latent representation obtained from the integration is added
   * ```MG_05-01-2023_adata_to_obs.ipynb```: Saves the Anndata annotation. These are needed later, for cl-TOME in the growth rate and driver gene correlation analysis.
</details>





### 1_Main_figure
Contains the notebooks to create all figures displayed in Figure 2 of the manuscript.

<details>
    <summary>0_Memory_and_runtime_benchmark </summary>
    &nbsp; 

Benchmarks the memory usage (on CPU) and runtime for moscot and WOT. For moscot low-rank we chose a rank 2,000 here.
    
   * 0_Slurm_jobs: contains the scripts to run moscot and WOT both on CPU and GPU (only for moscot). We scheduled the runs using using sbatch and I log the evaluation results to weights-and-biases.
   * 1_Results: Contains the benchmark results and notebooks to load and plot them.
</details>

<details>
    <summary>1_Cell_type_transition_analysis  </summary>
    Runs TOME and moscot and evaluates the predictions on a cell-type level
    &nbsp; 

   * ```0_Run_TOME.ipynb```: This notebook runs TOME on the Seurat-integrated embedding and saves the result.
   * ```1_Run_moscot.ipynb```: This notebook runs moscot on the Seurat-integrated embedding and saves the result.
   * ```2_Evaluate_and_plot.ipynb```: We use the germ-layer and curated transitions to evaluate the predicted transitions and create plots.
</details>

<details>
    <summary>2_Growth_rate_and_driver_gene_analysis  </summary>
    &nbsp; 
    
Again computes the transitions for both TOME and moscot, this time evalauting on a single cell level.
   

#### 0_cl-TOME:
  
Performs integration as done by Qiu et al. It contains the following notebooks:
  
   * ```0_TOME_Maps_for_growth_rate_and_driver_genes_analysis.ipynb```: Saves the identified neirest neighbors obtained while running TOME.
   * ```1_Transforming_Identified_Neigbors_to_Transport_Matrix.ipynb```: Takes the neirest neighbor files and shapes them into a sparse matrix.
   * ```2_TOME_transport_matrix_to_growth_rates.ipynb```: Uses the neirest neighbor matrix to calculate growth rates.
   * ```3_TOME_transport_matrix_to_pulls.ipynb```: Uses the neirest neighbor matrix to calculate pulls of selected cell types.
   
#### 1_moscot:

Runs moscot on all individual time points:

   * ```0_Run_moscot.ipynb```:  This saves the growth rates and also the pull of selected cell types for some time points.

#### 2_Evaluation:

Here we compute the scVI embedding to obtain an inferred count matrix. This scVI-matrix is then used to correlate the pull of selected cell populations to driver genes:

   * ```0_Run_scVI_for_gene_expression_inference.ipynb```:  This runs scVI on the four time points we are correlating pulls to driver gene expression.
   * ```1_Definitive_endoderm.ipynb```:  Loads the pull of E7.25 definitive endoderm for clTOME and msocot and correlates it to driver genes
   * ```2_Allantois.ipynb```:  Loads the pull of E8.0 allantois for clTOME and msocot and correlates it to driver genes
   * ```3_First_heart_field.ipynb```:  Loads the pull of E8.25 first heart field for clTOME and msocot and correlates it to driver genes
   * ```4_Pancreatic_epithelium.ipynb```:  Loads the pull of E11.5 pancreatic epithelium for clTOME and msocot and correlates it to driver genes

#### 3_Plots:

We create the UMAPs, growth rates and pull-correlation plots here:

   * ```0_E8_UMAPs.ipynb```:  Computes the UMAP for E8.0 to E8.25 and plots the growth rates and the driver gene correlations for the first heart field
   * ```1_Plotting_marker_gene_correlation.ipynb```:  This only creates the boxplot for driver gene to pull correlations

</details>



### 2_Supplementary_figures

This contains all the notebooks and scripts to create the supplementary figures that are connected to the mouse atlas data as well as the IPS reprogramming data


<details>
    <summary>0_Growth_rates </summary>
    &nbsp; 

Plots the growth rates for the mouse atlas data for all time pairs.
    
   * ```0_Create_obs.ipynb```: In order to not always have to relaod the full anndata-objects I only save the .obs-entry for convenience
   * ```1_Illustrate_growth_rates.ipynb```: Loads the growth rates for clTOME and moscot and plots them as histograms
</details>

<details>
    <summary>1_moscot_WOT_comparison </summary>
    &nbsp; 

For the reviews we made sure that the mapping obtained by moscot and WOT are comparable. This is done here:

   * ```0_Run_WOT.ipynb```: We run WOT on the mouse atlas up to E8.5 using the same initial growth rates, entropic regularization, and unbalancedness parameters as moscot.
   * ```1_Run_moscot.ipynb```: We rerun moscot using the cost normalized by the median (since WOT normalizes the cost matrix by the median, and msocot by the mean by default)
   * ```2_Single_cell_transitions_into_cell_type_transitions.ipynb```: Before we only saved the transition matrices. For evaluation we also need the cell type transitions, which is calculated here.
   * ```3_Evaluate_and_plot.ipynb```: Combine everything and make nice plots.
</details>


<details>
    <summary>2_PCA_and_scVI  </summary>
    &nbsp; 
    
In order to evaluate how resilient moscot and TOME are to batch effect we ran both models also on other embeddings than the Seurat integration.

#### 0_PCA:
  
Computes the PCA embedding (no integration) for all time points and runs TOME and moscot on it:
  
   * ```0_Comupte_PCA_for_TOME.ipynb```: Uses a standard scanpy procedure to cumpute a 30-dimensional PCA embedding.
   * ```1_Run_TOME_on_PCA.ipynb```: Runs TOME on this PCA embedding.
   * ```2_Run_moscot_on_PCA_embedding.ipynb```: Uses the same PCA embedding as used by TOME to run moscot on.
   
#### 1_scVI:
  
For a non-linear embedding we opted for scVI. We compute a scVI embedding for each of the 3 embrionic stages we defined:
  
   * ```0_Compute_stagewise_scVI.ipynb```: For each of the three stages (pre-gastrulation, gastrulation, and organogenesis) we concatenate all anndatas and compute a scVI embedding.
   * ```1_Run_TOME_on_scVI.ipynb```: Runs TOME using this scVI embedding.
   * ```2_Run_moscot_on_scVI_embedding.ipynb```: Runs moscot using this scVI embedding.

#### 2_Evaluate_and_plot:
  
Evaluates the mappings for the anchor- (=Seurat embedding), PCA- and scVI-embedding.
  
   * ```0_Evaluate_and_plot.ipynb```: Loads and evaluates the mappings for the three embeddings for TOME and moscot using germ-layer and curated transitions
   * ```1_Plot_scVI_3d_UMAP.ipynb```: For the E8.5 time point/pair we illustrate the evaluation results using a 3d UMAP (since TOME uses this 3d UMAP to calculate distances)

</details>

<details>
    <summary>3_Low_Rank </summary>
    &nbsp; 

We evaluated the runtime and performance of moscot for various choices of ranks for low-rank factorization.
    
   * ```0_Low_Rank_moscot.ipynb```: Calculates moscot mappings for all time points and different ranks and saves the evaluation results to weights-and-biases.
   * ```1_LR_wandb_plots.ipynb```: Loads the wandb-results and plots them.
</details>


<details>
    <summary>4_Growth_rates_for_IPS_reprogramming </summary>
    &nbsp; 

We had to switch to an IPS reprogramming dataset to show that growth rates that are predicted by moscot correlate with growth rates predicted using gene signatures since we did not have a working set of prolfieration/apoptosis genes spanning different cell types and developmental stages in mouse.

    
   * ```0_Compute_PCAs_and_UMAPs.ipynb```: Calculates PCA embedding for all analyzed time pairs (> day 8)
   * ```1_TOME_compute_growth_rates.ipynb```: Runs clTOME on the PCa embedding
   * ```2_TOME_Transforming_Identified_Neigbors_to_Transport_Matrix_to_growh_rates.ipynb```: Translate the custom TOME output to kNN growth rates
   * ```3_Run_moscot.ipynb```: Runs moscot on the IPS reprogramming PCA embeddings
   * ```4_Evaluate_and_plot.ipynb```: Combines all results and plots them
</details>


<details>
    <summary>5_Metacells </summary>
    &nbsp; 

Another way to speed up compuations is by combining cells into metacells. We show the shortcoming of this approach on the E9.5-E10.5-E11.5 - mouse atlas data:


   * ```0_Select_single_timepoint_anndatas.ipynb```: Saves the anndata objects for E9.5,E10.5, and E11.5
   * ```1_Metacells.ipynb```: I run Metacell-2[[2]] on all three time points.
   * ```2_Run_moscot_on_Metacells.ipynb```: I run moscot using the metacell aggregations.
   * ```3_Metacell_transition_matrix_to_single_cell_transition_matrix.ipynb```: I transform the metacell-transport matrix back to a single-cell transport matrix such that it can be compared to the conventional moscot results.
   * ```4_Cell_type_transition_accuracies.ipynb```: This accumulates the single-cell transport matrix into cell type transport matrix in order to compute germ-layer and curated transition accuracies.
   * ```5_Pancreatic_epithelium.ipynb```: We also want to evaluate how well the pull of pancreatic cells correlates to driver genes for the single-cell and the metacell cases.
   * ```6_Pancreatic_epithelium.ipynb```: This plots the result for the moscot maps on metacells but also illustrates on the E9.5 time point that rare cell types can easily be missed by the metacell aggregation.
</details>





[[1] Qiu et al.,  "Systematic reconstruction of cellular trajectories across mouse embryogenesis." Nat Genet (2022)](https://www.nature.com/articles/s41588-022-01018-x) <br> 
[[2] Schiebinger et al.,  "Optimal-Transport Analysis of Single-Cell Gene Expression Identifies Developmental Trajectories in Reprogramming." Cell (2019)](https://www.nature.com/articles/s41588-022-01018-x) <br> 
[[3]  Ben-Kiki et al.,  "Metacell-2: a divide-and-conquer metacell algorithm for scalable scRNA-seq analysis." Gen Biol (2022)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02667-1) 

