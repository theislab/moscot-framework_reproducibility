## Charting the spatiotemporal dynamics of mouse embryogenesis with moscot
a spatiotemporal atlas of mouse embryogenesis (MOSTA) over eight time points, from E9.5 to E16.5 [[1]](https://doi.org/10.1016/j.cell.2022.04.003). 

## Notebook folders


<details>
    <summary>0_Data_preparation  </summary>
    &nbsp; 
    
We use the data as provided in [[1]](https://doi.org/10.1016/j.cell.2022.04.003) downloaded from [MOSTA](https://db.cngb.org/stomics/mosta/) online db:
1. Full Embryo: [Mouse_embryo_all_stage.h5ad](https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/Mouse_embryo_all_stage.h5ad)
2. Brain: [E16.5_E1S3_cell_bin_whole_brain.h5ad](https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E1S3_cell_bin_whole_brain.h5ad)

#### Preprocessing notebooks:
  
   * ```ZP_2023-04-20_spatiotemporal_fullembryo-preprocess.ipynb```: Performs pre-processing for full embryo slides from all time points
   * ```ZP_2023-04-20_spatiotemporal_brain-preprocess.ipynb```: Combines and pre-processes brain cells from time points E9.5-E15.5 with the annotated cells at E16.5.

</details>


<details>
    <summary>1_Cell_type_transition_analysis  </summary>
    &nbsp; 

#### 0_grid_search:
  
  We compute the mapping between the time points using moscot _SpatioTemporalProblem_ and _TemporalProblem_. 
  For _SpatioTemporalProblem_ we run a grid search over the `alpha` values using SLURM.

  To run the grid search call within the directory:
    ```python3 run_mosta_st_map_grid.py```
  
   * ```run_st_map.sh```: Script to initialize sbatch runs. 
   * ```run_mosta_st_map_grid.py```: Main SLURM script to calculate the couplings between the time points.
   Calling the command `python3 run_mosta_st_map_grid.py` 
will instantiate sbatch calls to calculate all couplings over a range of `alpha` values. 
The accuracy of each mapping is saved as a `.csv` file under `data/space/spatiotemporal/output`.
   * ```mosta_st_map_accuracies.py```: Main function called by `run_mosta_st_map_grid.sh`. 
Evaluates the mapping for the give args and saves the accuracy of the mapping as a `.csv` file under `data/space/spatiotemporal/output`.
   * ```ZP_2023-04-20_spatiotemporal_fullembryo-accuracy.ipynb```: Imports the grid search run and visualizes the mapping accuracy

#### 1_mapping_across_timepoints:
  We compute the mapping between the time points using optimal `alpha` values found in `0_grid_search` and save them as `.pkl` files.
    To run the mapping call:
    ```python3 run_mosta_st_map_transitions.py```
  
   * ```run_st_map.sh```: Script to initialize sbatch runs. 
   * ```run_mosta_st_map_transitions.py```: Main SLURM script to calculate the couplings between the time points.
   Calling the command `python3 run_mosta_st_map_transitions.py` 
will instantiate sbatch calls to calculate couplings the optimal `alpha` values. 
The mapping as well as push forward of `Heart` cells will be saved as a `.pkl` files under `data/space/spatiotemporal/output`.
   * ```mosta_st_map_transitions.py```: Main function called by `run_mosta_st_map.sh`. 

</details>


<details>
    <summary>2_Full_embryo_CellRank_analysis  </summary>
    &nbsp;

* ```ZP_2023-04-20_spatiotemporal_fullembryo-heart.ipynb```: Imports the heart push forwards computed in 
`1_Cell_type_transition_analysis/1_mapping_across_timepoints` and uses them to visualize the heart cell mappings. 
* ```ZP_2023-04-20_spatiotemporal_fullembryo-cellrank.ipynb```: Imports the mappings computed in 
`1_Cell_type_transition_analysis/1_mapping_across_timepoints` and uses them to define a CellRank kernel. 

</details>


<details>
    <summary>3_Brain_analysis  </summary>
    &nbsp; 


#### 0_Brain_mapping:
  We compute the mapping between the brian cells across time points.
    To run the mapping call:
    ```python3 run_mosta_st_map_brain.py```
  
   * ```run_st_map.sh```: Script to initialize sbatch runs. 
   * ```run_mosta_st_map_brain.py```: Main SLURM script to calculate the couplings between the time points.
   Calling the command `python3 run_mosta_st_map_brain.py` 
The mapping will be saved as a `.pkl` files under `data/space/spatiotemporal/output_brain`.
   * ```mosta_st_map_brain.py```: Main function called by `run_mosta_st_map_brain.sh`. 
   * ```ZP_2023-04-20_spatiotemporal_brain-label.ipynb```: Imports the mappings and transfers the cell type 
annotation (labels cells at earlier time points). 

#### 1_CellRank_analysis:
* ```ZP_2023-04-20_spatiotemporal_brain-cellrank.ipynb```: Imports the mappings computed in 
`0_Brain_mapping` and uses them to define a CellRank kernel. 

</details>


[[1] Chen et al.,  "Spatiotemporal transcriptomic atlas of mouse organogenesis using DNA nanoball-patterned arrays." Cell (2022)](https://doi.org/10.1016/j.cell.2022.04.003) 