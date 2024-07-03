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
* ```ZP_2023-04-20_spatiotemporal_fullembryo-accuracy.ipynb```: Imports the grid search run and visualizes the mapping accuracy.

  ##### 0_moscot
  We compute the mapping between the time points using moscot _SpatioTemporalProblem_ and _TemporalProblem_. 
  For _SpatioTemporalProblem_ we run a grid search over the `alpha` values using SLURM.

  To run the grid search call within the directory:
    ```sbatch mosta_st_map_array.py```

The script requires: 
   * ```grid_config.txt```: Config file with tasks to launch 

Calling the command `sbatch mosta_st_map_array.py` 
will instantiate sbatch calls to calculate all couplings over a range of `alpha` values. 
The accuracy of each mapping is saved as a `.csv` file under `data/space/spatiotemporal/output`.
   * ```mosta_st_map_accuracies.py```: Main function called by `mosta_st_map_array.sh`. 
Evaluates the mapping for the give args and saves the accuracy of the mapping as a `.csv` file under `data/space/spatiotemporal/output`.

##### 1_TOME
 We compute the mapping using TOME[[2]](https://doi.org/10.1038/s41588-022-01018-x) the time points and save them as `.pkl` files.
    To run the mapping call:
    ```python3 run_mosta_st_map_tome.py```
  
   * ```mosta_st_map.sh```: Script to initialize sbatch runs. 
   * ```run_mosta_st_map_tome.py```: Main SLURM script to calculate the couplings between the time points.
   Calling the command `python3 run_mosta_st_map_tome.py` 
will instantiate sbatch calls to calculate couplings between all timepoints
The mapping accuracy will be saved under `data/space/spatiotemporal/output_tome`.
   * ```tome.py```: Main function called by `run_mosta_st_map.sh`. 
  
##### 2_PASTE2
 We compute the mapping using PASTE2[[2]](https://www.genome.org/cgi/doi/10.1101/gr.277670.123) the time points and save them as `.pkl` files.
    To run the mapping call:
    ```python3 run_mosta_st_map_paste2.py```
  
   * ```mosta_st_map.sh```: Script to initialize sbatch runs. 
   * ```run_mosta_st_map_paste2.py```: Main SLURM script to calculate the couplings between the time points.
   Calling the command `python3 run_mosta_st_map_tome.py` 
will instantiate sbatch calls to calculate couplings between all timepoints
The mapping accuracy will be saved under `data/space/spatiotemporal/output_paste`.
   * ```run_paste2.py```: Main function called by `run_mosta_st_map.sh`.

Files:
- requires `MOSTA_curated_transitions.csv` (Supplementary table 5 in Klein et al. (2023))

#### 1_mapping_across_timepoints:
  We compute the mapping between the time points and save them as `.pkl` files.
    To run the mapping call:
    ```python3 run_mosta_st_map_transitions.py```
  
   * ```run_st_map.sh```: Script to initialize sbatch runs. 
   * ```run_mosta_st_map_transitions.py```: Main SLURM script to calculate the couplings between the time points.
   Calling the command `python3 run_mosta_st_map_transitions.py` 
will instantiate sbatch calls to calculate couplings the given `alpha` values. 
The mapping as well as push forward of `Heart` cells will be saved as a `.pkl` files under `data/space/spatiotemporal/output_heart`.
   * ```mosta_st_map_transitions.py```: Main function called by `run_mosta_st_map.sh`. 
   * ```ZP_2023-04-20_spatiotemporal_fullembryo-heart.ipynb```: Imports the heart push forwards computed in 
`1_Cell_type_transition_analysis/1_mapping_across_timepoints` and uses them to visualize the heart cell mappings.

#### 2_AB_BC_CA_analysis:

</details>

<details>

<summary>2_Full_embryo_moscot_analysis  </summary>

We compute the mapping between the time points and save them as `.pkl` files.
    To run the mapping call:
    ```sbatch mosta_st_map_array_save.sh```
The script requires: 
   * ```save_config.txt```: Config file with tasks to launch 
   * ```mosta_st_map_accuracies_save.py```:  Main function called by `mosta_st_map_array_save.sh`.
 Calling the command ``sbatch mosta_st_map_array_save.sh```
will instantiate sbatch calls as defined in the config grid. 
The SpatioTemporal problems as well as push forward of `Liver` cells will be saved as a `.pkl` files under `data/space/spatiotemporal/output_liver`.
* ```ZP_2023-12-12_spatiotemporal_fullembryo-liver.ipynb``` Imports the SpatioTemporal problems and computes features correlation with Liver cells. 


</details>

<details>
    <summary>3_Full_embryo_CellRank_analysis  </summary>
    &nbsp;

* ```ZP_2023-04-20_spatiotemporal_fullembryo-cellrank.ipynb```: Imports the mappings computed in 
`1_Cell_type_transition_analysis/1_mapping_across_timepoints` and uses them to define a CellRank kernel. 

</details>


<details>
    <summary>4_Brain_analysis  </summary>
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

[[2] Qiu, C., Cao, J., Martin, B.K. et al. "Systematic reconstruction of cellular trajectories across mouse embryogenesis". Nat Genet 54, 328–341 (2022)](https://doi.org/10.1038/s41588-022-01018-x)

[[3] Liu, Xinhao, Ron Zeira, and Benjamin J. Raphael. "Partial alignment of multislice spatially resolved transcriptomics data." Genome Research 33.7 (2023)](https://www.genome.org/cgi/doi/10.1101/gr.277670.123)