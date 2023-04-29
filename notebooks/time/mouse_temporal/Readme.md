## Moscot enables atlas-scale temporal mapping of mouse embryogenesis
We applied moscot to a 1.7 M. cells mouse developental atlas. The data was preprocessed, annotated and partially created/ sequenced deeper by Qiu et al. (https://www.nature.com/articles/s41588-022-01018-x)
For reproducibility, the analysis was divided into 3 folders:


<details>
    <summary>0_Data_preparation:  </summary>
    The data was downloaded from http://tome.gs.washington.edu/ as .RDS files.
    
   ```MG_XXX_Integration```: Performs Seurat's anchor based batch correction analogous to Qiu et al, using using code obtianed from https://github.com/ChengxiangQiu/tome_code
  
   ```MG_XXX_Running_TOME```: Runs TOME using code obtianed from https://github.com/ChengxiangQiu/tome_code. In addittion saving the cell type transitions, we also save the identified neirest neighbors, which will be used in cell-level (cl)TOME
  
   ```MG_XXX_RDS_to_anndata```: Contains a notebooks used to transform .RDS objects to anndata objects using SeruatData/Disk. Some annotations are not carried over correctly, which is fixed by running in the notebooks ```MG_06-26-2022_Fix_anndata_annotations.ipynb```


    
</details>

<details>
    <summary>Data download</summary>
    The data was downloaded from http://tome.gs.washington.edu/
</details>

## Rules
1. Have an **empty line** after the `</summary>` tag or markdown/code blocks will not render.
1. Have an **empty line** after each `</details>` tag if you have multiple collapsible sections.
