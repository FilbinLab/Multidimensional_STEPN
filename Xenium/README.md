## Xenium 

- Raw data (fresh from Xenium machine) copied in Sara's folder, under `/n/scratch/users/s/sad167/EPN/Xenium/data/raw_data`
- Preliminary analysis from Carlos (clustering) copied in Sara's folder, under `/n/scratch/users/s/sad167/EPN/Xenium/data/processed_data_Carlos`

1 - Run scripts called `0_Xenium_rename_idents_` followed by date of run. The scripts contain the program annotations and creates preliminary graphs.
* *20231020: 0010652__Region_4
* *20231102: '0010575__Region_1', '0010575__Region_2', '0010575__Region_3', 
            '0010619__Region_1', '0010619__Region_2', '0010619__Region_3', '0010619__Region_4', '0010619__Region_5'
* *20231107: '0010501-Region_1', '0010501-Region_2'
* *20231109: '0010498__Region_1', '0010775__Region_1'
* *20231208 (samples used for postXenium-IF): '0010540__Region_1', '0010540__Region_2', '0010540__Region_3', '0010540__Region_4', 
                      '0010553__Region_1', '0010553__Region_2', '0010553__Region_3', '0010553__Region_4'

2 - Run scripts `1_Xenium_metaprogram_distribution.R` and `2_Xenium_images_distribution.R` to plot the barplot containing metaprogram frequency, respectively image plots of tumors colored by metaprogram annotation.

3- Run spatial calculations (clustering coefficients, neighborhood analyses) for each individual tumor running Python scripts called `3_Sara_neighbors_` followed by date of run.

4 - Put together results by averaging across tumor sections running Python script `4_Neighbors_average.ipynb`. To plot the neighborhood analysis heatmap, use R script `5_Neighbors_average_heatmap.R`



- Parameters for Jupyter session in O2:
  * *Modules to be preloaded*: `gcc/9.2.0 python/3.9.14 cmake/3.14.1 R/4.3.1`
  * *Jupyter Environment*: `source /home/cao385/envs/jupytervenv/bin/activate`
  * *Total Job Memory in GB for the job*: Usually `48` is enough
- To run niches analysis using Seurat package (click [here](https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#mouse-brain-10x-genomics-xenium-in-situ:) for more informations):
```