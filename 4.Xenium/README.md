# Xenium analysis

## Requirements

#### R information

```
R version 4.4.2 (2024-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Red Hat Enterprise Linux 9.5 (Plow)
```

#### Packages
```
 [1] RColorBrewer_1.1-3 paletteer_1.6.0    miceadds_3.17-44  
 [4] mice_3.16.0        progress_1.2.3     UCell_2.8.0       
 [7] yaml_2.3.10        crayon_1.5.3       patchwork_1.3.0   
[10] readxl_1.4.5       Seurat_5.3.0       SeuratObject_5.1.0
[13] sp_2.2-0           qs_0.27.3          glue_1.8.0        
[16] lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1     
[19] dplyr_1.1.4        purrr_1.0.4        readr_2.1.5       
[22] tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.2     
[25] tidyverse_2.0.0 
```

## Set Up

- In each script in the processing of Xenium data, we load the file in `resources/epn_functions.R`, and the functions in `resources/FLXenium`. User would need to change the paths to these files according to their filesystem in `0_preparation_scRNAseq_data.R`, `1_preprocessing.R` etc. so that they are loaded successfully.

- Secondly, the user would need to go to the file `resources/epn_functions.R` and change the default values of the arguments `home_dir` and `data_dir` in functions `SetUpEpendymomaGlobalVars` and `SetUpEpendymomaGlobalVarsGeneral`. home_dir point to `Multidimensional_STEPN/4.Xenium`, and the data_dir should point to the directory where the user wants the analyses results to be stored after processing. data_dir also indicates where the raw files are stored (`{data_dir}/raw_data/xenium_folders`, `{data_dir}/analysis/1_preparation/data/seurat_obj_ST_normal_malig_annotated.qs`, `{data_dir}/analysis/1_preparation/data/seurat_obj_malignant_annotated2.qs`)

## Processing 

### 0. Preparation of single-cell object from ZFTA-RELA patients for label transfer
Script `0_preparation_scRNAseq_data.R` subsets the sc/snRNA-seq object to ZFTA-RELA patients only (includes both malignant and non-malignant cells) and applies sctransform normalization.

### 1. Pre-processing of 10X Xenium output files
Script `1_preprocessing.R` preprocesses the output files from 10X Xenium by subsetting the files to cells with nCount>0, running SCTransform-based normalization, PCA and UMAP dimensionality reduction. 

### 2. Assign cell identities
Script `2_assign_programs.R` takes the sc/snRNA-seq object identities, and projects them onto the Xenium data. It then scores each cell for the identities used to build the Xenium panel, and assign a final cell identity to the Xenium cells.

### 3. Niche analysis
Script `3_NicheAnalysis.R` performs spatial niche analysis.

### 4. Coherence analysis
Script `4_coherence.R` performs the coherence analysis and saves the results.

### 5. Figures
Scripts to plot final results can be found in the `figures` folder.

### 6. CellCharter analysis 
All scripts to run the CellCharter analysis and plot results can be found in the `cellcharter` folder. 

CellCharter analysis is performed in three steps using python.
(1) Preprocessing by `scanpy`: `ST_EPN_scanpy.ipynb`. 
(2) Dimension reduction by `scVI`: `ST_EPN_scvi.ipynb`. 
(3) cellcharter analysis: `ST_EPN_cellcharter_k26.ipynb`.

Figures are generated using R.
(1) Figure-4C: Proportions of metaprograms of spatial clusters (`SpatialCluster_Metaprogram_Composition.Rmd`)
(2) Figure-4D: Heatmap of spatial clusters propotions of samples (`Sample_Spatial_Cluster_Heatmap.Rmd`)
(3) Figure-E6D: Proportions of spatial clusters of samples (`Sample_Spatial_Cluster_Heatmap.Rmd`)
(4) Figure-E6F: Pairwise correlation in metaprogram compositions between CellCharter and Meta-niche analysis (`SpatialCluster_MetaNiche_Correlation.Rmd`)
(5) Figure-E6G: Comparing entropies of sample composition between Tumor-enriched and TME-enriched spatial clusters (`Tumor_TME_Entropy_Comparison.Rmd`)



