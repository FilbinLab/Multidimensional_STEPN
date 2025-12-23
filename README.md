# Multidimensional_STEPN

Code repository for single-cell multidimensional profiling of supratentorial ependymomas, preprint: https://www.biorxiv.org/content/10.1101/2024.08.07.607066v1

Processed data for sc/nRNAseq and 10X Xenium are deposited on gene expression omnibus (GEO). 
scRNA-seq: GSE300150
Spatial: GSE300146

Raw sc/snRNAseq data are deposited in EGA (EGAS50000001513). 

Contact: djeong3 at mgh.harvard.edu, sara.g.danielli at gmail.com, carlosa_biagijunior at dfci.harvard.edu

# Overview
0. Software version requirements 
1. sc/snRNA-seq Preprocessing
2. Plots
3. Pseudobulk analysis (projections & plots)
4. Xenium (preprocessing, spatial coherence, niche, neighborhood enrichment, cellchat)
5. Co-culture (preprocessing, analysis, plots)
6. Preclinical_models (preprocessing, analysis, plots)
7. ShinyApp

## 0. Software version requirements 

All code was run in R version 4.3.1, Python v3.10.11, and Seurat v.5.0.2. Other packages include Hisat2 (v2.1.0), RSEM (v1.3.0), infercnv (v1.8.0), SingleR (v1.6.1), PAGODA2 (v0.1.4), NMF (v0.28), edgeR (v0.27), clusterProfiler (v4.6.2), org.Hs.eg.db (v3.18.0), RColorBrewer (v1.1-3), paletteer (v1.6.0), miceadds (v3.17-44), mice (v3.16.0), progress (v1.2.3), UCell (v2.8.0), yaml (v2.3.10), crayon (v1.5.3), patchwork (v1.3.0), readxl (v1.4.5), sp (v2.2-0), qs (v0.27.3), glue (v1.8.0)
lubridate (v1.9.4), forcats (v1.0.0), stringr (v1.5.1), dplyr (v1.1.4), purrr (v1.0.4), readr (v2.1.5), tidyr (v1.3.1), tibble (v3.2.1), ggplot2 (v3.5.2), tidyverse (v2.0.0), car (v3.1.3), scanpy (v1.11.2), scvi (1.3.1post1), squidpy (v1.6.5), cellcharter (v0.3.4), pandas (v2.2.3), numpy (v1.24.2), scVelo (v0.3.3), velocyto (v0.17.17)

## 1. sc/snRNA-seq Preprocessing
Codes to process sc/snRNA-seq data from ST-EPN patient samples. Codes are ordered logically, from QC-filtering (`1-filtering.R`) to inferCNV to determine malignant/non-malignant cells (the code is different for fresh and frozen samples - `3a-infercnv_frozen_step1_part1.R`, `3b-infercnv_fresh_step1_part1.R`, `4a-infercnv_frozen_step2.R`, `4b-infercnv_fresh_step2.R`, `5a-infercnv_frozen_step3.R`, `5b-infercnv_fresh_step3.R`), to preparing counts for NMF (`6-prepare_nmf_counts.R` and `6b-prepare_nmf_counts_fresh.R`), to running NMF (`7-NMF_rank-forO2.R`) and visualizing NMF results (`8-NMF.Rmd`)

Final datasets are cleaned up, annotated and reprocessed using `9_Cleanup_datasets.R` and `10_classification_malignant_normal.R`

Code to display oncoprint of Fig 1 is also in this folder (`Oncoprint.R`)

## 2. Plots
Codes to plot most figures for ST-EPN patient samples.

## 3. Pseudobulk analysis 
Codes to plot Extended Data Fig 3 and developmental projections of Fig 1f and Extended Data Fig 2a-b. To prepare the reference datasets, the codes stored in the folder `Developmental_dataset_processing` are used.

## 4. 10X Xenium

Code to perform the spatial data analysis. It contains following scripts:
- `0_preparation_scRNAseq_data.R`: prepare single cell objects for label transfer onto spatial data for annotations.
- `1_preprocessing.R`: load and preprocess Xenium data.
- `2_assign_programs.R`: perform annotation of the preprocessed Xenium data.
- `3_NicheAnalysis.R`: perform niche analysis on the data.
- `4_coherence.R`: perform coherence analysis on the data.

It also contains following folders:
- `resources`: contains metadata, and helper functions to process the data and generate results.
- `cellcharter`: contains scripts to run the CellCharter analysis and plot results.

## 5. Co-culture

Codes to pre-process and analyze patient models cocultured with rat E19 cortical cells and human iNs+iAs. 

## 6. Preclinical_models

Codes to pre-process patient models 

## 7. ShinyApp

ShinyApp under the ```epnApp``` folder contains all sc/sn-RNA-seq data, providing an accessible platform for data exploration including findings from the manuscript. Installation and execution instructions are conveniently located within the folder.
