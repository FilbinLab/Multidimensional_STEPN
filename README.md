# Multidimensional_STEPN


1. sc/snRNA-seq Preprocessing
2. Plots
3. Pseudobulk analysis (projections & plots)
4. Xenium (preprocessing, spatial coherence, niche, neighborhood enrichment, cellchat)
5. Co-culture (preprocessing, analysis, plots)  



## 1. sc/snRNA-seq Preprocessing
Codes to process sc/snRNA-seq data from ST-EPN patient samples. Codes are ordered logically, from QC-filtering (`1-filtering.R`) to inferCNV to determine malignant/non-malignant cells (the code is different for fresh and frozen samples - `3a-infercnv_frozen_step1_part1.R`, `3b-infercnv_fresh_step1_part1.R`, `4a-infercnv_frozen_step2.R`, `4b-infercnv_fresh_step2.R`, `5a-infercnv_frozen_step3.R`, `5b-infercnv_fresh_step3.R`), to preparing counts for NMF (`6-prepare_nmf_counts.R` and `6b-prepare_nmf_counts_fresh.R`), to running NMF (`7-NMF_rank-forO2.R`) and visualizing NMF results (`8-NMF.Rmd`)

Final datasets are cleaned up, annotated and reprocessed using `9_Cleanup_datasets.R` and `10_classification_malignant_normal.R`

Code to display oncoprint of figure 1 is also there (`Oncoprint.R`)


## 2. Plots
Codes to plot most paper figures from ST-EPN patient samples

## 3. Pseudobulk analysis 
Codes to plot Extended figure 3 and to plot developmental projections of Fig. 1f and Fig E2a. To prepare the reference datasets, the codes stored in the folder `Developmental_dataset_processing` are used.

## 4. Xenium
The folder contains the following subfolders:
`1-preprocessing`: to annotate Xenium cells with tumor cell states/types
`2-plots`: to display spatial maps
`3-neighborhood`: Python code to calculate cellular neighborhoods
`4-spatial_coherence`: code to display the spatial coherence results


## 5.Co-culture
Experiment to determine the cell state shifts between EP1NS cells cultures alone (monoculture) or with rat neurons and astrocytes (coculture).
In total, 3x4 plates coculture and 1x4 plates monoculture were sequenced. 1x run of coculture was bad QC, and therefore removed.
The remaining goodQC runs are stored under the following folders:

231207 - CoCulture, 4 plates `/labs/mfilbin/Demultiplexing/231207`
231116 - CoCulture, 4 plates `/labs/mfilbin/homes/biagi/Demultiplexing/231116` 
231103 - MonoCulture, 4 plates `/labs/mfilbin/homes/biagi/Demultiplexing/231103`

## 6.Preclinical_models
Codes to visualize FIg E10a/b