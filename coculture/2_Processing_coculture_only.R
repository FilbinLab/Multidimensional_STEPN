# Load packages -----------------------------------
rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratObject)
library(data.table)
library(tidyverse)
library(future)
library(clustree)
library(clusterProfiler)
library(SCpubr)
library(readxl)
library(qs)
library(SingleR)

# Organize environment  -----------------------------------
base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/9 - scRNAseq coculture"

analysis_dir <- file.path(base_dir, 'analysis')
qc_dir <- file.path(analysis_dir, 'qc')

resource_dir <- file.path(base_dir, 'scripts/resources')

plot_dir <- file.path(analysis_dir, 'plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}
data_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(data_dir)){dir.create(data_dir, recursive = T)}


source(file.path(resource_dir, 'Plot_style.R'))
source(file.path(resource_dir, 'single_cell_preprocessing_helper_functions.R'))

# Define color palette
col_model <- c( '#A7DBD8FF', '#F38630FF')


# load count matrices, create Seurat object and add annotation as metadata -----------------------------------

# load count matrices
cm <- do.call(cbind, lapply(list.files(qc_dir, pattern = '-cm_list.rds', full.names = T), function(x) {
  tmp <- readRDS(x)
  tmp$raw_data
}))

  # extract sample name from cell name
  CellName <- colnames(cm)
  samples <- sub('\\..*', '', CellName)

# Run Full Seurat 
seurat_obj <- RunFullSeurat_v5(cm, samples, doBatch = F, 
                 #var2batch = 'sample', batchMethod = 'harmony', 
                 resolution = c(0.3, 1), dims = 0.8, project = 'TEST', 
                 norm.type = 'RNA', verbose = T)


# rename name of mono and coculture in the metadata column metadata coculture
metadata_coculture <- seurat_obj@meta.data
  # rename name of mono and coculture in the metadata column metadata coculture
  metadata_coculture <- metadata_coculture %>% 
    mutate(model = ifelse(model == 'EP1NSMonoCulture', 'Monoculture', 'Coculture'))
  seurat_obj <- AddMetaData(seurat_obj, metadata_coculture)

## Reorder metadata info
seurat_obj$model <- factor(x = seurat_obj$model, 
                                  levels = c("Monoculture", "Coculture"))



## Score cells for EPN (ZFTA::RELA) scores identified in patient cohort 
# Load tumor signatures
ref.gene.EPN <- qread(file.path(resource_dir, "ZFTA_frozen_DE_list_signature1_MergedNPCs2.qs"))
# rename
names(ref.gene.EPN) <- c('NPC_like', 'Mesenchymal', 'Radial_glia', 'Ependymal', 'Neuroepithelial', 'Cycling')
# select top 100 genes per program
ref.gene.EPN <- lapply(ref.gene.EPN, select_top_100)



# score tumors
seurat_obj <- AddModuleScore(object = seurat_obj, assay = 'RNA', features = ref.gene.EPN, name = names(ref.gene.EPN))
# rename metadata names of scores
  # identify number of first column with metadata scores
  col_start <- length(colnames(seurat_obj@meta.data)) - length(names(ref.gene.EPN)) +1
  # identify number of last column with metadata scores
  col_end <- length(colnames(seurat_obj@meta.data))
  # rename columns with score name
  colnames(seurat_obj@meta.data)[col_start:col_end] <- names(ref.gene.EPN)

# Cell cycle scoring -----------------------------------
  ## Cell cycle scoring
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  ## Define high-cycling vs low-cycling cells
  Cycling <- WhichCells(seurat_obj, expression = S.Score > 0 | G2M.Score > 0)
  nonCycling <- WhichCells(seurat_obj, expression = S.Score <= 0 & G2M.Score <= 0)
  
  # Add cycling info to metadata column called 'Cycling_prop'
  seurat_obj$Cycling_prop <- ifelse(colnames(seurat_obj) %in% Cycling, "Cycling", "Non-cycling")
  
  # scale data
  seurat_obj <- ScaleData(seurat_obj)
  
  
  
# Use SingleR to transfer labels
    ## SingleR
    # transform into single cell object 
    seurat.sce <- as.SingleCellExperiment(seurat_obj)
    
    # load reference dataset 
    ref <- qread(file.path(resource_dir, "seurat_obj_ZFTA_Frozen_mergedNPCs.qs"))
    ref.sce <- as.SingleCellExperiment(ref)
    
    # use SingleR to transfer labels
    pred_seurat <- SingleR(test=seurat.sce, 
                           ref = ref.sce, 
                           labels=ref$Metaprogram, de.n=50)
    # visualize results
    table <- table(pred_seurat$labels)
    
    # add predicted labels to Seurat object
    seurat_obj[["SingleR.labels"]] <- pred_seurat$labels
    
# Save dataset
saveRDS(seurat_obj, file.path(data_dir, "Mono_co_culture_merged.rds"))
    
