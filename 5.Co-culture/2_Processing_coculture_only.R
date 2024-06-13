# Load packages -----------------------------------
rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratObject)
library(data.table)
library(tidyverse)
library(readxl)
library(qs)

# Organize environment  -----------------------------------
base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/9 - scRNAseq coculture"

analysis_dir <- file.path(base_dir, 'analysis')
qc_dir <- file.path(analysis_dir, 'qc')

resource_dir <- file.path(base_dir, 'scripts/resources')

plot_dir <- file.path(analysis_dir, 'plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}
data_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(data_dir)){dir.create(data_dir, recursive = T)}


source(file.path(resource_dir, 'Plot_style_v2.R'))
source(file.path(resource_dir, 'single_cell_preprocessing_helper_functions.R'))
source(file.path(resource_dir, 'NMF_helper_function.R'))


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
    mutate(sample = ifelse(sample == 'EP1NSMonoCulture', 'Monoculture', 'Coculture'))
  seurat_obj <- AddMetaData(seurat_obj, metadata_coculture)

## Reorder metadata info
seurat_obj$sample <- factor(x = seurat_obj$sample, 
                                  levels = c("Monoculture", "Coculture"))

## Score cells for Metaprogram scores identified in patient cohort  -----------------------------------
# Load tumor signatures
ref.gene.EPN <- readRDS(file.path(base_dir, 'data/patient/nmf_marker_genes_final_annotated.rds'))
names(ref.gene.EPN) <- c("Neuroepithelial-like", "MES-like","Ependymal-like","Radial-glia-like","Cycling", "Embryonic-neuronal-like", "Neuronal-like", "Embryonic-like")

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



# Use NMF scoring to transfer labels -----------------------------------

cm <- seurat_obj[["RNA"]]$counts
cm_norm <- as.matrix(log2(cm/10+1))
cm_mean <- log2(Matrix::rowMeans(cm)+1)
cm_center <- cm_norm - rowMeans(cm_norm)

nmf_marker_genes_final <- ref.gene.EPN

nmf_score_final <- scoreNmfGenes(cm_center, cm_mean, nmf_marker_genes_final)
nmf_score_final <- t(nmf_score_final)

num_meta_program <- length(nmf_marker_genes_final)
nmf_score_final_t <- data.frame(nmf_score_final)
for (i in 1:3){
  nmf_score_final_t <- metagene_score_signature(nmf_score_final_t, num_meta_program, i)
}

nmf_score_final_t$signature_1 <- gsub('\\.', '-', nmf_score_final_t$signature_1)
nmf_score_final_t$signature_2 <- gsub('\\.', '-', nmf_score_final_t$signature_2)
nmf_score_final_t$signature_3 <- gsub('\\.', '-', nmf_score_final_t$signature_3)

meta <- nmf_score_final_t
metadata <- seurat_obj@meta.data
metadata <- cbind(metadata, meta)

seurat_obj <- AddMetaData(seurat_obj, metadata)

# replace names metaprograms-----------------------------------
colnames(seurat_obj@meta.data)[29] <- 'Metaprogram'

# Save dataset -----------------------------------
saveRDS(seurat_obj, file.path(data_dir, "Mono_co_culture_merged.rds"))






# Combine with patient dataset  -----------------------------------
colnames(seurat_obj@meta.data)[4] <- 'model'

## Load seurat object of patients
seurat_patients <- qread(file.path(base_dir, 'data/patient/seurat_obj_ZFTARELA_combined_malig.qs'))
seurat_patients@meta.data <- seurat_patients@meta.data[ , c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'sample', 'Metaprogram_aggregate')]
DefaultAssay(seurat_patients) <- 'RNA'

# add metadata
seurat_patients$model <- 'Patient'


## Merge models and patient datasets
seurat_obj_merged  <- merge(seurat_patients, y = seurat_obj)
# merge layers
seurat_obj_merged <- JoinLayers(seurat_obj_merged)

## Reorder metadata info
seurat_obj_merged$model <- factor(x = seurat_obj_merged$model, levels = c("Patient", "Monoculture", "Coculture"))

## Extract metadata information
metadata_seurat_obj_merged <- seurat_obj_merged@meta.data


# Process datasets patient and model  -----------------------------------
## Process patient and model dataset with full seurat
cm <- seurat_obj_merged[["RNA"]]$counts
# extract sample name from cell name
CellName <- colnames(cm)
samples <- sub('\\..*', '', CellName)

# Run Full Seurat 
seurat_obj <- RunFullSeurat_v5(cm, samples, doBatch = F, 
                               resolution = c(0.3, 1), dims = 0.8, project = 'ZFTARELA', 
                               norm.type = 'RNA', verbose = T)
seurat_obj_merged <- seurat_obj
# Add original metadata information
seurat_obj_merged <- AddMetaData(seurat_obj_merged, metadata_seurat_obj_merged)  


## Reorder metadata info
seurat_obj_merged$model <- factor(x = seurat_obj_merged$model, 
                                  levels = c("Patient", "Monoculture", "Coculture"))





## Score cells for Metaprogram scores identified in patient cohort  -----------------------------------
# Load tumor signatures
ref.gene.EPN <- readRDS(file.path(base_dir, 'data/patient/nmf_marker_genes_final_annotated.rds'))
names(ref.gene.EPN) <- c("Neuroepithelial-like", "MES-like","Ependymal-like","Radial-glia-like","Cycling", "Embryonic-neuronal-like", "Neuronal-like",            "Embryonic-like")

# score tumors
seurat_obj_merged <- AddModuleScore(object = seurat_obj_merged, assay = 'RNA', features = ref.gene.EPN, name = names(ref.gene.EPN))
# rename metadata names of scores
  # identify number of first column with metadata scores
  col_start <- length(colnames(seurat_obj_merged@meta.data)) - length(names(ref.gene.EPN)) +1
  # identify number of last column with metadata scores
  col_end <- length(colnames(seurat_obj_merged@meta.data))
  # rename columns with score name
  colnames(seurat_obj_merged@meta.data)[col_start:col_end] <- names(ref.gene.EPN)

# Cell cycle scoring -----------------------------------
  ## Cell cycle scoring
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  seurat_obj_merged <- CellCycleScoring(seurat_obj_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  ## Define high-cycling vs low-cycling cells
  Cycling <- WhichCells(seurat_obj_merged, expression = S.Score > 0 | G2M.Score > 0)
  nonCycling <- WhichCells(seurat_obj_merged, expression = S.Score <= 0 & G2M.Score <= 0)
  
  # Add cycling info to metadata column called 'Cycling_prop'
  seurat_obj_merged$Cycling_prop <- ifelse(colnames(seurat_obj_merged) %in% Cycling, "Cycling", "Non-cycling")
  
  # scale data
  seurat_obj_merged <- ScaleData(seurat_obj_merged)
  
  
  
# Use NMF scoring to transfer labels -----------------------------------
  
  cm <- seurat_obj_merged[["RNA"]]$counts
  cm_norm <- as.matrix(log2(cm/10+1))
  cm_mean <- log2(Matrix::rowMeans(cm)+1)
  cm_center <- cm_norm - rowMeans(cm_norm)
  
  nmf_marker_genes_final <- ref.gene.EPN
  
  nmf_score_final <- scoreNmfGenes(cm_center, cm_mean, nmf_marker_genes_final)
  nmf_score_final <- t(nmf_score_final)
  
  num_meta_program <- length(nmf_marker_genes_final)
  nmf_score_final_t <- data.frame(nmf_score_final)
  for (i in 1:3){
    nmf_score_final_t <- metagene_score_signature(nmf_score_final_t, num_meta_program, i)
  }
  
  nmf_score_final_t$signature_1 <- gsub('\\.', '-', nmf_score_final_t$signature_1)
  nmf_score_final_t$signature_2 <- gsub('\\.', '-', nmf_score_final_t$signature_2)
  nmf_score_final_t$signature_3 <- gsub('\\.', '-', nmf_score_final_t$signature_3)
  
  meta <- nmf_score_final_t
  metadata <- seurat_obj_merged@meta.data
  metadata <- cbind(metadata, meta)
  
  seurat_obj_merged <- AddMetaData(seurat_obj_merged, metadata)
  
  # replace names metaprograms-----------------------------------
  colnames(seurat_obj_merged@meta.data)[44] <- 'Metaprogram-patient'

# Save dataset -----------------------------------
saveRDS(seurat_obj_merged, file.path(data_dir, "Mono_co_culture_patient_merged.rds"))
    
