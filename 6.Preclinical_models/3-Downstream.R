# Load packages -----------------------------------
rm(list = ls())

library(data.table)
library(tidyverse)
library(crayon)
library(ape)
library(readxl)
library(ComplexHeatmap)
library(qs)
library(cowplot)
library(SingleR)

# Organize environment  -----------------------------------
base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/6 - scRNAseq models"
analysis_dir <- file.path(base_dir, 'analysis')
resource_dir <- file.path(base_dir, 'scripts/resources')

plot_dir <- file.path(analysis_dir, 'plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

data_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(data_dir)){dir.create(data_dir, recursive = T)}

source(file.path(resource_dir, "single_cell_preprocessing_helper_functions.R"))
source(file.path(resource_dir, 'Plot_style.R'))
source(file.path(resource_dir, 'NMF_helper_function.R'))

metadata <- read_excel(file.path(base_dir, 'metadata.xlsx')) 


# load Seurat object -----------------------------------
## Load count matrix of EPN models
  cm_raw <- pbapply::pblapply(file.path(base_dir, 'analysis/qc/data/individual', paste0(metadata$FileName, '-cm_list.rds')), function(x) {
    readRDS(x)$raw_data
  })
  
  ## Create Seurat object of EPN models
    seurat_obj_models <- list()
    for (i in seq_along(cm_raw)) {
      seurat_obj_models[[i]] <- CreateSeuratObject(cm_raw[[i]])
    }
    
  ## Add metadata
    for (i in seq_along(seurat_obj_models)) {
      seurat_obj_models[[i]]$cell_line <- metadata$CellLine[[i]]
      seurat_obj_models[[i]]$model <- metadata$Model[[i]]
      seurat_obj_models[[i]]$type <- metadata$Type[[i]]
      seurat_obj_models[[i]]$FileName <- metadata$FileName[[i]]
    }
    
## Merge model datasets
    seurat_obj_models_merged  <- merge(seurat_obj_models[[1]], y = c(seurat_obj_models[2:length(seurat_obj_models)]))

## Load seurat object of patients
  seurat_patients <- qread(file.path(base_dir, 'data/patient/seurat_obj_ZFTARELA_combined_malig.qs'))
  seurat_patients@meta.data <- seurat_patients@meta.data[ , c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'sample', 'Metaprogram_aggregate')]
  DefaultAssay(seurat_patients) <- 'RNA'
  
  # add metadata
  seurat_patients$model <- 'Patient'
  seurat_patients$type <- 'frozen'
  seurat_patients$cell_line <- seurat_patients$sample
  seurat_patients$FileName <- seurat_patients$sample

## Merge models and patient datasets
  seurat_obj_merged  <- merge(seurat_patients, y = seurat_obj_models_merged)
  # merge layers
  seurat_obj_merged <- JoinLayers(seurat_obj_merged)

## Reorder metadata info
  seurat_obj_merged$model <- factor(x = seurat_obj_merged$model, levels = c("Patient", "PDX", "adherent", "spheres"))

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
 
    # Remove the "orig.ident" metadata column
  seurat_obj_merged$sample<- NULL
  
  ## Reorder metadata info
  seurat_obj_merged$model <- factor(x = seurat_obj_merged$model, 
                                    levels = c("Patient", "PDX", "adherent", "spheres"))
  
  
  
# Score cells for metamodules -----------------------------------

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
  
  
# SingleR to transfer labels -------------------------------------------
  # transform into single cell object 
  seurat_obj.sce <- as.SingleCellExperiment(seurat_obj_merged)
  
  # load reference dataset - Human Primary Cell Atlas (Mabbott et al. 2013)
  ref <- as.SingleCellExperiment(seurat_patients)
  
  
  # use SingleR to transfer labels
  pred.seurat_obj <- SingleR(test=seurat_obj.sce, 
                             ref = ref, 
                             #clusters = seurat_obj_merged$seurat_clusters, 
                             labels = ref$Metaprogram_aggregate, 
                             assay.type.ref = 1,
                             prune = TRUE)
  
  # add predicted labels to Seurat object
  seurat_obj_merged[["SingleR.labels"]] <- pred.seurat_obj$pruned.labels
  
  
  
  
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
  colnames(seurat_obj_merged@meta.data)[27] <- 'Metaprogram NMF scoring'
  
  
## Save dataset -----------------------------------
  saveRDS(seurat_obj_merged, file.path(data_dir, "Patient_models_merged.rds"))
  
  # Export details about number of cells
  seurat_obj_mal_table <- seurat_obj_merged@meta.data %>% 
    dplyr::group_by(model, cell_line) %>% 
    dplyr::summarise(nCount_RNA = mean(nCount_RNA),
              nFeature_RNA = mean(nFeature_RNA),
              n = n()) 
  
  seurat_obj_mal_table
  write.csv(seurat_obj_mal_table, file.path(data_dir, "0_scRNAseq_details_STEPN.csv"))
  
  
