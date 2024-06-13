# Load packages -----------------------------------
library(tidyverse)
library(ggpubr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(writexl)
library(paletteer)
library(readxl)
library(writexl)
library(qs)
library(ggrastr)
library(cowplot)
library(openxlsx)
library(SingleR)

# Organize environment  -----------------------------------
base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/11 - Plots scRNAseq"

data_dir <- file.path("/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/16 - scRNAseq processing/analysis")
analysis_dir <- file.path(base_dir, 'analysis')
plot_dir <- file.path(analysis_dir, "plots/STEPN")
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}
result_dir <- file.path(analysis_dir, "result")
if (!dir.exists(result_dir)){dir.create(result_dir, recursive = T)}

resource_dir <- file.path(base_dir, 'scripts/resources')
source(file.path(resource_dir, 'single_cell_preprocessing_helper_functions_CBJr.R'))

# load seurat objects  -----------------------------------

#frozen: malignant and normal
seurat_obj <- readRDS(file = file.path(data_dir, "infercnv/frozen/part_1/seurat_obj.rds"))
  # subset by removing reference
  Idents(seurat_obj) <- 'type'
  seurat_obj <- subset(seurat_obj, idents = "obs")

#frozen: malignant 
  seurat_obj_mal <- qread(file = file.path(data_dir, "NMF/data/seurat_obj_annotated.qs"))

#fresh: malignant
  seurat_obj_fresh <- readRDS(file = file.path(data_dir, "infercnv/fresh/seurat_obj_malignant_annotated.rds"))
  
  
# load malignant vs non-malignant annotation used for NMF
seurat_obj[['malignant']] <- ifelse(Cells(seurat_obj) %in% Cells(seurat_obj_mal), 'Malignant', 'Non-malignant')

# reprocess normal+malignant object with RUnFullseurat -----------------------------------
metadata <- seurat_obj@meta.data

cm <- seurat_obj[["RNA"]]$counts
cm_norm <- as.matrix(log2(cm/10+1))
cm_mean <- log2(Matrix::rowMeans(cm)+1)
cm_center <- cm_norm - rowMeans(cm_norm)
seurat_obj <- RunFullSeurat_v5(cm = cm, metadata = metadata,  doBatch = T, var2batch = 'sample', batchMethod = 'harmony', project = 'HOPE')


# load metadata
metadata <- read_xlsx(file.path(base_dir, 'data/metadata/metadata_ST_EPN_Patient.xlsx'))
metadata <- subset(metadata, metadata$`sc/snRNAseq` == 'snRNA-seq')
FileName <- unique(metadata$FileName)
index <- match(FileName, metadata$FileName)
metadata <- metadata[index, ]




# add metadata to object with frozen malignant cells only -----------------------------------
metadata <- read_xlsx(file.path(base_dir, 'data/metadata/metadata_ST_EPN_Patient.xlsx'))
metadata <- subset(metadata, metadata$`sc/snRNAseq` == 'snRNA-seq')
FileName <- unique(metadata$FileName)
index <- match(FileName, metadata$FileName)
metadata <- metadata[index, ]


subset_metadata_seurat <- list()
metadata_seurat <- seurat_obj_mal@meta.data
for (i in 1:length(FileName)) {
  subset_metadata_seurat[[i]] <- metadata_seurat[metadata_seurat$sample == FileName[[i]], ]
  subset_metadata_seurat[[i]]$Sample_deID <- metadata$Sample_deID[[i]]
  subset_metadata_seurat[[i]]$PatientID <- metadata$PatientID[[i]]
  subset_metadata_seurat[[i]]$Subtype <- metadata$Subtype[[i]]
  subset_metadata_seurat[[i]]$Age <- metadata$Age[[i]]
  subset_metadata_seurat[[i]]$Gender <- metadata$Gender[[i]]
  subset_metadata_seurat[[i]]$Sampling <- metadata$Sampling[[i]]
  subset_metadata_seurat[[i]]$Fusion <- metadata$Fusion[[i]]
  subset_metadata_seurat[[i]]$sequencing <- metadata$`sc/snRNAseq`[[i]]
}

# merge metadata information
merged_metadata_seurat <- subset_metadata_seurat[[1]]

for (i in 2:length(subset_metadata_seurat)) {
  merged_metadata_seurat <- rbind(merged_metadata_seurat, subset_metadata_seurat[[i]])
}

# Add metadata to original seurat object
seurat_obj_mal <- AddMetaData(seurat_obj_mal, merged_metadata_seurat)


# add metadata to object with malignant+nonmalignant -----------------------------------
metadata <- read_xlsx(file.path(base_dir, 'data/metadata/metadata_ST_EPN_Patient.xlsx'))
metadata <- subset(metadata, metadata$`sc/snRNAseq` == 'snRNA-seq')
FileName <- unique(metadata$FileName)
index <- match(FileName, metadata$FileName)
metadata <- metadata[index, ]

# first remove mg, od and BT

metadata_seurat <- seurat_obj@meta.data
for (i in 1:length(FileName)) {
  subset_metadata_seurat[[i]] <- metadata_seurat[metadata_seurat$sample == FileName[[i]], ]
  subset_metadata_seurat[[i]]$Sample_deID <- metadata$Sample_deID[[i]]
  subset_metadata_seurat[[i]]$PatientID <- metadata$PatientID[[i]]
  subset_metadata_seurat[[i]]$Subtype <- metadata$Subtype[[i]]
  subset_metadata_seurat[[i]]$Age <- metadata$Age[[i]]
  subset_metadata_seurat[[i]]$Gender <- metadata$Gender[[i]]
  subset_metadata_seurat[[i]]$Sampling <- metadata$Sampling[[i]]
  subset_metadata_seurat[[i]]$Fusion <- metadata$Fusion[[i]]
  subset_metadata_seurat[[i]]$sequencing <- metadata$`sc/snRNAseq`[[i]]
}

# merge metadata information
merged_metadata_seurat <- subset_metadata_seurat[[1]]

for (i in 2:length(subset_metadata_seurat)) {
  merged_metadata_seurat <- rbind(merged_metadata_seurat, subset_metadata_seurat[[i]])
}

# Add metadata to original seurat object
seurat_obj <- AddMetaData(seurat_obj, merged_metadata_seurat)




# add metadata to object with fresh malignant cells only -----------------------------------
metadata <- read_xlsx(file.path(base_dir, 'data/metadata/metadata_ST_EPN_Patient.xlsx'))
metadata <- subset(metadata, metadata$`sc/snRNAseq` == 'scRNA-seq')
FileName <- unique(metadata$FileName)
index <- match(FileName, metadata$FileName)
metadata <- metadata[index, ]

subset_metadata_seurat <- list()
metadata_seurat <- seurat_obj_fresh@meta.data

for (i in 1:length(FileName)) {
  subset_metadata_seurat[[i]] <- metadata_seurat[metadata_seurat$sample == FileName[[i]], ]
  subset_metadata_seurat[[i]]$Sample_deID <- metadata$Sample_deID[[i]]
  subset_metadata_seurat[[i]]$PatientID <- metadata$PatientID[[i]]
  subset_metadata_seurat[[i]]$Subtype <- metadata$Subtype[[i]]
  subset_metadata_seurat[[i]]$Age <- metadata$Age[[i]]
  subset_metadata_seurat[[i]]$Gender <- metadata$Gender[[i]]
  subset_metadata_seurat[[i]]$Sampling <- metadata$Sampling[[i]]
  subset_metadata_seurat[[i]]$Fusion <- metadata$Fusion[[i]]
  subset_metadata_seurat[[i]]$sequencing <- metadata$`sc/snRNAseq`[[i]]
}

# merge metadata information
merged_metadata_seurat <- subset_metadata_seurat[[1]]

for (i in 2:length(subset_metadata_seurat)) {
  merged_metadata_seurat <- rbind(merged_metadata_seurat, subset_metadata_seurat[[i]])
}

# Add metadata to original seurat object
seurat_obj_fresh <- AddMetaData(seurat_obj_fresh, merged_metadata_seurat)




# replace names metaprograms-----------------------------------
seurat_obj_mal@meta.data$Metaprogram <- str_replace_all(seurat_obj_mal@meta.data$Metaprogram, "NPC-like2", "Neuronal-like")
seurat_obj_mal@meta.data$Metaprogram <- str_replace_all(seurat_obj_mal@meta.data$Metaprogram, "Progenitor-like", "Embryonic-like")
seurat_obj_mal@meta.data$Metaprogram <- str_replace_all(seurat_obj_mal@meta.data$Metaprogram, "NPC-like", "Embryonic-neuronal-like")

# reorder metaprograms -----------------------------------
seurat_obj_mal$Metaprogram <- factor(x = seurat_obj_mal$Metaprogram, levels = c("Cycling", "Embryonic-like", "Neuroepithelial-like", "Radial-glia-like", 
                                                                                "Embryonic-neuronal-like", "Neuronal-like"  ,"Ependymal-like", "MES-like"))

seurat_obj_mal$Subtype <- factor(x = seurat_obj_mal$Subtype, levels = c("ZFTA-RELA" ,"ZFTA-Cluster 1",  "ZFTA-Cluster 2",  "ZFTA-Cluster 3",  "ZFTA-Cluster 4",  "ST-YAP1"))


colnames(seurat_obj_fresh@meta.data)[19] <- 'Metaprogram'


# save objects -----------------------------------
qsave(seurat_obj, file.path(base_dir, 'data/patient/seurat_obj_ST_normal_malig_annotated.qs'))
qsave(seurat_obj_mal, file.path(base_dir, 'data/patient/seurat_obj_ST_malig_annotated.qs'))
qsave(seurat_obj_fresh, file.path(base_dir, 'data/patient/seurat_obj_ST_fresh_malig_annotated.qs'))




# combine fresh and frozen ----------------------------------------------
rm(seurat_obj)
seurat_obj_combined_mal <- merge(seurat_obj_fresh, seurat_obj_mal)
seurat_obj_combined_mal <- JoinLayers(seurat_obj_combined_mal)

rm(seurat_obj_fresh, seurat_obj_mal)

# reprocess normal+malignant object with RUnFullseurat
metadata <- seurat_obj_combined_mal@meta.data

cm <- seurat_obj_combined_mal[["RNA"]]$counts
cm_norm <- as.matrix(log2(cm/10+1))
cm_mean <- log2(Matrix::rowMeans(cm)+1)
cm_center <- cm_norm - rowMeans(cm_norm)
seurat_obj_combined_mal <- RunFullSeurat_v5(cm = cm, metadata = metadata,  doBatch = T, var2batch = 'sample', batchMethod = 'harmony', project = 'HOPE')


# add metadata to object -----------------------------------
metadata <- read_xlsx(file.path(base_dir, 'data/metadata/metadata_ST_EPN_Patient.xlsx'))
metadata <- subset(metadata, metadata$`sc/snRNAseq` %in% c('scRNA-seq', 'snRNA-seq'))
FileName <- unique(metadata$FileName)
index <- match(FileName, metadata$FileName)
metadata <- metadata[index, ]

# first remove mg, od and BT

metadata_seurat <- seurat_obj_combined_mal@meta.data
for (i in 1:length(FileName)) {
  subset_metadata_seurat[[i]] <- metadata_seurat[metadata_seurat$sample == FileName[[i]], ]
  subset_metadata_seurat[[i]]$Sample_deID <- metadata$Sample_deID[[i]]
  subset_metadata_seurat[[i]]$PatientID <- metadata$PatientID[[i]]
  subset_metadata_seurat[[i]]$Subtype <- metadata$Subtype[[i]]
  subset_metadata_seurat[[i]]$Age <- metadata$Age[[i]]
  subset_metadata_seurat[[i]]$Gender <- metadata$Gender[[i]]
  subset_metadata_seurat[[i]]$Sampling <- metadata$Sampling[[i]]
  subset_metadata_seurat[[i]]$Fusion <- metadata$Fusion[[i]]
  subset_metadata_seurat[[i]]$sequencing <- metadata$`sc/snRNAseq`[[i]]
}

# merge metadata information
merged_metadata_seurat <- subset_metadata_seurat[[1]]

for (i in 2:length(subset_metadata_seurat)) {
  merged_metadata_seurat <- rbind(merged_metadata_seurat, subset_metadata_seurat[[i]])
}

# Add metadata to original seurat object
seurat_obj_combined_mal <- AddMetaData(seurat_obj_combined_mal, merged_metadata_seurat)




qsave(seurat_obj_combined_mal, file.path(base_dir, 'data/patient/seurat_obj_ST_combined_malig_annotated.qs'))


