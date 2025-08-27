rm(list = ls())

# Load packages -----------------------------------
library(data.table)
library(tidyverse)
library(crayon)
library(ape)
library(readxl)
library(qs)
library(cowplot)
library(Seurat)
library(ggpubr)
library(future)
library(SingleCellExperiment)
library(SingleR)
library(celldex)

# Set up environment  -------------------------------------
base_dir <- '/Users/sdaniell/Partners HealthCare Dropbox/Sara Danielli/Sara Danielli/Project/Ependymoma/12 - Xenium'
analysis_dir <- file.path(base_dir, 'analysis/1_preparation')
output_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
plot_dir <- file.path(analysis_dir, 'plot')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

data_dir <- "/Users/sdaniell/Partners HealthCare Dropbox/Sara Danielli/Sara Danielli/Project/Ependymoma/11 - Plots scRNAseq/data/patient"

resource_dir <- file.path('/Users/sdaniell/Partners HealthCare Dropbox/Sara Danielli/Sara Danielli/Project/Ependymoma/11 - Plots scRNAseq/scripts/resources')
source(file.path(resource_dir, 'single_cell_preprocessing_helper_functions_CBJr.R'))

# Load seurat object with both malignant and normal cells ---------------------
seurat_object <- qread(file = file.path(data_dir, "seurat_obj_ST_normal_malig_annotated.qs"))

# subset to ZFTA-RELA cells only
table(seurat_object$malignant)
Idents(seurat_object) <- 'Subtype'
seurat_object_ZR <- subset(seurat_object, idents = 'ZFTA-RELA')

# subset to normal cells only
table(seurat_object_ZR$malignant)
Idents(seurat_object_ZR) <- 'malignant'
seurat_object_ZR_normal <- subset(seurat_object_ZR, idents = 'Non-malignant')
sort(unique(seurat_object_ZR_normal$Sample_deID))

# load seurat object with malignant cells ------------------------------
seurat_obj_mal <- qread(file = file.path(data_dir, "seurat_obj_malignant_annotated2.qs"))

# subset to ZFTA-RELA cells only
table(seurat_obj_mal$Subtype)
Idents(seurat_obj_mal) <- 'Subtype'
seurat_object_ZR_mal <- subset(seurat_obj_mal, idents = 'ZFTA-RELA')
sort(unique(seurat_object_ZR_mal$Sample_deID))


# Create combined annotated object --------------------

# Add malignant cell state to object with all cells 
table(seurat_object_ZR_mal$Metaprogram_noCC)
table(seurat_object_ZR_mal$Metaprogram)
seurat_object_ZR_mal$cell_type <- seurat_object_ZR_mal$Metaprogram_noCC

# concatenate
seurat_object <- merge(seurat_object_ZR_normal, seurat_object_ZR_mal)
table(seurat_object$cell_type)
seurat_object <- JoinLayers(seurat_object)

# remove unclassified
Idents(seurat_object) <- 'cell_type'
seurat_object <- subset(seurat_object, idents = 'Unclassified', invert = T)

# reprocess normal+malignant object with RUnFullseurat
metadata <- seurat_object@meta.data

cm <- seurat_object[["RNA"]]$counts
cm_norm <- as.matrix(log2(cm/10+1))
cm_mean <- log2(Matrix::rowMeans(cm)+1)
cm_center <- cm_norm - rowMeans(cm_norm)
seurat_object <- RunFullSeurat_v5(cm = cm, metadata = metadata,  doBatch = F,  project = 'EPN')

# plot
colors_tumor_calling <- c('grey90', 'lightyellow3', '#3EBCB6FF',  '#BDA14DFF', '#09090CFF', '#D5114EFF')
names(colors_tumor_calling) <- c('Unclassified', 'Malignant', 'Myeloid', 'T-cells', 'Oligodendrocytes', 'Endothelial')
colors_metaprograms<- c("gray30","#F99E93FF","#9E5E9BFF","#74ADD1FF",'#0F4F8B', "#ACD39EFF","#96410EFF", 'mistyrose1')

names(colors_metaprograms) <- c("Cycling", "Neuroepithelial-like", "Radial-glia-like", 
                                "Embryonic-neuronal-like", "Neuronal-like" ,"Ependymal-like", "MES-like", "Embryonic-like")

colors_pop <- c(colors_tumor_calling, colors_metaprograms)

p1 <- DimPlot(seurat_object, reduction = 'umap', group.by = 'cell_type', cols = colors_pop, label.size = 4, label = T)
p2 <- DimPlot(seurat_object, reduction = 'umap', group.by = 'Sample_deID', label.size = 4, label = T)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 0.8))
ggsave(file.path(plot_dir, '1_UMAP.pdf'), width = 12, height = 4)

# run SCT -------------------------------------------------
seurat_object <- SCTransform(seurat_object,  verbose = T)

# save ----------------------------
qsave(seurat_object, file.path(output_dir, "ZFTARELA_Xenium_projection.qs"))
  
