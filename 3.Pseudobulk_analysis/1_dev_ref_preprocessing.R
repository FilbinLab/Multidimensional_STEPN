## This script imports two developmental references from the original publications, and prepares them for developmental mapping

# Load packages -----------------------------------
rm(list = ls())

library(data.table)
library(tidyverse)
library(qs)
library(Seurat)
library(data.table)
library(SeuratData)
library(SeuratDisk)

# Organize environment  -----------------------------------
base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Developmental datasets"

data_dir <- file.path(base_dir, 'data')

resource_dir <- file.path(base_dir, 'scripts/resources')
source(file.path(resource_dir, 'single_cell_preprocessing_helper_functions_SD_240131.R'))
source(file.path(resource_dir, 'NMF_helper_function.R'))



# Load and process Nowakowski et al. 2017 (cortex) -----------------------------------

analysis_dir <- file.path(base_dir, 'data/Nowakowski_2017')

# Load reference
mat <- fread(file.path(analysis_dir, "raw/exprMatrix.tsv.gz"))
meta <- read.table(file.path(analysis_dir, "raw/meta.tsv"), header=T, sep="\t", as.is=T, row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
ref <- CreateSeuratObject(counts = mat, project = "Nowakowski", meta.data=meta)

# normalize
ref <- NormalizeData(ref, normalization.method = "LogNormalize", scale.factor = 10000)
# HVGs
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000)

# scale
ref <- ScaleData(ref, features = rownames(ref))

# remove cluster with no annotation
ref <- subset(ref, subset = WGCNAcluster != "")

# save dataset
qsave(ref, file.path(analysis_dir, "Nowakowski_2017.qs"))


# Pseudobulk by cell_type
pseudobulk_mean <- AverageExpression(ref, return.seurat = F, slot = "counts", assays = "RNA", group.by = "WGCNAcluster")$RNA %>%
  as.data.frame()
colnames(pseudobulk_mean) <- sort(unique(ref$WGCNAcluster))
qsave(pseudobulk_mean, file.path(analysis_dir, "agg_cm_mean_WGCNAcluster.qs"))

# Identify markers
Idents(ref) <- ref$WGCNAcluster
markers <- RunPrestoAll(ref, only.pos = T, min.pct = 0.15) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T)
qsave(markers, file.path(analysis_dir, "markers_WGCNAcluster.qs"))
qsave(VariableFeatures(ref), file.path(analysis_dir,"hvgs_WGCNAcluster.qs"))


# group markers
cluster_interpetation <- read_xlsx(file.path(analysis_dir, "Cluster_interpretation.xlsx"))

# Use match() to get the indices of matching Cluster_names in df2
matching_indices <- match(ref@meta.data$WGCNAcluster, cluster_interpetation$`Cluster Name`)

ref@meta.data$cell_type_group <- cluster_interpetation$Annotation[matching_indices]

# remove cluster with annotation NA (either unknown or non-brain origin)
ref <- subset(ref, subset = cell_type_group != c("NA"))

# Pseudobulk by cell_type aggregate
pseudobulk_mean <- AverageExpression(ref, return.seurat = F, slot = "counts", assays = "RNA", group.by = "cell_type_group")$RNA %>%
  as.data.frame()
colnames(pseudobulk_mean) <- sort(unique(ref$cell_type_group))
qsave(pseudobulk_mean, file.path(analysis_dir, "agg_cm_mean_cell_type_group.qs"))

# Identify markers
Idents(ref) <- ref$cell_type_group
markers <- RunPrestoAll(ref, only.pos = T, min.pct = 0.15) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T)
qsave(markers, file.path(analysis_dir, "markers_cell_type_group.qs"))
qsave(VariableFeatures(ref), file.path(analysis_dir,"hvgs_cell_type_group.qs"))


# Pseudobulk by time point
pseudobulk_mean <- AverageExpression(ref, return.seurat = F, slot = "counts", assays = "RNA", group.by = "Age_in_Weeks")$RNA %>%
  as.data.frame()
colnames(pseudobulk_mean) <- sort(unique(ref$Age_in_Weeks))
qsave(pseudobulk_mean, file.path(analysis_dir, "agg_cm_mean_Age_in_Weeks.qs"))

# Identify markers
Idents(ref) <- 'Age_in_Weeks'
markers <- RunPrestoAll(ref, only.pos = T, min.pct = 0.15) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T)
qsave(markers, file.path(analysis_dir, "markers_Age_in_Weeks.qs"))
qsave(VariableFeatures(ref), file.path(analysis_dir,"hvgs_Age_in_Weeks.qs"))



# 2) Load and process Eze et al. 2021  -----------------------------------

analysis_dir <- file.path(base_dir, 'data/Eze_2021')

# Load reference
mat <- fread(file.path(analysis_dir, "raw/counts_exprMatrix.tsv.gz"))
meta <- read.table(file.path(analysis_dir, "raw/meta-2.tsv"), header=T, sep="\t", as.is=T, row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
ref <- CreateSeuratObject(counts = mat, project = "Nowakowski", meta.data=meta)
rm(mat)

# normalize
ref <- NormalizeData(ref, normalization.method = "LogNormalize", scale.factor = 10000)


# remove genes that ahve negative value (error with HVGs otherwise)
RowsNA<-names(which(rowSums(is.na(GetAssayData(ref)))>0))
'%!in%' <- function(x,y)!('%in%'(x,y)) #this is a NOT IN function
RowsKEEP<-rownames(ref)[rownames(ref) %!in% RowsNA]
ref <- subset(ref,features=RowsKEEP)


# HVGs
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000, layer = 'data')

# scale
ref <- ScaleData(ref)

# add cluster annotation
cluster_interpetation <- read_xlsx(file.path(analysis_dir, "raw/annotation.xlsx"), sheet = "5 Cortex Annotations")

# Use match() to get the indices of matching Cluster_names in df2
matching_indices <- match(ref@meta.data$Cluster, cluster_interpetation$`Cluster`)

ref@meta.data$cell_type <- cluster_interpetation$`Cell Type`[matching_indices]

# save dataset
qsave(ref, file.path(analysis_dir, "Eze_2021.qs"))

ref <- qread(file.path(analysis_dir, "Eze_2021.qs"))

# Pseudobulk by cell_type
pseudobulk_mean <- AverageExpression(ref, return.seurat = F, slot = "counts", assays = "RNA", group.by = "cell_type")$RNA %>%
  as.data.frame()
colnames(pseudobulk_mean) <- sort(unique(ref$cell_type))
qsave(pseudobulk_mean, file.path(analysis_dir, "agg_cm_mean.qs"))

# Identify markers
Idents(ref) <- ref$cell_type
markers <- RunPrestoAll(ref, only.pos = T, densify = T, max.cells.per.ident = 500) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T)
qsave(markers, file.path(analysis_dir, "markers.qs"))
qsave(VariableFeatures(ref), file.path(analysis_dir,"hvgs.qs"))



desired_order <-  rev(c("Neuroepithelial", "Radial Glial", "Neuronal", "Mesenchymal", "IPC" ))
ref$cell_type <- factor(x = ref$cell_type, levels = desired_order)

DotPlot(ref, group.by = 'cell_type', 
        features = c("PTGDS", 'FBLN5', 'CRABP2', 'GPC3', "HES3"),
        dot.scale = 6,
        #col.min = 0, 
        #cols = c("white", "red3"),
        scale = TRUE, assay = 'RNA') + 
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) + 
  RotatedAxis()
ggsave(file.path(analysis_dir, "1_Dotplot_markers_neuroepithelial.pdf"), width=6.5, height=4)





