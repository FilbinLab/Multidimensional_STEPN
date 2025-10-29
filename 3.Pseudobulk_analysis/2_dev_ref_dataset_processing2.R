## This script imports the two references preprocessed in the script before and combines them into a single dataset
## It also calculates pseudobulked cm, HVGs and DEGs for projection in the next script

# Load packages -----------------------------------
rm(list = ls())

library(data.table)
library(tidyverse)
library(qs)
library(Seurat)
library(data.table)
library(SeuratData)
library(SeuratDisk)
library(readxl)
library(SeuratWrappers)

# Organize environment  -----------------------------------
base_dir <- "/n/scratch/users/s/sad167/developmental_dataset"

data_dir <- file.path(base_dir, 'data')
analysis_dir <- file.path(base_dir, 'analysis/Nowakowski_Eze/data')
if (!dir.exists(analysis_dir)){dir.create(analysis_dir, recursive = T)}
plot_dir <- file.path(base_dir, 'analysis/Nowakowski_Eze/plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}
analysis_dir <- file.path(base_dir, 'analysis/Nowakowski_Eze')

resource_dir <- file.path(base_dir, 'scripts/resources')
source(file.path(resource_dir, 'single_cell_preprocessing_helper_functions_SD_240131.R'))
source(file.path(resource_dir, 'Plot_style_v2.R'))

# Read datasets -----------------------------------
Eze <- qread(file.path(data_dir, "Eze_2021/Eze_2021.qs"))
Nowakowski <- qread(file.path(data_dir, "Nowakowski_2017/Nowakowski_2017.qs"))


# Unify metadata -----------------------------------

# Eze - add times
cluster_interpetation_Eze <- read_xlsx(file.path(data_dir, "Eze_2021/Cluster_interpretation.xlsx"))
# Use match() to get the indices of matching elements
matching_indices <- match(Eze@meta.data$Individual, cluster_interpetation_Eze$`Individual`)
Eze@meta.data$dev_time_point <- cluster_interpetation_Eze$Trimester[matching_indices]
Eze@meta.data$Age_in_Weeks <- cluster_interpetation_Eze$Age_in_Weeks[matching_indices]

# rename cell type labels to be consistent with Nowakowski
Eze@meta.data <- Eze@meta.data %>%
  mutate(cell_type = ifelse(cell_type == "Radial Glial", "Radial glia", cell_type),
         cell_type = ifelse(cell_type == "IPC", "IPCs", cell_type))

# remove cluster with annotation Other
Eze <- subset(Eze, subset = cell_type != c("Other"))

# add origin
Eze[['origin']] <- 'Eze_2021'

# Nowakowski - add times
cluster_interpetation_Nowakowski <- read_xlsx(file.path(data_dir, "Nowakowski_2017/Cluster_interpretation_time.xlsx"))
# Use match() to get the indices of matching elements
matching_indices <- match(Nowakowski@meta.data$Age_in_Weeks, cluster_interpetation_Nowakowski$`Age_in_Weeks`)
Nowakowski@meta.data$dev_time_point <- cluster_interpetation_Nowakowski$Trimester[matching_indices]

# add grouped labels
# group markers
cluster_interpetation_Nowakowski <- read_xlsx(file.path(data_dir, "Nowakowski_2017/Cluster_interpretation.xlsx"))
# Use match() to get the indices of matching Cluster_names in df2
matching_indices <- match(Nowakowski@meta.data$WGCNAcluster, cluster_interpetation_Nowakowski$`Cluster Name`)
Nowakowski@meta.data$cell_type <- cluster_interpetation_Nowakowski$Annotation[matching_indices]
# remove cluster with annotation NA (either unknown or non-brain origin)
Nowakowski <- subset(Nowakowski, subset = cell_type != c("NA"))

# add origin
Nowakowski[['origin']] <- 'Nowakowski_2017'


# Merge datasets -----------------------------------
data <- merge(Eze, Nowakowski)
data <- JoinLayers(data)

rm(Nowakowski, Eze)

# reprocess with Full Seurat
data@meta.data$origin -> samples

cm <- data[["RNA"]]$counts
cm_norm <- as.matrix(log2(cm/10+1))
cm_mean <- log2(Matrix::rowMeans(cm)+1)
cm_center <- cm_norm - rowMeans(cm_norm)
data2 <- RunFullSeurat(cm = cm, samples = samples, RunHarmony = TRUE, project = 'PanCancer')

# Add metadata
data2 <- AddMetaData(data2, data@meta.data)
data <- data2

# save
qsave(data, file.path(analysis_dir, 'data/Nowakowski_Eze.qs'))


# Explore combined dataset -----------------------------------
data <- qread(file.path(analysis_dir, 'data/Nowakowski_Eze.qs'))

# color scheme
color_cell_type <- c("#F99E93FF","#9E5E9BFF",'#e01e37',
                     '#8be8d7', '#51ccd1', '#0377a8', '#023e7d',
                     "#96410EFF", '#e41b60', '#d0a03b', "#ACD39EFF")

names(color_cell_type) <- c( "Neuroepithelial", "Radial glia", "Early radial glia", 
                             "Newborn/maturing neuron",  "MGE inhibitory neuron lineage", "Neuronal", "IPCs" ,
                             "Mesenchymal",  "OPC",  "Astrocyte" , "Choroid")

color_dev_time_points <- c('#eeef20','#aacc00', '#55a630', '#51ccd1','#023e7d')
names(color_dev_time_points) <- c("First - Early",  "First - Middle", "First - Late", "Second" , "Third")
# reorder cell_type
data$cell_type <- factor(x = data$cell_type, 
                         levels = c("Neuroepithelial", "Early radial glia", "Radial glia", 
                                    "Newborn/maturing neuron",  "MGE inhibitory neuron lineage", "Neuronal", "IPCs" ,
                                    "Mesenchymal",  "OPC",  "Astrocyte" , "Choroid")) 

data$dev_time_point <- factor(x = data$dev_time_point, 
                         levels = c("First - Early",  "First - Middle", "First - Late", "Second" , "Third" )) 


# tSNE plots -----------------------------------
DimPlot(object = data, 
        reduction = 'umap', 
        group.by = "cell_type", 
        #split.by = 'dev_time_point',
        pt.size = 5, 
        shuffle = TRUE, 
        label = TRUE,
        label.size = 0,
        cols = color_cell_type,
        raster=TRUE, raster.dpi = c(1500, 1500)
) + labs(title = 'Nowakowsi et al. (2017) + Eze et al. (2021)', subtitle = 'n = 42,449 cells') + theme_tSNE
ggsave(file.path(plot_dir, '1_umap_cell_type.pdf'), width=7.5, height=5)

DimPlot(object = data, 
        reduction = 'umap', 
        group.by = "dev_time_point", 
        #split.by = 'dev_time_point',
        pt.size = 5, 
        shuffle = TRUE, 
        label = TRUE,
        label.size = 0,
        cols = color_dev_time_points,
        raster=TRUE, raster.dpi = c(1500, 1500)
) + labs(title = 'Nowakowsi et al. (2017) + Eze et al. (2021)', subtitle = 'n = 42,449 cells') + theme_tSNE
ggsave(file.path(plot_dir, '2_tumap_dev_time_point.pdf'), width=6, height=5)



# Find markers by cell type -----------------------------------
Idents(data) <- 'cell_type'
markers <- RunPrestoAll(data, only.pos = T, densify = T, max.cells.per.ident = 500, min.pct = 0.25) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T)
qsave(markers, file.path(analysis_dir, "data/markers_single_cell.qs"))
gene_list <- split(markers$gene, markers$cluster)
qsave(gene_list, file.path(analysis_dir, "data/markers_list_single_cell.qs"))



# Process for pseudobulk -----------------------------------
data@meta.data$cell_type <- Idents(data)

# Pseudobulk by cell_type
pseudobulk_mean <- AverageExpression(data, return.seurat = F, slot = "counts", assays = "RNA", group.by = "cell_type")$RNA %>%
  as.data.frame()
qsave(pseudobulk_mean, file.path(analysis_dir, "data/agg_cm_mean.qs"))

# Identify markers
Idents(data) <- data$cell_type
markers <- RunPrestoAll(data, only.pos = T, densify = T, max.cells.per.ident = 500) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T)
qsave(markers, file.path(analysis_dir, "data/markers.qs"))
qsave(VariableFeatures(data), file.path(analysis_dir,"data/hvgs.qs"))


# Pseudobulk by dev_time_point
Idents(data) <- data$dev_time_point

pseudobulk_mean <- AverageExpression(data, return.seurat = F, slot = "counts", assays = "RNA", group.by = "dev_time_point")$RNA %>%
  as.data.frame()
qsave(pseudobulk_mean, file.path(analysis_dir, "data/agg_cm_mean_dev_time.qs"))

# Identify markers
Idents(data) <- data$dev_time_point
markers <- RunPrestoAll(data, only.pos = T, densify = T, max.cells.per.ident = 500) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T)
qsave(markers, file.path(analysis_dir, "data/markers_dev_time.qs"))
qsave(VariableFeatures(data), file.path(analysis_dir,"data/hvgs.qs"))

