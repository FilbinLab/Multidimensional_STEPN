args <- commandArgs(trailingOnly = TRUE)
print(args)

SampleName_arg <- args[1]

# Load packages -----------------------------------
library(Seurat)
library(tidyverse)
library(glue)
library(gplots)
library(qs)


# Organize environment  -----------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

analysis_dir  <- file.path(base_dir, 'analysis/6_coherence')

seurat_obj_dir <- file.path(base_dir, 'analysis/3_program_annotation/data')

path_to_output <- file.path(analysis_dir, "data")
if (!dir.exists(path_to_output)){dir.create(path_to_output, recursive = T)}

plot_dir  <- file.path(base_dir, 'analysis/6_coherence/plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

resource_dir  <- file.path(base_dir, 'scripts_revisions/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))
source(file.path(resource_dir, "xenium_preprocessing_helper_functions_SD_v2.R"))
source(file.path(resource_dir, 'color_palette.R'))

metaprogram_colors <- colors_metaprograms_Xenium



## Run coherence ------------------------------------------------------

# Read annotated datasets
data <- qread(file.path(seurat_obj_dir, paste0('seurat_obj_', SampleName_arg, '.qs')))
print(paste0('Successful reading of ', SampleName_arg))

adata <- LayerData(data, 'counts', 'SCT')
metadata <- data.frame(cell_id = data@images$fov$centroids@cells,
                       x_centroid = data@images$fov$centroids@coords[,1],
                       y_centroid = data@images$fov$centroids@coords[,2])


## Pre-processing data
preprocessing <- prepareData(adata, metadata)

## Normalizing data
normalized_matrix <- norm_data(t(as.matrix(preprocessing$cell_feature_matrix)))
normalized_matrix <- sqrt(normalized_matrix) + sqrt(normalized_matrix + 1)
normalized_matrix <- norm_data(normalized_matrix)

metaprogram <- plyr::mapvalues(rownames(normalized_matrix), from = names(data$cell_type), to = as.character(unname(data$cell_type)), warn_missing = F)

## Gridding
gridding <- grid_spatial(norm_data = t(normalized_matrix), spatial = preprocessing$spatial, 
                         variable = metaprogram, nbins = round(sqrt(dim(data@assays$SCT)[2])))
print(paste0('nbins = ', round(sqrt(dim(data@assays$SCT)[2]))))

coherence <- coherence_score(grid_df = gridding, variable = 'Metaprogram')


# Save output  ----------------------
qsave(coherence, file.path(path_to_output, paste0('results_', SampleName_arg, '.qs')))



# Plot grids with density results ----------------------
gridding_coherence_density(gridding, plot_dir, color_coherence_density, colors_metaprograms_Xenium)

