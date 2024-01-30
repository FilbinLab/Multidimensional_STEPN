# Load packages -----------------------------------
#library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(data.table)
library(dplyr)
library(cowplot)
library(Seurat)
library(ggrastr)

# Organize environment  -----------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

plot_dir <- file.path(base_dir, 'analysis/plots/spatial_coherence')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

resource_dir  <- file.path('/n/scratch/users/s/sad167/EPN/Ependymoma2023/Xenium/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))


seurat_obj_dir <- c('analysis/20231020__200939__BT2126_BT1745/data/seurat_objects/Seurat_obj0010652-Region_4.qs',
                    
                    #'analysis/20231102__215055__7EP1_7EP41_3EP8/data/seurat_objects/Seurat_obj0010575-Region_1.qs',
                    #'analysis/20231102__215055__7EP1_7EP41_3EP8/data/seurat_objects/Seurat_obj0010575-Region_2.qs',
                    'analysis/20231102__215055__7EP1_7EP41_3EP8/data/seurat_objects/Seurat_obj0010575-Region_3.qs',
                    #'analysis/20231102__215055__7EP1_7EP41_3EP8/data/seurat_objects/Seurat_obj0010619-Region_1.qs',
                    #'analysis/20231102__215055__7EP1_7EP41_3EP8/data/seurat_objects/Seurat_obj0010619-Region_2.qs',
                    'analysis/20231102__215055__7EP1_7EP41_3EP8/data/seurat_objects/Seurat_obj0010619-Region_3.qs',
                    #'analysis/20231102__215055__7EP1_7EP41_3EP8/data/seurat_objects/Seurat_obj0010619-Region_4.qs',
                    'analysis/20231102__215055__7EP1_7EP41_3EP8/data/seurat_objects/Seurat_obj0010619-Region_5.qs',
                    
                    # 'analysis/20231107__203958__BT1717_BT775/data/seurat_objects/Seurat_obj0010501-Region_1.qs',
                    'analysis/20231107__203958__BT1717_BT775/data/seurat_objects/Seurat_obj0010501-Region_2.qs',
                    
                    'analysis/20231107__203958__BT1717_BT775/data/seurat_objects/Seurat_obj0010814-Region_2.qs', 
                    
                    'analysis/20231109__203408__BT1804_BT2169/data/seurat_objects/Seurat_obj0010498-Region_1.qs', 
                    'analysis/20231109__203408__BT1804_BT2169/data/seurat_objects/Seurat_obj0010775-Region_1.qs', 
                    
                    'analysis/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/seurat_objects/Seurat_obj0010540-Region_1.qs',
                    #'analysis/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/seurat_objects/Seurat_obj0010540-Region_2.qs',
                    #'analysis/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/seurat_objects/Seurat_obj0010540-Region_3.qs',
                    #'analysis/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/seurat_objects/Seurat_obj0010540-Region_4.qs',
                    #'analysis/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/seurat_objects/Seurat_obj0010553-Region_1.qs',
                    #'analysis/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/seurat_objects/Seurat_obj0010553-Region_2.qs',
                    'analysis/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/seurat_objects/Seurat_obj0010553-Region_3.qs'
                    #'analysis/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/seurat_objects/Seurat_obj0010553-Region_4.qs'
)

samples <- c('BT2126', '7EP41', '3EP8', '7EP1',
             'BT775', 'BT1717', 'BT2169', 'BT1804',
             '11EP22',  'BT1743')


for (i in seq_along(seurat_obj_dir)) {
  # read file
  seurat_obj <- qread(file.path(base_dir, seurat_obj_dir[i]))
  DefaultAssay(seurat_obj) <- 'Xenium'
  
  # slim down to reduce memory consumption
  seurat_obj <- DietSeurat(seurat_obj,  assays = 'Xenium')
  seurat_obj
  
  # print image
  plot <- ImageDimPlot(seurat_obj, fov = 'fov', group.by = 'group', cols = colors_groups, border.size = NA, size = 0.35, 
                                 dark.background = F, axes = TRUE) + NoLegend()
  plot

  # extract metaprogram
  metadata <- as.data.frame(seurat_obj@meta.data)
  metaprogram <- metadata$group
  
  # extract spatial coordinates
  seurat_obj_coordinates <- GetTissueCoordinates(seurat_obj)
  
  # add metaprogram info to spatial coordinates
  seurat_obj_coordinates$group <- metaprogram
  
  # create normalized coordinates (that intercept at 0 = min)
  min_x <- min(seurat_obj_coordinates$x)
  min_y <- min(seurat_obj_coordinates$y)
  
  # round data to have integers
  seurat_obj_coordinates$x_norm <- round(seurat_obj_coordinates$x - min_x)
  seurat_obj_coordinates$y_norm <- round(seurat_obj_coordinates$y - min_y)
  
  seurat_obj_coordinates_norm <- seurat_obj_coordinates[, c('x_norm', 'y_norm', 'cell')]
  
  # Create new spatial coordinates in seurat object
  cents <- CreateCentroids(seurat_obj_coordinates_norm)
  
  segmentations.data <- list(
    "centroids" = cents,
    "segmentation" = NULL
  )
  
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = NULL,
    assay = "Xenium"
  )
  
  seurat_obj[["norm_coordinates"]] <- coords
  
  # print image
  plot_normalized <- ImageDimPlot(seurat_obj, fov = "norm_coordinates", group.by = 'group', cols = colors_groups, border.size = NA, size = 0.35, 
                       dark.background = F, axes = TRUE) + NoLegend()
  plot + plot_normalized
  
  
  
  
  # creare matrix fatta di 0/1 per ogni color
  # Get the spatial data for the current FOV
  spatial_data <- seurat_obj$Xenium[["norm_coordinates"]]
  
  # Get the image matrix
  image_matrix <- spatial_data$image
  
  # aggiungere una riga e colonna
  # contare i punti che non hanno cambiato valore ma che non sono 0

  
}