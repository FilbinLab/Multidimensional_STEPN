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

plot_dir <- file.path(base_dir, 'analysis/plots')

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

# create empty lists where to store objects
seurat_obj_list <- list()
plot_list <- list()


for (i in seq_along(seurat_obj_dir)) {
  # read file
  seurat_obj_list <- qread(file.path(base_dir, seurat_obj_dir[i]))
  print(seurat_obj_dir[i])
  
  # Visualize  distribution clusters 
  plot_list[[i]] <- ImageDimPlot(seurat_obj_list, group.by = 'group', cols = colors_groups, border.size = NA, size = 0.5, 
                       dark.background = F) + NoLegend()
  rasterize(plot_list[[i]], dpi=400)
}


combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol=5) 
rasterize(combined_plot, dpi=400)
ggsave(file.path(plot_dir, 'Spatial_summary.pdf'), width=16, height=5)

