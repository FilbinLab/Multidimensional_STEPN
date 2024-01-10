# Load packages -----------------------------------
#library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
#library(ggpubr)
library(data.table)
#library(readxl)
library(dplyr)
#library(plyr)
library(cowplot)

# Organize environment  -----------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

resource_dir  <- file.path('/n/scratch/users/s/sad167/EPN/Ependymoma2023/Xenium/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))


seurat_obj_dir <- c('data/processed_data_Carlos/20231020__200939__BT2126_BT1745/0010652-Region_4/0010652-Region_4.qs',
                    'data/processed_data_Carlos/20231102__215055__7EP1_7EP41_3EP8/0010575-Region_3/0010575-Region_3.qs',
                    'data/processed_data_Carlos/20231102__215055__7EP1_7EP41_3EP8/0010619-Region_1/0010619-Region_1.qs',
                    'data/processed_data_Carlos/20231102__215055__7EP1_7EP41_3EP8/0010619-Region_5/0010619-Region_5.qs',
                    'data/processed_data_Carlos/20231107__203958__BT1717_BT775/0010501-Region_2/0010501-Region_2.qs',
                    'data/processed_data_Carlos/20231109__203408__BT1804_BT2169/0010498-Region_1/0010498-Region_1.qs', 
                    'data/processed_data_Carlos/20231109__203408__BT1804_BT2169/0010775-Region_1/0010775-Region_1.qs', 
                    'data/processed_data_Carlos/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/0010540-Region_1/0010540-Region_1.qs',
                    'data/processed_data_Carlos/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/0010553-Region_3/0010553-Region_3.qs')

data_files <- c('analysis/20231020__200939__BT2126_BT1745/data/individual/cell_ID_0010652-Region_4.csv',
                'analysis/20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010575-Region_3.csv',
                'analysis/20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010619-Region_1.csv',
                'analysis/20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010619-Region_5.csv',
                'analysis/20231107__203958__BT1717_BT775/data/individual/cell_ID_0010501-Region_2.csv',
                'analysis/20231109__203408__BT1804_BT2169/data/individual/cell_ID_0010498-Region_1.csv', 
                'analysis/20231109__203408__BT1804_BT2169/data/individual/cell_ID_0010775-Region_1.csv', 
                'analysis/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/individual/cell_ID_0010540-Region_1.csv',
                'analysis/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/individual/cell_ID_0010553-Region_3.csv')

metadata <- list()
proportions <- list()
for (i in seq_along(data_files)) {
  metadata[i] <- read.table(file.path(base_dir, data_files[i]), header = TRUE, sep = ',', row.names = "cell_id")
  proportions[[i]] <- prop.table(table(metadata[i]))
}

seurat_obj_list <- list()

for (i in seq_along(seurat_obj_dir)) {
  seurat_obj_list[[i]] <- qread(file.path(base_dir, seurat_obj_dir[i]))
  print(seurat_obj_dir[i])
}

# names(seurat_obj_list) <- c('BT2126_Region_4', '7EP41_Region_3', '3EP8_Region_1', '7EP1_Region_5',
#                             'BT775_Region_2', 'BT2169_Region_1', 'BT1804_Region_1',
#                             '11EP22_Region_1',  'BT1743_Region_3')

#' names(metadata) <- c('BT2126_Region_4', '7EP41_Region_1', '7EP41_Region_2', '7EP41_Region_3', 
#'                         '3EP8_Region_1.1', '3EP8_Region_2.1', '3EP8_Region_3', '7EP1_Region_4.1', '7EP1_Region_5', 
#'                         'BT775_Region_1', 'BT775_Region_2', 
#'                         #'BT1717_Region_1', 'BT1717_Region_2',
#'                         'BT2169_Region_1', 'BT1804_Region_1', 
#'                         '11EP22_Region_1', '11EP22_Region_2', '11EP22_Region_3', '7EP41_Region_4', 
#'                         '3EP8_Region_1.2', '3EP8_Region_2.2', 'BT1743_Region_3', '7EP1_Region_4.2')

# Concatenate seurat objects
seurat_obj_combined <- merge(seurat_obj_list[[1]], 
                             y = c(seurat_obj_list[c(2:length(seurat_obj_list))]), 
                             add.cell.ids = c('BT2126_Region_4', '7EP41_Region_3', '3EP8_Region_1', '7EP1_Region_5',
                                              'BT775_Region_2', 'BT2169_Region_1', 'BT1804_Region_1',
                                              '11EP22_Region_1',  'BT1743_Region_3'))


# remove single objects to empty memory
rm(seurat_obj_list)

# save combined object
qsave(seurat_obj_combined, file.path(base_dir, 'analysis/Seurat_obj_combined.qs'))

