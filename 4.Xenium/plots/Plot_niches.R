rm(list = ls())

# Load packages -----------------------------------
library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(data.table)
library(readxl)
library(dplyr)
library(ggrastr)
library(cowplot)
library(glue)
library(RColorBrewer)
library(xlsx)
library(ggrepel)
library(ComplexHeatmap)

options(future.globals.maxSize = 10000000000000000000)

# Organize environment  ----------------------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

resource_dir  <- file.path(base_dir, 'scripts_revisions/resources')
#source(file.path(resource_dir, 'Plotting_functions.R'))
#source(file.path(resource_dir, 'NMF_helper_function.R'))
#source(file.path(resource_dir, 'Plotting_helper_functions.R'))
source(file.path(resource_dir, 'xenium_preprocessing_helper_functions_SD_v2.R'))
source(file.path(resource_dir, 'color_palette.R'))

plot_dir  <- file.path(base_dir, 'analysis/4_niche/plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

data_dir <- file.path(base_dir, 'analysis/4_niche/data')


metadata <- read_xlsx(file.path(base_dir, 'scripts_revisions/SampleIdentifier.xlsx'))

# reorder by sampleID
metadata <- metadata %>% arrange(SampleID)



# Plot all tumor niches together --------------------------------------------
SampleName_vector <- metadata$SampleName
colors_metaprogram <- colors_metaprograms_Xenium
order_metaprograms <- c("Neuroepithelial-like", "Radial-glia-like", 
                        "Neuronal-like" ,"Ependymal-like", "MES-like", 
                        "T-cells", "Myeloid",  "Endothelial",  "Oligodendrocytes",  'Neurons', 'Unknown')


PlotNicheResults(data_dir, SampleName_vector, colors_metaprogram, colors_niches, niches_k = 4, plot_spatial_maps = F, max_col = 50, mid_col = 25) 
PlotNicheResults(data_dir, SampleName_vector, colors_metaprogram, colors_niches, niches_k = 5, plot_spatial_maps = F) 
PlotNicheResults(data_dir, SampleName_vector, colors_metaprogram, colors_niches, niches_k = 6, plot_spatial_maps = F) 



# Plot recurrent niches -------------------------------------------------------

# all samples niche = 4
HeatmapNicheCorrelation(data_dir, SampleName_vector, colors_metaprogram = colors_metaprograms_Xenium, 
                        colors_niches = paletteer::paletteer_d("beyonce::X80"), niches_k = 4, nClust = 6,
                        outputname = 'all')

# all samples niche = 5
HeatmapNicheCorrelation(data_dir, SampleName_vector, colors_metaprogram = colors_metaprograms_Xenium, 
                        colors_niches = paletteer::paletteer_d("beyonce::X80"), niches_k = 5, nClust = 5,
                        outputname = 'all')

# all samples niche = 6
HeatmapNicheCorrelation(data_dir, SampleName_vector, colors_metaprogram = colors_metaprograms_Xenium, 
                        colors_niches = paletteer::paletteer_d("beyonce::X80"), niches_k = 6, nClust = 6,
                        outputname = 'all')



# structured tumors
ordered <- c('STEPN10_Region_2',
             'STEPN15_Region_1',
             'STEPN01_Region_3',
             'STEPN18_Region_1',
             'STEPN12_Region_4',
             'STEPN12_Region_1',
             'STEPN19_Region_2',
             'STEPN19_Region_1',
             'STEPN12_Region_3',
             'STEPN12_Region_2')
  
HeatmapNicheCorrelation(data_dir, SampleName_vector = ordered, colors_metaprogram = colors_metaprograms_Xenium, 
                        colors_niches = paletteer::paletteer_d("beyonce::X80"), niches_k = 4, nClust = 5,
                        outputname = 'structured')

# disorganized tumors
disordered <- c('STEPN16_Region_1',
                'STEPN06_Region_5',
                'STEPN06_Region_4',
                'STEPN10_Region_3',
                'STEPN06_Region_1',
                'STEPN06_Region_3',
                'STEPN06_Region_2',
                'STEPN01_Region_2',
                'STEPN10_Region_1',
                'STEPN01_Region_1',
                'STEPN17_Region_1')

HeatmapNicheCorrelation(data_dir, SampleName_vector = disordered, colors_metaprogram = colors_metaprograms_Xenium, 
                        colors_niches = paletteer::paletteer_d("beyonce::X80"), niches_k = 4, nClust = 4,
                        outputname = 'disorganized')





# To plot zoomed images -----------------------------------------------------
PlotSpatialMapsZoom(SampleName = "STEPN12_Region_3", input_dir = data_dir, group_display = 'cell_type', 
                    colors = colors_metaprograms_Xenium, size_dot = 0.7, geom_rect = F,  
                    x_min = 1700, x_max = 2300, y_min = 2600, y_max = 3200,
                    ratio = c(2, 0.6))
ggsave(file.path(plot_dir, paste0('10_ImageDimPlot_malignant_STEPN12_Region_3.pdf')), width=12, height=10)

PlotSpatialMapsZoom(SampleName = "STEPN12_Region_3", input_dir = data_dir, group_display = 'niche_4', 
                    colors = colors_niches, size_dot = 0.7, geom_rect = F,  
                    x_min = 1700, x_max = 2300, y_min = 2600, y_max = 3200,
                    ratio = c(1, 1))
ggsave(file.path(plot_dir, paste0('10_ImageDimPlot_niche_STEPN12_Region_3.pdf')), width=12, height=10)

PlotSpatialMapsZoom(SampleName = "STEPN12_Region_3", input_dir = data_dir, group_display = 'cell_type', 
                    colors = colors_metaprograms_Xenium, size_dot = 0.7, geom_rect = F,  
                    x_min = 3700, x_max = 4800, y_min = 4400, y_max = 5400,
                    ratio = c(2, 0.6))
ggsave(file.path(plot_dir, paste0('11_ImageDimPlot_malignant_STEPN12_Region_3.pdf')), width=12, height=10)

PlotSpatialMapsZoom(SampleName = "STEPN12_Region_3", input_dir = data_dir, group_display = 'niche_4', 
                    colors = colors_niches, size_dot = 0.7, geom_rect = F,  
                    x_min = 3700, x_max = 4800, y_min = 4400, y_max = 5400,
                    ratio = c(1, 1))
ggsave(file.path(plot_dir, paste0('11_ImageDimPlot_niche_STEPN12_Region_3.pdf')), width=12, height=10)






PlotSpatialMapsZoom(SampleName = "STEPN16_Region_1", input_dir = data_dir, group_display = 'cell_type', 
                    colors = colors_metaprograms_Xenium, size_dot = 0.7, geom_rect = F,  
                    x_min = 500, x_max = 1500, y_min = 1000, y_max = 2000,
                    ratio = c(1, 1))
ggsave(file.path(plot_dir, paste0('12_ImageDimPlot_niche_STEPN16_Region_1.pdf')), width=12, height=10)



PlotSpatialMapsZoom(SampleName = "STEPN17_Region_1", input_dir = data_dir, group_display = 'cell_type', 
                    colors = colors_metaprograms_Xenium, size_dot = 0.7, geom_rect = F,  
                    x_min = 500, x_max = 1500, y_min = 1000, y_max = 2000,
                    ratio = c(1, 1))
ggsave(file.path(plot_dir, paste0('12_ImageDimPlot_niche_STEPN17_Region_1.pdf')), width=12, height=10)


