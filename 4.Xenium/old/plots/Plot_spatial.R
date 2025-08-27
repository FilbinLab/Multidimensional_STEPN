# Load packages -----------------------------------
library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(data.table)
library(readxl)
library(dplyr)
library(ggrastr)
library(ComplexHeatmap)
library(cowplot)


# Organize environment  -----------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

resource_dir  <- file.path(base_dir, 'scripts_revisions/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))
source(file.path(resource_dir, 'xenium_preprocessing_helper_functions_SD_v2.R'))
source(file.path(resource_dir, 'color_palette.R'))


plot_dir  <- file.path(base_dir, 'analysis/5_plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

input_dir <- file.path(base_dir, 'analysis/4_niche/data')
cellid_dir <- file.path(base_dir, 'analysis/3_program_annotation/data')

seurat_object_path <- list.files(file.path(base_dir, 'analysis/1_preparation/data'), pattern = "\\.qs$", full.names = TRUE)

metadata <- read_xlsx(file.path(base_dir, 'scripts_revisions/SampleIdentifier.xlsx'))

# remove these two samples because of low QC data
metadata <- metadata %>% filter(!SampleName %in% c('STEPN14_Region_1', 'STEPN14_Region_2'))
# reorder by sampleID
metadata <- metadata %>% arrange(SampleName)

# Plot distribution clusters  -------------------------------------
plots_cell_state <- list()
plots_malignant <- list()
plots_niche_4 <- list()
plots_niche_5 <- list()
plots_niche_6 <- list()



# Spatial maps ------------------------------------------------------

for (i in seq_along(metadata$SampleName)) {
  
  # read data
  data <- qread(file.path(input_dir, paste0('seurat_obj_', metadata$SampleName[i], '.qs')))

  
  # plot by cell state
  plots_cell_state[[i]] <- ImageDimPlot(data, group.by = 'cell_type', cols = colors_metaprograms_Xenium,
                                        border.size = NA, size = 0.5, dark.background = F)  +
    NoLegend() + labs(title = metadata$SampleName[i])

  plots_cell_state[[i]]
  ggsave(file.path(plot_dir, paste0('2_ImageDimPlot_cell_type_', metadata$SampleName[i], '.pdf')), width=6, height=5)



  # plot by malignant vs non-malignant
  plots_malignant[[i]] <- ImageDimPlot(data, group.by = 'malignant', cols = col_normal_malignant,
                                       border.size = NA, size = 0.5, dark.background = F) + NoLegend() +
    labs(title = metadata$SampleName[i])

  plots_malignant[[i]]
  ggsave(file.path(plot_dir, paste0('3_ImageDimPlot_malignant_', metadata$SampleName[i], '.pdf')), width=6, height=5)




  # plot by niche
  plots_niche_4[[i]] <- ImageDimPlot(data, group.by = 'niche_4', cols = colors_niches,
                                       border.size = NA, size = 0.5, dark.background = F) + NoLegend() +
    labs(title = metadata$SampleName[i])

  plots_niche_4[[i]]
  ggsave(file.path(plot_dir, paste0('4a_ImageDimPlot_niches_4_', metadata$SampleName[i], '.pdf')), width=6, height=5)


  plots_niche_5[[i]] <- ImageDimPlot(data, group.by = 'niche_5', cols = colors_niches,
                                     border.size = NA, size = 0.5, dark.background = F) + NoLegend() +
    labs(title = metadata$SampleName[i])

  plots_niche_5[[i]]
  ggsave(file.path(plot_dir, paste0('4b_ImageDimPlot_niches_5_', metadata$SampleName[i], '.pdf')), width=6, height=5)


  plots_niche_6[[i]] <- ImageDimPlot(data, group.by = 'niche_6', cols = colors_niches,
                                     border.size = NA, size = 0.5, dark.background = F) + NoLegend() +
    labs(title = metadata$SampleName[i])

  plots_niche_6[[i]]
  ggsave(file.path(plot_dir, paste0('4c_ImageDimPlot_niches_6_', metadata$SampleName[i], '.pdf')), width=6, height=5)


}



# Plot all tumors together
plot_grid(plotlist = plots_cell_state, ncol = 5)
ggsave(file.path(plot_dir, paste0('5_ALL_ImageDimPlot_cell_type.pdf')), width=20, height=16)

plot_grid(plotlist = plots_malignant, ncol = 5)
ggsave(file.path(plot_dir, paste0('6_ALL_ImageDimPlot_malignant.pdf')), width=20, height=16)

plot_grid(plotlist = plots_niche_4, ncol = 5)
ggsave(file.path(plot_dir, paste0('7_ALL_ImageDimPlot_niche_4.pdf')), width=20, height=16)

plot_grid(plotlist = plots_niche_5, ncol = 5)
ggsave(file.path(plot_dir, paste0('7b_ALL_ImageDimPlot_niche_5.pdf')), width=20, height=16)

plot_grid(plotlist = plots_niche_6, ncol = 5)
ggsave(file.path(plot_dir, paste0('7c_ALL_ImageDimPlot_niche_6.pdf')), width=20, height=16)



# Plot metaprogram distribution  ------------------------------------------------------

order_metaprograms <- c("Cycling",  "Embryonic-like", "Neuroepithelial-like", "Radial-glia-like", 
                        "Embryonic-neuronal-like", 
                        "Neuronal-like" ,"Ependymal-like", "MES-like", 
                        "T-cells", "Myeloid",  "Endothelial",  "Oligodendrocytes", 'Astrocyte', 'Neurons', 'Unknown')

SampleID_order <- c('STEPN-16','STEPN-06','STEPN-10','STEPN-01','STEPN-17','STEPN-15',
                    'STEPN-18','STEPN-19','STEPN-12')

SampleName_vector <- metadata$SampleName
SampleID_vector <- metadata$SampleID

# plot each technical replicate
PlotXeniumMetaprograms(cellid_dir, SampleName_vector, SampleID_vector, SampleID_order = SampleName_vector,
                       order_metaprograms, colors_metaprograms_Xenium, replicate = F)
ggsave(file.path(plot_dir, '8_Metaprogram_proportion.pdf'), width=8, height=5)


# plot averaging technical replicates
PlotXeniumMetaprograms(cellid_dir, SampleName_vector, SampleID_vector, SampleID_order = SampleID_order,
                       order_metaprograms, colors_metaprograms_Xenium, replicate = T)
ggsave(file.path(plot_dir, '9_Metaprogram_proportion_average.pdf'), width=6, height=5)





# To plot zoomed images -----------------------------------------------------
PlotSpatialMapsZoom(SampleName = "STEPN12_Region_3", input_dir, group_display = 'cell_type', 
                    colors = colors_metaprograms_Xenium, size_dot = 0.7, geom_rect = F,  
                    x_min = 1700, x_max = 2300, y_min = 2600, y_max = 3200)
  
ggsave(file.path(plot_dir, paste0('10_ImageDimPlot_malignant_STEPN12_Region_3.pdf')), width=6, height=5)

PlotSpatialMapsZoom(SampleName = "STEPN12_Region_3", input_dir, group_display = 'cell_type', 
                    colors = colors_metaprograms_Xenium, size_dot = 0.7, geom_rect = F,  
                    x_min = 4900, x_max = 5500, y_min = 5000, y_max = 5600)

ggsave(file.path(plot_dir, paste0('10b_ImageDimPlot_malignant_STEPN12_Region_3.pdf')), width=6, height=5)


PlotSpatialMapsZoom(SampleName = "STEPN12_Region_3", input_dir, group_display = 'cell_type', 
                    colors = colors_metaprograms_Xenium, size_dot = 0.7, geom_rect = F,  
                    x_min = 3500, x_max = 5000, y_min = 4400, y_max = 5640)

ggsave(file.path(plot_dir, paste0('10c_ImageDimPlot_malignant_STEPN12_Region_3.pdf')), width=6, height=5)




PlotSpatialMapsZoom(SampleName = "STEPN17_Region_1", input_dir, group_display = 'cell_type', 
                    colors = colors_metaprograms_Xenium, size_dot = 0.7, geom_rect = F,  
                    x_min = 2500, x_max = 3500, y_min = 1500, y_max = 2500)
ggsave(file.path(plot_dir, paste0('11_ImageDimPlot_malignant_STEPN17_Region_1.pdf')), width=6, height=5)




PlotSpatialMapsZoom(SampleName = "STEPN16_Region_1", input_dir, group_display = 'cell_type', 
                    colors = colors_metaprograms_Xenium, size_dot = 0.7, geom_rect = F,  
                    x_min = 1000, x_max = 2500, y_min = 1000, y_max = 2500)
ggsave(file.path(plot_dir, paste0('12_ImageDimPlot_malignant_STEPN16_Region_1.pdf')), width=6, height=5)




PlotSpatialMapsZoom(SampleName = "STEPN12_Region_4", input_dir, group_display = 'cell_type', 
                    colors = colors_metaprograms_Xenium, size_dot = 0.7, geom_rect = F,  
                    x_min = 1200, x_max = 2700, y_min = 1200, y_max = 2700)
ggsave(file.path(plot_dir, paste0('13_ImageDimPlot_malignant_STEPN12_Region_4.pdf')), width=6, height=5)





# Plot example to show MP assignment -----------------------------------------------------

# read data
SampleName = 'STEPN19_Region_2'
data <- qread(file.path(input_dir, paste0('seurat_obj_', SampleName, '.qs')))

# transform all colors in grey except for program of interest
colors <- colors_metaprograms_Xenium

colors <- ifelse(
  names(colors) == "MES-like",
  colors, # Keep original color for MES-like
  "grey80"                   # Set to grey80 for others
)
names(colors) <- names(colors_metaprograms_Xenium)

# plot by cell type
p1 <- ImageDimPlot(data, group.by = 'cell_type',  fov = "fov", cols = colors,  border.size = NA, size = 1, 
                   dark.background = F, axes = F) 
# plot MES like score
p2 <- ImageFeaturePlot(data, features = 'MES-like_UCell',  fov = "fov", cols = c('grey80', '#96410EFF'), 
                   border.size = NA, size = 1, dark.background = F, axes = F, min.cutoff = 'q50', max.cutoff = "q95") 

# create colors for Louvain clusters
colors2 <- c(rep('grey80', 13))
names(colors2) <- as.character(1:13)
colors <- ifelse(
  names(colors2) == "6",
  "#96410EFF", # Keep original color for MES-like
  "grey80"                   # Set to grey80 for others
)
names(colors) <- names(colors2)

# plot by Louvain cluster
p3 <- ImageDimPlot(data, group.by = 'SCT_snn_res.0.4',  fov = "fov", cols = colors,  
                   border.size = NA, size = 1, 
                   dark.background = F, axes = F)



# transform all colors in grey except for program of interest
colors <- colors_metaprograms_Xenium

colors <- ifelse(
  names(colors) == "Endothelial",
  colors, # Keep original color for MES-like
  "grey80"                   # Set to grey80 for others
)
names(colors) <- names(colors_metaprograms_Xenium)

# plot by cell type
p4 <- ImageDimPlot(data, group.by = 'cell_type',  fov = "fov", cols = colors,  border.size = NA, size = 1, 
                   dark.background = F, axes = F) 
# plot endothelial score
colnames(data@meta.data) <- gsub(" ","_", colnames(data@meta.data))
p5 <- ImageFeaturePlot(data, features = 'Endothelial_1_UCell',  fov = "fov", cols = c('grey80', "violetred3" ), 
                       border.size = NA, size = 1, dark.background = F, axes = F, min.cutoff = 'q50', max.cutoff = "q95") 

# create colors for Louvain clusters
colors2 <- c(rep('grey80', 13))
names(colors2) <- as.character(1:13)
colors <- ifelse(
  names(colors2) == "10",
  "violetred3" , # Keep original color for MES-like
  "grey80"                   # Set to grey80 for others
)
names(colors) <- names(colors2)

# plot by Louvain cluster
p6 <- ImageDimPlot(data, group.by = 'SCT_snn_res.0.4',  fov = "fov", cols = colors,  
                   border.size = NA, size = 1, 
                   dark.background = F, axes = F)




# plot combined
plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3,
          rel_widths = c(1, 1, 1))
ggsave(file.path(plot_dir, paste0('14_Example_MP_assignment_STEPN19_Region_2.pdf')), width=15, height=10)




# Plot correlation Xenium to sc/snRNA-seq -----------------------------------------------------
SampleName_vector <- metadata$SampleName
CalculateCorrelationMPXeniumSingleCell(seurat_object_path, SampleName_vector, cellid_dir, colors_metaprograms_Xenium) 
ggsave(file.path(plot_dir, "15_correlation.pdf"), width=5.5, height=3.5) 
