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

resource_dir  <- file.path(base_dir, 'scripts/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))

analysis_dir  <- file.path(base_dir, 'analysis/all_Xenium_runs')

plot_dir <- file.path(analysis_dir, 'plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

data_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(data_dir)){dir.create(data_dir, recursive = T)}

seurat_obj_dir <- file.path(data_dir, 'seurat_objects')
if (!dir.exists(seurat_obj_dir)){dir.create(seurat_obj_dir, recursive = T)}


metadata <- read_xlsx(file.path(base_dir, 'SampleIdentifier.xlsx'))

# reorder by sampleID
metadata <- metadata %>% arrange(SampleID)


# Plot distribution clusters  -------------------------------------
plots_cell_state <- list()
plots_malignant <- list()

# Read data
for (i in seq_along(metadata$Sample)) {
  data <- qread(file.path(seurat_obj_dir, paste0('seurat_obj_', metadata$SampleID[i], "_", metadata$SampleName[i], '.qs')))
   
  # # DotPlot markers   ------------------------------------------------------
  # DotPlot(data,
  #         features = c('ZFTA-RELA-Fusion1'),
  #         assay = 'SCT', scale = TRUE) +
  #   paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 12))
  # ggsave(file.path(plot_dir, paste0('1_ZFTARELA_DotPlot_', metadata$SampleID[i], "_", metadata$SampleName[i],'.pdf')), width=6, height=4)
  # 
  # # Export metadata with cell names and new annotation
  # cell_id <- data.frame(rownames(data@meta.data), data@meta.data$group)
  # # rename column names
  # colnames(cell_id)[colnames(cell_id) == "rownames.data.meta.data."] <- "cell_id"
  # colnames(cell_id)[colnames(cell_id) == "data.meta.data.group"] <- "group"
  # # save
  # write_csv(cell_id, file.path(data_dir, paste0('cell_ID_', metadata$SampleID[i], "_", metadata$SampleName[i],'.csv')))


  # Visualize  distribution clusters    ------------------------------------------------------
  plots_cell_state[[i]] <- ImageDimPlot(data, group.by = 'group', cols = colors_metaprograms_Xenium,
                                      border.size = NA, size = 0.5, dark.background = F)  + NoLegend() + 
    labs(title = paste0(metadata$SampleID[i],"_", metadata$SampleName[i]))
 # ggsave(file.path(plot_dir, paste0('2_ImageDimPlot_', metadata$SampleID[i], "_", metadata$SampleName[i],'.pdf')), width=6, height=5)

  # ImageDimPlot(data, group.by = 'malignant', cols = col_normal_malignant,
  #              border.size = NA, size = 0.5, dark.background = F) + NoLegend()
  # ggsave(file.path(plot_dir, paste0('3_ImageDimPlot_malignant_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)
  
  plots_malignant[[i]] <- ImageDimPlot(data, group.by = 'malignant', cols = col_normal_malignant, 
                                     border.size = NA, size = 0.5, dark.background = F) + NoLegend() + 
    labs(title = paste0(metadata$SampleID[i],"_", metadata$SampleName[i]))
}
  
  
# Plot all tumors together
plot_grid(plotlist = plots_cell_state, ncol = 5) 
ggsave(file.path(plot_dir, paste0('2_ImageDimPlot_all.pdf')), width=20, height=16)

plot_grid(plotlist = plots_malignant, ncol = 5) 
ggsave(file.path(plot_dir, paste0('2_ImageDimPlot_all_malignant.pdf')), width=20, height=16)




# Metaprogram distribution  ------------------------------------------------------

metaprogram_frequency <- list()


for (i in seq_along(metadata$Sample)) { 
  
  # read data
  data <- qread(file.path(seurat_obj_dir, paste0('seurat_obj_', metadata$SampleID[i], "_", metadata$SampleName[i], '.qs')))
  print(paste0("Reading ", metadata$Sample[i]))
 
  # Save data table with frequency
  metaprogram_frequency[[i]] <- data@meta.data %>%
    group_by(group) %>%
    summarise(n = n()) %>%
    mutate(freq = (n/sum(n)*100))
}

names(metaprogram_frequency) <- metadata$SampleName  

# bind rows together and transform into df
metaprogram_proportion <- bind_rows(metaprogram_frequency, .id = "Xenium_Region")
metaprogram_proportion <- as.data.frame(metaprogram_proportion)

# add sample name
matching_indices <- match(metaprogram_proportion$Xenium_Region, metadata$SampleName)

metaprogram_proportion$SampleID <- metadata$SampleID[matching_indices]
metaprogram_proportion$FullSampleID <- paste0(metaprogram_proportion$SampleID, "_", metaprogram_proportion$Xenium_Region)


# reorder so that it is the same as in the paper
metaprogram_proportion$group <- factor(metaprogram_proportion$group, 
                                levels = c('Neuroepithelial-like', "Radial glia-like", 'Neuronal-like',
                                           "Ependymal-like", "Mesenchymal", 'Myeloid', 'T-cell', 'Endothelial', 'VLMCs', 'Microglia',
                                           'Oligodendrocytes', 'Unassigned'))


# reorder by sampleID
metaprogram_proportion <- metaprogram_proportion %>% arrange(SampleID)

# Plot a stacked bar plot
ggplot(metaprogram_proportion, aes(x = FullSampleID, y = freq, fill = group)) +
  scale_fill_manual(values = colors_groups_barplot) +
  geom_bar(stat = "identity", position = "fill", color="black") +
  labs (y='Proportion', x='') + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=14, angle = 90, vjust = 0.5, hjust=1, colour="black"),
        axis.text.y = element_text(size=14, colour="black"),
        axis.title=element_text(size=14),
        legend.text = element_text(size = 14),
        legend.title = element_blank()) 
ggsave(file.path(plot_dir, '7_Metaprogram_proportion.pdf'), width=10, height=5)


# group by sample (average technical replicates)
metaprogram_proportion2 <- metaprogram_proportion %>% 
  group_by(group, SampleID) %>% 
  summarise(avg_frequency = mean(freq))

order <- c('STEPN-14','STEPN-10','STEPN-16','STEPN-06','STEPN-12','STEPN-15',
           'STEPN-18','STEPN-17','STEPN-01','STEPN-19')
metaprogram_proportion2$SampleID <- factor(metaprogram_proportion2$SampleID , levels = order)


ggplot(metaprogram_proportion2, aes(x = SampleID, y = avg_frequency, fill = group)) +
  scale_fill_manual(values = colors_groups_barplot) +
  geom_bar(stat = "identity", position = "fill", color="black") +
  labs (y='Proportion', x='') + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=14, angle = 90, vjust = 0.5, hjust=1, colour="black"),
        axis.text.y = element_text(size=14, colour="black"),
        axis.title=element_text(size=14),
        legend.text = element_text(size = 14),
        legend.title = element_blank()) 
ggsave(file.path(plot_dir, '7_Metaprogram_proportion_average.pdf'), width=6, height=4)

write.csv(metaprogram_proportion, file.path(data_dir, 'metaprogram_frequency.csv'))




# To plot zoomed images -----------------------------------------------------
i = 14
data <- qread(file.path(seurat_obj_dir, paste0('seurat_obj_', metadata$SampleID[i], "_", metadata$SampleName[i], '.qs')))

ImageDimPlot(data, group.by = 'group',  fov = "fov", cols = colors_metaprograms_Xenium,  border.size = NA, size = 0.5,
             dark.background = T, axes = TRUE, #border.color = "white", border.size = 0.1, coord.fixed = FALSE, nmols = 10000
) + 
  #scale_x_continuous(breaks = seq(0,5000,500)) + 
  #scale_y_continuous(breaks = seq(0,5000,500)) + 
  geom_rect(xmin = 3500, xmax = 5000, ymin = 4400, ymax = 5400, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 4900, xmax = 5500, ymin = 5000, ymax = 5600, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 1500, xmax = 2000, ymin = 2800, ymax = 3300, color = 'black', fill = NA, lwd = 0.2) + 
  seurat_theme()
ggsave(file.path(plot_dir, paste0('4_ImageDimPlot_malignant_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=6, height=5)


data[["zoom_1"]] <- Crop(data[["fov"]], x = c(1700, 2300), y = c(2600, 3200), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1", size = 1, border.color = "white", border.size = 0.1, cols = colors_metaprograms_Xenium) + seurat_theme()
ggsave(file.path(plot_dir, paste0('4_Zoom1_ImageDimPlot_celltype_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)

data[["zoom_1"]] <- Crop(data[["fov"]], x = c(4900, 5500), y = c(5000, 5600), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1",  size = 1, border.color = "white", border.size = 0.1, cols = colors_metaprograms_Xenium) + seurat_theme()
ggsave(file.path(plot_dir, paste0('4_Zoom2_ImageDimPlot_niche_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)

data[["zoom_1"]] <- Crop(data[["fov"]], x = c(3500, 5000), y = c(4400, 5400), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1",  size = 1, border.color = "white", border.size = 0.1, cols = colors_metaprograms_Xenium) + seurat_theme()
ggsave(file.path(plot_dir, paste0('4_Zoom3_ImageDimPlot_celltype_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)

data[["zoom_1"]] <- Crop(data[["fov"]], x = c(3500, 5000), y = c(4400, 5400), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1", group.by = 'niches', size = 1, border.color = "white", border.size = 0.1, cols = colors_to_use_niches) + seurat_theme()
ggsave(file.path(plot_dir, paste0('5_Zoom3_ImageDimPlot_niche_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)


i = 20
data <- qread(file.path(seurat_obj_dir, paste0('seurat_obj_', metadata$SampleID[i], "_", metadata$SampleName[i], '.qs')))

ImageDimPlot(data, group.by = 'group',  fov = "fov", cols = colors_metaprograms_Xenium,  border.size = NA, size = 0.5,
             dark.background = T, axes = TRUE, #border.color = "white", border.size = 0.1, coord.fixed = FALSE, nmols = 10000
) + 
  #scale_x_continuous(breaks = seq(0,5000,500)) + 
  #scale_y_continuous(breaks = seq(0,5000,500)) + 
  geom_rect(xmin = 2500, xmax = 3500, ymin = 1500, ymax = 2500, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 4900, xmax = 5500, ymin = 5000, ymax = 5600, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 1500, xmax = 2000, ymin = 2800, ymax = 3300, color = 'black', fill = NA, lwd = 0.2) + 
  seurat_theme()
ggsave(file.path(plot_dir, paste0('5_ImageDimPlot_malignant_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=6, height=5)


data[["zoom_1"]] <- Crop(data[["fov"]], x = c(2500, 3500), y = c(1500, 2500), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1", size = 1, border.color = "white", border.size = 0.1, cols = colors_metaprograms_Xenium) + seurat_theme()
ggsave(file.path(plot_dir, paste0('5_Zoom1_ImageDimPlot_celltype_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)

data[["zoom_1"]] <- Crop(data[["fov"]], x = c(2500, 3500), y = c(1500, 2500), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1", group.by = 'niches', size = 1, border.color = "white", border.size = 0.1, cols = colors_to_use_niches) + seurat_theme()
ggsave(file.path(plot_dir, paste0('5_Zoom1_ImageDimPlot_niche_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)



i = 19
data <- qread(file.path(seurat_obj_dir, paste0('seurat_obj_', metadata$SampleID[i], "_", metadata$SampleName[i], '.qs')))

ImageDimPlot(data, group.by = 'group',  fov = "fov", cols = colors_metaprograms_Xenium,  border.size = NA, size = 0.5,
             dark.background = T, axes = TRUE, #border.color = "white", border.size = 0.1, coord.fixed = FALSE, nmols = 10000
) + 
  #scale_x_continuous(breaks = seq(0,5000,500)) + 
  #scale_y_continuous(breaks = seq(0,5000,500)) + 
  geom_rect(xmin = 1000, xmax = 2500, ymin = 1000, ymax = 2500, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 4900, xmax = 5500, ymin = 5000, ymax = 5600, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 1500, xmax = 2000, ymin = 2800, ymax = 3300, color = 'black', fill = NA, lwd = 0.2) + 
  seurat_theme()
ggsave(file.path(plot_dir, paste0('6_ImageDimPlot_malignant_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=6, height=5)


data[["zoom_1"]] <- Crop(data[["fov"]], x = c(1000, 2500), y = c(1000, 2500), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1", size = 1, border.color = "white", border.size = 0.1, cols = colors_metaprograms_Xenium) + seurat_theme()
ggsave(file.path(plot_dir, paste0('6_Zoom1_ImageDimPlot_celltype_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)


data[["zoom_1"]] <- Crop(data[["fov"]], x = c(1000, 2500), y = c(1000, 2500), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1",group.by = 'niches', size = 1, border.color = "white", border.size = 0.1, cols = colors_to_use_niches) + seurat_theme()
ggsave(file.path(plot_dir, paste0('6_Zoom1_ImageDimPlot_niche_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)



i = 15
data <- qread(file.path(seurat_obj_dir, paste0('seurat_obj_', metadata$SampleID[i], "_", metadata$SampleName[i], '.qs')))

ImageDimPlot(data, group.by = 'group',  fov = "fov", cols = colors_metaprograms_Xenium,  border.size = NA, size = 0.5,
             dark.background = T, axes = TRUE, #border.color = "white", border.size = 0.1, coord.fixed = FALSE, nmols = 10000
) + 
  #scale_x_continuous(breaks = seq(0,5000,500)) + 
  #scale_y_continuous(breaks = seq(0,5000,500)) + 
  geom_rect(xmin = 1200, xmax = 2700, ymin = 1200, ymax = 2700, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 4900, xmax = 5500, ymin = 5000, ymax = 5600, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 1500, xmax = 2000, ymin = 2800, ymax = 3300, color = 'black', fill = NA, lwd = 0.2) + 
  seurat_theme()
ggsave(file.path(plot_dir, paste0('7_ImageDimPlot_malignant_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=6, height=5)


data[["zoom_1"]] <- Crop(data[["fov"]], x = c(1200, 2700), y = c(1200, 2700), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1", size = 1, border.color = "white", border.size = 0.1, cols = colors_metaprograms_Xenium) + seurat_theme()
ggsave(file.path(plot_dir, paste0('7_Zoom1_ImageDimPlot_celltype_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)


data[["zoom_1"]] <- Crop(data[["fov"]], x = c(1200, 2700), y = c(1200, 2700), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1",group.by = 'niches', size = 1, border.color = "white", border.size = 0.1, cols = colors_to_use_niches) + seurat_theme()
ggsave(file.path(plot_dir, paste0('7_Zoom1_ImageDimPlot_niche_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)
