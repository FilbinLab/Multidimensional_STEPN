# Load packages -----------------------------------
library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(data.table)
library(readxl)
library(dplyr)
library(ggrastr)

# Organize environment  -----------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

analysis_dir  <- file.path(base_dir, 'analysis/20231102__215055__7EP1_7EP41_3EP8')

plot_dir <- file.path(analysis_dir, 'plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

data_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(data_dir)){dir.create(data_dir, recursive = T)}

seurat_obj_dir <- file.path(analysis_dir, 'data/seurat_objects')
if (!dir.exists(seurat_obj_dir)){dir.create(seurat_obj_dir, recursive = T)}

resource_dir  <- file.path('/n/scratch/users/s/sad167/EPN/Ependymoma2023/Xenium/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))

## Create directories to store individual outputs
if (!dir.exists(file.path(plot_dir, 'individual'))){dir.create(file.path(plot_dir, 'individual'), recursive = T)}

samples <- c('0010575-Region_1', '0010575-Region_2', '0010575-Region_3', 
             '0010619-Region_1', '0010619-Region_2', '0010619-Region_3', '0010619-Region_4',
             '0010619-Region_5')
names(samples) <- c('0010575-Region_1', '0010575-Region_2', '0010575-Region_3', 
                    '0010619-Region_1', '0010619-Region_2', '0010619-Region_3', '0010619-Region_4',
                    '0010619-Region_5')

for (i in seq_along(samples)) {
  if (!dir.exists(file.path(plot_dir, paste0('individual/', names(samples)[i])))){dir.create(file.path(plot_dir, paste0('individual/', names(samples)[i])), recursive = T)}
}

for (i in seq_along(samples)) {
  if (!dir.exists(file.path(data_dir, paste0('individual/', names(samples)[i])))){dir.create(file.path(data_dir, paste0('individual/', names(samples)[i])), recursive = T)}
}



# Load gene list of Xenium panel -------------------------------------
Xenium_panel <- read_excel(file.path(base_dir, 'Xenium_panel.xlsx'))
unique(Xenium_panel$Program)

#EPN_programs <- c("Radial glia/NPC-like","Radial glia-like", "NPC-like-2" ,"NPC-like-1" ,"MES_Hypoxia","Glio_Angiogenesis", "Ependymal-like-2","Ependymal-like-1","Ependymal")

Xenium_panel_EPN <- Xenium_panel %>% filter(Type == 'Ependymoma' | Type == 'Intersection')

Xenium_panel_EPN <- Xenium_panel_EPN %>% group_by(Program) %>% 
  summarise(merged_values = paste(Gene, collapse = ",")) 

Xenium_panel_EPN_list <- as.list(Xenium_panel_EPN$merged_values)
names(Xenium_panel_EPN_list) <- Xenium_panel_EPN$Program


# split markers (which are now separated by ",") into vector
gene_list <- list()
for (i in seq_along(Xenium_panel_EPN_list)) {
  gene_list[[i]] <- unlist(strsplit(Xenium_panel_EPN_list[[i]], ","))
}
names(gene_list) <- names(Xenium_panel_EPN_list)

# combine list into new programs
gene_list$`Ependymal-like` <- c(gene_list$`Ependymal`, gene_list$`Ependymal-like-1`)
gene_list$`Neuroepithelial-like` <- c(gene_list$`Ependymal-like-2`, gene_list$Glio_Angiogenesis)
gene_list$`NPC-like` <- c(gene_list$`NPC-like-1`, gene_list$`NPC-like-2`)

gene_list$`Ependymal` <- NULL
gene_list$`Ependymal-like-1` <- NULL
gene_list$`Ependymal-like-2` <- NULL
gene_list$Glio_Angiogenesis <- NULL
gene_list$`NPC-like-1` <- NULL
gene_list$`NPC-like-2` <- NULL

# save Xenium gene list (for tumor subpopulations)
saveRDS(gene_list, file.path(data_dir, 'Xenium_tumor_gene_list.rds'))




# Define resolution to use and annotations for each tissue  -------------------------------------
resolutions_to_use <- c('SCT_snn_res.0.5', 'SCT_snn_res.0.5', 'SCT_snn_res.0.5', 
                        'SCT_snn_res.0.5', 'SCT_snn_res.0.5', 'SCT_snn_res.0.5', 'SCT_snn_res.0.5', 'SCT_snn_res.0.7')

annotation_clusters <- list (
  '0010575-Region_1' = c('0' = 'Ependymal-like',
                         '1' = 'Ependymal-like',
                         '2' = 'NPC-like',
                         '3' = 'Neuroepithelial-like', 
                         '4' = 'Mesenchymal',
                         '5' = 'Endothelial',
                         '6' = 'Myeloid',
                         '7' = 'NPC-like', 
                         '8' = 'NPC-like', 
                         '9' = 'Myeloid', 
                         '10' = 'T-cell', 
                         '11' = 'Myeloid', 
                         '12' = 'NPC-like'),
  '0010575-Region_2' = c('0' = 'Ependymal-like',
                         '1' = 'NPC-like',
                         '2' = 'Ependymal-like',
                         '3' = 'Neuroepithelial-like', 
                         '4' = 'Mesenchymal',
                         '5' = 'Myeloid',
                         '6' = 'Endothelial',
                         '7' = 'NPC-like', 
                         '8' = 'NPC-like', 
                         '9' = 'Myeloid', 
                         '10' = 'T-cell', 
                         '11' = 'Myeloid', 
                         '12' = 'NPC-like'),
  '0010575-Region_3' = c('0' = 'Ependymal-like',
                         '1' = 'NPC-like',
                         '2' = 'Ependymal-like',
                         '3' = 'Mesenchymal', 
                         '4' = 'Myeloid',
                         '5' = 'Ependymal-like',
                         '6' = 'Neuroepithelial-like',
                         '7' = 'Endothelial', 
                         '8' = 'NPC-like', 
                         '9' = 'NPC-like', 
                         '10' = 'Myeloid', 
                         '11' = 'T-cell', 
                         '12' = 'Myeloid',
                         '13' = 'NPC-like'),
  '0010619-Region_1' = c('0' = 'Ependymal-like',
                         '1' = 'NPC-like',
                         '2' = 'Ependymal-like',
                         '3' = 'Neuroepithelial-like', 
                         '4' = 'Myeloid',
                         '5' = 'VLMCs',
                         '6' = 'VLMCs',
                         '7' = 'Myeloid', 
                         '8' = 'VLMCs', 
                         '9' = 'Ependymal-like', 
                         '10' = 'NPC-like', 
                         '11' = 'Myeloid', 
                         '12' = 'Myeloid',
                         '13' = 'NPC-like'),
  '0010619-Region_2' = c('0' = 'Ependymal-like',
                         '1' = 'Ependymal-like',
                         '2' = 'NPC-like',
                         '3' = 'Myeloid', 
                         '4' = 'Neuroepithelial-like',
                         '5' = 'VLMCs',
                         '6' = 'VLMCs',
                         '7' = 'VLMCs', 
                         '8' = 'Ependymal-like', 
                         '9' = 'NPC-like',
                         '10' = 'NPC-like'),
  '0010619-Region_3' = c('0' = 'Ependymal-like',
                         '1' = 'NPC-like',
                         '2' = 'Ependymal-like',
                         '3' = 'VLMCs', 
                         '4' = 'Ependymal-like',
                         '5' = 'Myeloid',
                         '6' = 'Neuroepithelial-like',
                         '7' = 'VLMCs', 
                         '8' = 'NPC-like', 
                         '9' = 'Radial glia-like',
                         '10' = 'NPC-like'),
  '0010619-Region_4' = c('0' = 'Ependymal-like',
                         '1' = 'Endothelial',
                         '2' = 'NPC-like',
                         '3' = 'Myeloid', 
                         '4' = 'NPC-like',
                         '5' = 'Ependymal-like',
                         '6' = 'Ependymal-like',
                         '7' = 'VLMCs', 
                         '8' = 'T-cell', 
                         '9' = 'NPC-like',
                         '10' = 'NPC-like'),
  '0010619-Region_5' = c('0' = 'Ependymal-like',
                         '1' = 'Myeloid',
                         '2' = 'Ependymal-like',
                         '3' = 'NPC-like', 
                         '4' = 'NPC-like',
                         '5' = 'Ependymal-like',
                         '6' = 'NPC-like',
                         '7' = 'Endothelial', 
                         '8' = 'NPC-like', 
                         '9' = 'Endothelial',
                         '10' = 'T-cell',
                         '11' = 'Myeloid',
                         '12' = 'NPC-like',
                         '13' = 'VLMCs',
                         '14' = 'Endothelial')
)




# Process Xenium data  -------------------------------------

# Read data
for (i in seq_along(samples)) {
  data <- qread(file.path(base_dir, paste0('data/processed_data_Carlos/20231102__215055__7EP1_7EP41_3EP8/', names(samples)[i], '/', names(samples)[i], '.qs')))
  
  #DotPlot markers for different resolutions
  nres <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
  for (res in nres) {
    Idents(data) <- paste0("SCT_snn_res.", res)
    # DotPlot markers
    DotPlot(data, cols = c('lightgrey', 'red'), features = gene_list, assay = 'SCT', scale = TRUE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 6))
    ggsave(file.path(plot_dir, paste0('individual/', names(samples)[i], '/0_DotPlot_res', res ,'.pdf')), width=16, height=5)
  }

  # Change name identities
  Idents(data) <- resolutions_to_use[i]
  data <- RenameIdents(data, annotation_clusters[[i]])
  data[["group"]] <- Idents(data)
  
  # Add information about malignant or non-malignant
  # Add information about malignant or non-malignant
  data$malignant <- ifelse(data$group %in% c('NPC-like', 'Ependymal-like', 'Neuroepithelial-like',
                                             'Mesenchymal', 'Radia glia-like'), "Malignant", "Non-malignant")
  
  # DotPlot markers
  DotPlot(data, 
          features = c(gene_list, 'ZFTA-RELA-Fusion1'), 
          assay = 'SCT', scale = TRUE) + 
    paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 8))
  ggsave(file.path(plot_dir, paste0('individual/', names(samples)[i], '/1_DotPlot_final.pdf')), width=16, height=3)
  
  # Export metadata with cell names and new annotation
  cell_id <- data.frame(rownames(data@meta.data), data@meta.data$group)
  # rename column names
  colnames(cell_id)[colnames(cell_id) == "rownames.data.meta.data."] <- "cell_id"
  colnames(cell_id)[colnames(cell_id) == "data.meta.data.group"] <- "group"
  # save
  write_csv(cell_id, file.path(data_dir, paste0('individual/cell_ID_', names(samples)[i], '.csv')))
  
  
  
  # Visualize  distribution clusters 
  plot <- ImageDimPlot(data, group.by = 'group', cols = colors_groups, border.size = NA, size = 0.4, 
                       dark.background = F) + ggtitle(names(samples)[i]) + NoLegend()
  rasterize(plot, layers='Point', dpi=800)
  ggsave(file.path(plot_dir, paste0('individual/', names(samples)[i], '/2_ImageDimPlot.pdf')), width=6, height=5)
  
  plot <- ImageDimPlot(data, group.by = 'malignant', cols = col_normal_malignant, border.size = NA, size = 0.4, 
                       dark.background = F) + ggtitle(names(samples)[i])
  rasterize(plot, layers='Point', dpi=800)
  ggsave(file.path(plot_dir, paste0('individual/', names(samples)[i], '/3_ImageDimPlot_malignant.pdf')), width=8, height=5)
  

  # Perform niche analysis
  data <- BuildNicheAssay(object = data, fov = "fov", group.by = "group", niches.k = 3, neighbors.k = 30)

  # Plot niche images
  celltype.plot <- ImageDimPlot(data, group.by = 'group', fov = "fov",  cols = colors_groups, border.size = NA, size = 0.4,
                                dark.background = F) + ggtitle("Cell type")
  niche.plot <- ImageDimPlot(data, group.by = "niches", fov = "fov",  cols = col_niches, border.size = NA, size = 0.4,
                             dark.background = F) + ggtitle("Niches")
  celltype.plot | niche.plot
  ggsave(file.path(plot_dir, paste0('individual/', names(samples)[i], '/4_ImageDimPlot_niche.pdf')), width=16, height=5)

  # Save data table with frequency
  niche_freq <- as.data.frame(t(table(data$group, data$niches)))
  write_csv(niche_freq, file.path(data_dir, paste0('individual/', names(samples)[i], '/3_0010540-Region_1_niche_cell_number.csv')))

  # Plot niche frequency
  plot_bar(data, data$niches, data$group, colors_groups)
  ggsave(file.path(plot_dir, paste0('individual/', names(samples)[i], '/5_BarPlot_niches.pdf')), width=6, height=5)
  
  # save annotated object
  qsave(data, file.path(seurat_obj_dir, paste0('Seurat_obj', names(samples)[i], '.qs')))
  
}

