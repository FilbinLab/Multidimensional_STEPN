# Load packages -----------------------------------
library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(ggpubr)
library(data.table)
library(readxl)
library(dplyr)

# Organize environment  -----------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

analysis_dir  <- file.path(base_dir, 'analysis/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41')

plot_dir <- file.path(analysis_dir, 'plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

data_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(data_dir)){dir.create(data_dir, recursive = T)}

resource_dir  <- file.path(base_dir, 'scripts/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))

## Create directories to store individual outputs
if (!dir.exists(file.path(plot_dir, 'individual'))){dir.create(file.path(plot_dir, 'individual'), recursive = T)}

samples <- c('0010540-Region_1', '0010540-Region_2', '0010540-Region_3', '0010540-Region_4', 
             '0010553-Region_1', '0010553-Region_2', '0010553-Region_3', '0010553-Region_4')
names(samples) <- c('0010540-Region_1', '0010540-Region_2', '0010540-Region_3', '0010540-Region_4', 
                    '0010553-Region_1', '0010553-Region_2', '0010553-Region_3', '0010553-Region_4')

for (i in seq_along(samples)) {
  if (!dir.exists(file.path(plot_dir, paste0('individual/', names(samples)[i])))){dir.create(file.path(plot_dir, paste0('individual/', names(samples)[i])), recursive = T)}
}

for (i in seq_along(samples)) {
  if (!dir.exists(file.path(data_dir, paste0('individual/', names(samples)[i])))){dir.create(file.path(data_dir, paste0('individual/', names(samples)[i])), recursive = T)}
}

# Color palettes -------------------------------------
colors_groups <- c("gray30","#F99E93FF","#9E5E9BFF","#74ADD1FF","#ACD39EFF","#96410EFF", 'grey80',
                   '#ef3c2d',  '#ffba08', '#002855')
names(colors_groups) <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
                          "NPC-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
                          "Immune", "Endothelial", "Neurons")

col_niches <- c('#58148e', '#15a2a2', '#ea515f','#d0a03b','#0466c8')


col_normal_malignant <- c('#FFC72CFF', '#582C83FF')
names(col_normal_malignant) <- c('Non-malignant', 'Malignant')



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
resolutions_to_use <- c('SCT_snn_res.0.8', 'SCT_snn_res.0.9', 'SCT_snn_res.1', 'SCT_snn_res.0.3',
                        'SCT_snn_res.0.7', 'SCT_snn_res.0.7', 'SCT_snn_res.0.7', 'SCT_snn_res.0.3')

annotation_clusters <- list (
  '0010540-Region_1' = c('0' = 'NPC-like', 
                         '1' = 'Ependymal-like', 
                         '2' = 'Ependymal-like', 
                         '3' = 'NPC-like', 
                         '4' = 'NPC-like', 
                         '5' = 'NPC-like',
                         '6' = 'NPC-like', 
                         '7' = 'Ependymal-like', 
                         '8' = 'Neuroepithelial-like', 
                         '9' = 'Immune', 
                         '10' = 'NPC-like', 
                         '11' = 'Mesenchymal',
                         '12' = 'NPC-like', 
                         '13' = 'Neuroepithelial-like'),
  '0010540-Region_2' =  c('0' = 'NPC-like', 
                          '1' = 'NPC-like', 
                          '2' = 'Ependymal-like', 
                          '3' = 'Ependymal-like', 
                          '4' = 'NPC-like', 
                          '5' = 'NPC-like',
                          '6' = 'NPC-like', 
                          '7' = 'NPC-like', 
                          '8' = 'Ependymal-like', 
                          '9' = 'Neuroepithelial-like', 
                          '10' = 'Mesenchymal', 
                          '11' = 'Immune',
                          '12' = 'NPC-like', 
                          '13' = 'Neuroepithelial-like', 
                          '14' = 'Neuroepithelial-like', 
                          '15' = 'Immune'),
  '0010540-Region_3' =  c('0' = 'NPC-like', 
                          '1' = 'Ependymal-like', 
                          '2' = 'Mesenchymal', 
                          '3' = 'NPC-like', 
                          '4' = 'NPC-like', 
                          '5' = 'Ependymal-like',
                          '6' = 'NPC-like', 
                          '7' = 'NPC-like', 
                          '8' = 'NPC-like', 
                          '9' = 'Ependymal-like', 
                          '10' = 'Neuroepithelial-like', 
                          '11' = 'Immune',
                          '12' = 'Unassigned', 
                          '13' = 'Neuroepithelial-like', 
                          '14' = 'Unassigned', 
                          '15' = 'Mesenchymal',
                          '16' = 'NPC-like'),
  '0010540-Region_4' =  c('0' = 'Ependymal-like', 
                          '1' = 'Ependymal-like', 
                          '2' = 'NPC-like', 
                          '3' = 'Ependymal-like', 
                          '4' = 'Immune', 
                          '5' = 'Endothelial',
                          '6' = 'Mesenchymal', 
                          '7' = 'NPC-like', 
                          '8' = 'NPC-like', 
                          '9' = 'Immune', 
                          '10' = 'Immune'),
  '0010553-Region_1' = c('0' = 'Immune', 
                         '1' = 'Ependymal-like', 
                         '2' = 'Ependymal-like', 
                         '3' = 'NPC-like', 
                         '4' = 'Ependymal-like', 
                         '5' = 'Neuroepithelial-like',
                         '6' = 'Immune', 
                         '7' = 'NPC-like', 
                         '8' = 'Ependymal-like', 
                         '9' = 'Endothelial', 
                         '10' = 'Immune',
                         '11' = 'Immune',
                         '12' = 'Immune',
                         '13' = 'NPC-like'),
  '0010553-Region_2' = c('0' = 'Ependymal-like', 
                         '1' = 'NPC-like', 
                         '2' = 'Ependymal-like', 
                         '3' = 'Ependymal-like', 
                         '4' = 'Neuroepithelial-like', 
                         '5' = 'NPC-like',
                         '6' = 'NPC-like', 
                         '7' = 'Immune', 
                         '8' = 'NPC-like', 
                         '9' = 'Immune', 
                         '10' = 'Endothelial',
                         '11' = 'NPC-like',
                         '12' = 'Ependymal-like',
                         '13' = 'Immune',
                         '14' = 'NPC-like'),
  '0010553-Region_3' = c('0' = 'NPC-like', 
                         '1' = 'NPC-like', 
                         '2' = 'NPC-like', 
                         '3' = 'Ependymal-like', 
                         '4' = 'Mesenchymal', 
                         '5' = 'Mesenchymal',
                         '6' = 'NPC-like', 
                         '7' = 'Ependymal-like', 
                         '8' = 'NPC-like', 
                         '9' = 'Neuroepithelial-like', 
                         '10' = 'Neuroepithelial-like',
                         '11' = 'Immune',
                         '12' = 'NPC-like',
                         '13' = 'Neuroepithelial-like'),
  '0010553-Region_4' = c('0' = 'NPC-like', 
                         '1' = 'Ependymal-like', 
                         '2' = 'Ependymal-like', 
                         '3' = 'Endothelial', 
                         '4' = 'Immune', 
                         '5' = 'Neuroepithelial-like',
                         '6' = 'NPC-like', 
                         '7' = 'Immune', 
                         '8' = 'Immune', 
                         '9' = 'NPC-like')
)




# Process Xenium data  -------------------------------------

# Read data
for (i in seq_along(samples)) {
  data <- qread(file.path(base_dir, paste0('processed_data_Carlos/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/', names(samples)[i], '/', names(samples)[i], '.qs')))
  
  # DotPlot markers for different resolutions
  nres <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
  for (res in nres) {
    Idents(data) <- paste0("SCT_snn_res.", res)
    # DotPlot markers
    DotPlot(data, cols = c('lightgrey', 'red'), features = gene_list, assay = 'SCT', scale = TRUE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 6))
    ggsave(file.path(plot_dir, paste0('individual/', names(samples)[i], '/0_DotPlot_res', res ,'.pdf')), width=14, height=5)
  }
  
  # Change name identities
  Idents(data) <- resolutions_to_use[i]
  data <- RenameIdents(data, annotation_clusters[[i]])
  data[["group"]] <- Idents(data)
  
  # Add information about malignant or non-malignant
  data$malignant <- ifelse(data$group %in% c('NPC-like', 'Ependymal-like', 'Neuroepithelial-like',
                                             'Mesenchymal'), "Malignant", "Non-malignant")
  
  # DotPlot markers
  DotPlot(data, 
          features = c(gene_list, 'ZFTA-RELA-Fusion1'), 
          assay = 'SCT', scale = TRUE) + 
    paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 8))
  ggsave(file.path(plot_dir, paste0('individual/', names(samples)[i], '/1_DotPlot_final.pdf')), width=13, height=4)
  
  # Export metadata with cell names and new annotation
  cell_id <- data.frame(rownames(data@meta.data), data@meta.data$group)
  rownames(cell_id) <- cell_id$rownames.data.meta.data.
  cell_id$rownames.data.meta.data. <- NULL
  write_csv(cell_id, file.path(data_dir, paste0('individual/', names(samples)[i], '/cell_id.csv')))
  
  
  # Visualize  distribution clusters 
  ImageDimPlot(data, group.by = 'group', cols = colors_groups, border.size = NA, size = 0.4, 
               dark.background = F) + ggtitle("Clusters")
  ggsave(file.path(plot_dir, paste0('individual/', names(samples)[i], '/2_ImageDimPlot.pdf')), width=8, height=5)
  
  ImageDimPlot(data, group.by = 'malignant', cols = col_normal_malignant, border.size = NA, size = 0.4, 
               dark.background = F) + ggtitle("Malignancy")
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
  
}

