# Load packages -----------------------------------
library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(ggpubr)
library(data.table)
library(readxl)

# Organize environment  -----------------------------------
#base_dir <- file.path('/n/scratch/users/c/cao385/Immune/Xenium/results/20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/0010540-Region_1')
base_dir <- file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/12 - Xenium')

plot_dir <- file.path(base_dir, 'analysis/plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

data_dir <- file.path(base_dir, 'analysis/data')
if (!dir.exists(data_dir)){dir.create(data_dir, recursive = T)}

## Create directories to store individual outputs
if (!dir.exists(file.path(plot_dir, 'individual'))){dir.create(file.path(plot_dir, 'individual'), recursive = T)}

samples <- c('0010540-Region_1', '0010540-Region_2', '0010540-Region_3', '0010540-Region_4', 
             '0010553-Region_1', '0010553-Region_2', '0010553-Region_3', '0010553-Region_4')
names(samples) <- c('0010540-Region_1', '0010540-Region_2', '0010540-Region_3', '0010540-Region_4', 
             '0010553-Region_1', '0010553-Region_2', '0010553-Region_3', '0010553-Region_4')

for (i in seq_along(samples)) {
  if (!dir.exists(file.path(plot_dir, paste0('individual/', names(samples)[i])))){dir.create(file.path(plot_dir, paste0('individual/', names(samples)[i])), recursive = T)}
  
}
    


# Load gene list of Xenium panel -------------------------------------
Xenium_panel <- read_excel(file.path(base_dir, 'Xenium_panel.xlsx'))
unique(Xenium_panel$Program)

#EPN_programs <- c("Radial glia/NPC-like","Radial glia-like", "NPC-like-2" ,"NPC-like-1" ,"MES_Hypoxia","Glio_Angiogenesis", "Ependymal-like-2","Ependymal-like-1","Ependymal")

Xenium_panel_EPN <- Xenium_panel %>% filter(Type == 'Ependymoma')

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



# Load Xenium output 0010540-Region_1 (processed by Carlos) -------------------------------------

# Read data
data <- qread(file.path(base_dir, 'raw_data/0010540-Region_1.qs'))

# DotPlot markers for different resolutions
nres <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
for (res in nres) {
  Idents(data) <- paste0("SCT_snn_res.", res)
  # DotPlot markers
  DotPlot(data, cols = c('lightgrey', 'red'), features = gene_list, assay = 'SCT', scale = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 6))
  ggsave(file.path(plot_dir, paste0("individual/0010540-Region_1/1_DotPlot_res", res ,".pdf")), width=14, height=5)
}


# Change name identities
Idents(data) <- data$SCT_snn_res.0.8
data <- RenameIdents(data, 
                     '0' = 'NPC-like', 
                     '1' = 'Ependymal-like', 
                     '2' = 'Ependymal-like', 
                     '3' = 'NPC-like', 
                     '4' = 'Unassigned', 
                     '5' = 'Unassigned',
                     '6' = 'NPC-like', 
                     '7' = 'Ependymal-like', 
                     '8' = 'Neuroepithelial-like', 
                     '9' = 'Immune', 
                     '10' = 'Unassigned', 
                     '11' = 'Mesenchymal',
                     '12' = 'NPC-like', 
                     '13' = 'Neuroepithelial-like')

data[["group"]] <- Idents(data)

# DotPlot markers
DotPlot(data, 
        features = c(gene_list, 'ZFTA-RELA-Fusion1'), 
        assay = 'SCT', scale = TRUE) + 
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 8))
ggsave(file.path(plot_dir, "/1_DotPlot_0010540-Region_1.pdf"), width=13, height=4)

# Export metadata with cell names and new annotation
cell_id <- data@meta.data %>% select(group)
cell_id <-  cell_id %>%
  rownames_to_column(var = "cell_id")

write_csv(cell_id, file.path(data_dir, "0010540-Region_1_cell_id.csv"))

# Export table with total number of cells
table <- data@meta.data %>% 
  group_by(group) %>% summarize(n())
write_csv(table, file.path(data_dir, "2_cell_number_0010540-Region_1.csv"))


# Load Xenium output 0010540-Region_2 (processed by Carlos) -------------------------------------

# Read data
data <- qread(file.path(base_dir, 'raw_data/0010540-Region_2.qs'))

# DotPlot markers for different resolutions
nres <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
for (res in nres) {
  Idents(data) <- paste0("SCT_snn_res.", res)
  # DotPlot markers
  DotPlot(data, cols = c('lightgrey', 'red'), features = gene_list, assay = 'SCT', scale = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 6))
  ggsave(file.path(plot_dir, paste0("individual/0010540-Region_2/1_DotPlot_res", res ,".pdf")), width=14, height=5)
}

# Change name identities
Idents(data) <- data$SCT_snn_res.0.9
data <- RenameIdents(data, 
                     '0' = 'NPC-like', 
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
                     '15' = 'Immune')

data[["group"]] <- Idents(data)

# DotPlot markers
DotPlot(data, 
        features = c(gene_list, 'ZFTA-RELA-Fusion1'), 
        assay = 'SCT', scale = TRUE) + 
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 8))
ggsave(file.path(plot_dir, "1_DotPlot_0010540-Region_2.pdf"), width=13, height=4)

# Export metadata with cell names and new annotation
cell_id <- data@meta.data %>% select(group)
cell_id <-  cell_id %>%
  rownames_to_column(var = "cell_id")

write_csv(cell_id, file.path(data_dir, "0010540-Region_2_cell_id.csv"))

# Export table with total number of cells
table <- data@meta.data %>% 
  group_by(group) %>% summarize(n())
write_csv(table, file.path(data_dir, "2_cell_number_0010540-Region_2.csv"))



# Load Xenium output 0010540-Region_3 (processed by Carlos) -------------------------------------

# Read data
data <- qread(file.path(base_dir, 'raw_data/0010540-Region_3.qs'))

# DotPlot markers for different resolutions
nres <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
for (res in nres) {
  Idents(data) <- paste0("SCT_snn_res.", res)
  # DotPlot markers
  DotPlot(data, cols = c('lightgrey', 'red'), features = gene_list, assay = 'SCT', scale = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 6))
  ggsave(file.path(plot_dir, paste0("individual/0010540-Region_3/1_DotPlot_res", res ,".pdf")), width=14, height=5)
}


# Change name identities
Idents(data) <- data$SCT_snn_res.1
data <- RenameIdents(data, 
                     '0' = 'NPC-like', 
                     '1' = 'Ependymal-like', 
                     '2' = 'Mesenchymal', 
                     '3' = 'Tumor-c24', 
                     '4' = 'Tumor-c24', 
                     '5' = 'Ependymal-like',
                     '6' = 'NPC-like', 
                     '7' = 'Normal neurons', 
                     '8' = 'NPC-like', 
                     '9' = 'Ependymal-like', 
                     '10' = 'Neuroepithelial-like', 
                     '11' = 'Immune',
                     '12' = 'Unassigned', 
                     '13' = 'Neuroepithelial-like', 
                     '14' = 'Unassigned', 
                     '15' = 'Mesenchymal',
                     '16' = 'Tumor-c24')

data[["group"]] <- Idents(data)

# DotPlot markers
DotPlot(data, 
        features = c(gene_list, 'ZFTA-RELA-Fusion1'), 
        assay = 'SCT', scale = TRUE) + 
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 8))
ggsave(file.path(plot_dir, "1_DotPlot_0010540-Region_3.pdf"), width=13, height=4)

# Export metadata with cell names and new annotation
cell_id <- data@meta.data %>% select(group)
cell_id <-  cell_id %>%
  rownames_to_column(var = "cell_id")

write_csv(cell_id, file.path(data_dir, "0010540-Region_3_cell_id.csv"))

# Export table with total number of cells
table <- data@meta.data %>% 
  group_by(group) %>% summarize(n())
write_csv(table, file.path(data_dir, "2_cell_number_0010540-Region_3.csv"))



# Load Xenium output 0010540-Region_4 (processed by Carlos) -------------------------------------

# Read data
data <- qread(file.path(base_dir, 'raw_data/0010540-Region_4.qs'))

# DotPlot markers for different resolutions
nres <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
for (res in nres) {
  Idents(data) <- paste0("SCT_snn_res.", res)
  # DotPlot markers
  DotPlot(data, cols = c('lightgrey', 'red'), features = gene_list, assay = 'SCT', scale = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 6))
  ggsave(file.path(plot_dir, paste0("individual/0010540-Region_4/1_DotPlot_res", res ,".pdf")), width=14, height=5)
}


# Change name identities
Idents(data) <- data$SCT_snn_res.0.3
data <- RenameIdents(data, 
                     '0' = 'Ependymal-like', 
                     '1' = 'Ependymal-like', 
                     '2' = 'NPC-like', 
                     '3' = 'Ependymal-like', 
                     '4' = 'Immune', 
                     '5' = 'Endothelial',
                     '6' = 'Mesenchymal', 
                     '7' = 'NPC-like', 
                     '8' = 'NPC-like', 
                     '9' = 'Immune_c15', 
                     '10' = 'Immune')

data[["group"]] <- Idents(data)

# DotPlot markers
DotPlot(data, 
        features = c(gene_list, 'ZFTA-RELA-Fusion1'), 
        assay = 'SCT', scale = TRUE) + 
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 8))
ggsave(file.path(plot_dir, "1_DotPlot_0010540-Region_4.pdf"), width=13, height=4)

# Export metadata with cell names and new annotation
cell_id <- data@meta.data %>% select(group)
cell_id <-  cell_id %>%
  rownames_to_column(var = "cell_id")

write_csv(cell_id, file.path(data_dir, "0010540-Region_4_cell_id.csv"))

# Export table with total number of cells
table <- data@meta.data %>% 
  group_by(group) %>% summarize(n())
write_csv(table, file.path(data_dir, "2_cell_number_0010540-Region_4.csv"))


# Load Xenium output 0010553-Region_1 (processed by Carlos) -------------------------------------

# Read data
data <- qread(file.path(base_dir, 'raw_data/0010553-Region_1.qs'))

# DotPlot markers for different resolutions
nres <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
for (res in nres) {
  Idents(data) <- paste0("SCT_snn_res.", res)
  # DotPlot markers
  DotPlot(data, cols = c('lightgrey', 'red'), features = gene_list, assay = 'SCT', scale = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 6))
  ggsave(file.path(plot_dir, paste0("individual/0010553-Region_1/1_DotPlot_res", res ,".pdf")), width=14, height=5)
}


# Change name identities
Idents(data) <- data$SCT_snn_res.0.7
data <- RenameIdents(data, 
                     '0' = 'Immune', 
                     '1' = 'Ependymal-like', 
                     '2' = 'Ependymal-like', 
                     '3' = 'NPC-like', 
                     '4' = 'Ependymal-like', 
                     '5' = 'Neuroepithelial-like',
                     '6' = 'Immune', 
                     '7' = 'NPC-like', 
                     '8' = 'Unassigned', 
                     '9' = 'Endothelial', 
                     '10' = 'Immune',
                     '11' = 'Immune',
                     '12' = 'Immune',
                     '13' = 'NPC-like')

data[["group"]] <- Idents(data)

# DotPlot markers
DotPlot(data, 
        features = c(gene_list, 'ZFTA-RELA-Fusion1'), 
        assay = 'SCT', scale = TRUE) + 
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 8))
ggsave(file.path(plot_dir, "1_DotPlot_0010553-Region_1.pdf"), width=13, height=4)

# Export metadata with cell names and new annotation
cell_id <- data@meta.data %>% select(group)
cell_id <-  cell_id %>%
  rownames_to_column(var = "cell_id")

write_csv(cell_id, file.path(data_dir, "0010553-Region_1_cell_id.csv"))

# Export table with total number of cells
table <- data@meta.data %>% 
  group_by(group) %>% summarize(n())
write_csv(table, file.path(data_dir, "2_cell_number_0010553-Region_1.csv"))



# Load Xenium output 0010553-Region_2 (processed by Carlos) -------------------------------------

# Read data
data <- qread(file.path(base_dir, 'raw_data/0010553-Region_2.qs'))

# DotPlot markers for different resolutions
nres <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
for (res in nres) {
  Idents(data) <- paste0("SCT_snn_res.", res)
  # DotPlot markers
  DotPlot(data, cols = c('lightgrey', 'red'), features = gene_list, assay = 'SCT', scale = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 6))
  ggsave(file.path(plot_dir, paste0("individual/0010553-Region_2/1_DotPlot_res", res ,".pdf")), width=14, height=5)
}

# Change name identities
Idents(data) <- data$SCT_snn_res.0.7
data <- RenameIdents(data, 
                     '0' = 'Ependymal-like', 
                     '1' = 'Unassigned', 
                     '2' = 'Ependymal-like', 
                     '3' = 'Ependymal-like', 
                     '4' = 'Neuroepithelial-like', 
                     '5' = 'NPC-like',
                     '6' = 'NPC-like', 
                     '7' = 'Immune', 
                     '8' = 'Unassigned', 
                     '9' = 'Immune', 
                     '10' = 'Endothelial',
                     '11' = 'Unassigned',
                     '12' = 'Ependymal-like',
                     '13' = 'Immune',
                     '14' = 'NPC-like')

data[["group"]] <- Idents(data)

# DotPlot markers
DotPlot(data, 
        features = c(gene_list, 'ZFTA-RELA-Fusion1'), 
        assay = 'SCT', scale = TRUE) + 
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 8))
ggsave(file.path(plot_dir, "1_DotPlot_0010553-Region_2.pdf"), width=13, height=4)

# Export metadata with cell names and new annotation
cell_id <- data@meta.data %>% select(group)
cell_id <-  cell_id %>%
  rownames_to_column(var = "cell_id")

write_csv(cell_id, file.path(data_dir, "0010553-Region_2_cell_id.csv"))

# Export table with total number of cells
table <- data@meta.data %>% 
  group_by(group) %>% summarize(n())
write_csv(table, file.path(data_dir, "2_cell_number_0010553-Region_2.csv"))





# Load Xenium output 0010553-Region_3 (processed by Carlos) -------------------------------------

# Read data
data <- qread(file.path(base_dir, 'raw_data/0010553-Region_3.qs'))

# DotPlot markers for different resolutions
nres <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
for (res in nres) {
  Idents(data) <- paste0("SCT_snn_res.", res)
  # DotPlot markers
  DotPlot(data, cols = c('lightgrey', 'red'), features = gene_list, assay = 'SCT', scale = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 6))
  ggsave(file.path(plot_dir, paste0("individual/0010553-Region_3/1_DotPlot_res", res ,".pdf")), width=14, height=5)
}

# Change name identities
Idents(data) <- data$SCT_snn_res.0.7
data <- RenameIdents(data, 
                     '0' = 'Unassigned', 
                     '1' = 'NPC-like', 
                     '2' = 'Unassigned', 
                     '3' = 'Ependymal-like', 
                     '4' = 'Neurons', 
                     '5' = 'Mesenchymal',
                     '6' = 'Unassigned', 
                     '7' = 'Ependymal-like', 
                     '8' = 'Unassigned', 
                     '9' = 'Neuroepithelial-like', 
                     '10' = 'Neuroepithelial-like',
                     '11' = 'Immune',
                     '12' = 'NPC-like',
                     '13' = 'Neuroepithelial-like')

data[["group"]] <- Idents(data)

# DotPlot markers
DotPlot(data, 
        features = c(gene_list, 'ZFTA-RELA-Fusion1'), 
        assay = 'SCT', scale = TRUE) + 
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 8))
ggsave(file.path(plot_dir, "1_DotPlot_0010553-Region_3.pdf"), width=13, height=4)

# Export metadata with cell names and new annotation
cell_id <- data@meta.data %>% select(group)
cell_id <-  cell_id %>%
  rownames_to_column(var = "cell_id")

write_csv(cell_id, file.path(data_dir, "0010553-Region_3_cell_id.csv"))


# Export table with total number of cells
table <- data@meta.data %>% 
  group_by(group) %>% summarize(n())
write_csv(table, file.path(data_dir, "2_cell_number_0010553-Region_3.csv"))




# Load Xenium output 0010553-Region_4 (processed by Carlos) -------------------------------------

# Read data
data <- qread(file.path(base_dir, 'raw_data/0010553-Region_4.qs'))

# DotPlot markers for different resolutions
nres <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
for (res in nres) {
  Idents(data) <- paste0("SCT_snn_res.", res)
  # DotPlot markers
  DotPlot(data, cols = c('lightgrey', 'red'), features = gene_list, assay = 'SCT', scale = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 6))
  ggsave(file.path(plot_dir, paste0("individual/0010553-Region_4/1_DotPlot_res", res ,".pdf")), width=14, height=5)
}

# Change name identities
Idents(data) <- data$SCT_snn_res.0.3
data <- RenameIdents(data, 
                     '0' = 'NPC-like', 
                     '1' = 'Ependymal-like', 
                     '2' = 'Ependymal-like', 
                     '3' = 'Endothelial', 
                     '4' = 'Immune', 
                     '5' = 'Neuroepithelial-like',
                     '6' = 'NPC-like', 
                     '7' = 'Immune', 
                     '8' = 'Immune', 
                     '9' = 'NPC-like')

data[["group"]] <- Idents(data)

# DotPlot markers
DotPlot(data, 
        features = c(gene_list, 'ZFTA-RELA-Fusion1'), 
        assay = 'SCT', scale = TRUE) + 
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 8))
ggsave(file.path(plot_dir, "1_DotPlot_0010553-Region_4.pdf"), width=13, height=4)

# Export metadata with cell names and new annotation
cell_id <- data@meta.data %>% select(group)
cell_id <-  cell_id %>%
  rownames_to_column(var = "cell_id")

write_csv(cell_id, file.path(data_dir, "0010553-Region_4_cell_id.csv"))



# Export table with total number of cells
table <- data@meta.data %>% 
  group_by(group) %>% summarize(n())
write_csv(table, file.path(data_dir, "2_cell_number_0010553-Region_4.csv"))




# Combine tables with number of cells -------------------------------------
tables <- list()
for (i in seq_along(samples)) {
  tables[[i]] <- read_csv(file.path(data_dir, paste0("2_cell_number_", names(samples)[i], ".csv")))
}
names(tables) <- names(samples)

# bind dataframes into single
table_combined <- bind_col(tables)



