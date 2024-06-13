# Load packages -----------------------------------
rm(list = ls())

library(data.table)
library(tidyverse)
library(crayon)
library(ape)
library(readxl)
library(qs)
library(cowplot)
library(Seurat)

# Organize environment  -----------------------------------
base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/10 - Developmental mapping"
analysis_dir <- file.path(base_dir, 'analysis/data')

plot_dir <- file.path(base_dir, "analysis/plot")
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

data_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Developmental datasets/data"

resource_dir <- file.path(base_dir, 'scripts/resources')
#source(file.path(resource_dir, 'single_cell_preprocessing_helper_functions_SD_240131.R'))
source(file.path(resource_dir, 'single_cell_preprocessing_helper_functions_CBJr.R'))
source(file.path(resource_dir, 'NMF_helper_function.R'))


# Check if the file with tumor aggregated programs already exists
if (file.exists(file.path(analysis_dir, "agg_cm_mean_tumor.qs"))) {
  # If the file exists, print a message or perform any other action
  print("The file already exists. Skipping this part of the code.")
} else {
  
  # Process tumor dataset by sample -----------------------------------
  
  # Read tumor file (processed and annotated)
  data <- qread(file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/11 - Plots scRNAseq/data/patient/seurat_obj_ST_combined_malig_annotated.qs'))
  
  # Pseudobulk
  pseudobulk_mean <- AggregateExpression(data, return.seurat = F, slot = "counts", assays = "RNA", group.by = "Sample_deID")$RNA %>%
    as.data.frame()
  colnames(pseudobulk_mean) <- sort(unique(data$Sample_deID))
  qsave(pseudobulk_mean, file.path(analysis_dir, "agg_cm_mean_tumor.qs"))
  
  # Identify markers
  Idents(data) <- data$Sample_deID
  markers <- RunPrestoAll(data, 
                          only.pos = T, 
                          densify = T, 
                          max.cells.per.ident = 500,
                          #logfc.threshold = 0.25,
                          #min.pct = 0.15
                          ) %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(-avg_log2FC, .by_group = T)
  qsave(markers, file.path(analysis_dir, "markers_tumor.qs"))
  
  # HVGs
  qsave(VariableFeatures(data), file.path(analysis_dir, "hvgs.qs"))
  
  # Process tumor dataset by subtype -----------------------------------
  
  # Pseudobulk
  pseudobulk_mean <- AggregateExpression(data, return.seurat = F, slot = "counts", assays = "RNA", group.by = "Subtype")$RNA %>%
    as.data.frame()
  colnames(pseudobulk_mean) <- sort(unique(data$Subtype))
  qsave(pseudobulk_mean, file.path(analysis_dir, "agg_cm_mean_tumor_subtype.qs"))
  
  # Identify markers
  Idents(data) <- data$Subtype
  markers <- RunPrestoAll(data, 
                          only.pos = T, 
                          densify = T, 
                          max.cells.per.ident = 500,
                          #logfc.threshold = 0.25,
                          #min.pct = 0.15
                          )  %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(-avg_log2FC, .by_group = T)
  qsave(markers, file.path(analysis_dir, "markers_tumor_subtype.qs"))
  
  
  # Process tumor dataset by Metaprogram -----------------------------------
  
  # Pseudobulk
  pseudobulk_mean <- AggregateExpression(data, return.seurat = F, slot = "counts", assays = "RNA", group.by = "Metaprogram")$RNA %>%
    as.data.frame()
  colnames(pseudobulk_mean) <- sort(unique(data$Metaprogram))
  qsave(pseudobulk_mean, file.path(analysis_dir, "agg_cm_mean_tumor_metaprogram.qs"))
  
  # Identify markers
  Idents(data) <- data$Metaprogram
  markers <- RunPrestoAll(data, 
                          only.pos = T, 
                          densify = T, 
                          max.cells.per.ident = 500,
                          #logfc.threshold = 0.25,
                         #min.pct = 0.15
                         )  %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(-avg_log2FC, .by_group = T)
  qsave(markers, file.path(analysis_dir, "markers_tumor_metaprogram.qs"))
  
}



# Project tumor metaprograms onto reference Nowakowski 2017+Eze et al. 2021-----------------------------------

## Loading reference data 
ref_cm <- qread(file.path(data_dir, "Nowakowski_Eze/agg_cm_mean.qs"))
# change names with "_", otherwise I get an error
colnames(ref_cm) <- str_replace_all(colnames(ref_cm), "_", " ")

# remove some elements and reorder 
desired_order <-  c("Neuroepithelial", "Mesenchymal", "Early radial glia",  "Radial glia", 
                    "Newborn/maturing neuron",  "MGE inhibitory neuron lineage", "Neuronal", "IPCs" ,
                      "OPC",  "Astrocyte" , "Choroid")
ref_cm <- ref_cm[, desired_order]

ref_degs <- qread(file.path(data_dir, "Nowakowski_Eze/markers.qs"))
ref_degs <- ref_degs %>% group_by(cluster) %>% top_n(100, wt = avg_log2FC)
ref_degs <- lapply(split(ref_degs, f = ref_degs$cluster), function(x) x$gene)

# change names with "_", otherwise I get an error
names(ref_degs) <- str_replace_all(names(ref_degs), "_", " ")

# remove some elements and reorder 
ref_degs <- ref_degs[desired_order]

ref_hvgs <- qread(file.path(data_dir, "Nowakowski_Eze/hvgs.qs"))


## Loading query data
query_cm <- qread(file.path(analysis_dir, "agg_cm_mean_tumor.qs")) %>% as.data.frame()
query_degs <- qread(file.path(analysis_dir, 'markers_tumor.qs'))
query_degs <- split(query_degs, query_degs$cluster)
query_degs <- lapply(query_degs, function(x) unique(x$gene))
query_hvgs <- qread(file.path(analysis_dir, "hvgs.qs"))
query_metagene_order <-  rev(colnames(query_cm))

## project data
p1 <- projectDataNoScaling(query_cm, query_degs, query_hvgs, query_metagene_order,
                     ref_cm, ref_degs, ref_hvgs,
                     outFile = file.path(plot_dir, 'Nowakowski_2017_Eze_2021.pdf'),
                     desired_order, fig_width = 8, fig_height = 12)


## Loading query data (subtype)
query_cm <- qread(file.path(analysis_dir, "agg_cm_mean_tumor_subtype.qs")) %>% as.data.frame()
query_degs <- qread(file.path(analysis_dir, 'markers_tumor_subtype.qs'))
query_degs <- split(query_degs, query_degs$cluster)
query_degs <- lapply(query_degs, function(x) unique(x$gene))
query_metagene_order <- c("ST-YAP1", "ZFTA-Cluster 4", "ZFTA-Cluster 3", "ZFTA-Cluster 2", "ZFTA-Cluster 1", "ZFTA-RELA"    )

## project data
p2 <- projectDataNoScaling(query_cm, query_degs, query_hvgs, query_metagene_order,
                     ref_cm, ref_degs, ref_hvgs,
                     outFile = file.path(plot_dir, 'Nowakowski_2017_Eze_2021_subtype.pdf'),
                     desired_order, fig_width = 8, fig_height = 5)


## Loading query data (metaprograms)
query_cm <- qread(file.path(analysis_dir, "agg_cm_mean_tumor_metaprogram.qs")) %>% as.data.frame()
query_degs <- qread(file.path(analysis_dir, 'markers_tumor_metaprogram.qs'))
query_degs <- split(query_degs, query_degs$cluster)
query_degs <- lapply(query_degs, function(x) unique(x$gene))
query_metagene_order <- rev(c("Cycling",  "MES-like","Embryonic-like", "Neuroepithelial-like", "Radial-glia-like", 
                              "Embryonic-neuronal-like", "Neuronal-like" ,"Ependymal-like"))

## project data
p3 <- projectDataNoScaling(query_cm, query_degs, query_hvgs, query_metagene_order,
                           ref_cm, ref_degs, ref_hvgs,
                           outFile = file.path(plot_dir, 'Nowakowski_2017_Eze_2021_metaprogram.pdf'),
                           desired_order, fig_width = 8, fig_height = 6)



# Project tumor metaprograms onto reference Nowakowski 2017+Eze et al. 2021 (TIME) -----------------------------------

## Loading reference data 
ref_cm <- qread(file.path(data_dir, "Nowakowski_Eze/agg_cm_mean_dev_time.qs"))
# change names with "_", otherwise I get an error
colnames(ref_cm) <- str_replace_all(colnames(ref_cm), "_", " ")

# remove some elements and reorder 
desired_order <-  c("First - Early",  "First - Middle", "First - Late", "Second" , "Third")
ref_cm <- ref_cm[, desired_order]

ref_degs <- qread(file.path(data_dir, "Nowakowski_Eze/markers_dev_time.qs"))
ref_degs <- ref_degs %>% group_by(cluster) %>% top_n(100, wt = avg_log2FC)
ref_degs <- lapply(split(ref_degs, f = ref_degs$cluster), function(x) x$gene)

# change names with "_", otherwise I get an error
names(ref_degs) <- str_replace_all(names(ref_degs), "_", " ")

# remove some elements and reorder 
ref_degs <- ref_degs[desired_order]

ref_hvgs <- qread(file.path(data_dir, "Nowakowski_Eze/hvgs.qs"))


## Loading query data
query_cm <- qread(file.path(analysis_dir, "agg_cm_mean_tumor.qs")) %>% as.data.frame()
query_degs <- qread(file.path(analysis_dir, 'markers_tumor.qs'))
query_degs <- split(query_degs, query_degs$cluster)
query_degs <- lapply(query_degs, function(x) unique(x$gene))
query_hvgs <- qread(file.path(analysis_dir, "hvgs.qs"))
query_metagene_order <- rev(colnames(query_cm))
## project data
p4 <- projectDataNoScaling(query_cm, query_degs, query_hvgs, query_metagene_order,
                     ref_cm, ref_degs, ref_hvgs,
                     outFile = file.path(plot_dir, 'Nowakowski_2017_Eze_2021_times.pdf'),
                     desired_order, fig_width = 6, fig_height = 12)


## Loading query data (subtype)
query_cm <- qread(file.path(analysis_dir, "agg_cm_mean_tumor_subtype.qs")) %>% as.data.frame()
query_degs <- qread(file.path(analysis_dir, 'markers_tumor_subtype.qs'))
query_degs <- split(query_degs, query_degs$cluster)
query_degs <- lapply(query_degs, function(x) unique(x$gene))
query_metagene_order <- c("ST-YAP1", "ZFTA-Cluster 4", "ZFTA-Cluster 3", "ZFTA-Cluster 2", "ZFTA-Cluster 1", "ZFTA-RELA"    )

## project data
p5 <- projectDataNoScaling(query_cm, query_degs, query_hvgs, query_metagene_order,
                     ref_cm, ref_degs, ref_hvgs,
                     outFile = file.path(plot_dir, 'Nowakowski_2017_Eze_2021_subtype_times.pdf'),
                     desired_order, fig_width = 6, fig_height = 4)




## Loading query data (metaprograms)
query_cm <- qread(file.path(analysis_dir, "agg_cm_mean_tumor_metaprogram.qs")) %>% as.data.frame()
query_degs <- qread(file.path(analysis_dir, 'markers_tumor_metaprogram.qs'))
query_degs <- split(query_degs, query_degs$cluster)
query_degs <- lapply(query_degs, function(x) unique(x$gene))
query_metagene_order <- rev(c("Cycling",  "MES-like","Embryonic-like", "Neuroepithelial-like", "Radial-glia-like", 
                              "Embryonic-neuronal-like", "Neuronal-like" ,"Ependymal-like"))

## project data
p6 <- projectDataNoScaling(query_cm, query_degs, query_hvgs, query_metagene_order,
                           ref_cm, ref_degs, ref_hvgs,
                           outFile = file.path(plot_dir, 'Nowakowski_2017_Eze_2021_metaprograms_times.pdf'),
                           desired_order, fig_width = 6, fig_height = 4.5)


plot_grid(p1, p4, nrow = 1, rel_widths = c(1, 0.8),
          align = "h", axis = "b")
ggsave(file.path(plot_dir, 'Nowakowski_2017_Eze_2021_SUMMARY.pdf'), height = 12, width = 12)

plot_grid(p2, p5, nrow = 1, rel_widths = c(1, 0.8),
          align = "h", axis = "b")
ggsave(file.path(plot_dir, 'Nowakowski_2017_Eze_2021_SUMMARY.pdf'), height = 4.5, width = 12.5)

plot_grid(p3, p6, nrow = 1, rel_widths = c(1, 0.8),
          align = "h", axis = "b")
ggsave(file.path(plot_dir, 'Nowakowski_2017_Eze_2021_metaprogram_SUMMARY.pdf'), height = 5.5, width = 15)

