# Load packages -----------------------------------
rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(readxl)
library(data.table)
library(ComplexHeatmap)
#library(magick)
#library(RColorBrewer)
library(circlize)
library(qs)
library(writexl)

# Organize environment  -----------------------------------
#base_dir <- '/Volumes/Sara_PhD/scRNAseq_data'
base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/11 - Plots scRNAseq"

genelist_dir <- file.path(base_dir, "analysis/plots/STEPN/malignant")

plot_dir <- file.path(base_dir, 'analysis/plots/heatmap')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}



# -----------------------------------------------------------------------
# (4) Heatmap top marker genes 
# -----------------------------------------------------------------------
# load Seurat object -----------------------------------
seurat_obj_mal <- qread(file = file.path(base_dir, "data/patient/seurat_obj_ST_combined_malig_annotated.qs"))

seurat_obj_mal$Metaprogram <- factor(x = seurat_obj_mal$Metaprogram, levels = c("Cycling", "MES-like", "Embryonic-like", "Neuroepithelial-like", "Radial-glia-like", 
                                                                                "Embryonic-neuronal-like", "Neuronal-like"  ,"Ependymal-like"))


# save metadata as df
metadata <- data.frame(seurat_obj_mal@meta.data)

# aggregate cells from same cluster and identify top marker genes-----------------------------------
# aggregate cells from same cluster -----------------------------------
Idents(seurat_obj_mal) <- 'Metaprogram'

seurat_obj_mal <- AggregateExpression(seurat_obj_mal, return.seurat = TRUE)


# load gene lists -----------------------------------

markers <- read_xlsx(file.path(genelist_dir, '1_marker_genes.xlsx'))

marker_list <- markers %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T)
marker_list <- as.data.frame(marker_list)
# convert into list
marker_list <- split(marker_list[, -c(1:6)], f = marker_list$cluster)

# select top 200 genes per cluster
marker_list <- lapply(marker_list,head,200)

desired_order <-  c("Cycling", "MES-like", "Embryonic-like", "Neuroepithelial-like", "Radial-glia-like", 
                    "Embryonic-neuronal-like", "Neuronal-like"  ,"Ependymal-like")
marker_list <- marker_list[desired_order]



# Subset seurat object to genes of interest -----------------------------------
# create df with genes of interest
genes <- marker_list
genes_df <- data.frame(
  Values = unlist(marker_list),
  Group = rep(names(marker_list), sapply(marker_list, length))
)
rownames(genes_df) <- NULL

# subset Seurat object to genes of interest 
seurat_df <- seurat_obj_mal@assays$RNA$data[genes_df$Values, ]


# define annotations -----------------------------------

# scale dataset
seurat_df <- seurat_df - rowMeans(seurat_df)
seurat_df <- as.data.frame(seurat_df)


# define row annotation
row_ha = rowAnnotation(Group = genes_df$Group)


elements_to_find <- c('TOP2A','KIF18B','CDC20',
                      'HOXD3','GRHL2','GRAP',
                      'CDH6','FRZB','DLK1',
                      'PAX2','DLCK1','CSMD2',
                      'PADI4','NBPF6','PIWIL3',
                      'DCX','KIRREL3','KCNQ3',
                      'DNAH10','DNAH11','CCDC37',
                      'CD44','PDK1','VEGFA')


rows_genes_to_mark <- which(genes_df$Values %in% elements_to_find)
rows_genes_to_mark_name <- genes_df$Values[rows_genes_to_mark]

right_ha = rowAnnotation(foo = anno_mark(at = rows_genes_to_mark, 
                                         labels = rows_genes_to_mark_name))

# define colors heatmap
col_fun = colorRamp2(c(-1, 0, 1), c("#2E5A87FF", "#FCFDFEFF", "#A90C38FF"))


# plot heatmap 
Cairo::CairoPDF(file.path(plot_dir, "1_seurat_obj_mal_heatmap_aggregate_all_FindMarkers.pdf"), width=6, height=6)
ht <- Heatmap(seurat_df,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              #bottom_annotation = bottom_ha,
              border = TRUE,
              # row_split = factor(c(genes_df$Group), 
              #                    levels = c("Cycling", "Embryonic-like", "Neuroepithelial-like", "Radial-glia-like", 
              #                               "Embryonic-neuronal-like", "Neuronal-like"  ,"Ependymal-like", "MES-like")),
              #top_annotation = column_ha,
              right_annotation = right_ha,
              use_raster = TRUE,
              raster_quality = 10,
              cluster_rows = FALSE,
              col = col_fun,
              row_names_rot = 30
)
draw(ht)
dev.off()

