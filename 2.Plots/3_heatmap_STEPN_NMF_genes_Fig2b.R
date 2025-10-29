# Load packages -----------------------------------
rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(readxl)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(qs)
library(writexl)

# Organize environment  -----------------------------------
base_dir <- "/Users/sdaniell/Partners HealthCare Dropbox/Sara Danielli/Filbin lab - sample location/Filbin Lab Active/Sara/Project/Ependymoma/11 - Plots scRNAseq"

plot_dir <- file.path(base_dir, 'analysis/plots/heatmap')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}


# load Seurat object -----------------------------------
seurat_obj <- qread(file = file.path(base_dir, "data/patient/seurat_obj_ST_malig_annotated.qs"))

seurat_obj$Metaprogram <- factor(x = seurat_obj$Metaprogram, levels = c("Cycling", "MES/Hypoxia", "Embryonic-like", "Neuroepithelial-like", "Radial-glia-like", 
                                                                        "Embryonic-neuronal-like", "Neuronal-like"  ,"Ependymal-like"))



# -----------------------------------------------------------------------
# Heatmap top NMF marker genes 
# -----------------------------------------------------------------------

# load genes
markers_NMF <- readRDS(file.path(base_dir, 'data/signatures/nmf_marker_genes_final_annotated.rds'))
names(markers_NMF) <- c("Neuroepithelial-like", "MES/Hypoxia","Ependymal-like","Radial-glia-like","Cycling", "Embryonic-neuronal-like", "Neuronal-like", "Embryonic-like")
# select top 50
markers_NMF_30 <- lapply(markers_NMF, head, 30)


# save metadata as df
metadata <- data.frame(seurat_obj@meta.data)

# aggregate cells from same cluster and identify top marker genes-----------------------------------

# aggregate cells from same cluster 
Idents(seurat_obj) <- 'Metaprogram'
seurat_obj_agg <- AggregateExpression(seurat_obj, return.seurat = TRUE)



desired_order <-  c("Cycling", "MES/Hypoxia", "Embryonic-like", "Neuroepithelial-like", "Radial-glia-like", 
                    "Embryonic-neuronal-like", "Neuronal-like"  ,"Ependymal-like")
markers_NMF_30 <- markers_NMF_30[desired_order]


# Subset seurat object to genes of interest -----------------------------------
# create df with genes of interest
genes <- markers_NMF_30
genes_df <- data.frame(
  Values = unlist(markers_NMF_30),
  Group = rep(names(markers_NMF_30), sapply(markers_NMF_30, length))
)
rownames(genes_df) <- NULL

# subset Seurat object to genes of interest 
seurat_df <- seurat_obj_agg@assays$RNA$data[genes_df$Values, ]


# define annotations -----------------------------------

# scale dataset
seurat_df <- seurat_df - rowMeans(seurat_df)
seurat_df <- as.data.frame(seurat_df)


# define row annotation
row_ha = rowAnnotation(Group = genes_df$Group)


elements_to_find <- c('CENPF','HELLS','MKI67', # cycling
                      'ACTG1', 'CALM1', 'HMGN2', #Neuroepithelial
                      'LRRTM4', 'DCLK1', 'LRP1B', 'LRIG1', #radial glia
                      'RBFOX1', 'GRID2', # Embryonic-neuronal-like
                      'MEG3','NRXN3', # neuronal-like
                      'LRIG1', # RG like
                      'CD44', 'PGK1', 'ENO1',  #MES
                      'DNAH6', 'DNAH9', # ependymal
                      'ALDH1A2','CTBP2','ACTN4' # embryonic
)


rows_genes_to_mark <- which(genes_df$Values %in% elements_to_find)
rows_genes_to_mark_name <- genes_df$Values[rows_genes_to_mark]

right_ha = rowAnnotation(foo = anno_mark(at = rows_genes_to_mark, 
                                         labels = rows_genes_to_mark_name))

# define colors heatmap
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#2E5A87FF", "#FCFDFEFF", "#A90C38FF"))


# plot heatmap 
Cairo::CairoPDF(file.path(plot_dir, "1_heatmap_top30_NMF.pdf"), width=5, height=7)
ht <- Heatmap(seurat_df,
              cluster_columns = FALSE,
              show_row_names = F,
              show_column_names = T,
              border = TRUE,
              right_annotation = right_ha,
              use_raster = TRUE,
              raster_quality = 10,
              cluster_rows = FALSE,
              col = col_fun
)
draw(ht)
dev.off()


