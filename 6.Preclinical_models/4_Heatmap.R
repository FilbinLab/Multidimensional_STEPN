# Load packages -----------------------------------
rm(list = ls())

library(data.table)
library(tidyverse)
library(R.utils)
library(ggpubr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(writexl)
library(tidyverse)
library(paletteer)
library(readxl)
library(cowplot)
#library(scCustomize)
library(ComplexHeatmap)
library(qs)


# Organize environment  -----------------------------------
base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/6 - scRNAseq models"
analysis_dir <- file.path(base_dir, 'analysis')
resource_dir <- file.path(base_dir, 'scripts/resources')

plot_dir <- file.path(analysis_dir, 'plots')
data_dir <- file.path(analysis_dir, 'data')

source(file.path(resource_dir, 'Plot_style.R'))
source(file.path(resource_dir, "single_cell_preprocessing_helper_functions.R"))


# load Seurat object -----------------------------------
seurat_obj_merged <- readRDS(file.path(analysis_dir, "data/Patient_models_merged.rds"))

# export metadata as df
metadata <- data.frame(seurat_obj_merged@meta.data)

# load gene lists
NMFmarkers <- readRDS(file.path(base_dir, 'data/patient/nmf_marker_genes_final_annotated_200.rds'))
names(NMFmarkers) <- c("Neuroepithelial-like", "MES-like","Ependymal-like","Radial-glia-like","Cycling" ,            
                         "Embryonic-neuronal-like", "Neuronal-like", "Embryonic-like")

NMFmarkers <- lapply(NMFmarkers,head,50)


# create df with genes of interest
genes <- c(NMFmarkers$`Ependymal-like`, NMFmarkers$`Neuroepithelial-like`, NMFmarkers$`Neuronal-like`)
genes_df <- data.frame(
  Values = unlist(NMFmarkers),
    Group = rep(names(NMFmarkers), sapply(NMFmarkers, length))
    )
rownames(genes_df) <- NULL

genes_df_subset <- genes_df %>%
  filter(Group == 'Ependymal-like' |  Group == 'Neuronal-like' | Group == 'Neuroepithelial-like' )


# subset Seurat object to genes of interest -----------------------------------
seurat_df <- seurat_obj_merged@assays$RNA$data[genes_df_subset$Values, ]


# define annotations
order <- data.frame(metadata %>% 
                      dplyr::arrange(metadata$Ependymal.like - metadata$Neuronal.like))
ordered_rows <- rownames(order)
ordered_scores <- order$Ependymal.like - order$Neuronal.like


# scale dataset
seurat_df <- seurat_df - rowMeans(seurat_df)
seurat_df <- as.data.frame(seurat_df)

# order based on neuronal score
seurat_df <- seurat_df[, ordered_rows]

# define annotations (important: order based on new order of cells)
metadata2 <- metadata[colnames(seurat_df), ]
metadata <- metadata2
annotation_info <- metadata[, c("model", "Cycling_prop", "signature_1")]

# define colors
col_model = c('grey80', paletteer::paletteer_d("beyonce::X18")[3:5])
names(col_model) = unique(annotation_info$model)



col_cycling_prop = c('grey80', 'grey10')
names(col_cycling_prop) = unique(annotation_info$Cycling_prop)

colors_metaprograms<- c("gray30","#F99E93FF","#9E5E9BFF","#74ADD1FF",'#0F4F8B', "#ACD39EFF","#96410EFF", 'mistyrose1')

names(colors_metaprograms) <- c("Cycling", "Neuroepithelial-like", "Radial-glia-like", 
                                "Embryonic-neuronal-like", "Neuronal-like" ,"Ependymal-like", "MES-like", "Embryonic-like")


col_column_ha <- list(model = col_model,
                      Cycling_prop = col_cycling_prop,
                      signature_1 = colors_metaprograms
                      )


# define column annotation
column_ha = HeatmapAnnotation(df = annotation_info,
                              col = col_column_ha)

bottom_ha = HeatmapAnnotation(scores = anno_barplot(ordered_scores))

# define row annotation
row_ha = rowAnnotation(Group = genes_df_subset$Group)

# genes to mark
elements_to_find <- c('VIM', 'TAGLN2', 'DCX', 'MEG3', 'NRXN3', 'DNAH9', 'DNAH7')
rows_genes_to_mark <- which(genes_df_subset$Values %in% elements_to_find)
rows_genes_to_mark_name <- genes_df_subset$Values[rows_genes_to_mark]

right_ha = rowAnnotation(foo = anno_mark(at = rows_genes_to_mark, 
                                         labels = rows_genes_to_mark_name))

# generate heatmaplibrary(circlize)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c('#2E5A87FF' , '#FCFDFEFF', '#A90C38FF'))
col_fun(seq(-3, 3))


# plot heatmap 
Cairo::CairoPDF(file.path(plot_dir, "1_heatmap_NPC_Ependymal.pdf"), width=12, height=6)
Heatmap(seurat_df,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              bottom_annotation = bottom_ha,
              #col = col_fun,
              #col = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", n=10)),
              column_split =c(annotation_info$model), border = TRUE,
              row_split = factor(c(genes_df_subset$Group), levels = c("Neuroepithelial-like", "Neuronal-like", "Ependymal-like")),
              top_annotation = column_ha,
              right_annotation = right_ha,
              #left_annotation = row_ha,
              use_raster = TRUE,
              raster_quality = 30
              )
dev.off()
 


