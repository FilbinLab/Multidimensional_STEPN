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
library(scCustomize)
library(ComplexHeatmap)


# Organize environment  -----------------------------------
base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/9 - scRNAseq coculture"
analysis_dir <- file.path(base_dir, 'analysis')
resource_dir <- file.path(base_dir, 'scripts/resources')

plot_dir <- file.path(base_dir, 'analysis/plots/synaptic')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

data_dir <- file.path(base_dir, 'analysis/data')


# load Seurat object -----------------------------------
seurat_obj_merged <- readRDS(file.path(data_dir, "Mono_co_culture_merged.rds"))


# load synaptic gene lists -----------------------------------
synaptic_genes <- read_excel(file.path(resource_dir, 'Neurotransmitter_Genes.xlsx')) 
synaptic_genes <- as.data.frame(synaptic_genes)

# transform df into a list
synaptic_list <- split(synaptic_genes$Approved.symbol, synaptic_genes$Group)


# score Seurat object -----------------------------------
seurat_obj_merged <- AddModuleScore(object = seurat_obj_merged, assay = 'RNA', features = synaptic_list, name = names(synaptic_list))

# rename metadata names of scores
  # identify number of first column with metadata scores
  col_start <- length(colnames(seurat_obj_merged@meta.data)) - length(names(synaptic_list)) +1
  # identify number of last column with metadata scores
  col_end <- length(colnames(seurat_obj_merged@meta.data))
  # rename columns with score name
  colnames(seurat_obj_merged@meta.data)[col_start:col_end] <- names(synaptic_list)


# Dotplot expression synaptic scores -----------------------------------
# Dotplot of scores
Idents(seurat_obj_merged) = 'model'
DotPlot(seurat_obj_merged, features = names(synaptic_list), assay = 'RNA', scale = FALSE)  + 
  #scale_colour_distiller(palette="RdBu") +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        axis.title=element_blank(),
        legend.text = element_text(size = 12),
        #legend.title = element_blank()
  ) 
ggsave(file.path(plot_dir, "1_Dotplot.pdf"), width=8, height=3)


# Dotplot of genes
DotPlot(seurat_obj_merged, features = synaptic_genes$Approved.symbol, assay = 'RNA', scale = FALSE)  + 
  #scale_colour_distiller(palette="RdBu") +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        axis.title=element_blank(),
        legend.text = element_text(size = 12),
        #legend.title = element_blank()
  ) 
ggsave(file.path(plot_dir, "1_Dotplot_genes.pdf"), width=35, height=3)

