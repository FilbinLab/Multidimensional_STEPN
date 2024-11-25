# Load packages -----------------------------------
library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(data.table)
library(readxl)
library(dplyr)
library(ggrastr)
library(UCell)
library(cowplot)
library(SingleCellExperiment)
library(SingleR)
library(celldex)

options(future.globals.maxSize = 8000 * 1024^2)

# Organize environment  -----------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

# folder path ----------------------------------------------

# folder with reference single cell object for label transfer (has to contain a .qs object with annotation)
seurat_object_path <- list.files(file.path(base_dir, 'analysis/1_preparation/data'), pattern = "\\.qs$", full.names = TRUE)

# folder with Xenium objects from 1_preprocessing
Xenium_seurat_obj_dir <- file.path(base_dir, 'analysis/2_preprocessing')

# folders to store results
analysis_dir  <- file.path(base_dir, 'analysis/3_program_annotation')
plot_dir  <- file.path(analysis_dir, 'plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}
data_dir  <- file.path(analysis_dir, 'data')
if (!dir.exists(data_dir)){dir.create(data_dir, recursive = T)}

resource_dir  <- file.path(base_dir, 'scripts_revisions/resources')
source(file.path(resource_dir, 'color_palette.R'))
source(file.path(resource_dir, "xenium_preprocessing_helper_functions_SD_v2.R"))



# read metadata ---------------------------------------------------------
metadata <- read_xlsx(file.path(base_dir, 'scripts_revisions/SampleIdentifier.xlsx'))
metadata <- metadata %>% filter(SampleName == SampleName_arg)

Xenium_dir <- file.path(Xenium_seurat_obj_dir, metadata$SampleName)


# Read marker genes of Xenium panel genes -------------------------------------
marker_genes <- read_xlsx(file.path(base_dir, 'Xenium_panel.xlsx'))
marker_genes <- marker_genes[, c('Gene', 'Annotation_Sara')]
marker_genes_list <- split(marker_genes$Gene, marker_genes$Annotation_Sara)


# select cell types with priority ---------------

# select cell types with priority from Xenium panel 
labels_priority_panel <- c('Endothelial', 'Astrocytes', 'Neurons', 'Oligodendrocytes', 'T-cells', 'Myeloid')

# select cell types with priority from single cell object (all tumor cells 
# except for embryonic-like (they all had low ZFTA-RELA expression and neuron markers, 
# most likely normal cells)
labels_priority_scObject <-c( "Neuronal-like" , "Radial-glia-like", "MES-like", "Ependymal-like",
                             'Oligodendrocytes', 'T-cells', 'Myeloid', 'Endothelial')

malignant_cells <- c("Neuroepithelial-like", "Radial-glia-like", 
                     "Embryonic-neuronal-like", "Neuronal-like" ,"Ependymal-like", "MES-like", "Embryonic-like")



colors <- colors_metaprograms_Xenium


# Assign labels
XeniumCellAssignment(seurat_object_path, Xenium_dir, marker_genes_list, 
                                 labels_priority_panel, labels_priority_scObject,
                                 malignant_cells,
                                 colors)
  

