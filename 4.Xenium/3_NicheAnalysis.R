args <- commandArgs(trailingOnly = TRUE)
print(args)

SampleName_arg <- args[1]

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

# folder with Xenium objects from 2_assign_programs
Xenium_seurat_obj_dir <- file.path(base_dir, 'analysis/3_program_annotation/data')

# folders to store results
analysis_dir  <- file.path(base_dir, 'analysis/4_niche')
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

# Perform niche analysis --------------------------------------------
NicheAnalysis(Xenium_dir, data_dir, niches_min = 4, niches_max = 6, n_neighbors = 20)
  


