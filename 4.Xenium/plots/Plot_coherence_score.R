# Load packages -----------------------------------
library(Seurat)
library(tidyverse)
library(glue)
library(gplots)
library(qs)
library(reshape)


# Organize environment  -----------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

coherence_dir <- file.path(base_dir, 'analysis/6_coherence/data')
cellid_dir <- file.path(base_dir, 'analysis/3_program_annotation/data')

plot_dir  <- file.path(base_dir, 'analysis/6_coherence/plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}


resource_dir  <- file.path(base_dir, 'scripts_revisions/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))
source(file.path(resource_dir, "xenium_preprocessing_helper_functions_SD_v2.R"))
source(file.path(resource_dir, 'color_palette.R'))




## Read metadata ------------------------------------------------------
metadata <- read_xlsx(file.path(base_dir, 'scripts_revisions/SampleIdentifier.xlsx'))
metadata <- metadata %>% filter(!SampleName %in% c('STEPN14_Region_1', 'STEPN14_Region_2'))

# reorder by sampleID
metadata <- metadata %>% arrange(SampleName)

FileName <- metadata$SampleName

order_metaprograms <- c("Cycling",  "Embryonic-like", "Neuroepithelial-like", "Radial-glia-like", 
                        "Embryonic-neuronal-like", 
                        "Neuronal-like" ,"Ependymal-like", "MES-like", 
                        "T-cells", "Myeloid",  "Endothelial",  "Oligodendrocytes", 'Astrocyte', 'Neurons', 'Unknown')


# Box plots by MP and sample ----------------------
plot_coherence_score(FileName, cellid_dir, coherence_dir, plot_dir, colors_metaprograms_Xenium, col_sampling) 



# Linear regression  ----------------------
CalculateLinearRegression(FileName, cellid_dir, coherence_dir, plot_dir)
