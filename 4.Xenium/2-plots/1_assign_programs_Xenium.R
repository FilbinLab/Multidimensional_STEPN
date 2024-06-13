# Load packages -----------------------------------
library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(data.table)
library(readxl)
library(dplyr)
library(ggrastr)


# Organize environment  -----------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

resource_dir  <- file.path(base_dir, 'scripts/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))

analysis_dir  <- file.path(base_dir, 'analysis/all_Xenium_runs')

plot_dir <- file.path(analysis_dir, 'plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

data_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(data_dir)){dir.create(data_dir, recursive = T)}

seurat_obj_dir <- file.path(data_dir, 'seurat_objects')
if (!dir.exists(seurat_obj_dir)){dir.create(seurat_obj_dir, recursive = T)}


metadata <- read_xlsx(file.path(base_dir, 'SampleIdentifier.xlsx'))


# Define resolution to use and annotations for each tissue  -------------------------------------

annotation_clusters <- list (
  '0010652-Region_4' = c('0' = 'Neuroepithelial-like',
                         '1' = 'Neuronal-like',
                         '2' = 'Myeloid',
                         '3' = 'Neuronal-like', 
                         '4' = 'Endothelial',
                         '5' = 'Ependymal-like',
                         '6' = 'Neuroepithelial-like',
                         '7' = 'Endothelial', 
                         '8' = 'Neuronal-like', 
                         '9' = 'Ependymal-like', 
                         '10' = 'T-cell'),
  
  '0010575-Region_1' = c('0' = 'Ependymal-like',
                         '1' = 'Ependymal-like',
                         '2' = 'Neuronal-like',
                         '3' = 'Neuroepithelial-like', 
                         '4' = 'Mesenchymal',
                         '5' = 'Endothelial',
                         '6' = 'Myeloid',
                         '7' = 'Neuronal-like', 
                         '8' = 'Neuronal-like', 
                         '9' = 'Myeloid', 
                         '10' = 'T-cell', 
                         '11' = 'Myeloid', 
                         '12' = 'Neuronal-like'),
  
  '0010575-Region_2' = c('0' = 'Ependymal-like',
                         '1' = 'Neuronal-like',
                         '2' = 'Ependymal-like',
                         '3' = 'Neuroepithelial-like', 
                         '4' = 'Mesenchymal',
                         '5' = 'Myeloid',
                         '6' = 'Endothelial',
                         '7' = 'Neuronal-like', 
                         '8' = 'Neuronal-like', 
                         '9' = 'Myeloid', 
                         '10' = 'T-cell', 
                         '11' = 'Myeloid', 
                         '12' = 'Neuronal-like'),
  
  '0010575-Region_3' = c('0' = 'Ependymal-like',
                         '1' = 'Neuronal-like',
                         '2' = 'Ependymal-like',
                         '3' = 'Mesenchymal', 
                         '4' = 'Myeloid',
                         '5' = 'Ependymal-like',
                         '6' = 'Neuroepithelial-like',
                         '7' = 'Endothelial', 
                         '8' = 'Neuronal-like', 
                         '9' = 'Neuronal-like', 
                         '10' = 'Myeloid', 
                         '11' = 'T-cell', 
                         '12' = 'Myeloid',
                         '13' = 'Neuronal-like'),
  
  '0010619-Region_1' = c('0' = 'Ependymal-like',
                         '1' = 'Neuronal-like',
                         '2' = 'Ependymal-like',
                         '3' = 'Neuroepithelial-like', 
                         '4' = 'Myeloid',
                         '5' = 'VLMCs',
                         '6' = 'VLMCs',
                         '7' = 'Myeloid', 
                         '8' = 'VLMCs', 
                         '9' = 'Ependymal-like', 
                         '10' = 'Neuronal-like', 
                         '11' = 'Myeloid', 
                         '12' = 'Myeloid',
                         '13' = 'Neuronal-like'),
  
  '0010619-Region_2' = c('0' = 'Ependymal-like',
                         '1' = 'Ependymal-like',
                         '2' = 'Neuronal-like',
                         '3' = 'Myeloid', 
                         '4' = 'Neuroepithelial-like',
                         '5' = 'VLMCs',
                         '6' = 'VLMCs',
                         '7' = 'VLMCs', 
                         '8' = 'Ependymal-like', 
                         '9' = 'Neuronal-like',
                         '10' = 'Neuronal-like'),
  
  '0010619-Region_3' = c('0' = 'Ependymal-like',
                         '1' = 'Neuronal-like',
                         '2' = 'Ependymal-like',
                         '3' = 'VLMCs', 
                         '4' = 'Ependymal-like',
                         '5' = 'Myeloid',
                         '6' = 'Neuroepithelial-like',
                         '7' = 'VLMCs', 
                         '8' = 'Neuronal-like', 
                         '9' = 'Radial glia-like',
                         '10' = 'Neuronal-like'),
  
  '0010619-Region_4' = c('0' = 'Ependymal-like',
                         '1' = 'Endothelial',
                         '2' = 'Neuronal-like',
                         '3' = 'Myeloid', 
                         '4' = 'Neuronal-like',
                         '5' = 'Ependymal-like',
                         '6' = 'Ependymal-like',
                         '7' = 'VLMCs', 
                         '8' = 'T-cell', 
                         '9' = 'Neuronal-like',
                         '10' = 'Neuronal-like'),
  
  '0010619-Region_5' = c('0' = 'Ependymal-like',
                         '1' = 'Myeloid',
                         '2' = 'Ependymal-like',
                         '3' = 'Neuronal-like', 
                         '4' = 'Neuronal-like',
                         '5' = 'Ependymal-like',
                         '6' = 'Neuronal-like',
                         '7' = 'Endothelial', 
                         '8' = 'Neuronal-like', 
                         '9' = 'Endothelial',
                         '10' = 'T-cell',
                         '11' = 'Myeloid',
                         '12' = 'Neuronal-like',
                         '13' = 'VLMCs',
                         '14' = 'Endothelial'),
  
  '0010501-Region_1' = c('0' = 'Ependymal-like',
                         '1' = 'Neuronal-like',
                         '2' = 'Mesenchymal',
                         '3' = 'Myeloid',
                         '4' = 'Neuroepithelial-like',
                         '5' = 'Neuronal-like',
                         '6' = 'Endothelial',
                         '7' = 'Neuronal-like',
                         '8' = 'Neuronal-like',
                         '9' = 'T-cell'),
  
  '0010501-Region_2' = c('0' = 'Ependymal-like',
                         '1' = 'Neuroepithelial-like',
                         '2' = 'Mesenchymal',
                         '3' = 'Neuronal-like',
                         '4' = 'Myeloid',
                         '5' = 'Endothelial',
                         '6' = 'Myeloid',
                         '7' = 'Neuronal-like',
                         '8' = 'Ependymal-like',
                         '9' = 'Neuronal-like',
                         '10' = 'Ependymal-like'),
  
  '0010814-Region_1' = c('0' = 'Ependymal-like',
                         '1' = 'Neuronal-like',
                         '2' = 'Neuronal-like',
                         '3' = 'Endothelial',
                         '4' = 'Neuronal-like',
                         '5' = 'Myeloid',
                         '6' = 'Radial glia-like',
                         '7' = 'Neuroepithelial-like',
                         '8' = 'Ependymal-like',
                         '9' = 'Unassigned',
                         '10' = 'Radial glia-like',
                         '11' = 'Ependymal-like',
                         '12' = 'Unassigned',
                         '13' = 'Unassigned'), 
  
  '0010814-Region_2' = c('0' = 'Ependymal-like',
                         '1' = 'Neuronal-like',
                         '2' = 'Endothelial',
                         '3' = 'Neuroepithelial-like',
                         '4' = 'Ependymal-like',
                         '5' = 'Neuronal-like',
                         '6' = 'Myeloid',
                         '7' = 'Mesenchymal',
                         '8' = 'Endothelial',
                         '9' = 'Ependymal-like',
                         '10' = 'Unassigned',
                         '11' = 'Ependymal-like',
                         '12' = 'Neuronal-like',
                         '13' = 'Unassigned',
                         '14' = 'Unassigned',
                         '15' = 'Unassigned',
                         '16' = 'Unassigned',
                         '17' = 'Unassigned',
                         '18' = 'Unassigned',
                         '19' = 'Unassigned',
                         '20' = 'Unassigned'),
  
  '0010498-Region_1' = c('0' = 'Ependymal-like', 
                         '1' = 'Neuroepithelial-like', 
                         '2' = 'Mesenchymal', 
                         '3' = 'Neuronal-like', 
                         '4' = 'Myeloid', 
                         '5' = 'Microglia',
                         '6' = 'Microglia', 
                         '7' = 'Microglia', 
                         '8' = 'Microglia', 
                         '9' = 'Mesenchymal', 
                         '10' = 'Microglia', 
                         '11' = 'Neuroepithelial-like',
                         '12' = 'Microglia', 
                         '13' = 'Ependymal-like', 
                         '14' = 'Neuroepithelial-like', 
                         '15' = 'Neuronal-like',
                         '16' = 'Microglia',
                         '17' = 'Neuronal-like',
                         '18' = 'Myeloid',
                         '19' = 'Oligodendrocytes'),
  
  '0010775-Region_1' =  c('0' = 'Neuronal-like', 
                          '1' = 'Myeloid', 
                          '2' = 'Neuronal-like', 
                          '3' = 'Ependymal-like', 
                          '4' = 'VLMCs', 
                          '5' = 'Neuronal-like',
                          '6' = 'Mesenchymal', 
                          '7' = 'Neuronal-like', 
                          '8' = 'Neuronal-like', 
                          '9' = 'Neuronal-like', 
                          '10' = 'T-cell', 
                          '11' = 'Neuronal-like',
                          '12' = 'Neuronal-like', 
                          '13' = 'Neuronal-like'),
  
  '0010540-Region_1' = c('0' = 'Neuronal-like', 
                         '1' = 'Ependymal-like', 
                         '2' = 'Ependymal-like', 
                         '3' = 'Neuronal-like', 
                         '4' = 'Neuronal-like', 
                         '5' = 'Neuronal-like',
                         '6' = 'Neuronal-like', 
                         '7' = 'Ependymal-like', 
                         '8' = 'Neuroepithelial-like', 
                         '9' = 'Myeloid', 
                         '10' = 'Neuronal-like', 
                         '11' = 'Mesenchymal',
                         '12' = 'Neuronal-like', 
                         '13' = 'Neuroepithelial-like'),
  
  '0010540-Region_2' =  c('0' = 'Neuronal-like', 
                          '1' = 'Neuronal-like', 
                          '2' = 'Ependymal-like', 
                          '3' = 'Ependymal-like', 
                          '4' = 'Neuronal-like', 
                          '5' = 'Neuronal-like',
                          '6' = 'Neuronal-like', 
                          '7' = 'Neuronal-like', 
                          '8' = 'Ependymal-like', 
                          '9' = 'Neuroepithelial-like', 
                          '10' = 'Mesenchymal', 
                          '11' = 'Myeloid',
                          '12' = 'Neuronal-like', 
                          '13' = 'Neuroepithelial-like', 
                          '14' = 'Neuroepithelial-like', 
                          '15' = 'Myeloid'),
  
  '0010540-Region_3' =  c('0' = 'Neuronal-like', 
                          '1' = 'Ependymal-like', 
                          '2' = 'Mesenchymal', 
                          '3' = 'Neuronal-like', 
                          '4' = 'Neuronal-like', 
                          '5' = 'Ependymal-like',
                          '6' = 'Neuronal-like', 
                          '7' = 'Neuronal-like', 
                          '8' = 'Neuronal-like', 
                          '9' = 'Ependymal-like', 
                          '10' = 'Neuroepithelial-like', 
                          '11' = 'Myeloid',
                          '12' = 'Unassigned', 
                          '13' = 'Neuroepithelial-like', 
                          '14' = 'Unassigned', 
                          '15' = 'Mesenchymal',
                          '16' = 'Neuronal-like'),
  
  '0010540-Region_4' =  c('0' = 'Ependymal-like', 
                          '1' = 'Ependymal-like', 
                          '2' = 'Neuronal-like', 
                          '3' = 'Ependymal-like', 
                          '4' = 'Myeloid', 
                          '5' = 'Endothelial',
                          '6' = 'Mesenchymal', 
                          '7' = 'Neuronal-like', 
                          '8' = 'Neuronal-like', 
                          '9' = 'T-cell', 
                          '10' = 'Myeloid'),
  
  '0010553-Region_1' = c('0' = 'Myeloid', 
                         '1' = 'Ependymal-like', 
                         '2' = 'Ependymal-like', 
                         '3' = 'Neuronal-like', 
                         '4' = 'Ependymal-like', 
                         '5' = 'Neuroepithelial-like',
                         '6' = 'Myeloid', 
                         '7' = 'Neuronal-like', 
                         '8' = 'Ependymal-like', 
                         '9' = 'Endothelial', 
                         '10' = 'Myeloid',
                         '11' = 'Myeloid',
                         '12' = 'Myeloid',
                         '13' = 'Neuronal-like'),
  
  '0010553-Region_2' = c('0' = 'Ependymal-like', 
                         '1' = 'Neuronal-like', 
                         '2' = 'Ependymal-like', 
                         '3' = 'Ependymal-like', 
                         '4' = 'Neuroepithelial-like', 
                         '5' = 'Neuronal-like',
                         '6' = 'Neuronal-like', 
                         '7' = 'Myeloid', 
                         '8' = 'Neuronal-like', 
                         '9' = 'Myeloid', 
                         '10' = 'Endothelial',
                         '11' = 'Neuronal-like',
                         '12' = 'Ependymal-like',
                         '13' = 'Myeloid',
                         '14' = 'Neuronal-like'),
  
  '0010553-Region_3' = c('0' = 'Neuronal-like', 
                         '1' = 'Neuronal-like', 
                         '2' = 'Neuronal-like', 
                         '3' = 'Ependymal-like', 
                         '4' = 'Mesenchymal', 
                         '5' = 'Mesenchymal',
                         '6' = 'Neuronal-like', 
                         '7' = 'Ependymal-like', 
                         '8' = 'Neuronal-like', 
                         '9' = 'Neuroepithelial-like', 
                         '10' = 'Neuroepithelial-like',
                         '11' = 'Myeloid',
                         '12' = 'Neuronal-like',
                         '13' = 'Neuroepithelial-like'),
  
  '0010553-Region_4' = c('0' = 'Neuronal-like', 
                         '1' = 'Ependymal-like', 
                         '2' = 'Ependymal-like', 
                         '3' = 'Endothelial', 
                         '4' = 'Myeloid', 
                         '5' = 'Neuroepithelial-like',
                         '6' = 'Neuronal-like', 
                         '7' = 'T-cell', 
                         '8' = 'Myeloid', 
                         '9' = 'Neuronal-like')
)





# Process Xenium data  -------------------------------------

# Read data
for (i in seq_along(metadata$Sample)) {
  data <- qread(file.path(base_dir, paste0('data/processed_data_Carlos/',  metadata$SamplePath[i], "/", metadata$SampleName[i], '.qs')))
  print(paste0('Analyzing: ', metadata$SampleID[i],' ' , metadata$SampleName[i]))
  
  # Change name identities
  Idents(data) <- metadata$Resolution[i]
  data <- RenameIdents(data, annotation_clusters[[i]])
  data[["group"]] <- Idents(data)
  
  # Add information about malignant or non-malignant
  data$malignant <- ifelse(data$group %in% c('Neuronal-like', 'Ependymal-like', 'Neuroepithelial-like',
                                             'Mesenchymal', 'Radia glia-like'), "Malignant", "Non-malignant")
  
  # Perform niche analysis
  data <- BuildNicheAssay(object = data, fov = "fov", group.by = "group", niches.k = 4, neighbors.k = 20)

  
  # save annotated object
  qsave(data, file.path(seurat_obj_dir, paste0('seurat_obj_', metadata$SampleID[i], "_", metadata$SampleName[i], '.qs')))
}

