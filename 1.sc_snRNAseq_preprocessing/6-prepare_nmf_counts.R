library(Seurat)
library(tidyverse)
library(glue)
library(readxl)


base_dir <- "/n/scratch/users/s/sad167/EPN/scRNAseq"
resources_dir <- file.path(base_dir, 'scripts/resources')

source(file.path(resources_dir, 'single_cell_preprocessing_helper_functions.R'))

## Dealing with counts for each sample
metadata <- read_excel(file.path(base_dir, 'metadata.xlsx')) 

data_frozen <- readRDS("/n/scratch/users/s/sad167/EPN/scRNAseq/analysis/infercnv/frozen/part_1/seurat_obj_malignant.rds")
data_frozen <- AddMetaData(data_frozen, metadata = 'frozen', col.name = 'type')
data <- data_frozen

nmf_dir <- file.path(base_dir, 'analysis/NMF/counts')
if (!dir.exists(nmf_dir)){dir.create(nmf_dir, recursive = T)}

tmp <- data

for (smp in unique(tmp$sample)) {
  tmp2 <- subset(tmp, subset = (sample == smp))
  counts <- tmp2[["RNA"]]$counts
  saveRDS(counts, file.path(nmf_dir, glue('{smp}.rds')))
}
