args <- commandArgs(trailingOnly = TRUE)

current_sample <- args[1]


## Loading packages
library(Seurat)
library(NMF)
library(tidyverse)
library(pagoda2)
library(glue)
library(readxl)



## Sourcing necessary functions
hvgPagoda <- function(cm, n.OdGenes=3000, gam.k=10, plot_var=TRUE, n.cores=6, verbose=TRUE){
  pagoda_obj <- Pagoda2$new(x = Matrix(as.matrix(cm), sparse=TRUE), n.cores=n.cores)
  ## Adjust the variance
  pagoda_obj$adjustVariance(plot=plot_var, gam.k=gam.k)
  pagoda_hvg <- pagoda_obj$misc$varinfo
  pagoda_hvg <- pagoda_hvg[order(pagoda_hvg$lpa), ]
  if (verbose){
    cat(paste0("The number of highly variable genes identified by Pagoda2: ", length(pagoda_obj$misc$odgenes), '\n'))
  }
  return(pagoda_hvg)
}

nmf_df_preprocessing <- function(cm){
  ## convert negative values to zero
  cm <- ifelse(cm < 0, 0, cm)
  ## remove genes with zeros in all cells 
  cm <- cm[Matrix::rowSums(cm) != 0,]
}


## Defining data directory to save outputs
base_dir <- "/n/scratch/users/s/sad167/EPN/scRNAseq"


## Reading data
metadata <- read_excel(file.path(base_dir, 'metadata.xlsx')) 

# subset to frozen
metadata <- metadata %>%
  filter(`sc/snRNAseq` == 'snRNA-seq') 

FileName <- unique(metadata$SampleName)

idx <- which(metadata$SampleName == current_sample)

data_dir <- file.path(base_dir, 'analysis/NMF/ranks')
if (!dir.exists(data_dir)){dir.create(data_dir, recursive = T)}

cm <- readRDS(file.path(base_dir, 'analysis/NMF/counts', glue('{current_sample}.rds')))
cm <- cm[-grep(pattern = "^MT-|^RPS|^RPL", x = rownames(cm)), ]

## Finding High Variable Genes (HGVs)
ode <- hvgPagoda(cm, plot_var = F)
genes_to_keep <- rownames(ode)[1:10000]


## Keep only tumor cells and top 10000 highly variable genes for raw and log transformed data
cm <- cm[genes_to_keep, ]
cm_norm <- as.matrix(log2(cm/10+1))
cm_mean <- log2(Matrix::rowMeans(cm)+1)
cm_center <- cm_norm - rowMeans(cm_norm)


## Checking which HVGs are in current sample
cm_sample <- nmf_df_preprocessing(cm_center)
rm(cm, cm_center, cm_norm, cm_mean, metadata, genes_to_keep, ode); gc()


## NMF with optimal rank
start_time <- Sys.time()
message(paste0("Current sample: ", current_sample, " Starting time: ", start_time))
nmf_rank <- nmf(cm_sample, method = 'snmf/r', rank = 6, nrun = 100, seed = 123456, .options = "tvP7")
saveRDS(nmf_rank, file = file.path(data_dir, glue('{current_sample}.rds')))
end_time <- Sys.time()
message(paste0("Ending time: ", end_time - start_time))