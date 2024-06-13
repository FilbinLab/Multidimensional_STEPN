library(Seurat)
library(tidyverse)
library(cluster)
#library(factoextra)
library(dendextend)
library(weights)
library(ggpubr)
library(matrixStats)
library(readxl)

base_dir <- "/n/scratch/users/s/sad167/EPN/scRNAseq"
resources_dir <- file.path(base_dir, 'scripts/resources')

source(file.path(resources_dir, "single_cell_preprocessing_helper_functions_CBJr.R"))
source(file.path(resources_dir, 'inferCNV_helper_functions.R'))
source(file.path(resources_dir, 'Plotting_helper_functions.R'))

metadata <- read_excel(file.path(base_dir, 'metadata.xlsx')) 


tp <- 'fresh'
metadata_subset <- metadata %>%
  filter(Type == tp)


qc_data <- file.path(base_dir, "analysis/qc/data")
infercnv_dir <- file.path(base_dir, 'analysis/infercnv/fresh')
infercnv_plots <- file.path(base_dir, 'analysis/infercnv/fresh/plots')
if (!dir.exists(infercnv_dir)){dir.create(infercnv_dir, recursive = T)}
if (!dir.exists(infercnv_plots)){dir.create(infercnv_plots, recursive = T)}


cm <- readRDS(file.path(infercnv_dir, "cm_exp_ctrl.rds"))
samples <- readRDS(file.path(infercnv_dir, "samples_exp_ctrl.rds"))
all_cnv_values <- readRDS(file.path(infercnv_dir, "all_cnv_values.rds"))
obs_ref <- readRDS(file.path(infercnv_dir, "obs_ref.rds"))


## Hierarchical clustering
cnv_hc_each_sample <- list()
cnv_cutree_each_sample <- list()
weights_list <- list()

for (sample in sort(unique(samples))){
  ## Subset each sample + (mg and od)
  to_subset <- (samples == sample | samples == "mg" | samples == 'od')
  message("======================================")
  message("Processing sample: ", sample)
  message("Start time: ", Sys.time())
  ## Compute weights
  weights <- rowVars(as.matrix(all_cnv_values[,to_subset]), useNames = T)
  weights_list[[sample]] <- weights
  ## Compute pairwise correlations and distances
  pairwise_cor <- cov.wt(all_cnv_values[,to_subset], wt = weights, cor = TRUE)$cor
  pairwise_dist <- 1 - pairwise_cor
  ## Hierarchical clustering and cutree
  hc <- hclust(as.dist(pairwise_dist), method="ward.D2")
  cnv_hc_each_sample[[sample]] <- hc
  cnv_cutree_each_sample[[sample]] <- cutree(hc, 4)
  message("End time: ", Sys.time())
}
saveRDS(cnv_hc_each_sample, file.path(infercnv_dir, "cnv_hc_raw.rds"))

## Print out cutree results
cnv_cutree_table <- list()
for (sample in sort(unique(samples))){
  cnv_cutree_table[[sample]] <- table(samples[names(cnv_cutree_each_sample[[sample]])], cnv_cutree_each_sample[[sample]])
}
names(cnv_cutree_table) <- sort(unique(samples))

## Store normal clusters in a list
normal_clusters <- lapply(cnv_cutree_table, function(x) which((x["mg",] + x["od",])>=10))
try (if(any(names(cnv_cutree_each_sample) != names(normal_clusters))) stop("cutree results list and normal clusters list should have identical names for elements"))

## Classify each cell as Normal (F) vs Tumor (T) based on normal_clusters list
call_cnv_res <- NULL
for (i in seq(length(normal_clusters))){
  sample <- names(normal_clusters)[i]
  #print(sample)
  normal_clust <- normal_clusters[[i]]
  #print(normal_clust)
  hc_clust_res <- cnv_cutree_each_sample[[i]]
  #print(head(hc_clust_res))
  sample_names <- sapply(names(hc_clust_res), function(x) unlist(strsplit(x, split = ".", fixed = T))[1])
  #print(head(sample_names))
  hc_clust_res <- hc_clust_res[sample_names == sample]
  tmp <- ifelse(hc_clust_res %in% normal_clust, F, T)
  names(tmp) <- names(hc_clust_res)[sample_names == sample]
  call_cnv_res <- c(call_cnv_res, tmp)
}

## Generate a named vector for mg and od, and label them as normal
mg_res <- structure(rep(F, length(grep('mg', samples))), names = names(samples)[samples == "mg"])
od_res <- structure(rep(F, length(grep('od', samples))), names = names(samples)[samples == "od"])

## Concatenate CNV calling for all samples and normal control
call_cnv_res_complete <- c(call_cnv_res, mg_res, od_res)
call_cnv_res_complete <- call_cnv_res_complete[names(samples)]
call_cnv_res <- call_cnv_res_complete[names(samples)[obs_ref == "obs"]]

## Store CNV calling for all samples alone and all samples + normal control
saveRDS(call_cnv_res, file.path(infercnv_dir, "call_cnv.rds"))
saveRDS(call_cnv_res_complete, file.path(infercnv_dir, "call_cnv_w_ctrl.rds"))
