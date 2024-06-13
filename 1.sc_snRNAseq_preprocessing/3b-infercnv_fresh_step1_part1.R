library(data.table)
library(infercnv)
library(tidyverse)
library(crayon)
library(ape)
library(readxl)

options(scipen = 100)

base_dir <- "/n/scratch/users/s/sad167/EPN/scRNAseq"
resources_dir <- file.path(base_dir, 'scripts/resources')
ext_ctrl <- file.path(resources_dir, 'ctrl_data')

source(file.path(resources_dir, "single_cell_preprocessing_helper_functions_CBJr.R"))
source(file.path(resources_dir, "inferCNV_helper_functions.R"))

metadata <- read_excel(file.path(base_dir, 'metadata.xlsx')) 

type <- c('frozen', 'fresh')

tp <- 'fresh'
metadata_subset <- metadata %>%
  filter(QC_pas == 'TRUE')  %>%
  filter(Type == tp) 

infercnv_dir <- file.path(base_dir, 'analysis/infercnv', tp)
if (!dir.exists(infercnv_dir)){dir.create(infercnv_dir, recursive = T)}

if (tp == 'frozen') {
  gene_order_file <- file.path(resources_dir, 'hg19_clean_gen_pos.txt') ## frozen
} else if (tp == 'fresh') {
  gene_order_file <- file.path(resources_dir, 'gencode_v19_gene_pos.txt') ## fresh
}

## Load sample names and count matrix
cm_raw <- pbapply::pblapply(file.path(base_dir, 'analysis/qc/data/individual', paste0(metadata_subset$SampleName, '-cm_list.rds')), function(x) {
  readRDS(x)$raw_data
})
cm_raw <- do.call(cbind, cm_raw)

samples <- lapply(file.path(base_dir, 'analysis/qc/data/individual', paste0(metadata_subset$SampleName, '-samples.rds')), function(x) {
  readRDS(x)
})
samples <- unlist(samples)

## Name original samples vector
orig_samples <- samples

## Load chromosome orders of genes
gene_order <- read.table(gene_order_file, header = F, sep = "\t", row.names = 1)

## Load and name normal control data, merge them with all samples and save the complete df
if (tp == 'frozen') {
  addNormalControl(ext_ctrl, infercnv_dir, cm_raw, orig_samples, type = 'frozen') ## frozen
} else if (tp == 'fresh') {
  addNormalControl(ext_ctrl, infercnv_dir, cm_raw, orig_samples, type = 'fresh') ## fresh
}
rm(cm_raw); gc()


## Load complete df for cm and samples
cm <- readRDS(file.path(infercnv_dir, "cm_exp_ctrl.rds"))
samples <- readRDS(file.path(infercnv_dir, "samples_exp_ctrl.rds"))

## Prepare input data (cm and annotations) for inferCNV
## Construct a vector of initial characterization of malignancy for each cell
## All normal controls are normal (either mg or od)
## All cells from patient samples are treated as malignant at this moment
malig_status <- ifelse(samples == "mg" | samples == "od", samples, paste("malignant", samples, sep="_"))

## Write out cm and malignancy_status as inputs for inferCNV
fwrite(cm, file = file.path(infercnv_dir, 'counts.matrix'), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

## Write tab deliminated cell annotations
write.table(malig_status, file = file.path(infercnv_dir, "cellAnnotations.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

## Run inferCNV
set.seed(1234)
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = file.path(infercnv_dir, "counts.matrix"),
                                     annotations_file = file.path(infercnv_dir, "cellAnnotations.txt"),
                                     delim = "\t",
                                     gene_order_file = gene_order_file,
                                     ref_group_names = c("mg", "od"))

# perform infercnv operations to reveal cnv signal
out_dir <- file.path(infercnv_dir, "out_dir")
set.seed(1234)
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir = out_dir,  # dir is auto-created for storing outputs
                              cluster_by_groups = T,   # cluster
                              denoise = T,
                              HMM = T,
                              plot_steps = F,
                              num_threads = 10,
                              analysis_mode = 'subclusters',
                              tumor_subcluster_partition_method = 'random_trees',
                              tumor_subcluster_pval = 0.05,
                              write_expr_matrix = F,
                              no_plot = T)

## Saving observations and references predicted tumor cells files
save_cnv_files(infercnv_obj, out_dir = out_dir)

## Normalized gene expressions of predicted tumor cells
observations <- fread(file.path(out_dir, 'infercnv.observations.txt.gz'), sep = ",", data.table = F) %>%
  column_to_rownames('V1')

## Normalized gene expressions of predicted normal cells
references <- fread(file.path(out_dir, 'infercnv.references.txt.gz'), sep = ",", data.table = F) %>%
  column_to_rownames('V1')

## Concatnate
all_cnv_values <- cbind(observations, references)

## Labels for obs vs ref
obs_ref <- c(rep("obs", dim(observations)[2]), rep("ref", dim(references)[2]))
names(obs_ref) <- c(colnames(observations), colnames(references))

## CRITICAL: sort the order of cells in cnv_values and obs_ref to the same as cm
all_cnv_values <- all_cnv_values[,names(malig_status)]
obs_ref <- obs_ref[names(malig_status)]
all_cnv_values <-  all_cnv_values - 1

saveRDS(all_cnv_values, file.path(infercnv_dir, "all_cnv_values.rds"))
saveRDS(obs_ref, file.path(infercnv_dir, "obs_ref.rds"))

## Compute CNV correlations
mean_cnv_values <- calculateAverageCnvValue(all_cnv_values, samples, malig_status)
cor_coef <- calculateCnvCor(all_cnv_values, mean_cnv_values, samples)
saveRDS(cor_coef, file.path(infercnv_dir, "cor_coef.rds"))
