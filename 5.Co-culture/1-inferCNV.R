library(biomaRt)
library(Seurat)
library(tidyverse)
library(qs)
library(glue)
library(data.table)
library(infercnv)
library(matrixStats)

source('~/Projects/General-Codes/Resources/Plotting_helper_functions.R')
source('~/Projects/General-Codes/Resources/single_cell_preprocessing_helper_functions.R')
source('~/Projects/General-Codes/Resources/inferCNV_helper_functions.R')


################################################################################
## Retrieve orthologs genes (rat x human)
################################################################################
rat_mart <- useEnsembl("ensembl", dataset = "rnorvegicus_gene_ensembl", mirror = "useast")

orthologs_rat_human <- getBM(
  attributes = c(
    "ensembl_gene_id", "external_gene_name", "description",
    "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name", 
    "hsapiens_homolog_orthology_type", "hsapiens_homolog_perc_id"
  ),
  mart = rat_mart
)

# Optional: keep only ortholog_one2one relationships (high confidence)
orthologs_rat_human_filtered <- orthologs_rat_human %>% 
  filter(hsapiens_homolog_orthology_type == "ortholog_one2one", external_gene_name != '', hsapiens_homolog_associated_gene_name != '') %>% 
  dplyr::select(c('external_gene_name', 'hsapiens_homolog_associated_gene_name')) %>% 
  rename_all(~c('rat_genes', 'human_genes'))



################################################################################
################################################################################
#> Case 1: patient models using 7LP as reference
################################################################################
################################################################################
## Loading samples
data <- qread('/n/scratch/users/c/cao385/Ependymoma/10x/infercnv/patient_models.qs')
data$sample <- as.character(data$sample)
orig_samples <- data$sample


## Loading reference
ref <- qread('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-7LP.qs') %>%
  RenameCells(add.cell.id = '7LP')
ref$Annotation <- 'normal'
rownames(ref) <- plyr::mapvalues(rownames(ref), from = orthologs_rat_human_filtered$rat_genes, to = orthologs_rat_human_filtered$human_genes, warn_missing = F)


## Adjusting orthologs genes for query and ref
genes_to_use <- intersect(rownames(data), rownames(ref))

data <- subset(data, features = genes_to_use)
ref <- subset(ref, features = genes_to_use)


cm_raw <- LayerData(data, 'counts')[intersect(rownames(data), rownames(ref)), ]
ext_ctrl_cm <- LayerData(ref, 'counts')[intersect(rownames(data), rownames(ref)), ]

normal_samples <- structure(rep('normal', ncol(ref)), names = colnames(ref)) ## Generate control sample list

ext_ctrl_cm <- ext_ctrl_cm[,names(normal_samples)] ## Reorder the ext_ctrl cm and subset genes
ext_ctrl_cm <- ext_ctrl_cm[rownames(cm_raw),] ## Reorder the ext_ctrl cm and subset genes

cm <- cbind(cm_raw, ext_ctrl_cm) ## Concatenate samples and normal controls
samples <- c(structure(unname(orig_samples), names = names(orig_samples)), normal_samples) ## Concatenate samples and normal controls

#malig_status <- ifelse(samples == "normal", samples, paste("malignant", samples, sep="_"))
malig_status <- ifelse(samples == "normal", samples, samples)

## Write out cm and malignancy_status as inputs for inferCNV
fwrite(as.data.frame(cm), file = "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv/counts.matrix.gz", sep = '\t', quote = F, row.names = T, col.names = T, nThread = 4, compress = 'auto')
## Write tab deliminated cell annotations
fwrite(as.data.frame(malig_status), file = "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv/cellAnnotations.txt.gz", sep = "\t", quote = FALSE, row.names = T, col.names = FALSE, compress = 'auto')


set.seed(1234)
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv/counts.matrix.gz",
                                     annotations_file = "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv/cellAnnotations.txt.gz",
                                     delim = "\t",
                                     gene_order_file = "/n/scratch/users/c/cao385/Ependymoma/10x/gene_order.txt",
                                     ref_group_names = 'normal', min_max_counts_per_cell = c(-Inf,+Inf))

# perform infercnv operations to reveal cnv signal
set.seed(1234)
out_dir <- "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv/out_dir2"
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir = out_dir,  # dir is auto-created for storing outputs
                              cluster_by_groups = T,   # cluster
                              denoise = T,
                              HMM = T,
                              plot_steps = F,
                              num_threads = 12,
                              analysis_mode = 'samples',
                              tumor_subcluster_partition_method = 'leiden',
                              tumor_subcluster_pval = 0.05,
                              leiden_resolution = 0.01,
                              png_res = 600,
                              write_expr_matrix = T)


## Normalized gene expressions of predicted tumor cells
observations <- fread(file.path(out_dir, 'infercnv.observations.txt'), data.table = F) %>%
  column_to_rownames('V1')

## Normalized gene expressions of predicted normal cells
references <- fread(file.path(out_dir, 'infercnv.references.txt'), data.table = F) %>%
  column_to_rownames('V1')

## Concatnate
all_cnv_values <- cbind(observations, references)

## Labels for obs vs ref
obs_ref <- c(rep("obs", dim(observations)[2]), rep("ref", dim(references)[2]))
names(obs_ref) <- c(colnames(observations), colnames(references))

## CRITICAL: sort the order of cells in cnv_values and obs_ref to the same as cm
all_cnv_values <- all_cnv_values[,names(malig_status)[names(malig_status) %in% colnames(all_cnv_values)]]
obs_ref <- obs_ref[names(malig_status)[names(malig_status) %in% colnames(all_cnv_values)]]
all_cnv_values <-  all_cnv_values - 1

qsave(all_cnv_values, "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv/all_cnv_values.qs")
qsave(obs_ref, "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv/obs_ref.qs")

counts <- cbind(cm_raw, ext_ctrl_cm)
data <- RunFullSeurat_v5(cm = counts, metadata = malig_status, doBatch = T, var2batch = 'sample', batchMethod = 'harmony', verbose = T, project = 'Daeun')

p1 <- DimPlot(data, reduction = 'umap.unintegrated', group.by = 'sample', cols = colors_to_use) + ggtitle('Unintegrated') + seurat_theme()
p2 <- DimPlot(data, reduction = 'umap.harmony', group.by = 'sample', cols = colors_to_use) + ggtitle('Harmony Integration') + seurat_theme()
p1 + p2



################################################################################
################################################################################
#> Case 2: coculture rat only samples using 7LP as reference
################################################################################
################################################################################
## Loading samples
data_tmp <- qread('/n/scratch/users/c/cao385/Ependymoma/10x/infercnv/patient_models.qs')
data_tmp$sample <- as.character(data_tmp$sample)
data_tmp <- subset(data_tmp, subset = sample %in% c("BT165_Coculture_rats", "EP1NS_Coculture_rats", "VBT242_Coculture_rats"))


tmp_BT2214Rat <- qread(glue('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-BT2214Rat.qs'))
tmp_BT2214Rat <- subset(tmp_BT2214Rat, subset = species == 'Human')
tmp_BT2214Rat <- subset(tmp_BT2214Rat, features = grep('^GRCh38----', rownames(tmp_BT2214Rat), value = T))
tmp_BT2214Rat$sample <- 'BT2214Rat'
rownames(tmp_BT2214Rat) <- gsub("GRCh38----", "", rownames(tmp_BT2214Rat))
tmp_BT2214Rat <- subset(tmp_BT2214Rat, features = grep("^ENSG", rownames(tmp_BT2214Rat), value = TRUE, invert = T))


tmp_EP1NSminusTTXRat <- qread(glue('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-EP1NSminusTTXRat.qs'))
tmp_EP1NSminusTTXRat <- subset(tmp_EP1NSminusTTXRat, subset = species == 'Human')
tmp_EP1NSminusTTXRat <- subset(tmp_EP1NSminusTTXRat, features = grep('^GRCh38----', rownames(tmp_EP1NSminusTTXRat), value = T))
tmp_EP1NSminusTTXRat$sample <- 'EP1NSminusTTXRat'
rownames(tmp_EP1NSminusTTXRat) <- gsub("GRCh38----", "", rownames(tmp_EP1NSminusTTXRat))
tmp_EP1NSminusTTXRat <- subset(tmp_EP1NSminusTTXRat, features = grep("^ENSG", rownames(tmp_EP1NSminusTTXRat), value = TRUE, invert = T))


tmp_EP1NSplusTTXRat <- qread(glue('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-EP1NSplusTTXRat.qs'))
tmp_EP1NSplusTTXRat <- subset(tmp_EP1NSplusTTXRat, subset = species == 'Human')
tmp_EP1NSplusTTXRat <- subset(tmp_EP1NSplusTTXRat, features = grep('^GRCh38----', rownames(tmp_EP1NSplusTTXRat), value = T))
tmp_EP1NSplusTTXRat$sample <- 'EP1NSplusTTXRat'
rownames(tmp_EP1NSplusTTXRat) <- gsub("GRCh38----", "", rownames(tmp_EP1NSplusTTXRat))
tmp_EP1NSplusTTXRat <- subset(tmp_EP1NSplusTTXRat, features = grep("^ENSG", rownames(tmp_EP1NSplusTTXRat), value = TRUE, invert = T))


## Merging datasets to use
data <- merge(data_tmp, list(tmp_BT2214Rat, tmp_EP1NSminusTTXRat, tmp_EP1NSplusTTXRat))
data <- JoinLayers(data)

orig_samples <- data$sample


## Loading reference
ref <- qread('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-7LP.qs') %>%
  RenameCells(add.cell.id = '7LP')
ref$Annotation <- 'normal'
rownames(ref) <- plyr::mapvalues(rownames(ref), from = orthologs_rat_human_filtered$rat_genes, to = orthologs_rat_human_filtered$human_genes, warn_missing = F)


## Adjusting orthologs genes for query and ref
genes_to_use <- intersect(rownames(data), rownames(ref))

data <- subset(data, features = genes_to_use)
ref <- subset(ref, features = genes_to_use)


cm_raw <- LayerData(data, 'counts')[intersect(rownames(data), rownames(ref)), ]
ext_ctrl_cm <- LayerData(ref, 'counts')[intersect(rownames(data), rownames(ref)), ]

normal_samples <- structure(rep('normal', ncol(ref)), names = colnames(ref)) ## Generate control sample list

ext_ctrl_cm <- ext_ctrl_cm[,names(normal_samples)] ## Reorder the ext_ctrl cm and subset genes
ext_ctrl_cm <- ext_ctrl_cm[rownames(cm_raw),] ## Reorder the ext_ctrl cm and subset genes

cm <- cbind(cm_raw, ext_ctrl_cm) ## Concatenate samples and normal controls
samples <- c(structure(unname(orig_samples), names = names(orig_samples)), normal_samples) ## Concatenate samples and normal controls

#malig_status <- ifelse(samples == "normal", samples, paste("malignant", samples, sep="_"))
malig_status <- ifelse(samples == "normal", samples, samples)

## Write out cm and malignancy_status as inputs for inferCNV
fwrite(as.data.frame(cm), file = "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv3/counts.matrix.gz", sep = '\t', quote = F, row.names = T, col.names = T, nThread = 4, compress = 'auto')
## Write tab deliminated cell annotations
fwrite(as.data.frame(malig_status), file = "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv3/cellAnnotations.txt.gz", sep = "\t", quote = FALSE, row.names = T, col.names = FALSE, compress = 'auto')


set.seed(1234)
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv3/counts.matrix.gz",
                                     annotations_file = "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv3/cellAnnotations.txt.gz",
                                     delim = "\t",
                                     gene_order_file = "/n/scratch/users/c/cao385/Ependymoma/10x/gene_order.txt",
                                     ref_group_names = 'normal', min_max_counts_per_cell = c(-Inf,+Inf))

# perform infercnv operations to reveal cnv signal
set.seed(1234)
out_dir <- "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv3/out_dir"
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir = out_dir,  # dir is auto-created for storing outputs
                              cluster_by_groups = T,   # cluster
                              denoise = T,
                              HMM = T,
                              plot_steps = F,
                              num_threads = 12,
                              analysis_mode = 'samples',
                              tumor_subcluster_partition_method = 'leiden',
                              tumor_subcluster_pval = 0.05,
                              leiden_resolution = 0.01,
                              png_res = 600,
                              write_expr_matrix = T)


## Normalized gene expressions of predicted tumor cells
observations <- fread(file.path(out_dir, 'infercnv.observations.txt'), data.table = F) %>%
  column_to_rownames('V1')

## Normalized gene expressions of predicted normal cells
references <- fread(file.path(out_dir, 'infercnv.references.txt'), data.table = F) %>%
  column_to_rownames('V1')

## Concatnate
all_cnv_values <- cbind(observations, references)

## Labels for obs vs ref
obs_ref <- c(rep("obs", dim(observations)[2]), rep("ref", dim(references)[2]))
names(obs_ref) <- c(colnames(observations), colnames(references))

## CRITICAL: sort the order of cells in cnv_values and obs_ref to the same as cm
all_cnv_values <- all_cnv_values[,names(malig_status)[names(malig_status) %in% colnames(all_cnv_values)]]
obs_ref <- obs_ref[names(malig_status)[names(malig_status) %in% colnames(all_cnv_values)]]
all_cnv_values <-  all_cnv_values - 1

qsave(all_cnv_values, "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv3/all_cnv_values.qs")
qsave(obs_ref, "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv3/obs_ref.qs")

counts <- cbind(cm_raw, ext_ctrl_cm)
data <- RunFullSeurat_v5(cm = counts, metadata = malig_status, doBatch = T, var2batch = 'sample', batchMethod = 'harmony', verbose = T, project = 'Daeun')
#data <- AddMetaData(data, metadata = glue('cnv_{cnv_cutree_each_sample$`2LP`}'), col.name = 'cnv_clusters')
data$tumor <- ifelse(data$seurat_clusters %in% c(1,3,4,5,7,8,9,10,11,12,14,16,19,20), 'normal', 'tumor')
data$tdTomato <- ifelse(FetchData(data, 'tdTomato-1') %>% pull('tdTomato-1') > 0, 'yes', 'no')

p1 <- DimPlot(data, reduction = 'umap.unintegrated', group.by = 'sample', cols = colors_to_use, pt.size = 1) + ggtitle('Unintegrated') + seurat_theme()
p2 <- DimPlot(data, reduction = 'umap.unintegrated', group.by = 'tumor', cols = colors_to_use, pt.size = 1) + ggtitle('Unintegrated') + seurat_theme()
p3 <- DimPlot(data, reduction = 'umap.unintegrated', group.by = 'tdTomato', cols = colors_to_use, pt.size = 0.5, order = T) + ggtitle('Unintegrated') + seurat_theme()
p1 + p2 + p3

qsave(data, '/n/scratch/users/c/cao385/Ependymoma/10x/output/CocultureRat-post-inferCNV.qs')



################################################################################
################################################################################
#> Case 3: 2LP and EP1NSmono using 1LP as reference
################################################################################
################################################################################
## Loading samples
data_1 <- qread('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-2LP.qs')
data_1 <- RenameCells(data_1, add.cell.id = '2LP')
data_1$sample <- '2LP'

data_2 <- qread('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-EP1NSmono.qs')
data_2 <- RenameCells(data_2, add.cell.id = '12LP')
data_2$sample <- '12LP'

data <- merge(data_1, data_2)
data <- JoinLayers(data)

orig_samples <- data$sample


ref <- qread('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-1LP_filtered.qs') %>%
  RenameCells(add.cell.id = '1LP')
ref$Annotation <- 'normal'


cm_raw <- LayerData(data, 'counts')[intersect(rownames(data), rownames(ref)), ]
ext_ctrl_cm <- LayerData(ref, 'counts')[intersect(rownames(data), rownames(ref)), ]

normal_samples <- structure(rep('normal', ncol(ref)), names = colnames(ref)) ## Generate control sample list

ext_ctrl_cm <- ext_ctrl_cm[,names(normal_samples)] ## Reorder the ext_ctrl cm and subset genes
ext_ctrl_cm <- ext_ctrl_cm[rownames(cm_raw),] ## Reorder the ext_ctrl cm and subset genes

cm <- cbind(cm_raw, ext_ctrl_cm) ## Concatenate samples and normal controls
samples <- c(structure(unname(orig_samples), names = names(orig_samples)), normal_samples) ## Concatenate samples and normal controls

#malig_status <- ifelse(samples == "normal", samples, paste("malignant", samples, sep="_"))
malig_status <- ifelse(samples == "normal", samples, samples)

## Write out cm and malignancy_status as inputs for inferCNV
fwrite(as.data.frame(cm), file = "/n/scratch/users/c/cao385/Ependymoma/10x/inferCNV-2LP_12LP_ref-1LP/counts.matrix.gz", sep = '\t', quote = F, row.names = T, col.names = T, nThread = 4, compress = 'auto')
## Write tab deliminated cell annotations
fwrite(as.data.frame(malig_status), file = "/n/scratch/users/c/cao385/Ependymoma/10x/inferCNV-2LP_12LP_ref-1LP/cellAnnotations.txt.gz", sep = "\t", quote = FALSE, row.names = T, col.names = FALSE, compress = 'auto')


set.seed(1234)
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = "/n/scratch/users/c/cao385/Ependymoma/10x/inferCNV-2LP_12LP_ref-1LP/counts.matrix.gz",
                                     annotations_file = "/n/scratch/users/c/cao385/Ependymoma/10x/inferCNV-2LP_12LP_ref-1LP/cellAnnotations.txt.gz",
                                     delim = "\t",
                                     gene_order_file = "/n/scratch/users/c/cao385/Ependymoma/10x/gene_order.txt",
                                     ref_group_names = 'normal', min_max_counts_per_cell = c(-Inf,+Inf))

# perform infercnv operations to reveal cnv signal
set.seed(1234)
out_dir <- "/n/scratch/users/c/cao385/Ependymoma/10x/inferCNV-2LP_12LP_ref-1LP/out_dir2"
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir = out_dir,  # dir is auto-created for storing outputs
                              cluster_by_groups = T,   # cluster
                              denoise = T,
                              HMM = T,
                              plot_steps = F,
                              num_threads = 12,
                              analysis_mode = 'samples',
                              tumor_subcluster_partition_method = 'leiden',
                              tumor_subcluster_pval = 0.05,
                              leiden_resolution = 0.01,
                              png_res = 600,
                              write_expr_matrix = T)


## Normalized gene expressions of predicted tumor cells
observations <- fread(file.path(out_dir, 'infercnv.observations.txt'), data.table = F) %>%
  column_to_rownames('V1')

## Normalized gene expressions of predicted normal cells
references <- fread(file.path(out_dir, 'infercnv.references.txt'), data.table = F) %>%
  column_to_rownames('V1')

## Concatnate
all_cnv_values <- cbind(observations, references)

## Labels for obs vs ref
obs_ref <- c(rep("obs", dim(observations)[2]), rep("ref", dim(references)[2]))
names(obs_ref) <- c(colnames(observations), colnames(references))

## CRITICAL: sort the order of cells in cnv_values and obs_ref to the same as cm
all_cnv_values <- all_cnv_values[,names(malig_status)[names(malig_status) %in% colnames(all_cnv_values)]]
obs_ref <- obs_ref[names(malig_status)[names(malig_status) %in% colnames(all_cnv_values)]]
all_cnv_values <-  all_cnv_values - 1

qsave(all_cnv_values, "/n/scratch/users/c/cao385/Ependymoma/10x/inferCNV-2LP_12LP_ref-1LP/all_cnv_values.qs")
qsave(obs_ref, "/n/scratch/users/c/cao385/Ependymoma/10x/inferCNV-2LP_12LP_ref-1LP/obs_ref.qs")


counts <- cbind(cm_raw, ext_ctrl_cm)
data <- RunFullSeurat_v5(cm = counts, metadata = malig_status, doBatch = T, var2batch = 'sample', batchMethod = 'harmony', verbose = T, project = 'Daeun')

DimPlot(data, reduction = 'umap.unintegrated', group.by = 'unintegrated_clusters', cols = colors_to_use, pt.size = 1, label = T, label.box = T) + ggtitle('Unintegrated') + seurat_theme()

data$tumor <- ifelse(data$unintegrated_clusters %in% c(9,10,11,12,13,15,16,18,19), 'normal', 'tumor')
data$tdTomato <- ifelse(FetchData(data, 'tdTomato-1') %>% pull('tdTomato-1') > 0, 'yes', 'no')

p1 <- DimPlot(data, reduction = 'umap.unintegrated', group.by = 'sample', cols = colors_to_use, pt.size = 1) + ggtitle('Unintegrated') + seurat_theme()
p2 <- DimPlot(data, reduction = 'umap.unintegrated', group.by = 'tumor', cols = colors_to_use, pt.size = 1) + ggtitle('Unintegrated') + seurat_theme()
p3 <- DimPlot(data, reduction = 'umap.unintegrated', group.by = 'tdTomato', cols = colors_to_use, pt.size = 0.5, order = T) + ggtitle('Unintegrated') + seurat_theme()
p1 + p2 + p3

qsave(data, '/n/scratch/users/c/cao385/Ependymoma/10x/output/2LP_12LP-post-inferCNV.qs')


## 12LP
data_12LP <- subset(data, subset = sample == '12LP')
data_12LP <- RunFullSeurat_v5(cm = LayerData(data_12LP, 'counts'), metadata = data_12LP@meta.data, doBatch = F, verbose = T, project = '12LP')
qsave(data_12LP, '/n/scratch/users/c/cao385/Ependymoma/10x/output/12LP-post-inferCNV.qs')


## 2LP
data_2LP <- subset(data, subset = sample == '2LP')
data_2LP <- RunFullSeurat_v5(cm = LayerData(data_2LP, 'counts'), metadata = data_2LP@meta.data, doBatch = F, verbose = T, project = '2LP')
qsave(data_2LP, '/n/scratch/users/c/cao385/Ependymoma/10x/output/2LP-post-inferCNV-latest.qs')


data_2LP <- qread('/n/scratch/users/c/cao385/Ependymoma/10x/output/2LP-post-inferCNV-latest.qs')

set.seed(1234)
k = 50
knn.norm <- get.knn(as.matrix(Embeddings(data_2LP, 'umap')), k = k)
knn.norm <- data.frame(from = rep(1:nrow(knn.norm$nn.index), k), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
nw.norm <- graph_from_data_frame(knn.norm, directed = FALSE)
nw.norm <- simplify(nw.norm)

lc.norm <- cluster_louvain(nw.norm, resolution = 2)
lc.norm <- as.factor(membership(lc.norm))
names(lc.norm) <- rownames(Embeddings(data_2LP, 'umap'))
data_2LP[['louvain_clusters']] <- lc.norm

p1 <- DimPlot(data_2LP, reduction = 'umap', group.by = 'tumor', cols = colors_to_use, pt.size = 0.8) + seurat_theme()
p2 <- DimPlot(data_2LP, reduction = 'umap', group.by = 'louvain_clusters', cols = colors_to_use, pt.size = 0.8, label = T, label.box = T, repel = T) + seurat_theme() + NoLegend()
p1 + p2

data_2LP$tumor <- ifelse(data_2LP$louvain_clusters %in% c(5,9,10,15,22), 'normal', 'tumor')
DimPlot(data_2LP, reduction = 'umap', group.by = 'tumor', cols = colors_to_use, pt.size = 1) + ggtitle('2LP') + seurat_theme()
qsave(data_2LP, '/n/scratch/users/c/cao385/Ependymoma/10x/output/2LP-post-inferCNV-latest.qs')


p1 <- DimPlot(data_2LP, reduction = 'umap', group.by = 'tumor', cols = colors_to_use, pt.size = 1) + ggtitle('2LP') + seurat_theme()
p2 <- DimPlot(data_12LP, reduction = 'umap', group.by = 'tumor', cols = colors_to_use[2], pt.size = 1) + ggtitle('12LP') + seurat_theme()
pt <- p1 + p2
ggsave(plot = pt, glue('{base_dir}/2LP_12LP-UMAP.pdf'), width = 12, height = 5)



################################################################################
################################################################################
#> Case 4: 2LP using 1LP as reference
################################################################################
################################################################################
## Loading samples
data <- qread('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-2LP.qs')
data$sample <- '2LP'
orig_samples <- data$sample


## Loading reference
ref <- qread('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-1LP.qs') %>%
  RenameCells(add.cell.id = '1LP')
ref$Annotation <- 'normal'


cm_raw <- LayerData(data, 'counts')[intersect(rownames(data), rownames(ref)), ]
ext_ctrl_cm <- LayerData(ref, 'counts')[intersect(rownames(data), rownames(ref)), ]

normal_samples <- structure(rep('normal', ncol(ref)), names = colnames(ref)) ## Generate control sample list

ext_ctrl_cm <- ext_ctrl_cm[,names(normal_samples)] ## Reorder the ext_ctrl cm and subset genes
ext_ctrl_cm <- ext_ctrl_cm[rownames(cm_raw),] ## Reorder the ext_ctrl cm and subset genes

cm <- cbind(cm_raw, ext_ctrl_cm) ## Concatenate samples and normal controls
samples <- c(structure(unname(orig_samples), names = names(orig_samples)), normal_samples) ## Concatenate samples and normal controls

#malig_status <- ifelse(samples == "normal", samples, paste("malignant", samples, sep="_"))
malig_status <- ifelse(samples == "normal", samples, samples)

## Write out cm and malignancy_status as inputs for inferCNV
fwrite(as.data.frame(cm), file = "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv2/counts.matrix.gz", sep = '\t', quote = F, row.names = T, col.names = T, nThread = 4, compress = 'auto')
## Write tab deliminated cell annotations
fwrite(as.data.frame(malig_status), file = "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv2/cellAnnotations.txt.gz", sep = "\t", quote = FALSE, row.names = T, col.names = FALSE, compress = 'auto')


set.seed(1234)
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv2/counts.matrix.gz",
                                     annotations_file = "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv2/cellAnnotations.txt.gz",
                                     delim = "\t",
                                     gene_order_file = "/n/scratch/users/c/cao385/Ependymoma/10x/gene_order.txt",
                                     ref_group_names = 'normal', min_max_counts_per_cell = c(-Inf,+Inf))

# perform infercnv operations to reveal cnv signal
set.seed(1234)
out_dir <- "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv2/out_dir2"
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir = out_dir,  # dir is auto-created for storing outputs
                              cluster_by_groups = T,   # cluster
                              denoise = T,
                              HMM = T,
                              plot_steps = F,
                              num_threads = 12,
                              analysis_mode = 'samples',
                              tumor_subcluster_partition_method = 'leiden',
                              tumor_subcluster_pval = 0.05,
                              leiden_resolution = 0.01,
                              png_res = 600,
                              write_expr_matrix = T)


## Normalized gene expressions of predicted tumor cells
observations <- fread(file.path(out_dir, 'infercnv.observations.txt'), data.table = F) %>%
  column_to_rownames('V1')

## Normalized gene expressions of predicted normal cells
references <- fread(file.path(out_dir, 'infercnv.references.txt'), data.table = F) %>%
  column_to_rownames('V1')

## Concatnate
all_cnv_values <- cbind(observations, references)

## Labels for obs vs ref
obs_ref <- c(rep("obs", dim(observations)[2]), rep("ref", dim(references)[2]))
names(obs_ref) <- c(colnames(observations), colnames(references))

## CRITICAL: sort the order of cells in cnv_values and obs_ref to the same as cm
all_cnv_values <- all_cnv_values[,names(malig_status)[names(malig_status) %in% colnames(all_cnv_values)]]
obs_ref <- obs_ref[names(malig_status)[names(malig_status) %in% colnames(all_cnv_values)]]
all_cnv_values <-  all_cnv_values - 1

qsave(all_cnv_values, "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv2/all_cnv_values.qs")
qsave(obs_ref, "/n/scratch/users/c/cao385/Ependymoma/10x/infercnv2/obs_ref.qs")

counts <- cbind(cm_raw, ext_ctrl_cm)
data <- RunFullSeurat_v5(cm = counts, metadata = malig_status, doBatch = T, var2batch = 'sample', batchMethod = 'harmony', verbose = T, project = 'Daeun')
data <- AddMetaData(data, metadata = glue('cnv_{cnv_cutree_each_sample$`2LP`}'), col.name = 'cnv_clusters')
data$tumor <- ifelse(data$seurat_clusters %in% c(3,6,7,10,11,15,17,18), 'normal', 'tumor')
data$tdTomato <- ifelse(FetchData(data, 'tdTomato-1') %>% pull('tdTomato-1') > 0, 'yes', 'no')

p1 <- DimPlot(data, reduction = 'umap.unintegrated', group.by = 'sample', cols = colors_to_use, pt.size = 1) + ggtitle('Unintegrated') + seurat_theme()
p2 <- DimPlot(data, reduction = 'umap.unintegrated', group.by = 'cnv_clusters', cols = colors_to_use, pt.size = 1) + ggtitle('Unintegrated') + seurat_theme()
p3 <- DimPlot(data, reduction = 'umap.unintegrated', group.by = 'tumor', cols = colors_to_use, pt.size = 1) + ggtitle('Unintegrated') + seurat_theme()
p4 <- DimPlot(data, reduction = 'umap.unintegrated', group.by = 'tdTomato', cols = colors_to_use, pt.size = 0.5, order = T) + ggtitle('Unintegrated') + seurat_theme()
p1 + p2 + p3 + p4

qsave(data, '/n/scratch/users/c/cao385/Ependymoma/10x/output/2LP-post-inferCNV.qs')