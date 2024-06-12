library(tidyverse)
library(data.table)
library(glue)


################################################################################
################################################################################
## Organize counts and qc files
################################################################################
################################################################################

## BT1804
smp <- 'BT1804'
cm <- lapply(list.files('~/Downloads/Ependymoma/', pattern = glue('cm_counts_{smp}'), full.names = T), function(x) {
  fread(x) %>% 
    column_to_rownames('V1')
})
cm <- do.call(cbind, cm)
colnames(cm) <- gsub('Tumor', '', colnames(cm))
colnames(cm) <- gsub('-', '.', colnames(cm))

qc <- lapply(list.files('~/Downloads/Ependymoma/', pattern = glue('qc_{smp}'), full.names = T), function(x) {
  fread(x) %>% 
    column_to_rownames('name')
})
qc <- do.call(rbind, qc)
rownames(qc) <- gsub('Tumor', '', rownames(qc))
qc$sample <- gsub('Tumor', '', qc$sample)
qc <- qc[colnames(cm),]
qc <- qc %>% rownames_to_column('name')

fwrite(cm, glue('~/Downloads/Ependymoma/files/cm_{smp}.csv'), sep = ',', quote = F, row.names = T)
fwrite(qc, glue('~/Downloads/Ependymoma/files/qc_{smp}.csv'), sep = ',', quote = F, row.names = F)


## BT2126
smp <- 'BT2126'
cm <- lapply(list.files('~/Downloads/Ependymoma/', pattern = glue('cm_counts_{smp}'), full.names = T), function(x) {
  fread(x) %>% 
    column_to_rownames('V1')
})
cm <- do.call(cbind, cm)
colnames(cm) <- gsub('-', '.', colnames(cm))

qc <- lapply(list.files('~/Downloads/Ependymoma/', pattern = glue('qc_{smp}'), full.names = T), function(x) {
  fread(x) %>% 
    column_to_rownames('name')
})
qc <- do.call(rbind, qc)
qc <- qc[colnames(cm),]
qc <- qc %>% rownames_to_column('name')

fwrite(cm, glue('~/Downloads/Ependymoma/files/cm_{smp}.csv'), sep = ',', quote = F, row.names = T)
fwrite(qc, glue('~/Downloads/Ependymoma/files/qc_{smp}.csv'), sep = ',', quote = F, row.names = F)


## BT2169
smp <- 'BT2169'
cm <- lapply(list.files('~/Downloads/Ependymoma/', pattern = glue('cm_counts_{smp}'), full.names = T), function(x) {
  fread(x) %>% 
    column_to_rownames('V1')
})
cm <- do.call(cbind, cm)
colnames(cm) <- gsub('-', '.', colnames(cm))

qc <- lapply(list.files('~/Downloads/Ependymoma/', pattern = glue('qc_{smp}'), full.names = T), function(x) {
  fread(x) %>% 
    column_to_rownames('name')
})
qc <- do.call(rbind, qc)
qc <- qc[colnames(cm),]
qc <- qc %>% rownames_to_column('name')

fwrite(cm, glue('~/Downloads/Ependymoma/files/cm_{smp}.csv'), sep = ',', quote = F, row.names = T)
fwrite(qc, glue('~/Downloads/Ependymoma/files/qc_{smp}.csv'), sep = ',', quote = F, row.names = F)


## BT1717
smp <- 'BT1717'
cm <- lapply(list.files('~/Downloads/Ependymoma/', pattern = glue('cm_counts_{smp}'), full.names = T), function(x) {
  fread(x) %>% 
    column_to_rownames('V1')
})
cm <- do.call(cbind, cm)
colnames(cm) <- gsub('-', '.', colnames(cm))

qc <- lapply(list.files('~/Downloads/Ependymoma/', pattern = glue('qc_{smp}'), full.names = T), function(x) {
  fread(x) %>% 
    column_to_rownames('name')
})
qc <- do.call(rbind, qc)
qc$sample <- glue('{qc$sample}-{qc$plate}')
qc <- qc[colnames(cm),]
qc <- qc %>% rownames_to_column('name')

fwrite(cm, glue('~/Downloads/Ependymoma/files/cm_{smp}.csv'), sep = ',', quote = F, row.names = T)
fwrite(qc, glue('~/Downloads/Ependymoma/files/qc_{smp}.csv'), sep = ',', quote = F, row.names = F)


## MUV56
smp <- 'MUV56'
cm <- lapply(list.files('~/Downloads/Ependymoma/', pattern = glue('cm_counts_{smp}'), full.names = T), function(x) {
  fread(x) %>% 
    column_to_rownames('V1')
})
cm <- do.call(cbind, cm)
colnames(cm) <- gsub('-', '.', colnames(cm))

qc <- lapply(list.files('~/Downloads/Ependymoma/', pattern = glue('qc_{smp}'), full.names = T), function(x) {
  fread(x) %>% 
    column_to_rownames('name')
})
qc <- do.call(rbind, qc)
qc$sample <- glue('{qc$sample}-{qc$plate}')
qc <- qc[colnames(cm),]
qc <- qc %>% rownames_to_column('name')

fwrite(cm, glue('~/Downloads/Ependymoma/files/cm_{smp}.csv'), sep = ',', quote = F, row.names = T)
fwrite(qc, glue('~/Downloads/Ependymoma/files/qc_{smp}.csv'), sep = ',', quote = F, row.names = F)


## BT775 (Peds4)
smp <- 'Peds4'
cm <- lapply(list.files('~/Downloads/Ependymoma/', pattern = glue('cm_counts_{smp}'), full.names = T), function(x) {
  fread(x) %>% 
    column_to_rownames('V1')
})
cm <- do.call(cbind, cm)
colnames(cm) <- gsub('Peds4', 'BT775', colnames(cm))
colnames(cm) <- gsub('-', '.', colnames(cm))

qc <- lapply(list.files('~/Downloads/Ependymoma/', pattern = glue('qc_{smp}'), full.names = T), function(x) {
  fread(x) %>% 
    column_to_rownames('name')
})
qc <- do.call(rbind, qc)
rownames(qc) <- gsub('Peds4', 'BT775', rownames(qc))
qc$sample <- gsub('Peds4', 'BT775', qc$sample)
qc$sample <- glue('{qc$sample}-{qc$plate}')
qc <- qc[colnames(cm),]
qc <- qc %>% rownames_to_column('name')

fwrite(cm, glue('~/Downloads/Ependymoma/files/cm_BT775.csv'), sep = ',', quote = F, row.names = T)
fwrite(qc, glue('~/Downloads/Ependymoma/files/qc_BT775.csv'), sep = ',', quote = F, row.names = F)



################################################################################
################################################################################
## QC
################################################################################
################################################################################
base_dir <- '/n/scratch/users/c/cao385/Ependymoma/dj83'

source(file.path(base_dir, 'resources/single_cell_preprocessing_helper_functions.R'))

analysis_dir <- file.path(base_dir, "analysis/ZFTA_fresh/qc")
if (!dir.exists(analysis_dir)){dir.create(analysis_dir, recursive = T)}

## Load house keeping genes
hkgenes <- read.csv(file.path(base_dir, 'resources/tirosh_house_keeping.txt'), skip = 1) %>%
  pull(1)

sample_names <- list.files(file.path(base_dir, 'data/ZFTA_fresh/counts'), pattern = 'cm_*')
sample_names <- gsub('cm_|.csv.gz', '', sample_names)

for (sample in sample_names) {
  # sample <- sample_names[1]
  
  dir.create(file.path(analysis_dir, sample))
  
  ## Load count matrix
  cm <- fread(list.files(file.path(base_dir, "data/ZFTA_fresh/counts"), pattern = glue('cm_{sample}.csv'), full.names = T), data.table = F) %>% 
    column_to_rownames("V1")
  
  ## Load alignment qc
  qc_align <- read.csv(list.files(file.path(base_dir, "data/ZFTA_fresh/counts"), pattern = glue('qc_{sample}.csv'), full.names = T)) %>%
    column_to_rownames('name')
  
  ## Match cell names if necessary and check cell names match
  ## Subset qc_align to keep only those within cm
  qc_align <- qc_align[colnames(cm), ]
  
  ## Make a vector of sample names
  samples <- structure(qc_align[,"sample"], names = rownames(qc_align))
  
  ## Compute QC metrics on all cells
  qc_raw_data <- calculateQcMetrics(cm, hkgenes)
  ## Add total aligned reads
  qc_raw_data$nAligned <- qc_align$nAligned
  ## Compute alignment rate
  qc_raw_data$align_rate <- qc_align$nAligned/qc_align$nTotal
  ## Add sample-plate labels
  qc_raw_data$sample_plate <- paste(qc_raw_data$sample, qc_raw_data$plate, sep = "_")
  
  
  ## Histogram of total number of genes
  p1 <- ggplot(qc_raw_data, aes(x=gene)) + geom_histogram(bins=50) +
    #scale_x_continuous(breaks=seq(0,15000,1000)) +
    geom_vline(xintercept = 1E3, color="red", linetype="dashed", size=1) +
    ggtitle("Total number of gene") +
    xlab("Total number of genes") + ylab("Number of cells") +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=14),
                            axis.title = element_text(face="bold", size=14),
                            axis.text = element_text(face="bold", size=14),
                            legend.title = element_text(face="bold", size=16),
                            legend.text = element_text(size=14))
  
  ## Histogram of total number of genes (log10 scale)
  p2 <- ggplot(qc_raw_data, aes(x=gene)) + geom_histogram(bins=50) +
    scale_x_log10(breaks = c(10, 100, 1E3, 1E4)) +
    geom_vline(xintercept = 1E3, color="red", linetype="dashed", size=1) +
    ggtitle("Total number of gene (log10)") +
    xlab("Log10 total number of genes") + ylab("Number of cells") +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=14),
                            axis.title = element_text(face="bold", size=14),
                            axis.text = element_text(face="bold", size=14),
                            legend.title = element_text(face="bold", size=16),
                            legend.text = element_text(size=14))
  
  ## Histogram of total number of aligned reads (log10 scale)
  p3 <- ggplot(qc_raw_data, aes(x=nAligned)) + geom_histogram(bins=50) +
    scale_x_log10(breaks = c(10, 100, 1E3, 1E4, 1E5, 1E6)) +
    geom_vline(xintercept = 1E4, color="red", linetype="dashed", size=1) +
    ggtitle("Total number of aligned reads (log10)") +
    xlab("Log10 total number of aligned reads") + ylab("Number of cells") +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=14),
                            axis.title = element_text(face="bold", size=14),
                            axis.text = element_text(face="bold", size=14),
                            legend.title = element_text(face="bold", size=16),
                            legend.text = element_text(size=14))
  
  ## Histogram of alignment rate
  p4 <- ggplot(qc_raw_data, aes(x=align_rate)) + geom_histogram(bins=50) +
    ggtitle("Alignment rate") +
    geom_vline(xintercept = 0.1, color="red", linetype="dashed", size=1) +
    scale_x_continuous(breaks = seq(0,0.8,0.1)) +
    xlab("Alignment rate") + ylab("Number of cells") +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=14),
                            axis.title = element_text(face="bold", size=14),
                            axis.text = element_text(face="bold", size=14),
                            legend.title = element_text(face="bold", size=16),
                            legend.text = element_text(size=14))
  (p1+p2)/(p3+p4)
  ggsave(glue("{sample}-QC_histogram_nGene_nAligned.pdf"), path = file.path(analysis_dir, sample), width = 12, height = 8)
  
  
  ## Determine cutoff
  sample_plate <- qc_raw_data$sample_plate
  
  ## Cutoff per plate (mean - 2*sd)
  tmp <- log10(qc_raw_data$gene+1)
  cutoff_per_plate <- sapply(unique(sample_plate), function(x){
    tmp2 <- tmp[sample_plate == x]
    return(mean(tmp2) - 2*sd(tmp2))
  })
  qc_raw_data$log_gene <- tmp
  
  ## Try a single and/or plate-specific cutoff 
  nGene_cutoff <- 1000
  qc_raw_data$pass_qc_1 <- qc_raw_data$gene >= nGene_cutoff
  qc_raw_data$pass_qc_2 <- sapply(rownames(qc_raw_data), 
                                  function(x) qc_raw_data[x, "log_gene"] >= cutoff_per_plate[qc_raw_data[x, "sample_plate"]]) 
  
  ## Filter based on different cutoffs (nGene, nAligned, align_rate)
  align_rate_cutoff <- 0.1
  qc_raw_data$pass_qc <- qc_raw_data$pass_qc_1 & 
    qc_raw_data$pass_qc_2 &
    qc_raw_data$align_rate >= align_rate_cutoff
  
  
  ggplot(data = qc_raw_data, aes(x = gene, y = align_rate, color = pass_qc)) + geom_point() +
    geom_hline(yintercept = align_rate_cutoff, color="red", linetype="dashed", size=1) +
    geom_vline(xintercept = nGene_cutoff, color="red", linetype="dashed", size=1) +
    xlab("Total number of genes\n") + ylab("Alignment rate\n") +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16),
                            axis.title = element_text(face="bold", size=16),
                            axis.text = element_text(face="bold", size=14),
                            legend.title = element_text(face="bold", size=16),
                            legend.text = element_text(size=14))
  ggsave(glue("{sample}-genes_vs_hk.png"), path = file.path(analysis_dir, sample), width = 8, height=6)
  
  
  ## Preprocessing: filtering on cells
  ## Filter and normalize cm without filtering genes
  cm <- cm[,qc_raw_data$pass_qc]
  cm_norm <- log2(cm/10+1)
  cm_center <- cm_norm - rowMeans(cm_norm)
  cm_list <- list("raw_data" = cm, "norm_data" = cm_norm, "center_data" = cm_center)
  
  ## Identify sufficiently expressed genes
  gene_cutoff <- 4
  gene_cell_cutoff <- 20
  genes_to_keep <- rowSums(cm_list$raw_data >= 2^gene_cutoff) >= gene_cell_cutoff
  genes_to_keep <- rownames(cm_list$raw_data)[genes_to_keep]
  
  ## Subset sample
  samples <- samples[colnames(cm_list$raw_data)]
  ## Calculate and store qc metrics for filtered and unfiltered data
  qc_filtered_data <- calculateQcMetrics(cm_list$raw_data, hkgenes)
  
  rm(cm); rm(cm_norm); rm(cm_center); gc()
  
  ## Save preprocessed results
  saveRDS(cm_list, file = glue("{analysis_dir}/{sample}/{sample}-cm_list.rds"))
  saveRDS(qc_raw_data, file = glue("{analysis_dir}/{sample}/{sample}-qc_raw_data.rds"))
  saveRDS(qc_filtered_data, file = glue("{analysis_dir}/{sample}/{sample}-qc_filtered_data.rds"))
  saveRDS(samples, file = glue("{analysis_dir}/{sample}/{sample}-samples.rds"))
  saveRDS(genes_to_keep, file = glue("{analysis_dir}/{sample}/{sample}-genes_to_keep.rds"))
}

cm <- do.call(cbind, lapply(list.files(analysis_dir, pattern = '-cm_list.rds', recursive = T, full.names = T), function(x) {
  tmp <- readRDS(x)
  tmp$raw_data
}))

qc_raw_data <- do.call(rbind, lapply(list.files(analysis_dir, pattern = '-qc_raw_data.rds', recursive = T, full.names = T), function(x) {
  readRDS(x)
}))

cm_norm <- log2(cm/10+1)
cm_center <- cm_norm - rowMeans(cm_norm)
cm_list <- list("raw_data" = cm, "norm_data" = cm_norm, "center_data" = cm_center)

samples <- structure(unlist(lapply(strsplit(colnames(cm), '\\.'), `[[`, 1)), names = colnames(cm))

## Identify sufficiently expressed genes
gene_cutoff <- 4
gene_cell_cutoff <- 20
genes_to_keep <- rowSums(cm_list$raw_data >= 2^gene_cutoff) >= gene_cell_cutoff
genes_to_keep <- rownames(cm_list$raw_data)[genes_to_keep]

## Subset sample
samples <- samples[colnames(cm_list$raw_data)]
## Calculate and store qc metrics for filtered and unfiltered data
qc_filtered_data <- calculateQcMetrics(cm_list$raw_data, hkgenes)

rm(cm); rm(cm_norm); rm(cm_center); gc()

## Save preprocessed results
saveRDS(cm_list, file = file.path(analysis_dir, "cm_list.rds"))
saveRDS(qc_raw_data, file = file.path(analysis_dir, "qc_raw_data.rds"))
saveRDS(qc_filtered_data, file = file.path(analysis_dir, "qc_filtered_data.rds"))
saveRDS(samples, file = file.path(analysis_dir, "samples.rds"))
saveRDS(genes_to_keep, file = file.path(analysis_dir, "genes_to_keep.rds"))


## Barplot of good/poor quality cells for each sample
## # of cells that pass QC in each sample
crosstab <- table(qc_raw_data$sample, qc_raw_data$pass_qc)
crosstab <- data.frame(crosstab)
colnames(crosstab) <- c("Sample", "Good_quality", "Frequency")
crosstab$Sample = factor(crosstab$Sample, levels = crosstab%>% filter(Good_quality == "TRUE") %>% arrange(Frequency)%>% pull(Sample))
ggplot(crosstab, aes_string(x="Sample", y="Frequency", fill="Good_quality")) +
  geom_bar(stat="identity", color="black") + ggtitle("QC of each sample\n") +
  geom_hline(yintercept = 100, color = "black", linetype = "dashed", size = 1) +
  xlab("") + ylab("Frequency\n") +
  labs(fill="Good quality") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=32),
                          axis.title = element_text(face="bold", size=16),
                          axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1),
                          axis.text.y = element_text(face = "bold", size = 14),
                          legend.title = element_text(size=14, face="bold"),
                          legend.text = element_text(size=18),
                          plot.margin = margin(10,10,10,60))
ggsave("num_cells_pass_qc.pdf.pdf", path = analysis_dir, width = 8, height = 6)

## % of cells that pass QC in each sample
crosstab <- table(qc_raw_data$sample, qc_raw_data$pass_qc)
crosstab <- crosstab/rowSums(crosstab)
crosstab <- data.frame(crosstab)
colnames(crosstab) <- c("Sample", "Good_quality", "Percentage")
crosstab$Sample = factor(crosstab$Sample, levels = crosstab%>% filter(Good_quality == "TRUE") %>% arrange(Percentage)%>% pull(Sample))
ggplot(crosstab, aes(x = Sample, y = Percentage, fill = Good_quality)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", size = 1) +
  ggtitle("Percentage of good quality cells in each sample") +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
ggsave("percentage_cells_pass_qc.pdf", path = analysis_dir, width = 8, height = 6)





library(Seurat)
library(qs)
library(stringr)

source('~/Projects/General-Codes/Resources/single_cell_preprocessing_helper_functions.R')
source('~/Projects/General-Codes/Resources/NMF_helper_function.R')
source('~/Projects/General-Codes/Resources/Plotting_helper_functions.R')

base_dir <- '/n/scratch/users/c/cao385/Ependymoma/dj83'
analysis_dir <- file.path(base_dir, "analysis/ZFTA_fresh/qc")

cm_list <- readRDS(file.path(analysis_dir, "cm_list.rds"))
samples <- readRDS(file.path(analysis_dir, "samples.rds"))

#data <- RunFullSeurat_v5(cm = cm_list$raw_data, metadata = samples, doBatch = T, var2batch = 'sample', batchMethod = 'harmony', dims = 0.8, project = 'EPN_Fresh', norm.type = 'RNA', verbose = T)
#qsave(data, '/n/scratch/users/c/cao385/Ependymoma/dj83/New/seurat_fresh.qs')
#DimPlot(data, group.by = 'sample', cols = clusterExperiment::bigPalette, reduction = 'tsne.unintegrated')
#DimPlot(data, group.by = 'sample', cols = clusterExperiment::bigPalette, reduction = 'tsne.harmony') + seurat_theme()
#markers_frozen <- qread('/n/scratch/users/c/cao385/Ependymoma/dj83/data/ZFTA_frozen/gene_signatures/ZFTA_frozen_DE_list_signature1_mergedNPCs_top200.qs')

data <- readRDS('/n/scratch/users/c/cao385/Ependymoma/seurat_obj_frozen_NPC_merged_clean.rds')
Idents(data) <- data$Metaprogram
pb_degs <- FindAllMarkers(data, logfc.threshold = 1, only.pos = T) %>% filter(p_val_adj < 0.05)
pb_degs <- pb_degs[!str_detect(pb_degs$gene, "RP11"), ]
pb_degs <- pb_degs %>% group_by(cluster) %>% arrange(-avg_log2FC, .by_group = T) %>% top_n(50, wt = avg_log2FC)
pb_degs <- lapply(split(pb_degs, f = pb_degs$cluster), function(x) x$gene)


data <- qread('/n/scratch/users/c/cao385/Ependymoma/dj83/New/seurat_fresh.qs')

## Normalize and center cm (still only keep malignant cells and all genes that pass the filter)
cm_norm <- as.matrix(log2(cm_list$raw_data/10+1))
cm_mean <- log2(Matrix::rowMeans(cm_list$raw_data)+1)
cm_center <- cm_norm - rowMeans(cm_norm)

nmf_score <- t(scoreNmfGenes(cm_center, cm_mean, pb_degs, verbose = F))

nmf_score_final_t <- data.frame(nmf_score)
nmf_score_final_t <- metagene_score_signature(nmf_score_final_t, 6, 1)
data <- AddMetaData(data, nmf_score_final_t)

DimPlot(data, group.by = 'signature_1', cols = clusterExperiment::bigPalette, reduction = 'tsne.harmony') + seurat_theme()


plotProportion(data$sample, data$signature_1,
               unique(data$sample), 
               unique(data$signature_1),
               c("Sample", "Annotation", "Proportion"),
               "Sample", "Proportion", "Annotation", 
               clusterExperiment::bigPalette,
               y_title = "Proportion", 
               x_text_size = 12, y_text_size = 12, axis_title_size = 14, legend_title_size = 14, legend_text_size = 12)

plotProportion(x = data$signature_1, y = data$sample, 
               x_order = unique(data$signature_1), 
               y_order = unique(data$sample), 
               col_names = c("Annotation", "Sample", "Proportion"),
               x_var = "Annotation", y_var = "Proportion", fill_var = "Sample", 
               colors = paletteer::paletteer_c("scico::roma", n = length(unique(data$sample))),
               y_title = "Proportion", 
               x_text_size = 12, y_text_size = 12, axis_title_size = 14, legend_title_size = 14, legend_text_size = 12)

qsave(data, '/n/scratch/users/c/cao385/Ependymoma/dj83/New/seurat_fresh_annotated.qs')