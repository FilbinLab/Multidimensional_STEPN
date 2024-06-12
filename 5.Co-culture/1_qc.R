rm(list = ls())

library(data.table)
library(tidyverse)
library(glue)
library(ggpubr)

base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/9 - scRNAseq coculture"
resource_dir <- file.path(base_dir, 'scripts/resources')

source("~/Dropbox (Partners HealthCare)/Shared/Scripts/single_cell_preprocessing_helper_functions.R")

analysis_dir <- file.path(base_dir, "analysis/qc")
if (!dir.exists(analysis_dir)){dir.create(analysis_dir, recursive = T)}

## Load house keeping genes
hkgenes <- read.csv(file.path(resource_dir, 'tirosh_house_keeping.txt'), skip = 1) %>%
  pull(1)

sample_names <- list.files(file.path(base_dir, 'data/coculture'), pattern = 'cm_counts_*')
sample_names <- unique(gsub("\\-.*", "", gsub('cm_counts_|.csv.gz', '', sample_names)))

for (sample in sample_names) {
  # sample <- sample_names[1]
  
  ## Load count matrix
  cm <- lapply(list.files(file.path(base_dir, "data/coculture"), pattern = glue('cm_counts_{sample}'), full.names = T), function(x) {
    fread(x, data.table = F) %>% 
      column_to_rownames("V1") %>% 
      rename_with(~ gsub('-', '\\.', .x))
  })
  cm <- do.call(cbind, cm)
  
  ## Load alignment qc
  qc_align <- lapply(list.files(file.path(base_dir, "data/coculture"), pattern = glue('qc_{sample}'), full.names = T), function(x) {
    fread(x, data.table = F) %>% 
      column_to_rownames('name')
  })
  qc_align <- do.call(rbind, qc_align)
  
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
  ggsave(glue("{sample}-QC_histogram_nGene_nAligned.pdf"), path = analysis_dir, width = 12, height=8)
  
  
  ## Determine cutoff
  sample_plate <- qc_raw_data$sample_plate
  
  ## Cutoff per plate (mean - 2*sd)
  tmp <- log10(qc_raw_data$gene+1)
  cutoff_per_plate <- sapply(unique(sample_plate), function(x){
    tmp2 <- tmp[sample_plate == x]
    return(mean(tmp2) - 2*sd(tmp2))
  })
  qc_raw_data$log_gene <- tmp
  
  
  ## Fresh samples
  nGene_cutoff <- 1000
  qc_raw_data$pass_qc_1 <- qc_raw_data$gene >= nGene_cutoff
  qc_raw_data$pass_qc_2 <- sapply(rownames(qc_raw_data), 
                                  function(x) qc_raw_data[x, "log_gene"] >= cutoff_per_plate[qc_raw_data[x, "sample_plate"]]) 
  
  align_rate_cutoff = 0.1
  qc_raw_data$pass_qc = qc_raw_data$pass_qc_1 & 
    qc_raw_data$pass_qc_2 &
    qc_raw_data$align_rate >= align_rate_cutoff
  
  
  #hk_mean_expr_cutoff <- 2
  #qc_raw_data$pass_qc = ifelse(qc_raw_data$gene >= nGene_cutoff & qc_raw_data$hk >= hk_mean_expr_cutoff, T, F)
  #ggplot(data=qc_raw_data, aes(x=gene, y=hk, color=pass_qc)) + geom_point() +
  #  geom_hline(yintercept = hk_mean_expr_cutoff, color="red", linetype="dashed", size=1) +
  #  geom_vline(xintercept = nGene_cutoff, color="red", linetype="dashed", size=1) +
  #  xlab("Total number of genes\n") + ylab("Mean expression of house-keeping genes\n") +
  #  theme_classic() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16),
  #                          axis.title = element_text(face="bold", size=16),
  #                          axis.text = element_text(face="bold", size=14),
  #                          legend.title = element_text(face="bold", size=16),
  #                          legend.text = element_text(size=14))
  #ggsave(glue("{sample}-genes_vs_hk.png"), path = analysis_dir, width = 8, height=6)
  
  ggplot(data = qc_raw_data, aes(x = gene, y = align_rate, color = pass_qc)) + geom_point() +
    geom_hline(yintercept = align_rate_cutoff, color="red", linetype="dashed", size=1) +
    geom_vline(xintercept = nGene_cutoff, color="red", linetype="dashed", size=1) +
    xlab("Total number of genes\n") + ylab("Alignment rate\n") +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16),
                            axis.title = element_text(face="bold", size=16),
                            axis.text = element_text(face="bold", size=14),
                            legend.title = element_text(face="bold", size=16),
                            legend.text = element_text(size=14))
  ggsave(glue("{sample}-genes_vs_alignrate.png"), path = analysis_dir, width = 8, height=6)
  
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
  saveRDS(cm_list, file = glue("{analysis_dir}/{sample}-cm_list.rds"))
  saveRDS(qc_raw_data, file = glue("{analysis_dir}/{sample}-qc_raw_data.rds"))
  saveRDS(qc_filtered_data, file = glue("{analysis_dir}/{sample}-qc_filtered_data.rds"))
  saveRDS(samples, file = glue("{analysis_dir}/{sample}-samples.rds"))
  saveRDS(genes_to_keep, file = glue("{analysis_dir}/{sample}-genes_to_keep.rds"))
}


cm <- do.call(cbind, lapply(list.files(analysis_dir, pattern = '-cm_list.rds', full.names = T), function(x) {
  tmp <- readRDS(x)
  tmp$raw_data
}))

qc_raw_data <- do.call(rbind, lapply(list.files(analysis_dir, pattern = '-qc_raw_data.rds', full.names = T), function(x) {
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


## % of cells that pass QC in each sample/plate
crosstab <- table(qc_raw_data$sample_plate, qc_raw_data$pass_qc)
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
ggsave("percentage_cells_pass_qc_plate.pdf", path = analysis_dir, width = 8, height = 6)