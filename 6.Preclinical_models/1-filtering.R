rm(list = ls())

library(data.table)
library(tidyverse)
library(readxl)
library(crayon)


base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/6 - scRNAseq models"
resources_dir <- file.path(base_dir, 'scripts/resources')

source(file.path(resources_dir, "single_cell_preprocessing_helper_functions.R"))

analysis_dir <- file.path(base_dir, "analysis")
if (!dir.exists(analysis_dir)){dir.create(analysis_dir, recursive = T)}

## Load house keeping genes
hkgenes <- read.csv(file.path(resources_dir, 'tirosh_house_keeping.txt'), skip = 1) %>%
  pull(1)

metadata <- read_excel(file.path(base_dir, 'metadata.xlsx')) 

FileName <- unique(metadata$FileName)
Project <- unique(metadata$Project)
Type <- unique(metadata$Type)

## Create directories if needed
qc_dir <- file.path(base_dir, "analysis", 'qc')
qc_plots <- file.path(base_dir, "analysis", 'qc/plots')
qc_data <- file.path(base_dir, "analysis", 'qc/data')
if (!dir.exists(qc_dir)){dir.create(qc_dir, recursive = T)}
if (!dir.exists(qc_plots)){dir.create(qc_plots, recursive = T)}
if (!dir.exists(file.path(qc_plots, 'individual'))){dir.create(file.path(qc_plots, 'individual'), recursive = T)}
if (!dir.exists(file.path(qc_data, 'individual'))){dir.create(file.path(qc_data, 'individual'), recursive = T)}


for (i in seq_along(metadata$FileName)) {
  if (!file.exists(file.path(qc_data, 'individual', sprintf("%s-cm_list.rds", metadata$FileName[i]))) &
      !file.exists(file.path(qc_data, 'individual', sprintf("%s-qc_raw_data.rds", metadata$FileName[i]))) &
      !file.exists(file.path(qc_data, 'individual', sprintf("%s-qc_filtered_data.rds", metadata$FileName[i]))) &
      !file.exists(file.path(qc_data, 'individual', sprintf("%s-samples.rds", metadata$FileName[i]))) &
      !file.exists(file.path(qc_data, 'individual', sprintf("%s-genes_to_keep.rds", metadata$FileName[i])))) {
    cat(cyan(sprintf("Sample %s\n", metadata$FileName[i])))
    
    ## Load count matrix
    cm <- fread(list.files(file.path(base_dir, "data/models/counts"), pattern = paste0('_', metadata$FileName[i], '.csv'), full.names = T), data.table = F) %>%
      column_to_rownames("gene_id") %>%
      rename_with(~ gsub('-', '\\.', .x))
    
    ## Load alignment qc
    qc_align <- read.csv(list.files(file.path(base_dir, "data/models/qc"), pattern = paste0('_', metadata$FileName[i], '.csv'), full.names = T)) %>%
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
    ggsave(sprintf("%s-QC_histogram_nGene_nAligned.pdf", metadata$FileName[i]), path = file.path(qc_plots, 'individual'), width = 12, height=8)
    
    
    ## Determine cutoff
    sample_plate <- qc_raw_data$sample_plate
    
    ## Cutoff per plate (mean - 2*sd)
    tmp <- log10(qc_raw_data$gene+1)
    cutoff_per_plate <- sapply(unique(sample_plate), function(x){
      tmp2 <- tmp[sample_plate == x]
      return(mean(tmp2) - 2*sd(tmp2))
    })
    qc_raw_data$log_gene <- tmp
    
    if (metadata$Type[i] == 'frozen') {
      nGene_cutoff <- 1000
      align_rate_cutoff <- 0.3
      
      ## Filter by nGene
      pass_qc_1 <- qc_raw_data$gene >= nGene_cutoff
      
      ## Filter by plate-based cutoff
      pass_qc_2 <- sapply(rownames(qc_raw_data), function(x) qc_raw_data[x, "log_gene"] >= cutoff_per_plate[qc_raw_data[x, "sample_plate"]])
      ## Filter through single cutoff and plate-based cutoff
      qc_raw_data$pass_qc <- pass_qc_1 &
        pass_qc_2 &
        qc_raw_data$align_rate >= align_rate_cutoff
      
      ## Scatter plot of # of genes and align rate, colored by QC
      ggplot(data = qc_raw_data, aes(x=gene, y=align_rate, color=pass_qc)) + geom_point() +
        geom_hline(yintercept = align_rate_cutoff, color="red", linetype="dashed", size=1) +
        geom_vline(xintercept = nGene_cutoff, color="red", linetype="dashed", size=1) +
        xlab("Total number of genes\n") + ylab("Align Rate") +
        theme_classic() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16),
                                axis.title = element_text(face="bold", size=16),
                                axis.text = element_text(face="bold", size=14),
                                legend.title = element_text(face="bold", size=16),
                                legend.text = element_text(size=14))
      ggsave(sprintf("%s-genes_vs_alignrate.pdf", metadata$FileName[i]), path = file.path(qc_plots, 'individual'), width = 8, height=6)
    } else if (metadata$Type[i] == 'fresh') {
      nGene_cutoff <- 2000
      hk_mean_expr_cutoff <- 2
      
      qc_raw_data$pass_qc = ifelse(qc_raw_data$gene >= nGene_cutoff & qc_raw_data$hk >= hk_mean_expr_cutoff, T, F)
      
      ggplot(data=qc_raw_data, aes(x=gene, y=hk, color=pass_qc)) + geom_point() +
        geom_hline(yintercept = hk_mean_expr_cutoff, color="red", linetype="dashed", size=1) +
        geom_vline(xintercept = nGene_cutoff, color="red", linetype="dashed", size=1) +
        xlab("Total number of genes\n") + ylab("Mean expression of house-keeping genes\n") +
        theme_classic() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16),
                                axis.title = element_text(face="bold", size=16),
                                axis.text = element_text(face="bold", size=14),
                                legend.title = element_text(face="bold", size=16),
                                legend.text = element_text(size=14))
      ggsave(sprintf("%s-genes_vs_hk.pdf", metadata$FileName[i]), path = file.path(qc_plots, 'individual'), width = 8, height=6)
    }
    
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
    saveRDS(cm_list, file = file.path(qc_data, 'individual', sprintf("%s-cm_list.rds", metadata$FileName[i])))
    saveRDS(qc_raw_data, file = file.path(qc_data, 'individual', sprintf("%s-qc_raw_data.rds", metadata$FileName[i])))
    saveRDS(qc_filtered_data, file = file.path(qc_data, 'individual', sprintf("%s-qc_filtered_data.rds", metadata$FileName[i])))
    saveRDS(samples, file = file.path(qc_data, 'individual', sprintf("%s-samples.rds", metadata$FileName[i])))
    saveRDS(genes_to_keep, file = file.path(qc_data, 'individual', sprintf("%s-genes_to_keep.rds", metadata$FileName[i])))
  }
}


cm <- do.call(cbind, lapply(list.files(file.path(qc_data, 'individual'), pattern = '-cm_list.rds', full.names = T), function(x) {
  tmp <- readRDS(x)
  tmp$raw_data
}))

qc_raw_data <- do.call(rbind, lapply(list.files(file.path(qc_data, 'individual'), pattern = '-qc_raw_data.rds', full.names = T), function(x) {
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
saveRDS(cm_list, file = file.path(qc_data, 'individual', "combined-cm_list.rds"))
saveRDS(qc_raw_data, file = file.path(qc_data, 'individual', "combined-qc_raw_data.rds"))
saveRDS(qc_filtered_data, file = file.path(qc_data, 'individual', "combined-qc_filtered_data.rds"))
saveRDS(samples, file = file.path(qc_data, 'individual', "combined-samples.rds"))
saveRDS(genes_to_keep, file = file.path(qc_data, 'individual', "combined-genes_to_keep.rds"))


## Barplot of good/poor quality cells for each sample
## # of cells that pass QC in each sample
crosstab <- table(qc_raw_data$sample, qc_raw_data$pass_qc)
crosstab <- data.frame(crosstab)
colnames(crosstab) <- c("Sample", "Good_quality", "Frequency")
ggplot(crosstab, aes_string(x="Sample", y="Frequency", fill="Good_quality")) +
  geom_bar(stat="identity", color="black") + ggtitle("QC of each sample\n") +
  xlab("") + ylab("Frequency\n") +
  labs(fill="Good quality") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5, face="bold", size=32),
                          axis.title = element_text(face="bold", size=16),
                          axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1),
                          axis.text.y = element_text(face = "bold", size = 14),
                          legend.title = element_text(size=14, face="bold"),
                          legend.text = element_text(size=18),
                          plot.margin = margin(10,10,10,60))
figWidth <- ifelse(length(unique(samples)) > 20, 15, 8)
ggsave("num_cells_pass_qc.pdf.pdf", path = qc_plots, width = figWidth, height = 6)

## % of cells that pass QC in each sample
crosstab <- table(qc_raw_data$sample, qc_raw_data$pass_qc)
crosstab <- crosstab/rowSums(crosstab)
crosstab <- data.frame(crosstab)
colnames(crosstab) <- c("Sample", "Good_quality", "Percentage")
ggplot(crosstab, aes_string(x = "Sample", y = "Percentage", fill = "Good_quality")) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", size = 1) +
  ggtitle("Percentage of good quality cells in each sample") +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
figHeight <- ifelse(length(unique(samples)) > 20, 10, 6)
ggsave("percentage_cells_pass_qc.pdf", path = qc_plots, width = 8, height = figHeight)

