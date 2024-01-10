library(Seurat)
library(tidyverse)
library(enrichR)
library(patchwork)
library(glue)
library(Matrix)
library(DropletUtils)
library(SeuratWrappers)
library(dorothea)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)

#################################################################################
############################## Programs enrichment ##############################
#################################################################################
data <- readRDS('/n/scratch/users/c/cao385/Ependymoma/seurat_obj_frozen_NPC_merged_clean.rds')
Idents(data) <- data$Metaprogram

plotList <- list()
programs <- c("Cycling", "Neuroepithelial-like", "NPC-like", "Ependymal-like", "MES-like")
#programs <- levels(data$Metaprogram)
for (program in programs) {
  markers <- FindMarkers(data, ident.1 = program, logfc.threshold = 0.5, only.pos = F) %>% filter(p_val_adj < 0.05)
  #markers <- markers[!str_detect(rownames(markers), "RP11"), ]
  markers <- markers %>% arrange(-avg_log2FC) %>% rownames_to_column('gene')
  markers <- markers %>% pull(avg_log2FC, gene)
  
  gse <- gseGO(geneList = markers, 
               ont = "ALL", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "BH")
  gse@result <- gse@result %>% filter(NES > 0)
  plotList[[program]] <- dotplot(gse, x ='NES', showCategory = 20) + ggtitle(program) + theme(plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1))
}
pt <- wrap_plots(plotList$Cycling, plotList$`Neuroepithelial-like`, plotList$`NPC-like`, plotList$`Ependymal-like`, plotList$`MES-like`)
ggsave(plot = pt, '/n/scratch/users/c/cao385/enrich_GSEA.pdf', width = 30, height = 15)




plotList <- list()
programs <- levels(data$Metaprogram)
for (program in programs) {
  markers <- FindMarkers(data, ident.1 = program, logfc.threshold = 0.5, only.pos = F) %>% filter(p_val_adj < 0.05)
  markers <- markers %>% arrange(-avg_log2FC) %>% rownames_to_column('gene')
  markers <- markers %>% top_n(500, wt = avg_log2FC) %>% pull(gene)
  
  gse <- enrichGO(gene = markers, 
                  ont = "ALL", 
                  keyType = "SYMBOL", 
                  minGSSize = 3, 
                  maxGSSize = 800,
                  pvalueCutoff = 0.05, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")
  gse@result <- gse@result %>% arrange(p.adjust)
  plotList[[program]] <- dotplot(gse, showCategory = 20) + ggtitle(program) + theme(plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1))
}
pt <- wrap_plots(plotList$Cycling, plotList$`Neuroepithelial-like`, plotList$`Radial glia-like`, plotList$`NPC-like`, plotList$`Ependymal-like`, plotList$`MES-like`)
ggsave(plot = pt, '/n/scratch/users/c/cao385/enrich_GO.pdf', width = 30, height = 15)





#################################################################################
################################## cellphonedb ##################################
#################################################################################
show_condition <- function(code) {
  tryCatch(code,
           error = function(c) "error",
           warning = function(c) "warning",
           message = function(c) "message"
  )
}

source('/home/cao385/Projects/General-Codes/Resources/Plotting_helper_functions.R')


data <- readRDS('/n/scratch/users/c/cao385/Ependymoma/seurat_obj_frozen_NPC_merged_clean.rds')

metadata <- data.frame(Cell = Cells(data),
                       cell_type = data$Metaprogram)
write.table(metadata, '/n/scratch/users/c/cao385/Ependymoma/cpdb/data/Frozen-metadata.tsv', sep = '\t', quote = F, row.names = F)


count_raw <- LayerData(data, layer = 'counts', assay = 'SCT')[,Cells(data)]
write10xCounts('/n/scratch/users/c/cao385/Ependymoma/cpdb/data/Frozen-counts', count_raw, version = '2')
file.rename('/n/scratch/users/c/cao385/Ependymoma/cpdb/data/Frozen-counts/genes.tsv', '/n/scratch/users/c/cao385/Ependymoma/cpdb/data/Frozen-counts/features.tsv')


Idents(data) <- data$Metaprogram
markers <- RunPrestoAll(data, only.pos = T, logfc.threshold = 1) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T)
#top_n(50, wt = avg_log2FC)
markers <- data.frame(cell_type = markers$cluster,
                      gene = markers$gene,
                      logFC = markers$avg_log2FC,
                      P.Value = markers$p_val,
                      adj.P.Val = markers$p_val_adj)
write.table(markers, '/n/scratch/users/c/cao385/Ependymoma/cpdb/data/Frozen-degs.tsv', sep = '\t', quote = F, row.names = F)


data("entire_database", package = "dorothea")
TFs <- unique(entire_database$tf)[unique(entire_database$tf) %in% rownames(data)]

tf_results <- NULL
for (program in unique(as.character(data$Metaprogram))) {
  tf_results <- cbind(tf_results, head(names(sort(rowSums(LayerData(data, 'data')[TFs, WhichCells(data, idents = program)]), decreasing = T)), 20))
  message(program)
}
colnames(tf_results) <- unique(as.character(data$Metaprogram))

tf_results <- stack(as.list(as.data.frame(tf_results))) %>% rename_all(~c('TF', 'cell_type')) %>% relocate(cell_type, .before = TF)
write.table(tf_results, '/n/scratch/users/c/cao385/Ependymoma/cpdb/data/Frozen-TFs.tsv', sep = '\t', quote = F, row.names = F)


################################################################################
################################################################################
library(ktplots)
library(qs)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)

data <- readRDS('/n/scratch/users/c/cao385/Ependymoma/seurat_obj_frozen_NPC_merged_clean.rds')

sce <- SingleCellExperiment(list(counts = LayerData(data, 'counts', 'SCT')), colData = data@meta.data)

pvals <- read.delim("/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Frozen-method3_withScore/degs_analysis_relevant_interactions_Frozen.txt", check.names = FALSE)
pvals <- pvals[!duplicated(pvals$interacting_pair),]

means <- read.delim("/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Frozen-method3_withScore/degs_analysis_significant_means_Frozen.txt", check.names = FALSE)
means <- means[!duplicated(means$interacting_pair),]

decon <- read.delim("/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Frozen-method3_withScore/degs_analysis_deconvoluted_Frozen.txt", check.names = FALSE)

relevant_interactions_v5 <- read.delim("/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Frozen-method3_withScore/degs_analysis_relevant_interactions_Frozen.txt", check.names = FALSE)

pdf('/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Frozen/significant_interactions.pdf', width = 12, height = 10)
plot_cpdb_heatmap(pvals = relevant_interactions_v5, degs_analysis = TRUE, title = "Sum of significant interactions")
dev.off()


means_v5 <- read.delim("/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Frozen-method3_withScore/degs_analysis_means_Frozen.txt", check.names = FALSE)
means_v5 <- means_v5[!duplicated(means_v5$interacting_pair),]
means_v5[is.na(means_v5)] <- 0

interaction_scores_v5 <- read.delim("/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Frozen-method3_withScore/degs_analysis_interaction_scores_Frozen.txt", check.names = FALSE)
cellsign_v5 <- read.delim("/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Frozen-method3_withScore/degs_analysis_CellSign_active_interactions_Frozen.txt", check.names = FALSE)






#What the job is then of the Ependymal cells - do they secrete a ligand that has a receptor on the other cell populations?

plot_cpdb(
  scdata = data,
  cell_type1 = ".",
  cell_type2 = '.',
  means = means_v5 %>% select(1:13, starts_with("Ependymal-like")) %>% filter(gene_a != '') %>% filter(gene_b != ''),
  pvals = relevant_interactions_v5 %>% select(1:13, starts_with("Ependymal-like")) %>% filter(gene_a != '') %>% filter(gene_b != ''),
  celltype_key = "Metaprogram",
  max_size = 5,
  highlight_size = 0.5,
  degs_analysis = TRUE,
  standard_scale = TRUE,
  #genes = c(subset_pairs$gene_a, subset_pairs$gene_b)[c(subset_pairs$gene_a, subset_pairs$gene_b) != ''],
  interaction_scores = interaction_scores_v5 %>% select(1:13, starts_with("Ependymal-like")) %>% filter(gene_a != '') %>% filter(gene_b != ''),
  cellsign = cellsign_v5 %>% select(1:13, starts_with("Ependymal-like")) %>% filter(gene_a != '') %>% filter(gene_b != ''),
  filter_by_cellsign = T)
ggsave('/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Frozen/sign_interactions.pdf', width = 7, height = 17)


plot_cpdb(
  scdata = data,
  cell_type1 = ".",
  cell_type2 = '.',
  means = means_v5 %>% select(1:13, contains("|Ependymal-like")) %>% filter(gene_a != '') %>% filter(gene_b != ''),
  pvals = relevant_interactions_v5 %>% select(1:13, contains("|Ependymal-like")) %>% filter(gene_a != '') %>% filter(gene_b != ''),
  celltype_key = "Metaprogram",
  max_size = 5,
  highlight_size = 0.5,
  degs_analysis = TRUE,
  standard_scale = TRUE,
  #genes = c(subset_pairs$gene_a, subset_pairs$gene_b)[c(subset_pairs$gene_a, subset_pairs$gene_b) != ''],
  interaction_scores = interaction_scores_v5 %>% select(1:13, contains("|Ependymal-like")) %>% filter(gene_a != '') %>% filter(gene_b != ''),
  cellsign = cellsign_v5 %>% select(1:13, contains("|Ependymal-like")) %>% filter(gene_a != '') %>% filter(gene_b != ''),
  filter_by_cellsign = T)
ggsave('/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Frozen/sign_interactions2.pdf', width = 7, height = 17)





rel_interations <- unique(relevant_interactions_v5$classification)[unique(relevant_interactions_v5$classification) != '']
for (int in rel_interations) {
  subset_pairs <- relevant_interactions_v5 %>% filter(classification == int)
  
  checkpoint <- show_condition(
    plot_cpdb(
      scdata = data,
      cell_type1 = "Microglia",
      cell_type2 = '.',
      means = means_v5,
      pvals = relevant_interactions_v5,
      celltype_key = "Annotation_v2",
      max_size = 5,
      highlight_size = 0.5,
      degs_analysis = TRUE,
      standard_scale = TRUE,
      genes = c(subset_pairs$gene_a, subset_pairs$gene_b)[c(subset_pairs$gene_a, subset_pairs$gene_b) != ''],
      interaction_scores = interaction_scores_v5,
      cellsign = cellsign_v5,
      filter_by_cellsign = T)
  )
  
  if (!is.character(checkpoint)) {
    plot_cpdb(
      scdata = data,
      cell_type1 = "Microglia",
      cell_type2 = '.',
      means = means_v5,
      pvals = relevant_interactions_v5,
      celltype_key = "Annotation_v2",
      max_size = 5,
      highlight_size = 0.5,
      degs_analysis = TRUE,
      standard_scale = TRUE,
      genes = c(subset_pairs$gene_a, subset_pairs$gene_b)[c(subset_pairs$gene_a, subset_pairs$gene_b) != ''],
      interaction_scores = interaction_scores_v5,
      cellsign = cellsign_v5,
      filter_by_cellsign = T)
    tmp <- gsub(' |\\/', '_', int)
    ggsave(glue("/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Approach-1/dot_significant/{tmp}.pdf"), width = 7, height = 2, units = 'in', scale = 3)
  }
}


rel_interations <- unique(relevant_interactions_v5$classification)[unique(relevant_interactions_v5$classification) != '']
subset_pairs <- relevant_interactions_v5 %>% filter(classification %in% c("Signaling by Glycine", "Signaling by Serotonin", "Signaling by Brain-derived neurotrophic factor", "Signaling by Neurotrophin", "Signaling by GABA", "Signaling by Glutamate", "Signaling by Histamine"))
plot_cpdb(
  scdata = data,
  cell_type1 = "Microglia",
  cell_type2 = '.',
  means = means_v5,
  pvals = relevant_interactions_v5,
  celltype_key = "Annotation_v2",
  max_size = 5,
  highlight_size = 0.5,
  degs_analysis = TRUE,
  standard_scale = TRUE,
  genes = c(subset_pairs$gene_a, subset_pairs$gene_b)[c(subset_pairs$gene_a, subset_pairs$gene_b) != ''],
  interaction_scores = interaction_scores_v5,
  cellsign = cellsign_v5,
  filter_by_cellsign = T)
ggsave("/n/scratch/users/c/cao385/Ependymoma/cpdb_approach1_synapse.pdf", width = 15, height = 4)



means1 <- read.delim("/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Approach-1_method2_withScore/statistical_analysis_means_Frozen.txt", check.names = FALSE)
means1 <- means1[!duplicated(means1$interacting_pair),]

pvals1 <- read.delim("/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Approach-1_method2_withScore/statistical_analysis_pvalues_Frozen.txt", check.names = FALSE)
pvals1 <- pvals1[!duplicated(pvals1$interacting_pair),]

decon1 <- read.delim("/n/scratch/users/c/cao385/Ependymoma/cpdb/results/Approach-1_method2_withScore/statistical_analysis_deconvoluted_Frozen.txt", check.names = FALSE)
#data(cpdb_output2) # legacy reasons
meta <- data@meta.data
plot_cpdb2(
  scdata = sce,
  cell_type1 = "Microglia",
  cell_type2 = "Mature neurons",
  celltype_key = "Annotation_v2", # column name where the cell ids are located in the metadata
  means = means1,
  pvals = pvals1,
  deconvoluted = decon1, # new options from here on specific to plot_cpdb2
)


plot_cpdb3(
  scdata = sce,
  cell_type1 = "Microglia",
  cell_type2 = "Mature neurons",
  celltype_key = "Annotation_v2", # column name where the cell ids are located in the metadata
  means = means1,
  pvals = pvals1,
  deconvoluted = decon1, # new options from here on specific to plot_cpdb3
  
)