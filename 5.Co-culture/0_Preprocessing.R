library(Seurat)
library(tidyverse)
library(glue)
library(SeuratWrappers)
library(ggpubr)
library(SingleCellExperiment)
library(scDblFinder)
library(SoupX)
library(qs)

source('~/Projects/General-Codes/Resources/Plotting_helper_functions.R')
source('/home/cao385/Projects/General-Codes/Resources/single_cell_preprocessing_helper_functions.R')

base_dir <- '/n/scratch/users/c/cao385/Ependymoma/10x'


################################################################################
## Monoculture samples (Human only)
################################################################################
files <- glue("{base_dir}/data/{c('1LP', '2LP', '8LP', '9LP')}")
fnames <- basename(files)

metrics <- NULL
plot_list_tdTomato <- list()

for (i in seq_along(files)) {
  set.seed(1234)
  sc <- load10X(dataDir = files[i])
  tmp <- try(autoEstCont(sc, verbose = F, doPlot = F), silent = TRUE)
  if(!"try-error" %in% class(tmp)) {
    sc <- autoEstCont(sc)
  } else {
    sc <- setContaminationFraction(sc, 0.2)
  }
  out <- adjustCounts(sc, roundToInt = T)
  
  data <- CreateSeuratObject(counts = out, project = fnames[i])
  
  sce <- SingleCellExperiment(list(counts = LayerData(data, 'counts')))
  set.seed(1234)
  results <- scDblFinder(sce, returnType = 'table', clusters = F) %>% 
    as.data.frame() %>% 
    dplyr::filter(type == 'real')
  doublets <- results %>% dplyr::filter(class == 'doublet') %>% rownames()
  data$doublets_score <- results$score
  data$doublets <- ifelse(Cells(data) %in% doublets, 'Doublet', 'Singlet')
  
  
  data$log10GenesPerUMI <- log10(data$nFeature_RNA) / log10(data$nCount_RNA)
  data[["percent.mt"]] <- PercentageFeatureSet(object = data, pattern = "^MT-")
  data[["percent.ribo"]] <- PercentageFeatureSet(object = data, pattern = "^RPL|^RPS")
  data[["percent.hsprotein"]] <- PercentageFeatureSet(object = data, pattern = "^HSP|^DNAJ")
  
  
  qc_plot <- qc_plots(data)
  ggsave(plot = qc_plot, glue('{base_dir}/output/QC_Plot-{fnames[i]}.png'), width = 15, height = 10)
  
  
  data_filtered <- subset(data, subset = percent.mt < 10 & nCount_RNA > 500 & nCount_RNA < 8000 & nFeature_RNA > 500 & nFeature_RNA < 7000)
  data_filtered <- subset(data_filtered, subset = doublets == 'Singlet')
  
  
  tdTomato_Expression_Before <- FetchData(data, 'tdTomato-1') %>% pull(`tdTomato-1`)
  tdTomato_Expression_After <- FetchData(data_filtered, 'tdTomato-1') %>% pull(`tdTomato-1`)
  
  metrics <- rbind(metrics, data.frame(nCells_Before = ncol(data), 
                                       nCells_Normal_Before = length(which(tdTomato_Expression_Before == 0)), 
                                       nCells_Tumor_Before = sum(tdTomato_Expression_Before > 0),
                                       nDoublets = length(doublets),
                                       nCells_After = ncol(data_filtered), 
                                       nCells_Normal_After = length(which(tdTomato_Expression_After == 0)), 
                                       nCells_Tumor_After = sum(tdTomato_Expression_After > 0), 
                                       row.names = fnames[i]))
  
  
  data <- RunFullSeurat_v5(cm = LayerData(data_filtered), metadata = data_filtered@meta.data, doBatch = F, dims = 0.8, verbose = T, project = fnames[i])
  data$tdTomato <- ifelse(FetchData(data, 'tdTomato-1') %>% pull('tdTomato-1') > 0, 'yes', 'no')
  qsave(data, glue('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-{fnames[i]}.qs'))
  
  
  #pt <- DimPlot(data, group.by = 'tdTomato', cols = c("#E31A1C", "#1F78B4"), order = T) + seurat_theme() + ggtitle(fnames[i])
  plot_list_tdTomato[[fnames[i]]] <- FeaturePlot(data, 'tdTomato-1', order = T, cols = c('lightgrey', 'red'), max.cutoff = 'q95') + seurat_theme() + ggtitle(glue('tdTomato - {fnames[i]}'))
  #plot_list_tdTomato[[fnames[i]]] <- annotate_figure(ggarrange(p1, p2), top = text_grob(fnames[i], color = "black", face = "bold", size = 14))
  
  
  message(glue('{fnames[i]} \n'))
}

ggarrange(plotlist = plot_list_tdTomato)
#ggview::ggview(width = 15, height = 12)
ggsave(glue('{base_dir}/output/HumanOnly-tdTomato.png'), width = 15, height = 12)




################################################################################
## Human+Rat samples
################################################################################
files <- glue("{base_dir}/data/{c('BT165rat', 'EP1NSrat', 'VBT242rat')}")
fnames <- basename(files)

metrics <- NULL
seurat_list <- list()
plot_list_violin <- plot_list_tdTomato <- list()

for (i in seq_along(files)) {
  sc <- load10X(dataDir = files[i])
  
  tmp <- try(autoEstCont(sc, verbose = F, doPlot = F), silent = TRUE)
  if(!"try-error" %in% class(tmp)) {
    sc <- autoEstCont(sc)
  } else {
    sc <- setContaminationFraction(sc, 0.2)
  }
  out <- adjustCounts(sc, roundToInt = T)
  
  data <- CreateSeuratObject(counts = out, project = fnames[i])
  
  
  human_genes <- human_genes <- grep("^GRCh38----", rownames(data), value = TRUE)
  rat_genes <- grep("^mRatBN7-2-", rownames(data), value = TRUE)
  
  # Calculate human and rat counts for each cell
  human_counts <- Matrix::colSums(LayerData(data, 'counts')[human_genes, ])
  rat_counts <- Matrix::colSums(LayerData(data, 'counts')[rat_genes, ])
  
  # Add these counts as metadata to the Seurat object
  data <- AddMetaData(data, metadata = data.frame(human_counts = human_counts, rat_counts = rat_counts))
  
  # Add metadata for species assignment based on the proportion of reads
  data <- data %>%
    AddMetaData(metadata = data.frame(
      human_ratio = data$human_counts / (data$human_counts + data$rat_counts),
      rat_ratio = data$rat_counts / (data$human_counts + data$rat_counts)
    ))
  
  # Assign species based on a threshold (e.g., 80%)
  species_threshold <- 0.8
  data$species <- ifelse(data$human_ratio > species_threshold, "Human",
                         ifelse(data$rat_ratio > species_threshold, "Rat", "Doublet"))
  
  pt <- VlnPlot(data, features = c("human_ratio", "rat_ratio"), group.by = "species", combine = F)
  pt <- lapply(pt, function(x) x + seurat_theme())
  plot_list_violin[[fnames[i]]] <- annotate_figure(ggarrange(plotlist = pt, common.legend = T, legend = 'right'), top = text_grob(fnames[i], color = "black", face = "bold", size = 14))
  
  # Filter out doublets
  data <- subset(data, species != "Doublet")
  
  data <- RunFullSeurat_v5(cm = LayerData(data), metadata = data@meta.data, doBatch = F, dims = 0.8, verbose = T, project = fnames[i])
  data$tdTomato <- ifelse(FetchData(data, 'GRCh38----tdTomato-1') %>% pull('GRCh38----tdTomato-1') > 0, 'yes', 'no')
  qsave(data, glue('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-{fnames[i]}.qs'))
  
  
  p1 <- FeaturePlot(data, 'GRCh38----tdTomato-1', order = T, cols = c('lightgrey', 'red'), max.cutoff = 'q95') + seurat_theme() + ggtitle('tdTomato-1')
  #p2 <- DimPlot(data, group.by = 'species', cols = c("#E31A1C", "#1F78B4"), order = T) + seurat_theme() + ggtitle('')
  p2 <- DimPlot(data, group.by = 'species', cols = c("#E31A1C", "#1F78B4"), order = T, split.by = 'species', combine = F)
  p2 <- lapply(p2, function(x) x + seurat_theme() + NoLegend('') + ggtitle(''))
  p2 <- ggarrange(plotlist = p2)
  plot_list_tdTomato[[fnames[i]]] <- annotate_figure(ggarrange(p1, p2), top = text_grob(fnames[i], color = "black", face = "bold", size = 14))
  
}

ggarrange(plotlist = plot_list_violin, ncol = 1)
#ggview::ggview(width = 8, height = 25)
ggsave(glue('{base_dir}/output/HumanRat-qc.png'), width = 8, height = 25)

ggarrange(plotlist = plot_list_tdTomato, ncol = 1)
#ggview::ggview(width = 8, height = 20)
ggsave(glue('{base_dir}/output/HumanRat-tdTomato.png'), width = 8, height = 20)





################################################################################
## Rat samples
################################################################################
files <- glue("{base_dir}/data/{c('7LP')}")
fnames <- basename(files)

metrics <- NULL
seurat_list <- list()
plot_list_violin <- plot_list_tdTomato <- list()

for (i in seq_along(files)) {
  set.seed(1234)
  sc <- load10X(dataDir = files[i])
  tmp <- try(autoEstCont(sc, verbose = F, doPlot = F), silent = TRUE)
  if(!"try-error" %in% class(tmp)) {
    sc <- autoEstCont(sc)
  } else {
    sc <- setContaminationFraction(sc, 0.2)
  }
  out <- adjustCounts(sc, roundToInt = T)
  
  data <- CreateSeuratObject(counts = out, project = fnames[i])
  
  sce <- SingleCellExperiment(list(counts = LayerData(data, 'counts')))
  set.seed(1234)
  results <- scDblFinder(sce, returnType = 'table', clusters = F) %>% 
    as.data.frame() %>% 
    dplyr::filter(type == 'real')
  doublets <- results %>% dplyr::filter(class == 'doublet') %>% rownames()
  data$doublets_score <- results$score
  data$doublets <- ifelse(Cells(data) %in% doublets, 'Doublet', 'Singlet')
  
  
  data$log10GenesPerUMI <- log10(data$nFeature_RNA) / log10(data$nCount_RNA)
  data[["percent.mt"]] <- PercentageFeatureSet(object = data, pattern = "^Mt-")
  data[["percent.ribo"]] <- PercentageFeatureSet(object = data, pattern = "^Rpl|^Rps")
  data[["percent.hsprotein"]] <- PercentageFeatureSet(object = data, pattern = "^Hsp|^Dnaj")
  
  
  qc_plot <- qc_plots(data)
  ggsave(plot = qc_plot, glue('{base_dir}/output/QC_Plot-{fnames[i]}.png'), width = 15, height = 10)
  
  
  data_filtered <- subset(data, subset = percent.mt < 10 & nCount_RNA > 500 & nCount_RNA < 8000 & nFeature_RNA > 500 & nFeature_RNA < 7000)
  data_filtered <- subset(data_filtered, subset = doublets == 'Singlet')
  data_filtered <- subset(data_filtered, features = grep("^ENSRNOG", rownames(data_filtered), value = TRUE, invert = T))
  
  tdTomato_Expression_Before <- FetchData(data, 'tdTomato-1') %>% pull(`tdTomato-1`)
  tdTomato_Expression_After <- FetchData(data_filtered, 'tdTomato-1') %>% pull(`tdTomato-1`)
  
  metrics <- rbind(metrics, data.frame(nCells_Before = ncol(data), 
                                       nCells_Normal_Before = length(which(tdTomato_Expression_Before == 0)), 
                                       nCells_Tumor_Before = sum(tdTomato_Expression_Before > 0),
                                       nDoublets = length(doublets),
                                       nCells_After = ncol(data_filtered), 
                                       nCells_Normal_After = length(which(tdTomato_Expression_After == 0)), 
                                       nCells_Tumor_After = sum(tdTomato_Expression_After > 0), 
                                       row.names = fnames[i]))
  
  
  data <- RunFullSeurat_v5(cm = LayerData(data_filtered), metadata = data_filtered@meta.data, doBatch = F, dims = 0.8, verbose = T, project = fnames[i])
  data$tdTomato <- ifelse(FetchData(data, 'tdTomato-1') %>% pull('tdTomato-1') > 0, 'yes', 'no')
  qsave(data, glue('/n/scratch/users/c/cao385/Ependymoma/10x/output/SeuratObject-{fnames[i]}.qs'))
  
  
  #pt <- DimPlot(data, group.by = 'tdTomato', cols = c("#E31A1C", "#1F78B4"), order = T) + seurat_theme() + ggtitle(fnames[i])
  #plot_list_tdTomato[[fnames[i]]] <- FeaturePlot(data, 'tdTomato-1', order = T, cols = c('lightgrey', 'red'), max.cutoff = 'q95') + seurat_theme() + ggtitle(glue('tdTomato - {fnames[i]}'))
  #plot_list_tdTomato[[fnames[i]]] <- annotate_figure(ggarrange(p1, p2), top = text_grob(fnames[i], color = "black", face = "bold", size = 14))
  
  
  message(glue('{fnames[i]} \n'))
  
}















# Proceed with downstream analysis separately for human and rat cells
data_Human <- subset(data, species == "Human")
data_Human <- subset(data_Human, features = grep("^GRCh38----", rownames(data_Human), value = TRUE))
rownames(data_Human) <- gsub("GRCh38----", "", rownames(data_Human))
data_Human <- subset(data_Human, features = grep("^ENSG", rownames(data_Human), value = TRUE, invert = T))

data_Rat <- subset(data, species == "Rat")
data_Rat <- subset(data_Rat, features = grep("^mRatBN7-2-|tdTomato", rownames(data_Rat), value = TRUE))
rownames(data_Rat) <- gsub("mRatBN7-2-|GRCh38----", "", rownames(data_Rat))
data_Rat <- subset(data_Rat, features = grep("^ENSRNOG", rownames(data_Rat), value = TRUE, invert = T))





sce_Human <- SingleCellExperiment(list(counts = LayerData(data_Human, 'counts')))
set.seed(1234)
results <- scDblFinder(sce_Human, returnType = 'table', clusters = F) %>% 
  as.data.frame() %>% 
  dplyr::filter(type == 'real')
doublets <- results %>% dplyr::filter(class == 'doublet') %>% rownames()
data_Human$doublets_score <- results$score
data_Human$doublets <- ifelse(Cells(data_Human) %in% doublets, 'Doublet', 'Singlet')

data_Human[["percent.mt"]] <- PercentageFeatureSet(object = data_Human, pattern = "^MT-")
data_Human[["percent.ribo"]] <- PercentageFeatureSet(object = data_Human, pattern = "^RPL|^RPS")
data_Human[["percent.hsprotein"]] <- PercentageFeatureSet(object = data_Human, pattern = "^HSP|^DNAJ")

plot1 <- FeatureScatter(data_Human, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = 'lightgrey', group.by = 'orig.ident') + 
  geom_point(data = data_Human@meta.data %>% filter(percent.mt < 10, nCount_RNA > 500, nCount_RNA < 8000), aes(x = nCount_RNA, y = percent.mt), color = '#E31A1C') + 
  geom_hline(yintercept = 10, linetype = 'dashed') + 
  geom_vline(xintercept = 8000, linetype = 'dashed') + 
  seurat_theme() + 
  NoLegend() + 
  labs(title = '')

plot2 <- FeatureScatter(data_Human, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 'lightgrey', group.by = 'orig.ident') + 
  geom_point(data = data_Human@meta.data %>% filter(nFeature_RNA > 500, nFeature_RNA < 7000, nCount_RNA > 500, nCount_RNA < 8000), aes(x = nCount_RNA, y = nFeature_RNA), color = '#1F78B4') + 
  geom_hline(yintercept = c(500,7000), linetype = 'dashed') + 
  geom_vline(xintercept = c(500,8000), linetype = 'dashed') + 
  seurat_theme() + 
  NoLegend() + 
  labs(title = '')

plot_list_qc[[fnames[i]]] <- annotate_figure(ggarrange(plot1, plot2), top = text_grob(fnames[i], color = "black", face = "bold", size = 14))


data_filtered_Human <- subset(data_Human, subset = percent.mt < 10 & nFeature_RNA > 500 & nFeature_RNA < 8000 & nFeature_RNA > 500 & nFeature_RNA < 7000)
data_filtered_Human <- subset(data_filtered_Human, subset = doublets == 'Singlet')


seurat_list[[fnames[i]]] <- data_filtered_Human

tdTomato_Expression_Before <- FetchData(data_Human, 'tdTomato-1', layer = 'counts') %>% pull(`tdTomato-1`)
tdTomato_Expression_After <- FetchData(data_filtered_Human, 'tdTomato-1', layer = 'counts') %>% pull(`tdTomato-1`)

metrics <- rbind(metrics, data.frame(nCells_Before = ncol(data_Human), 
                                     nCells_Normal_Before = length(which(tdTomato_Expression_Before == 0)), 
                                     nCells_Tumor_Before = sum(tdTomato_Expression_Before > 0),
                                     nDoublets = length(doublets),
                                     nCells_After = ncol(data_filtered_Human), 
                                     nCells_Normal_After = length(which(tdTomato_Expression_After == 0)), 
                                     nCells_Tumor_After = sum(tdTomato_Expression_After > 0), 
                                     row.names = fnames[i]))


data_Human <- RunFullSeurat_v5(cm = LayerData(data_filtered_Human), metadata = data_filtered_Human@meta.data, doBatch = F, dims = 0.8, verbose = T, project = fnames[i])
data_Human$tdTomato <- ifelse(FetchData(data_Human, 'tdTomato-1') %>% pull('tdTomato-1') > 0, 'tumor', 'normal')

p1 <- DimPlot(data_Human, group.by = 'tdTomato', cols = c("#E31A1C", "#1F78B4"), order = T) + seurat_theme() + ggtitle('')
p2 <- FeaturePlot(data_Human, 'tdTomato-1', order = T, cols = c('lightgrey', 'red'), max.cutoff = 'q95') + seurat_theme()

plot_list_tdTomato[[fnames[i]]] <- annotate_figure(ggarrange(p1, p2), top = text_grob(fnames[i], color = "black", face = "bold", size = 14))

message(glue('{fnames[i]} \n'))
