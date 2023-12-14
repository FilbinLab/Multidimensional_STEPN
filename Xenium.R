library(Seurat)
library(ggplot2)
library(glue)
library(tidyverse)
library(qs)
library(future)
library(ggpubr)
library(data.table)
library(ComplexHeatmap)
library(SeuratWrappers)
library(patchwork)

resources_dir <- "/home/cao385/Projects/General-Codes/Resources"
source(file.path(resources_dir, "single_cell_preprocessing_helper_functions.R"))
source(file.path(resources_dir, "NMF_helper_function.R"))
source(file.path(resources_dir, "Plotting_helper_functions.R"))

#plan("multisession", workers = 5)
options(future.globals.maxSize = +Inf)

optimizePCA <- function(sobj, csum){
  dp <- Stdev(sobj)^2
  for (z in 1:length(dp)) {
    soma <- sum(dp[1:z])/sum(dp)
    if (soma >= csum) {
      best_pc <- z
      break()
    }
  }
  return(best_pc)
}


base_dir <- "/n/scratch/users/c/cao385/Immune"
run_dir <- '20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41'
samples <- list.files(glue('{base_dir}/runs/{run_dir}'))
Strings2remove <- 'output-XETG00083__|__20231208__193752'

for (smp in samples) {
  dir.create(glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}"), recursive = T)

  xenium.obj <- LoadXenium(glue('{base_dir}/runs/{run_dir}/{smp}'), fov = "fov")
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

  xenium.obj <- SCTransform(xenium.obj, assay = "Xenium", vars.to.regress = c('nFeature_Xenium', 'nCount_Xenium'), vst.flavor = "v2")
  xenium.obj <- RunPCA(xenium.obj, npcs = 50, features = rownames(xenium.obj))
  xenium.obj <- RunUMAP(xenium.obj, dims = 1:optimizePCA(xenium.obj, 0.8))
  xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:optimizePCA(xenium.obj, 0.8))
  xenium.obj <- FindClusters(xenium.obj, resolution = seq(0.1, 1, 0.1))

  qsave(xenium.obj, glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/{gsub('__', '-', gsub(Strings2remove, '', smp))}.qs"))
}




for (smp in samples) {
  xenium.obj <- qread(glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/{gsub('__', '-', gsub(Strings2remove, '', smp))}.qs"))

  clusters <- grep('SCT_snn_res', colnames(xenium.obj@meta.data), value = T)

  pt1 <- mapply(function(x,y) {
    DimPlot(xenium.obj, group.by = x, cols = clusterExperiment::bigPalette, pt.size = 0.5, raster = F) +
      theme_bw() +
      labs(title = y, x = '', y = '') +
      theme(panel.background = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
            axis.ticks.length = unit(0, "cm"), axis.text = element_text(size = 0),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.caption = element_text(hjust = 0.5, size = 16))
  }, x = clusters, y = glue('Resolution {seq(0.1,1,0.1)}'))
  pt1 <- wrap_plots(pt1, ncol = 4, nrow = 3)
  ggsave(plot = pt1, filename = glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/UMAP.png"), width = 18, height = 12)


  pt2 <- mapply(function(x,y) {
    ImageDimPlot(xenium.obj, cols = "polychrome", size = 0.75, group.by = x) + ggtitle(y)
    ggsave(filename = glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/ImagePlot/{gsub(' ', '', y)}.png"), width = 8, height = 6)
  }, x = clusters, y = glue('Resolution {seq(0.1,1,0.1)}'))


  markers <- lapply(clusters, function(x,y) {
    Idents(xenium.obj) <- xenium.obj[[x]][[1]]
    tmp <- RunPrestoAll(xenium.obj, only.pos = T, verbose = T) %>%
      filter(p_val_adj < 0.05) %>%
      group_by(cluster) %>%
      arrange(-avg_log2FC, .by_group = T)
    return(tmp)
  })
  names(markers) <- seq(0.1,1,0.1)
  markers <- bind_rows(markers, .id = 'resolution')
  write.csv(markers, glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/markers.csv"), quote = F, row.names = F)
}





base_panel <- fread(glue('{base_dir}/Xenium/Base_Panel.tsv'))
base_panel <- base_panel %>% filter(!Annotation %in% c('Astrocyte', 'OPC'))
base_panel <- lapply(split(base_panel, f = base_panel$Annotation), function(x) x$Gene)

custom_panel <- fread(glue('{base_dir}/Xenium/Customize_Panel.tsv'))
custom_panel <- custom_panel %>% filter(!Type %in% c('DMG', 'GBM', 'GBM/Normal', 'Intersection'))
custom_panel <- lapply(split(custom_panel, f = custom_panel$Program), function(x) x$Gene)

geneList <- c(base_panel, custom_panel)

samples <- list.files(glue('{base_dir}/runs/{run_dir}'))[8]

for (smp in samples) {
  dir.create(glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/corr"))

  xenium.obj <- qread(glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/{gsub('__', '-', gsub(Strings2remove, '', smp))}.qs"))
  markers <- read.csv(glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/markers.csv"))

  geneList2 <- lapply(geneList, function(x) x[x %in% rownames(xenium.obj)])
  geneList2 <- geneList2[unlist(lapply(geneList2, length)) > 2]

  cm_norm <- as.matrix(log2(LayerData(xenium.obj, 'counts')/10+1))
  cm_mean <- log2(Matrix::rowMeans(LayerData(xenium.obj, 'counts'))+1)
  cm_center <- cm_norm - rowMeans(cm_norm)

  #nres <- seq(0.1, 1, 0.1)
  nres <- c(0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1)
  for (res in nres) {
    mks <- markers %>% filter(resolution == res) %>% group_by(cluster) %>% top_n(30, wt = avg_log2FC)
    mks <- lapply(split(mks, f = mks$cluster), function(x) x$gene)
    names(mks) <- glue('Cluster_{names(mks)}')

    nmf_score <- t(scoreNmfGenes(cm_center, cm_mean, c(mks, geneList2), verbose = F))

    nmf_factor_hc <- clusterNmfFactors(nmf_score)
    nmf_factor_cor <- nmf_factor_hc$cor_coef[nmf_factor_hc$hc_obj$order, nmf_factor_hc$hc_obj$order]

    hm_colors <- rev((brewer.pal(n=9, name="RdBu")))
    hm_colors <- colorRampPalette(colors = hm_colors)

    fontcolors <- ifelse(grepl('Cluster_', rownames(nmf_factor_cor)), 'red', 'darkgreen')
    fontfaces <- ifelse(grepl('Cluster_', rownames(nmf_factor_cor)), 'bold', 'plain')
    rowAnno <- rowAnnotation(rows = anno_text(rownames(nmf_factor_cor), gp = gpar(fontface = fontfaces, col = fontcolors, fontsize = 8)))

    pdf(glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/corr/Correlation_Signatures_Res{res}.pdf"), width = 15, height = 10)
    draw(Heatmap(nmf_factor_cor, col = hm_colors(100), cluster_rows = F, cluster_columns = F, show_column_dend = F, right_annotation = rowAnno,
                 show_row_dend = F, show_row_names = F, show_column_names = F, name = 'cor', row_names_gp = gpar(fontsize = 8),
                 layer_fun = function(j, i, x, y, width, height, fill) {
                   grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", lwd = 0.5, fill = NA))}))
    dev.off()
  }
}

samples <- list.files(glue('{base_dir}/runs/{run_dir}'))
smp <- samples[8]

xenium.obj <- qread(glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/{gsub('__', '-', gsub(Strings2remove, '', smp))}.qs"))
markers <- read.csv(glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/markers.csv"))
geneList2 <- lapply(geneList, function(x) x[x %in% rownames(xenium.obj)])
geneList2 <- geneList2[unlist(lapply(geneList2, length)) > 2]
cm_norm <- as.matrix(log2(LayerData(xenium.obj, 'counts')/10+1))
cm_mean <- log2(Matrix::rowMeans(LayerData(xenium.obj, 'counts'))+1)
cm_center <- cm_norm - rowMeans(cm_norm)

res <- 0.3
mks <- markers %>% filter(resolution == res) %>% group_by(cluster) %>% top_n(30, wt = avg_log2FC)
mks <- lapply(split(mks, f = mks$cluster), function(x) x$gene)
names(mks) <- glue('Cluster_{names(mks)}')
nmf_score <- t(scoreNmfGenes(cm_center, cm_mean, c(mks, geneList2), verbose = F))
nmf_factor_hc <- clusterNmfFactors(nmf_score)
nmf_factor_cor <- nmf_factor_hc$cor_coef[nmf_factor_hc$hc_obj$order, nmf_factor_hc$hc_obj$order]
hm_colors <- rev((brewer.pal(n=9, name="RdBu")))
hm_colors <- colorRampPalette(colors = hm_colors)
fontcolors <- ifelse(grepl('Cluster_', rownames(nmf_factor_cor)), 'red', 'darkgreen')
fontfaces <- ifelse(grepl('Cluster_', rownames(nmf_factor_cor)), 'bold', 'plain')
rowAnno <- rowAnnotation(rows = anno_text(rownames(nmf_factor_cor), gp = gpar(fontface = fontfaces, col = fontcolors, fontsize = 8)))
pdf(glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/corr/Correlation_Signatures_Res{res}.pdf"), width = 15, height = 10)
draw(Heatmap(nmf_factor_cor, col = hm_colors(100), cluster_rows = F, cluster_columns = F, show_column_dend = F, right_annotation = rowAnno,
             show_row_dend = F, show_row_names = F, show_column_names = F, name = 'cor', row_names_gp = gpar(fontsize = 8),
             layer_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", lwd = 0.5, fill = NA))}))
dev.off()


res <- 0.7
mks <- markers %>% filter(resolution == res) %>% group_by(cluster) %>% top_n(30, wt = avg_log2FC)
mks <- lapply(split(mks, f = mks$cluster), function(x) x$gene)
names(mks) <- glue('Cluster_{names(mks)}')
nmf_score <- t(scoreNmfGenes(cm_center, cm_mean, c(mks, geneList2), verbose = F))
nmf_factor_hc <- clusterNmfFactors(nmf_score)
nmf_factor_cor <- nmf_factor_hc$cor_coef[nmf_factor_hc$hc_obj$order, nmf_factor_hc$hc_obj$order]
hm_colors <- rev((brewer.pal(n=9, name="RdBu")))
hm_colors <- colorRampPalette(colors = hm_colors)
fontcolors <- ifelse(grepl('Cluster_', rownames(nmf_factor_cor)), 'red', 'darkgreen')
fontfaces <- ifelse(grepl('Cluster_', rownames(nmf_factor_cor)), 'bold', 'plain')
rowAnno <- rowAnnotation(rows = anno_text(rownames(nmf_factor_cor), gp = gpar(fontface = fontfaces, col = fontcolors, fontsize = 8)))
pdf(glue("{base_dir}/Xenium/results/{run_dir}/{gsub('__', '-', gsub(Strings2remove, '', smp))}/corr/Correlation_Signatures_Res{res}.pdf"), width = 15, height = 10)
draw(Heatmap(nmf_factor_cor, col = hm_colors(100), cluster_rows = F, cluster_columns = F, show_column_dend = F, right_annotation = rowAnno,
             show_row_dend = F, show_row_names = F, show_column_names = F, name = 'cor', row_names_gp = gpar(fontsize = 8),
             layer_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", lwd = 0.5, fill = NA))}))
dev.off()


















################################################################################
## BT1873 (K27M): 20231115__203753__BT1873_BT1733 - 0011155-Region_1.qs
################################################################################
data <- qread('/n/scratch3/users/c/cao385/Xenium/results/20231115__203753__BT1873_BT1733/0011155-Region_1/0011155-Region_1.qs')

DimPlot(data, group.by = 'SCT_snn_res.0.3', cols = clusterExperiment::bigPalette, pt.size = 0.3, raster = FALSE) +
  theme_bw() +
  labs(title = '', x = '', y = '') +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
        axis.ticks.length = unit(0, "cm"), axis.text = element_text(size = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = 16))



base_panel <- fread(glue('{base_dir}/Base_Panel.tsv'))
base_panel <- base_panel %>% filter(!Annotation %in% c('Astrocyte', 'OPC'))
base_panel <- lapply(split(base_panel, f = base_panel$Annotation), function(x) x$Gene)
base_panel <- base_panel[c(grep('Immune', names(base_panel)), grep('NormalBrain', names(base_panel)))]

custom_panel <- fread(glue('{base_dir}/Customize_Panel.tsv'))
custom_panel <- custom_panel %>% filter(!Type %in% c('DMG', 'Ependymoma'))
custom_panel <- lapply(split(custom_panel, f = custom_panel$Program), function(x) x$Gene)

geneList <- c(custom_panel, base_panel)



markers <- read.csv('/n/scratch3/users/c/cao385/Xenium/results/20231115__203753__BT1873_BT1733/0011155-Region_1/markers.csv')
geneList2 <- lapply(geneList, function(x) x[x %in% rownames(data)])
geneList2 <- geneList2[unlist(lapply(geneList2, length)) > 2]
cm_norm <- as.matrix(log2(LayerData(data, 'counts')/10+1))
cm_mean <- log2(Matrix::rowMeans(LayerData(data, 'counts'))+1)
cm_center <- cm_norm - rowMeans(cm_norm)

mks <- markers %>% filter(resolution == 0.3) %>% group_by(cluster) %>% top_n(30, wt = avg_log2FC)
mks <- lapply(split(mks, f = mks$cluster), function(x) x$gene)
names(mks) <- glue('Cluster_{names(mks)}')
nmf_score <- t(scoreNmfGenes(cm_center, cm_mean, c(mks, geneList2), verbose = F))
nmf_factor_hc <- clusterNmfFactors(nmf_score)
nmf_factor_cor <- nmf_factor_hc$cor_coef[nmf_factor_hc$hc_obj$order, nmf_factor_hc$hc_obj$order]
hm_colors <- rev((brewer.pal(n=9, name="RdBu")))
hm_colors <- colorRampPalette(colors = hm_colors)
fontcolors <- ifelse(grepl('Cluster_', rownames(nmf_factor_cor)), 'red', 'darkgreen')
fontfaces <- ifelse(grepl('Cluster_', rownames(nmf_factor_cor)), 'bold', 'plain')
rowAnno <- rowAnnotation(rows = anno_text(rownames(nmf_factor_cor), gp = gpar(fontface = fontfaces, col = fontcolors, fontsize = 8)))
Heatmap(nmf_factor_cor, col = hm_colors(100), cluster_rows = F, cluster_columns = F, show_column_dend = F, right_annotation = rowAnno,
        show_row_dend = F, show_row_names = F, show_column_names = F, name = 'cor', row_names_gp = gpar(fontsize = 8),
        layer_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", lwd = 0.5, fill = NA))})


Idents(data) <- data$SCT_snn_res.0.3
data <- RenameIdents(data, '12' = 'AC-like', '7' = 'AC-like', '0' = 'AC-like', '6' = 'AC-like', '4' = 'AC-like', '10' = 'AC-like', '8' = 'Neurons', '9' = 'MES2-like', '1' = 'OPC-like', '3' = 'OPC-like', '11' = 'OPC-like', '5' = 'MES1-like', '2' = 'Tcells', '13' = 'Tcells')

DimPlot(data, cols = colors_to_use, pt.size = 0.5, raster = FALSE) +
  theme_bw() +
  labs(title = '', x = '', y = '') +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
        axis.ticks.length = unit(0, "cm"), axis.text = element_text(size = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = 16))


DotPlot(data, cols = c('lightgrey', 'red'),
        c('CCL4', ## Inflammatory
          'CD48', 'CD52', 'S100A4', ## Monocytes
          'IFITM3', ## IFN TAM
          'TGFBI', 'GPNMB', ## Hypoxic TAM
          'PDCD1', 'CTLA4'
        )) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size = 8))


base_panel <- fread(glue('{base_dir}/Base_Panel.tsv'))
custom_panel <- fread(glue('{base_dir}/Customize_Panel.tsv'))

DotPlot(data, cols = c('lightgrey', 'red'),
        base_panel %>% filter(grepl('Immune', Annotation)) %>% pull(Gene)) +
  theme(axis.text.x = element_text(angle = 90), axis.text=element_text(size = 8))

#custom_panel %>% filter(Type == 'Normal') %>% pull(Gene)
#base_panel %>% filter(grepl('Immune', Annotation)) %>% pull(Gene)

data <- qread('/n/scratch3/users/c/cao385/Xenium/results/20231115__203753__BT1873_BT1733/0011155-Region_1/0011155-Region_1.qs')
Idents(data) <- data$SCT_snn_res.0.3
data <- RenameIdents(data, '12' = 'AC-like', '7' = 'AC-like', '0' = 'AC-like', '6' = 'AC-like', '4' = 'AC-like', '10' = 'AC-like', '8' = 'Neurons', '9' = 'MES2-like', '1' = 'OPC-like', '3' = 'OPC-like', '11' = 'OPC-like', '5' = 'MES1-like', '2' = 'Tcells', '13' = 'Tcells')
data$Annotation <- ifelse(Idents(data) %in% 'Tcells', 'Immune', 'Malignant')

data <- BuildNicheAssay(object = data, fov = "fov", group.by = "Annotation", niches.k = 8, neighbors.k = 30)
celltype.plot <- ImageDimPlot(data, group.by = "Annotation", size = 0.5, cols = colors_to_use[1:2], dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(data, group.by = "niches", size = 0.5, dark.background = F) + ggtitle("Niches") + scale_fill_manual(values = colors_to_use)
celltype.plot | niche.plot

t(table(data$Annotation, data$niches))