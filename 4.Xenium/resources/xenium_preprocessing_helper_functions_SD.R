## Load libraries
library(tidyverse)
library(glue)
library(qs)
library(SeuratWrappers)
library(Seurat)
library(ggpubr)
library(crayon)
library(parallel)
library(readxl)

## Function to preprocess the xenium sequencing (creates seurat object, find markers for different resolutions and score the signatures)
## @para: rawDir (sample sequencing raw dir)
## @para: outDir (output directory)
## @para: base_panel (brain base panel provided by 10x)
## @para: custom_panel (customized panel)
## @returns: save files in output directory

RunFullXenium <- function(rawDir, outDir) {
  
  # Check if the seurat file already exists
  if (file.exists(file.path(glue('{outDir}/seurat.qs')))) {
    print(glue("The file", {outDir}, " already exists. Skipping this part."))
    
    } else {
    

  tmp <- read.csv(glue('{rawDir}/metrics_summary.csv'))
  smp <- glue('{tmp$cassette_name}_{tmp$region_name}')
  
  outDir <- glue('{outDir}/{smp}')
  if (!dir.exists(outDir)) {dir.create(outDir, recursive = T)}
  
  ## Preprocessing seurat objects
  message(green('-------------------------------------------------------'))
  message(blue(glue('Starting pre-processing for sample: {yellow$underline$bold(smp)}!!')))
  message(blue('[1/8] Loading Xenium raw files...'))
  
  suppressWarnings(suppressMessages(xenium.obj <- LoadXenium(rawDir, fov = "fov")))
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
  message(glue("{ncol(xenium.obj)} cells."))
  
  message(blue('[2/8] Running SCTransform-based normalization...'))
  xenium.obj <- SCTransform(xenium.obj, assay = "Xenium", vars.to.regress = c('nFeature_Xenium', 'nCount_Xenium'), vst.flavor = "v2", verbose = FALSE)
  
  message(blue('[3/8] Running PCA and UMAP reductions...'))
  xenium.obj <- RunPCA(xenium.obj, npcs = 50, features = rownames(xenium.obj), verbose = FALSE)
  xenium.obj <- RunUMAP(xenium.obj, dims = 1:optimizePCA(xenium.obj, 0.8), verbose = FALSE)
  
  # message(blue('[4/8] Clustering...'))
  # xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:optimizePCA(xenium.obj, 0.8), verbose = FALSE)
  # xenium.obj <- FindClusters(xenium.obj, resolution = seq(0.1, 1, 0.1), verbose = FALSE)
  # xenium.obj <- changeClusterNumbers(xenium.obj)
  qsave(xenium.obj, glue('{outDir}/seurat.qs'))
  }
  
}
  


scoreNmfGenes <- function(cm_center, cm_mean, nmf_gene_list, cores = 8, verbose = FALSE, simple = FALSE){
  #message("Scoring signatures... \n")
  
  #results <- pbapply::pblapply(nmf_gene_list, function(x) {
  results <- lapply(nmf_gene_list, function(x) {
    scores <- scoreSignature(cm_center, cm_mean, x, verbose = verbose, simple = simple, cores = cores)
  })
  results <- do.call('rbind', results)
  
  return(results)
}


scoreSignature <- function(X.center, X.mean, s, n = 100, cores, simple = FALSE, verbose = FALSE) {
  if(verbose) {
    message("cells: ", ncol(X.center))
    message("genes: ", nrow(X.center))
    message("genes in signature: ", length(s))
    message("Using simple average?", simple)
    message("processing...")
  }
  
  s <- intersect(rownames(X.center), s)
  if (verbose) {message("genes in signature, and also in this dataset: ", length(s))}
  ##message("These genes are: ", s)
  
  if (simple){
    s.score <- colMeans(X.center[s,])
  }else{
    if (length(s) > 100) {
      s.score <- Matrix::colMeans(do.call(rbind, mclapply(s, function(g) {
        g.n <- names(sort(abs(X.mean[g] - X.mean))[2:(n+1)])
        X.center[g, ] - Matrix::colMeans(X.center[g.n, ])
      }, mc.cores = cores)))
    } else {
      s.score <- colMeans(do.call(rbind, lapply(s, function(g) {
        g.n <- names(sort(abs(X.mean[g] - X.mean))[2:(n+1)])
        X.center[g, ] - colMeans(X.center[g.n, ])
      })))
    }
  }
  if(verbose) message(" done")
  return(s.score)
}



clusterNmfFactors <- function(nmf_basis, cor_method="pearson"){
  nmf_dist = 1-cor(nmf_basis, method=cor_method)
  hc_nmf = hclust(as.dist(nmf_dist), method="ward.D2")
  result = list()
  result[["cor_coef"]] = 1-nmf_dist
  result[["dist"]] = nmf_dist
  result[["hc_obj"]] = hc_nmf
  return(result)
}


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


## Customized ggplot2 theme for plotting
seurat_theme <- function(){
  theme_bw() +
    theme(panel.background = element_rect(colour = "black", linewidth = 0.1),
          plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
          axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}


## Make a bar graph to summarize proportion of clusters/programs/etc in each sample
## @para: x, x-axis data
## @para: y, y-axis data
## @para: x_order: order of x-axis variables
## @para: y_order: order of y-axis variables
## @para: col_names: column names for df for plotting
## @para: x_var: name of x-axis variable to plot
## @para: y_var: name of y-axis variable to plot
## @para: fill_var: name of variable for filling the bar graph
## @para: plot_title: title of the plot
## @para: x_title: title of x-axis
## @para: y_title: title of y-axis
## @para: bar_position: stack vs dodge
## The rest parameters are sizes of different labels
plotProportion <- function(x, y, x_order, y_order, col_names,
                           x_var, y_var, fill_var, colors,
                           plot_title="", x_title="", y_title="",
                           bar_position = "stack",
                           title_size=32, axis_title_size=28,
                           x_text_size=20, x_text_angle=45, x_text_hjust=1,
                           y_text_size=24, legend_title_size=28, legend_text_size=24){
  crosstab = table(x, y)
  crosstab = crosstab/rowSums(crosstab)
  crosstab = crosstab[x_order, y_order]
  crosstab = data.frame(crosstab)
  
  colnames(crosstab) = col_names
  ggplot(crosstab, aes_string(x=x_var, y=y_var, fill=fill_var)) +
    geom_bar(stat="identity", color="black", position = bar_position) +
    scale_fill_manual(values = colors) +
    ggtitle(plot_title) + xlab(x_title) + ylab(y_title) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, angle = 0, size = title_size, face = "bold", vjust = 1),
          axis.title = element_text(face="bold", size=axis_title_size),
          axis.text.x = element_text(angle=x_text_angle, hjust=x_text_hjust, size=x_text_size),
          axis.text.y = element_text(size=y_text_size),
          legend.title = element_text(face="bold", size=legend_title_size),
          legend.text = element_text(size=legend_text_size), 
          axis.ticks.length = unit(.2, "cm")) + 
    scale_y_continuous(expand = c(0.01,0)) +
    scale_x_discrete(expand = c(0,0.5))
}


changeClusterNumbers <- function(seurat_obj){ 
  seurat_obj$SCT_snn_res.0.1 <- factor(as.numeric(seurat_obj$SCT_snn_res.0.1), levels = 1:length(unique(as.numeric(seurat_obj$SCT_snn_res.0.1))))
  seurat_obj$SCT_snn_res.0.2 <- factor(as.numeric(seurat_obj$SCT_snn_res.0.2), levels = 1:length(unique(as.numeric(seurat_obj$SCT_snn_res.0.2))))
  seurat_obj$SCT_snn_res.0.3 <- factor(as.numeric(seurat_obj$SCT_snn_res.0.3), levels = 1:length(unique(as.numeric(seurat_obj$SCT_snn_res.0.3))))
  seurat_obj$SCT_snn_res.0.4 <- factor(as.numeric(seurat_obj$SCT_snn_res.0.4), levels = 1:length(unique(as.numeric(seurat_obj$SCT_snn_res.0.4))))
  seurat_obj$SCT_snn_res.0.5 <- factor(as.numeric(seurat_obj$SCT_snn_res.0.5), levels = 1:length(unique(as.numeric(seurat_obj$SCT_snn_res.0.5))))
  seurat_obj$SCT_snn_res.0.6 <- factor(as.numeric(seurat_obj$SCT_snn_res.0.6), levels = 1:length(unique(as.numeric(seurat_obj$SCT_snn_res.0.6))))
  seurat_obj$SCT_snn_res.0.7 <- factor(as.numeric(seurat_obj$SCT_snn_res.0.7), levels = 1:length(unique(as.numeric(seurat_obj$SCT_snn_res.0.7))))
  seurat_obj$SCT_snn_res.0.8 <- factor(as.numeric(seurat_obj$SCT_snn_res.0.8), levels = 1:length(unique(as.numeric(seurat_obj$SCT_snn_res.0.8))))
  seurat_obj$SCT_snn_res.0.9 <- factor(as.numeric(seurat_obj$SCT_snn_res.0.9), levels = 1:length(unique(as.numeric(seurat_obj$SCT_snn_res.0.9))))
  seurat_obj$SCT_snn_res.1 <- factor(as.numeric(seurat_obj$SCT_snn_res.1), levels = 1:length(unique(as.numeric(seurat_obj$SCT_snn_res.1))))
  return(seurat_obj)
}


  
  
# AssignLabels <- function(seurat_obj_dir, outDir, panel) { 
#   ## Assign labels
#   # Process Xenium data  -------------------------------------
#   
#   # Read data
#   data <- qread(file.path(seurat_obj_dir_ref, paste0(metadata$SampleName[i], '/seurat.qs')))
#   print(paste0('Analyzing: ', metadata$SampleName[i]))
#   
#   # 1) transfer anchors from single cell data
#   Idents(seurat_object) <- 'celltype'
#   anchors <- FindTransferAnchors(reference = seurat_object, query = data, normalization.method = "SCT", npcs = 20, reference.assay = "SCT", query.assay = "SCT")
#   # transfer labels
#   predictions <- TransferData(
#     anchorset = anchors,
#     refdata = seurat_object$celltype
#   )
#   data <- AddMetaData(object = data, metadata = predictions$predicted.id, col.name = 'predicted_label_snRNAseq')
#   data <- AddMetaData(object = data, metadata = predictions$prediction.score.max, col.name = 'predicted_label_snRNAseq_score')
#   
#   # If cell score < 0.3 --> unassigned
#   data$predicted_label_snRNAseq <- ifelse(data$predicted_label_snRNAseq_score < 0.3, 'Unknown', data$predicted_label_snRNAseq)
#   
#   
#   # 2) Add scores for Xenium panel marker genes
#   DefaultAssay(data) <- 'SCT'
#   data <- AddModuleScore_UCell(data, features = marker_genes_list)
#   U_Cell_signatures <- paste0(names(marker_genes_list), '_UCell')
#   
#   # retrieve highest UCell score
#   UCell_scores <- data@meta.data[, U_Cell_signatures]
#   
#   # Assign cell identity based on highest value (if cell score < 0.5 --> unassigned)
#   UCell_scores$cell_type_UCell <- apply(UCell_scores, 1, function(row) {
#     max_score <- max(row)
#     if (max_score < 0.5) {
#       return("Unknown")
#     } else {
#       return(names(UCell_scores)[which.max(row)])
#     }
#   })
#   
#   # Remove "_UCell" from each element
#   UCell_scores$cell_type_UCell <- gsub("_UCell", "", UCell_scores$cell_type_UCell)
#   UCell_scores$cell_type_average <- UCell_scores$cell_type_UCell
#   
#   # group into broad categories
#   UCell_scores$cell_type_average <- gsub("Neurons .*", "Neurons", UCell_scores$cell_type_average)
#   UCell_scores$cell_type_average <- gsub("Microglia", "Myeloid", UCell_scores$cell_type_average)
#   
#   data <- AddMetaData(object = data, metadata = UCell_scores$cell_type_average, col.name = 'predicted_label_UCell')
#   data <- AddMetaData(object = data, metadata = UCell_scores$cell_type_UCell, col.name = 'predicted_label_UCell_precise')
#   
#   
#   # 3) Check for congruence
#   data$cell_type <- ifelse(
#     data$predicted_label_snRNAseq == 'Unknown', 
#     data$predicted_label_UCell, 
#     ifelse(
#       data$predicted_label_UCell == data$predicted_label_snRNAseq, 
#       data$predicted_label_UCell, 
#       ifelse(
#         (data$predicted_label_UCell == 'RG-like' & data$predicted_label_snRNAseq == 'NPC-like') |
#           (data$predicted_label_UCell == 'NPC-like' & data$predicted_label_snRNAseq == 'RG-like'), 
#         'NPC/RG-like', 
#         ifelse(
#           data$predicted_label_UCell == 'RG-like' & data$predicted_label_snRNAseq == 'RG-like', 
#           'RG-like', 
#           ifelse(
#             data$predicted_label_UCell == 'NPC-like' & data$predicted_label_snRNAseq == 'NPC-like', 
#             'NPC-like', 
#             ifelse(
#               data$predicted_label_UCell %in% c('Blood vessel', 'Neurons', 'Oligodendrocyte', 'Astrocyte', 'T-cells', 'Myeloid'), 
#               data$predicted_label_UCell, 
#               data$predicted_label_snRNAseq
#             )
#           )
#         )
#       )
#     )
#   )
#   
#   
#   
#   
#   
#   
#   
#   # 4) Plot labels
#   colors_metaprograms_Xenium2 <- c(colors_metaprograms_Xenium, "#8be8d7", "#EFC0C5")
#   names(colors_metaprograms_Xenium2) <- c(names(colors_metaprograms_Xenium), 'immune', "NPC/RG-like" )
#   
#   p1 <- ImageDimPlot(data, group.by = 'predicted_label_UCell',  size = 0.5, border.size = NA,  dark.background = F, cols = colors_metaprograms_Xenium)
#   p2 <- ImageDimPlot(data, group.by = 'predicted_label_snRNAseq', size = 0.5, border.size = NA,  dark.background = F, cols = colors_metaprograms_Xenium2)
#   p3 <- ImageDimPlot(data, group.by = 'cell_type', size = 0.5, border.size = NA,  dark.background = F, cols = colors_metaprograms_Xenium2) 
#   
#   plot_grid(p1, p2, p3, ncol = 3)
#   ggsave(file.path(plot_dir, paste0('1_Initial_cell_assignment_', metadata$SampleName[i], '.pdf')), width = 15, height = 5)
#   
#   # Add information about malignant or non-malignant
#   data$malignant <- ifelse(data$cell_type %in% c('AC-like', 'RG-like', 'MES-like','Cycling', 'NPC-like', 'INPC-like',  'eIN-like', 'NPC/RG-like'),
#                            "Malignant", "Non-malignant")  
#   number_cells <- as.data.frame(table(data$cell_type))
#   write_csv(number_cells, file.path(plot_dir, '0_number_cells_final_classification.csv'))
#   
#   # Perform niche analysis
#   data <- BuildNicheAssay(object = data, fov = "fov", group.by = "cell_type", niches.k = 4, neighbors.k = 20)
#   
#   
#   # save annotated object
#   qsave(data, file.path(seurat_obj_dir, paste0('seurat_obj_', metadata$SampleName[i], '.qs')))
# 
# 
# 
# }
# 
#   
#  