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
library(cowplot)
library(circlize)

## Function to preprocess the xenium sequencing (creates seurat object, find markers for different resolutions and score the signatures)
## @para: rawDir (sample sequencing raw dir)
## @para: outDir (output directory)
## @para: base_panel (brain base panel provided by 10x)
## @para: custom_panel (customized panel)
## @returns: save files in output directory

RunFullXenium <- function(smp, rawDir, outDir) {
  
  
  ## Preprocessing seurat objects
  message(green('-------------------------------------------------------'))
  message(blue(glue('Starting pre-processing for sample: {yellow$underline$bold(smp)}!!')))
  message(blue('[1/5] Loading Xenium raw files...'))
  
  suppressWarnings(suppressMessages(xenium.obj <- LoadXenium(rawDir, fov = "fov")))
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
  message(glue("{ncol(xenium.obj)} cells."))
  
  message(blue('[2/5] Running SCTransform-based normalization...'))
  xenium.obj <- SCTransform(xenium.obj, assay = "Xenium", vars.to.regress = c('nFeature_Xenium', 'nCount_Xenium'), vst.flavor = "v2", verbose = FALSE)
  
  message(blue('[3/5] Running PCA and UMAP reductions...'))
  xenium.obj <- RunPCA(xenium.obj, npcs = 50, features = rownames(xenium.obj), verbose = FALSE)
  xenium.obj <- RunUMAP(xenium.obj, dims = 1:optimizePCA(xenium.obj, 0.8), verbose = FALSE)
  
  message(blue('[4/5] Clustering...'))
  xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:optimizePCA(xenium.obj, 0.8), verbose = FALSE)
  xenium.obj <- FindClusters(xenium.obj, resolution = seq(0.1, 1, 0.1), verbose = FALSE)

  xenium.obj <- changeClusterNumbers(xenium.obj)
  
  message(blue('[5/5] Saving'))
  
  qsave(xenium.obj, glue('{outDir}/seurat.qs'))
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




XeniumCellAssignment <- function(seurat_object_path, Xenium_dir, marker_genes_list, 
                                 labels_priority_panel, labels_priority_scObject,
                                 malignant_cells,
                                 colors) {
  
  
  ## 
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[1/8] Loading Xenium seurat object ', basename(Xenium_dir))))
  
  data <- qread(file.path(Xenium_dir, '/seurat.qs'))
  
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[2/8] Loading reference seurat object')))
  
  seurat_object <- qread(seurat_object_path)
  
  
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[3/8] Transfer anchors from single cell data')))
  
  # same approach as in Seurat vignette
  Idents(seurat_object) <- 'cell_type'
  
  anchors <- FindTransferAnchors(reference = seurat_object, query = data, 
                                 normalization.method = "SCT", npcs = 50)
  
  # transfer labels
  predictions <- TransferData(
    anchorset = anchors,
    refdata = seurat_object$cell_type,
    #prediction.assay = TRUE,
    weight.reduction = data[["pca"]], dims = 1:50
  )
  
  
  # add metadata labels
  data <- AddMetaData(object = data, metadata = predictions$predicted.id, col.name = 'predicted_label_snRNAseq')
  data <- AddMetaData(object = data, metadata = predictions$prediction.score.max, col.name = 'predicted_label_snRNAseq_score')
  
  # If cell score < 0.4 --> unassigned
  print('Mean max prediction score is ')
  mean <- mean(predictions$prediction.score.max)
  mean
  data$predicted_label_snRNAseq <- ifelse(data$predicted_label_snRNAseq_score < 0.4, 'Unknown', data$predicted_label_snRNAseq)
  
  
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[4/8] Add scores for Xenium panel marker genes')))
  
  DefaultAssay(data) <- 'SCT'
  data <- AddModuleScore_UCell(data, features = marker_genes_list)
  U_Cell_signatures <- paste0(names(marker_genes_list), '_UCell')
  
  # retrieve highest UCell score
  UCell_scores <- data@meta.data[, U_Cell_signatures]
  
  # Assign cell identity based on highest value (if cell score < 0.6 --> unassigned)
  UCell_scores$cell_type_UCell <- apply(UCell_scores, 1, function(row) {
    max_score <- max(row)
    if (max_score < 0.6) {
      return("Unknown")
    } else {
      return(names(UCell_scores)[which.max(row)])
    }
  })
  
  
  # Remove "_UCell" from each element
  UCell_scores$cell_type_UCell <- gsub("_UCell", "", UCell_scores$cell_type_UCell)
  UCell_scores$cell_type_average <- UCell_scores$cell_type_UCell
  
  # group into broad categories
  UCell_scores$cell_type_average <- gsub("Neurons .*", "Neurons", UCell_scores$cell_type_average)
  UCell_scores$cell_type_average <- gsub("Microglia", "Myeloid", UCell_scores$cell_type_average)
  UCell_scores$cell_type_average <- gsub("Endothelial .*", "Endothelial", UCell_scores$cell_type_average)
  UCell_scores$cell_type_average <- gsub("Dendritic cell", "Myeloid", UCell_scores$cell_type_average)
  
  # for EPN (astro are neuronal-like cells)
  UCell_scores$cell_type_average <- gsub("Astrocyte", "Neuronal-like", UCell_scores$cell_type_average)
  
  
  data <- AddMetaData(object = data, metadata = UCell_scores$cell_type_average, col.name = 'predicted_label_UCell')
  data <- AddMetaData(object = data, metadata = UCell_scores$cell_type_UCell, col.name = 'predicted_label_UCell_precise')
  
  
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[5/8] Assign labels')))
  
  data$cell_type <- ifelse(
    data$predicted_label_UCell %in% labels_priority_panel, 
    data$predicted_label_UCell,
    ifelse(
      data$predicted_label_snRNAseq %in% labels_priority_scObject,
      data$predicted_label_snRNAseq,
      ifelse(
        data$predicted_label_snRNAseq == "Unknown",
        data$predicted_label_UCell, 
        data$predicted_label_UCell
      )
    )
  )
  
  print('Classified cell numbers')
  print(table(data$cell_type))
  
  
  
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[6/8] Check for frequency of assigned cell type within each Seurat cluster')))

  # Create a contingency table ------------------------

  # resolution chosen in metadata excel file
  print('Chosen resolution of ')
  print(metadata$Resolution)

  contingency_table <- table(data$cell_type, data@meta.data[, metadata$Resolution])
  print('Contingency table')
  contingency_table
  write_csv(as.data.frame(contingency_table), file.path(data_dir, '1_cell_type_per_cluster.csv'))

  column_sums <- colSums(contingency_table)
  percentage_table <- sweep(contingency_table, 2, column_sums, FUN = "/") * 100
  rounded_percentage_table <- round(percentage_table, 2)

  print('Contingency table in percentage')
  print(rounded_percentage_table)
  write_csv(as.data.frame(rounded_percentage_table), file.path(data_dir, '1_cell_type_per_cluster_percentage.csv'))


  # # reclassify cells based on most likely cluster (only if it is composed of more than 70%)
  # # get cell type with highest frequency for each cluster
  # max_row_names <- apply(contingency_table, 2, function(col) {
  #   rowname <- rownames(contingency_table)[which.max(col)]
  #   return(rowname)
  # })
  # 
  # print(max_row_names)
  # 
  # # Identify the cell types that dominate (>70%) in each cluster
  # proportion_table <- sweep(contingency_table, 2, column_sums, FUN = "/")
  # dominated_clusters <- apply(proportion_table, 2, function(col) {
  #   max_prop <- max(col)
  #   return(max_prop > 0.7)  # Check if any cell type has >70% proportion
  # })
  # 
  # # reassign cells
  # Idents(data) <- metadata$Resolution
  # new.cluster.ids <- max_row_names
  # names(new.cluster.ids) <- levels(data)
  # data <- RenameIdents(data, new.cluster.ids)
  # data[["cell_type_fin"]] <- Idents(object = data)
  # data$cell_type_fin2 <- ifelse(data@meta.data[, metadata$Resolution] %in% dominated_clusters, data$cell_type_fin, data$cell_type)

  
  # Display the final results
  number_cells <- as.data.frame(table(data$cell_type))
  print(number_cells)
  write_csv(number_cells, file.path(data_dir, '2_number_cells_cell_type.csv'))

  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[7/8] Plot spatial maps')))
  
  p1 <- ImageDimPlot(data, group.by = 'predicted_label_UCell',  size = 0.5, border.size = NA,  dark.background = F, cols = colors)
  p2 <- ImageDimPlot(data, group.by = 'predicted_label_snRNAseq', size = 0.5, border.size = NA,  dark.background = F, cols = colors)
  p3 <- ImageDimPlot(data, group.by = 'cell_type', size = 0.5, border.size = NA,  dark.background = F, cols = colors) 
  p4 <- ImageDimPlot(data, group.by = metadata$Resolution, size = 0.5, border.size = NA,  dark.background = F) 
  
  plot_grid(p1, p2, p3, p4, ncol = 2)
  ggsave(file.path(plot_dir, paste0('1_SpatialMaps_', basename(Xenium_dir), '.pdf')), width = 10, height = 10)
  
  
  # Plot dotplot ZFTA-RELA expression
  DotPlot(data, assay = 'SCT', features = "ZFTA-RELA-Fusion1", group.by = 'cell_type')
  ggsave(file.path(plot_dir, paste0('2_DotPlot_ZR_', basename(Xenium_dir), '.pdf')), width = 4, height = 5)
  
  
  # Add information about malignant or non-malignant
  data$malignant <- ifelse(data$cell_type %in% malignant_cells, "Malignant", "Non-malignant")  
  
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[8/8] Save annotated object')))
  
  qsave(data, file.path(data_dir, paste0('seurat_obj_', basename(Xenium_dir), '.qs')))
  
  
  # export cell ID
  cell_id <- data.frame(rownames(data@meta.data), data@meta.data$cell_type)
  
  # rename column names
  colnames(cell_id)[colnames(cell_id) == "rownames.data.meta.data."] <- "cell_id"
  colnames(cell_id)[colnames(cell_id) == "data.meta.data.cell_type"] <- "group"
  # save
  write_csv(cell_id, file.path(data_dir, paste0('cell_ID_', basename(Xenium_dir), ".csv")))
  
  
}




NicheAnalysis <- function(input_dir, output_dir, 
                          niches_min, niches_max,
                          n_neighbors) {
  
  
  message(green('-------------------------------------------------------'))
  message(blue('[1/3] Loading Xenium seurat object with annotation'))
  
  data <- qread(file.path(input_dir, paste0('seurat_obj_', basename(Xenium_dir), '.qs')))
  
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[2/3] Performing niche analysis')))
  
  for (i in seq(niches_min, niches_max, 1)){
    print(glue('Calculating niche {i}')) 
    data <- BuildNicheAssay(object = data, fov = "fov", group.by = "cell_type", niches.k = i, 
                            neighbors.k = n_neighbors, cluster.name = glue('niche_{i}'))
    print(glue('Finished calculating niche {i}')) 
  }
  
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[3/3] Saving object')))
  
  qsave(data, file.path(output_dir, paste0('seurat_obj_', basename(Xenium_dir), '.qs')))
  
  
}






PlotXeniumMetaprograms <- function(cellid_dir, SampleName_vector, SampleID_vector, SampleID_order,
                                   order_metaprograms, colors_metaprogram, replicate = FALSE) {
  
  message(green('-------------------------------------------------------'))
  message(blue('[1/4] Loading cellIDs '))
  
  metaprogram_frequency <- list()
  
  for (i in seq_along(SampleName_vector)) { 
    # read data
    cellID <- suppressWarnings(suppressMessages(read_csv(file.path(cellid_dir, paste0('cell_ID_', SampleName_vector[i], '.csv')))))
    print(paste0("Reading ", SampleName_vector[i]))
    
    # transform in data table with frequency
    metaprogram_frequency[[i]] <- cellID %>%
      group_by(group) %>%
      summarise(n = n()) %>%
      mutate(freq = (n/sum(n)*100))
  }
  
  names(metaprogram_frequency) <- SampleName_vector  
  
  message(green('-------------------------------------------------------'))
  message(blue('[2/4] Creating df'))
  
  # bind rows together and transform into df
  metaprogram_proportion <- bind_rows(metaprogram_frequency, .id = "Xenium_Region")
  metaprogram_proportion <- as.data.frame(metaprogram_proportion)
  
  message(green('-------------------------------------------------------'))
  message(blue('[3/4] Adding sample names and reordering'))
  
  # add sample name
  matching_indices <- match(metaprogram_proportion$Xenium_Region, SampleName_vector)
  metaprogram_proportion$SampleID <- SampleID_vector[matching_indices]
  
  # reorder so that it is the same as in the paper
  metaprogram_proportion$group <- factor(metaprogram_proportion$group, 
                                         levels = order_metaprograms)
  
  # reorder by SampleID
  if (!is.null(SampleID_order)) {
    metaprogram_proportion$SampleID <- factor(metaprogram_proportion$SampleID, levels = SampleID_order)
    metaprogram_proportion <- metaprogram_proportion %>% arrange(SampleID)
  } else {
    metaprogram_proportion <- metaprogram_proportion %>% arrange(SampleID)
  }
  
  
  # Check if averaging by SampleID is required
  if (replicate) {
    metaprogram_proportion <- metaprogram_proportion %>%
      group_by(group, SampleID) %>%
      summarise(avg_frequency = mean(freq), .groups = 'drop')
    
    # Modify plot input to match averaged data
    plot_data <- metaprogram_proportion
    x_axis <- "SampleID"  # Update x-axis to reflect unique SampleIDs
    y_var <- "avg_frequency"
  } else {
    plot_data <- metaprogram_proportion
    x_axis <- "Xenium_Region"  # Default x-axis for individual replicates
    y_var <- "freq"
  }
  
  message(green('-------------------------------------------------------'))
  message(blue('[4/4] Plotting'))
  
  # Plot a stacked bar plot
  ggplot(plot_data, aes_string(x = x_axis, y = y_var, fill = "group")) +
    scale_fill_manual(values = colors_metaprogram) +
    geom_bar(stat = "identity", position = "fill", color = "black") +
    labs(y = 'Proportion', x = '') + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_blank())
}





PlotSpatialMapsZoom <- function(SampleName, input_dir, group_display, colors, size_dot, geom_rect = F,  x_min, x_max, y_min, y_max, ratio) {
  
  message(green('-------------------------------------------------------'))
  message(blue('[1/4] Loading seurat object'))
  
  data <- qread(file.path(input_dir, paste0('seurat_obj_', SampleName, '.qs')))
  
  
  message(green('-------------------------------------------------------'))
  message(blue('[2/4] Plotting spatial map'))
  
  p1 <- ImageDimPlot(data, group.by = group_display,  fov = "fov", cols = colors,  border.size = NA, size = size_dot,
                     dark.background = T, axes = TRUE) + 
    geom_rect(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, color = 'black', fill = NA, lwd = 0.5) + seurat_theme()
  
  message(green('-------------------------------------------------------'))
  message(blue('[3/4] Plotting zoomed map'))
  
  data[["zoom_1"]] <- Crop(data[["fov"]], x = c(x_min, x_max), y = c(y_min, y_max), coords = "plot")
  DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
  p2 <- ImageDimPlot(data, fov = "zoom_1", size = 1, group.by = group_display, border.color = "white", border.size = 0.1, 
                     cols = colors) + seurat_theme() + theme(legend.position = "none")
  
  message(green('-------------------------------------------------------'))
  message(blue('[4/4] Plotting and saving both images'))
  plot_grid(p1, p2, ncol = 2, rel_widths = ratio)
  
}






# spatial coherence
prepareData <- function(cell_feature_matrix, cells_info, scale = 1, filter_genes = 10, filter_cells = 10) {
  spatial <- cells_info %>% dplyr::select(c("cell_id", "x_centroid", "y_centroid")) %>%
    column_to_rownames("cell_id") %>%
    rename_all(~c("imagecol", "imagerow"))
  
  spatial[["imagecol"]] = spatial[["imagecol"]] * scale
  spatial[["imagerow"]] = spatial[["imagerow"]] * scale
  
  cell_feature_matrix <- cell_feature_matrix[rowSums(cell_feature_matrix) >= filter_genes, ]
  
  cell_feature_matrix <- cell_feature_matrix[, colSums(cell_feature_matrix) >= filter_cells]
  
  spatial <- spatial[colnames(cell_feature_matrix), ]
  
  return(list(spatial = spatial, cell_feature_matrix = cell_feature_matrix))
}


norm_data <- function(cell_feature_matrix) {
  counts_per_cell <- rowSums(cell_feature_matrix)
  counts_greater_than_zero <- counts_per_cell[counts_per_cell > 0]
  after <- median(counts_greater_than_zero)
  counts_per_cell <- counts_per_cell / after
  mat <- cell_feature_matrix/counts_per_cell
  return(mat)
}


grid_spatial <- function(norm_data, spatial, variable, nbins) {
  #grid_spatial <- function(norm_data, spatial, variable, nbins = 200) {
  
  message(glue('{nbins} by {nbins} has this many spots: {nbins*nbins}'))
  
  n_squares = nbins * nbins
  cell_bcs = colnames(norm_data)
  xs <- as.integer(unname(spatial$imagecol))
  ys <- as.integer(unname(spatial$imagerow))
  
  h2 <- hist2d(xs, ys, nbins = nbins, show = F)
  grid_counts <- h2$counts
  xedges <- h2$x.breaks
  yedges <- h2$y.breaks
  
  grid_expr <- matrix(data = 0, nrow = n_squares, ncol = nrow(norm_data))
  grid_coords <- matrix(data = 0, nrow = n_squares, ncol = 2)
  grid_cell_counts <- rep(0, n_squares)
  
  cell_labels <- as.character(variable)
  cell_set = sort(unique(cell_labels))
  cell_info <- matrix(data = 0, nrow = n_squares, ncol = length(cell_set))
  
  pb <- txtProgressBar(min = 0, max = nbins, style = 3, width = 50, char = "=")
  for (i in 1:nbins) {
    x_left <- xedges[i]
    x_right <- xedges[i + 1]
    for (j in 1:nbins) {
      n <- ((i-1) * nbins) + j
      
      y_down <- yedges[j]
      y_up <- yedges[j + 1]
      
      grid_coords[n, 1] = (x_right + x_left) / 2
      grid_coords[n, 2] = (y_up + y_down) / 2
      
      
      # Now determining the cells within the gridded area #
      #if ((i != nbins-1) & (j == nbins-1)) { # top left corner
      if ((i != nbins-1) & (j == nbins)) { # top left corner
        x_true = (xs >= x_left) & (xs < x_right)
        y_true = (ys <= y_up) & (ys > y_down)
        #} else if ((i == nbins - 1) & (j != nbins)) { # bottom right corner
      } else if ((i == nbins - 1) & (j != nbins-1)) { # bottom right corner
        x_true = (xs > x_left) & (xs <= x_right)
        y_true = (ys < y_up) & (ys >= y_down)
      } else {   # average case
        x_true = (xs >= x_left) & (xs < x_right)
        y_true = (ys < y_up) & (ys >= y_down)
      }
      
      cell_bool = x_true & y_true
      grid_cells = cell_bcs[cell_bool]
      
      grid_cell_counts[n] = length(grid_cells)
      
      # Summing the expression across these cells to get the grid expression #
      if (length(grid_cells) > 0) {
        if (length(grid_cells) != 1) {
          grid_expr[n,] = rowSums(norm_data[, cell_bool])
        } else {
          grid_expr[n,] = norm_data[, cell_bool]
        }
        
      }
      
      # If we have cell type information, will record #
      if (length(grid_cells) > 0) {
        grid_cell_types <- cell_labels[cell_bool]
        
        tmp <- c()
        for (ct in cell_set) {
          tmp <- c(tmp, length(which(grid_cell_types == ct)) / length(grid_cell_types))
        }
        cell_info[n, ] <- tmp
      }
      
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  
  grid_expr <- grid_expr %>% as.data.frame()
  rownames(grid_expr) <- glue('grid_{1:n_squares}')
  colnames(grid_expr) <- rownames(norm_data)
  
  grid_data <- list(expr = grid_expr,
                    image_col = grid_coords[,1],
                    image_row = grid_coords[,2],
                    n_cells = grid_cell_counts,
                    grid_coords = grid_coords)
  
  cell_info <- cell_info %>% as.data.frame()
  rownames(cell_info) <- rownames(grid_expr)
  colnames(cell_info) <- cell_set
  
  max_indices <- apply(cell_info, 1, function(x) names(which.max(x)))
  max_indices <- max_indices[grid_data$n_cells > 0]
  max_indices <- as.character(max_indices)
  
  grid_data$expr <- grid_data$expr[grid_data$n_cells > 0, ]
  grid_data$image_col <- grid_data$image_col[grid_data$n_cells > 0]
  grid_data$image_row <- grid_data$image_row[grid_data$n_cells > 0]
  
  dt <- data.frame(x = grid_data$image_col, y = grid_data$image_row, Metaprogram = max_indices)
  
  return(dt)
}


coherence_score <- function(grid_df, variable) {
  mat <- grid_df %>% mutate(x = as.character(x)) %>% mutate(y = as.character(y))
  mat <- reshape2::acast(mat, y~x, value.var = variable) %>% as.data.frame()
  rownames(mat) <- 1:nrow(mat)
  colnames(mat) <- 1:ncol(mat)
  
  programs <- unique(grid_df[[variable]])
  results <- list()
  for (program in programs) {
    tmp <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = list(rownames(mat), colnames(mat)))
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        
        if (!is.na(mat[i,j])) {
          
          if (mat[i,j] == program) {
            
            if (length(mat[i,j+1]) != 0) { ## right
              if (!is.na(mat[i,j+1])) {
                if (mat[i,j+1] == program) {
                  tmp[i,j+1] <- tmp[i,j+1] + 1
                }
              }
            }
            
            if (length(mat[i,j-1]) != 0) { ## left
              if (!is.na(mat[i,j-1])) {
                if (mat[i,j-1] == program) {
                  tmp[i,j-1] <- tmp[i,j-1] + 1
                }
              }
            }
            
            if (length(mat[i-1,j]) != 0) { ## top
              if (!is.na(mat[i-1,j])) {
                if (mat[i-1,j] == program) {
                  tmp[i-1,j] <- tmp[i-1,j] + 1
                }
              }
            }
            
            if (length(mat[i+1,j]) != 0) { ## bottom
              if (!is.na(mat[i+1,j])) {
                if (mat[i+1,j] == program) {
                  tmp[i+1,j] <- tmp[i+1,j] + 1
                }
              }
            }
            
            if (length(mat[i-1,j+1]) != 0) { ## top righ diagonal
              if (!is.na(mat[i-1,j+1])) {
                if (mat[i-1,j+1] == program) {
                  tmp[i-1,j+1] <- tmp[i-1,j+1] + 1
                }
              }
            }
            
            if (length(mat[i+1,j+1]) != 0) { ## bottom righ diagonal
              if (!is.na(mat[i+1,j+1])) {
                if (mat[i+1,j+1] == program) {
                  tmp[i+1,j+1] <- tmp[i+1,j+1] + 1
                }
              }
            }
            
            if (length(mat[i-1,j-1]) != 0) { ## top left diagonal
              if (!is.na(mat[i-1,j-1])) {
                if (mat[i-1,j-1] == program) {
                  tmp[i-1,j-1] <- tmp[i-1,j-1] + 1
                }
              }
            }
            
            if (length(mat[i+1,j-1]) != 0) { ## bottom left diagonal
              if (!is.na(mat[i+1,j-1])) {
                if (mat[i+1,j-1] == program) {
                  tmp[i+1,j-1] <- tmp[i+1,j-1] + 1
                }
              }
            }
            
          } else { next }
        }
      }
      #message(i)
    }
    results[[program]] <- tmp
    #message(glue('{program} Done!'))
  }
  
  suppressWarnings(counts <- grid_df %>% count_(variable))
  
  tmp <- data.frame()
  for (metaprogram in counts[[variable]]) {
    if (sum(results[[metaprogram]]) == 0) {
      tmp <- rbind(tmp, data.frame(metaprogram = metaprogram, average = 0))
    } else {
      tmp <- rbind(tmp, data.frame(metaprogram = metaprogram, average = (sum(results[[metaprogram]]))/counts$n[counts[[variable]] == metaprogram]))
    }
  }
  
  return(list(coherence_score = mean(tmp$average), coherence_score_program = tmp, results_df = results))
}



gridding_coherence_density <- function(gridding, plot_dir, color_coherence_density, colors_metaprograms_Xenium) {
  p1 <- ggplot(gridding, aes(x = y, y = x, color = Metaprogram)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = colors_metaprograms_Xenium) +
  theme_void() +
  guides(colour = guide_legend(override.aes = list(size = 4)))
  p1
  ggsave(file.path(plot_dir, paste0('1_Gridding_', SampleName_arg, '.pdf')), width = 9, height = 8)
  
  
  mat <- gridding %>% mutate(x = as.character(x)) %>% mutate(y = as.character(y))
  mat <- reshape2::acast(mat, y~x, value.var = 'Metaprogram') %>% as.data.frame()
  
  # Sum coherence score across all MPs for each position (regions with high coherence will score high)
  summed_matrix <- Reduce(`+`, coherence$results_df)
  
  #change x/y values with actual positions
  rownames(summed_matrix) <- rownames(mat) 
  colnames(summed_matrix) <- colnames(mat) 
  
  summed_matrix <- reshape2::melt(summed_matrix)
  
  df <- summed_matrix #
  colnames(df) <- c("x", "y", "value")  # Rename columns
  
  # Convert x and y to numeric if necessary
  df$x <- as.numeric(as.character(df$x))
  df$y <- as.numeric(as.character(df$y))
  df$value <- as.numeric(as.character(df$value))
  
  
  
  
  # Plot density of spatial cohernece score
  p2 <- ggplot(df, aes(x = x, y = y, color = factor(value))) +
    geom_point(size = 0.2) +
    scale_color_manual(values = color_coherence_density) +
    theme_void() +
    guides(colour = guide_legend(override.aes = list(size = 4)))
  p2
  ggsave(file.path(plot_dir, paste0('2_Gridding_density_coherence_score_', SampleName_arg, '.pdf')), width = 9, height = 8)
  
  # plot combined
  plot_grid(p1, p2, ncol = 2)
  ggsave(file.path(plot_dir, paste0('3_Gridding_results_', SampleName_arg, '.pdf')), width = 18, height = 8)

}






plot_coherence_score <- function(FileName, cellid_dir, coherence_dir, plot_dir, 
                                 colors_metaprograms_Xenium, col_sampling) {
  # read  proportion of each cell type in each tumor ----------------------
  
  message(green('-------------------------------------------------------'))
  message(blue('[1/4] Reading metaprogram proportions'))
  
  metaprogram_frequency <- list()
  
  for (i in seq_along(FileName)) { 
    # read data
    cellID <- suppressWarnings(suppressMessages(read_csv(file.path(cellid_dir, paste0('cell_ID_', FileName[i], '.csv')))))
    print(paste0("Reading ", FileName[i]))
    
    # transform in data table with frequency
    metaprogram_frequency[[i]] <- cellID %>%
      group_by(group) %>%
      summarise(n = n()) %>%
      mutate(freq = (n/sum(n)*100))
  }
  
  names(metaprogram_frequency) <- FileName  
  
  
  # bind rows together and transform into df
  metaprogram_proportion <- bind_rows(metaprogram_frequency, .id = "Xenium_Region")
  metaprogram_proportion <- as.data.frame(metaprogram_proportion)
  
  
  # reorganize dataframe with metaprogram proportion
  metaprogram_proportion <- metaprogram_proportion[ , c("Xenium_Region", "group", "freq")]
  
  proportion2 <- metaprogram_proportion %>% 
    complete(nesting(Xenium_Region, group), fill = list(freq = 0))  
  
  proportion3 <- cast(proportion2, Xenium_Region~group, mean)
  
  # export dataframe
  write_csv(as.data.frame(proportion3), file.path(coherence_dir, '4_Program_frequency.csv'))
  
  
  message(green('-------------------------------------------------------'))
  message(blue('[2/4] Read spatial coherence score'))
  
  spatial_coherence <- list()
  average_spatial_coherence_df <- NULL
  average_spatial_coherence_program <- list()
  
  for (i in seq_along(FileName)) { 
    spatial_coherence[[i]] <- qread(file.path(coherence_dir, paste0('results_', FileName[i], '.qs')))
    print(paste0("Reading spatial coherence score ", FileName[i]))
    
    average_spatial_coherence_df[i] <-spatial_coherence[[i]]$coherence_score
    average_spatial_coherence_program[[i]] <- spatial_coherence[[i]]$coherence_score_program
    print(paste0("Extracting global spatial coherence score ", FileName[i]))
  }
  
  names(spatial_coherence) <- names(average_spatial_coherence_df) <- names(average_spatial_coherence_program) <- FileName  
  
  
  message(green('-------------------------------------------------------'))
  message(blue('[3/4] Reorganize average spatial coherence per sample'))
  
  
  
  # reorganize and export average_spatial_coherence_df ------------------------------------------
  average_spatial_coherence_df <- as.data.frame(average_spatial_coherence_df)
  colnames(average_spatial_coherence_df)[1] <- 'average_spatial_coherence'
  
  # calculate scaled score
  min <- min(average_spatial_coherence_df$average_spatial_coherence)
  max <- max(average_spatial_coherence_df$average_spatial_coherence)
  range <- max - min
  average_spatial_coherence_df$scaled_spatial_coherence <- (average_spatial_coherence_df$average_spatial_coherence - min)/range
  
  # calculate average score per sample
  average_spatial_coherence_df$SampleName <- FileName
  average_spatial_coherence_df$SampleID <- metadata$SampleID
  average_spatial_coherence_df$Source <- metadata$Source
  
  average_spatial_coherence_average <- average_spatial_coherence_df %>% 
    dplyr::group_by(Source, SampleID) %>%
    dplyr::summarise(mean = mean(scaled_spatial_coherence))
  
  # re-scale
  min <- min(average_spatial_coherence_average$mean)
  max <- max(average_spatial_coherence_average$mean)
  range <- max - min
  average_spatial_coherence_average$scaled_spatial_coherence <- (average_spatial_coherence_average$mean - min)/range
  

  # export dataframe
  print('Saving spatial coherence score per sample')
  write_csv(as.data.frame(average_spatial_coherence_df), file.path(coherence_dir, '5_Spatial_coherence_score_by_sample.csv'))

  # export dataframe
  write_csv(as.data.frame(average_spatial_coherence_average), file.path(coherence_dir, '6_Spatial_coherence_score_average_by_sample.csv'))
  

  
  
  
  
  message(green('-------------------------------------------------------'))
  message(blue('[4/4] Reorganize average spatial coherence per program'))
  
  
  
  # reorganize and export average_spatial_coherence_program (per program) ------------------------------------------
  average_spatial_coherence_program
  
  # add column with sample name
  for (i in seq_along(average_spatial_coherence_program)) {
    average_spatial_coherence_program[[i]]$SampleName <- names(average_spatial_coherence_program)[i]
  }
  
  # merge information
  average_spatial_coherence_program_df <- average_spatial_coherence_program[[1]]
  
  for (i in 2:length(average_spatial_coherence_program)) {
    average_spatial_coherence_program_df <- rbind(average_spatial_coherence_program_df, average_spatial_coherence_program[[i]])
  }
  
  # remove unassigend cells 
  average_spatial_coherence_program_df <- average_spatial_coherence_program_df[average_spatial_coherence_program_df$metaprogram != "Unknown", ] 
  
  # calculate scaled score
  min <- min(average_spatial_coherence_program_df$average)
  max <- max(average_spatial_coherence_program_df$average)
  range <- max - min
  average_spatial_coherence_program_df$scaled_spatial_coherence <- (average_spatial_coherence_program_df$average - min)/range
  
  # add information on sample ID
  matching_indices <- match(average_spatial_coherence_program_df$SampleName, metadata$SampleName)
  average_spatial_coherence_program_df$SampleID <- metadata$SampleID[matching_indices]
  average_spatial_coherence_program_df$Source <- metadata$Source[matching_indices]
  
  
  
  # export dataframe
  write_csv(as.data.frame(average_spatial_coherence_program_df), file.path(coherence_dir, '7_Spatial_coherence_score_averge_by_MP.csv'))
  
  
  
  
  
  message(green('-------------------------------------------------------'))
  message(blue('[5/5] Plot scores'))
  
  
  # reorder
  average_spatial_coherence_program_df$metaprogram <- factor(average_spatial_coherence_program_df$metaprogram, 
                                                             levels = order_metaprograms)
  
  # Plot scores by metaprogram   
  ggplot(average_spatial_coherence_program_df, aes(x = reorder(metaprogram, scaled_spatial_coherence), y = scaled_spatial_coherence)) +
    geom_boxplot(fatten = NULL, outlier.shape = NA, width = 0.5, color = 'black') +
    geom_point(aes(group = SampleName, fill = metaprogram,  color = metaprogram), size = 2, shape = 21, stroke = 0, position = position_dodge(0.2)) +
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.4, size = 1, linetype = "solid") + 
    scale_fill_manual(values = colors_metaprograms_Xenium) +
    scale_color_manual(values = colors_metaprograms_Xenium) +
    labs(x = "Metaprogram", y = "Spatial coherence score") +
    theme_minimal() + theme(panel.border = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            plot.background = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                            axis.text.y = element_text(size=12, colour="black"),
                            axis.title=element_text(size=12),
                            plot.title = element_text(size=10, face="bold")) 
  ggsave(file.path(plot_dir, "8_coherence_by_MP.pdf"), width=6, height=4)
  
  
  
  
  # Plot scores by sample   
  ggplot(average_spatial_coherence_df, aes(x = reorder(SampleID, scaled_spatial_coherence), y = scaled_spatial_coherence)) +
    geom_boxplot(fatten = NULL, outlier.shape = NA, width = 0.5, color = 'black') +
    geom_point(aes(group = SampleName, fill = Source, color = Source), size = 4, shape = 21, stroke = 0, position = position_dodge(0.2)) +
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.4, size = 1, linetype = "solid") + 
    scale_color_manual(values = col_sampling) +
    scale_fill_manual(values = col_sampling) +
    labs(x = "Sample", y = "Spatial coherence score") +
    theme_minimal() + theme(panel.border = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            plot.background = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                            axis.text.y = element_text(size=12, colour="black"),
                            axis.title=element_text(size=12),
                            plot.title = element_text(size=10, face="bold")) 
  ggsave(file.path(plot_dir, "9_coherence_by_sample_color_source.pdf"), width=6, height=4)
  
  
  # Plot scores by sample   
  ggplot(average_spatial_coherence_df, aes(x = reorder(SampleID, scaled_spatial_coherence), y = scaled_spatial_coherence)) +
    geom_boxplot(fatten = NULL, outlier.shape = NA, width = 0.5, color = 'black') +
    geom_point(aes(group = SampleName, fill = SampleID, color = SampleID), size = 4, shape = 21, stroke = 0, position = position_dodge(0.2)) +
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.4, size = 1, linetype = "solid") + 
    #scale_color_manual(values = col_sampling) +
    #scale_fill_manual(values = col_sampling) +
    labs(x = "Sample", y = "Spatial coherence score") +
    theme_minimal() + theme(panel.border = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            plot.background = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                            axis.text.y = element_text(size=12, colour="black"),
                            axis.title=element_text(size=12),
                            plot.title = element_text(size=10, face="bold")) 
  ggsave(file.path(plot_dir, "10_coherence_by_sample.pdf"), width=6, height=4)
  
  
  # Plot scores by source   
  average_spatial_coherence3 <- average_spatial_coherence_df %>% 
    dplyr::group_by( Source, SampleID) %>%
    dplyr::summarise(mean = mean(scaled_spatial_coherence),
                     sem = sd(scaled_spatial_coherence) / sqrt(n()))
  
  
  ggplot(average_spatial_coherence3, aes(x = reorder(Source, mean), y = mean)) +
    #geom_violin(width = 1, color = 'black') + 
    geom_boxplot(fatten = NULL, outlier.shape = NA, width = 0.5, color = 'black') +
    geom_point(aes(group = SampleID, fill = Source, color = Source), size = 4, shape = 21, stroke = 0, position = position_dodge(0.2)) +
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.4, size = 1, linetype = "solid") + 
    scale_fill_manual(values = col_sampling) +
    scale_color_manual(values = col_sampling) +
    labs(x = "Metaprogram", y = "Spatial coherence score") +
    theme_minimal() + theme(panel.border = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            plot.background = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                            axis.text.y = element_text(size=12, colour="black"),
                            axis.title=element_text(size=12),
                            plot.title = element_text(size=10, face="bold")) +
    ggpubr::stat_compare_means(method = "t.test", size = 4)
  ggsave(file.path(plot_dir, "11_coherence_by_source_average_technical.pdf"), width=4, height=5)
  
  
  ggplot(average_spatial_coherence_df, aes(x = reorder(Source, scaled_spatial_coherence), y = scaled_spatial_coherence)) +
    geom_boxplot(fatten = NULL, outlier.shape = NA, width = 0.5, color = 'black') +
    geom_point(aes(group = SampleID, fill = Source, color = Source), size = 4, shape = 21, stroke = 0, position = position_dodge(0.2)) +
    stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.4, size = 1, linetype = "solid") + 
    scale_fill_manual(values = col_sampling) +
    scale_color_manual(values = col_sampling) +
    labs(x = "Metaprogram", y = "Spatial coherence score") +
    theme_minimal() + theme(panel.border = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            plot.background = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                            axis.text.y = element_text(size=12, colour="black"),
                            axis.title=element_text(size=12),
                            plot.title = element_text(size=10, face="bold")) +
    ggpubr::stat_compare_means(method = "t.test", size = 4)
  ggsave(file.path(plot_dir, "12_coherence_by_source.pdf"), width=4, height=5)

}







CalculateLinearRegression <- function(FileName, cellid_dir, coherence_dir, plot_dir) {
 
   # read  proportion of each cell type in each tumor ----------------------
  
  message(green('-------------------------------------------------------'))
  message(blue('[1/6] Reading metaprogram proportions'))
  
  metaprogram_frequency <- list()
  
  for (i in seq_along(FileName)) { 
    # read data
    cellID <- suppressWarnings(suppressMessages(read_csv(file.path(cellid_dir, paste0('cell_ID_', FileName[i], '.csv')))))
    print(paste0("Reading ", FileName[i]))
    
    # transform in data table with frequency
    metaprogram_frequency[[i]] <- cellID %>%
      group_by(group) %>%
      summarise(n = n()) %>%
      mutate(freq = (n/sum(n)*100))
  }
  
  names(metaprogram_frequency) <- FileName  
  
  
  # bind rows together and transform into df
  metaprogram_proportion <- bind_rows(metaprogram_frequency, .id = "Xenium_Region")
  metaprogram_proportion <- as.data.frame(metaprogram_proportion)
  
  
  # reorganize dataframe with metaprogram proportion
  metaprogram_proportion <- metaprogram_proportion[ , c("Xenium_Region", "group", "freq")]
  
  proportion2 <- metaprogram_proportion %>% 
    complete(nesting(Xenium_Region, group), fill = list(freq = 0))  
  
  proportion3 <- cast(proportion2, Xenium_Region~group, mean)

  
  message(green('-------------------------------------------------------'))
  message(blue('[2/4] Read spatial coherence score'))
  
  spatial_coherence <- list()
  average_spatial_coherence_df <- NULL
  average_spatial_coherence_program <- list()
  
  for (i in seq_along(FileName)) { 
    spatial_coherence[[i]] <- qread(file.path(coherence_dir, paste0('results_', FileName[i], '.qs')))
    print(paste0("Reading spatial coherence score ", FileName[i]))
    
    average_spatial_coherence_df[i] <-spatial_coherence[[i]]$coherence_score
    average_spatial_coherence_program[[i]] <- spatial_coherence[[i]]$coherence_score_program
    print(paste0("Extracting global spatial coherence score ", FileName[i]))
  }
  
  names(spatial_coherence) <- names(average_spatial_coherence_df) <- names(average_spatial_coherence_program) <- FileName  
  

  
  message(green('-------------------------------------------------------'))
  message(blue('[3/4] Reorganize average spatial coherence per sample'))
  
  
  
  # reorganize and export average_spatial_coherence_df ------------------------------------------
  average_spatial_coherence_df <- as.data.frame(average_spatial_coherence_df)
  colnames(average_spatial_coherence_df)[1] <- 'average_spatial_coherence'
  
  # calculate scaled score
  min <- min(average_spatial_coherence_df$average_spatial_coherence)
  max <- max(average_spatial_coherence_df$average_spatial_coherence)
  range <- max - min
  average_spatial_coherence_df$scaled_spatial_coherence <- (average_spatial_coherence_df$average_spatial_coherence - min)/range
  
  # calculate average score per sample
  average_spatial_coherence_df$SampleName <- FileName
  average_spatial_coherence_df$SampleID <- metadata$SampleID
  average_spatial_coherence_df$Source <- metadata$Source
  
 
  head(average_spatial_coherence_df)
  

  message(green('-------------------------------------------------------'))
  message(blue('[5/5] Calculate linear regression'))
  
  linear_regression_df <- average_spatial_coherence_df
 
   # add information on frequency MP
  matching_indices <- match(linear_regression_df$SampleName, proportion3$Xenium_Region)

  for (col_name in colnames(proportion3)[2:length(colnames(proportion3))]) {
    linear_regression_df[[col_name]] <- proportion3[[col_name]][matching_indices]
  }
  

  linear_regression_df[is.na(linear_regression_df)] <- 0
  
  linear_regression_df$Unknown <- NULL
  
  head(linear_regression_df)
  write.csv(as.data.frame(linear_regression_df), file.path(plot_dir, '13_Linear_regression_input.csv'))
  
  
  
  message(green('-------------------------------------------------------'))
  message(blue('[5/5] Build linear regression'))
  
  colnames(linear_regression_df) <- gsub("-", "", colnames(linear_regression_df) )
  
  predictors <- colnames(linear_regression_df)[-c(1:5)]

  formula <- as.formula(paste("scaled_spatial_coherence ~", paste(predictors, collapse = " + ")))

  # Fit the linear model
  model <- lm(formula, data = linear_regression_df)

  summary(model)
  summary(model)$coefficient
  write.csv(as.data.frame(summary(model)$coefficient), file.path(plot_dir, '14_Summary_linear_regression.csv'))
  
  linear_regression_df[is.na(linear_regression_df)] <- 0
  
  
  # Fit linear regression for each MP
  message(green('-------------------------------------------------------'))
  message(blue('[5/6] Build linear regression for each MP'))
  
  
  # Initialize a dataframe to store regression results
  results_df <- data.frame(
    Predictor = character(),
    Estimate = numeric(),
    StdError = numeric(),
    tValue = numeric(),
    pValue = numeric(),
    stringsAsFactors = FALSE
  )
  
  plot_list <- list()
  
  # Loop over each predictor
  for (predictor in predictors) {
    # Dynamically create the formula
    formula <- as.formula(paste("scaled_spatial_coherence ~", predictor))
    
    # Fit the linear model
    model <- lm(formula, data = linear_regression_df)
    
    # Extract coefficients summary
    summary_model <- summary(model)
    coefficients <- summary_model$coefficients
    
    # Extract R-squared and p-value
    r_squared <- summary_model$r.squared
    p_value <- coefficients[2, "Pr(>|t|)"]  # p-value of the predictor term
    
    # Append results to dataframe
    results_df <- rbind(
      results_df,
      data.frame(
        Predictor = predictor,
        Estimate = coefficients[2, "Estimate"],
        StdError = coefficients[2, "Std. Error"],
        tValue = coefficients[2, "t value"],
        pValue = p_value
      )
    )
    
    # Create and save the plot
    plot_list[[predictor]] <- ggplot(linear_regression_df, aes_string(x = "scaled_spatial_coherence", y = predictor)) +
      labs(x = "Scaled spatial coherence score", y = predictor) +
      geom_smooth(method = 'lm', se = TRUE, color = 'black') +
      geom_point(size = 4) +
      theme_minimal() +
      labs(subtitle = paste0('R² = ', round(r_squared, 3), ', p-value: ', signif(p_value, 3)))

  }
  
  # Export results to a CSV file
  write.csv(results_df, file.path(plot_dir, "15_Simple_linear_regression_results.csv"), row.names = FALSE)
  
  
  # plot 
  combined_plot <- plot_grid(plotlist = plot_list, ncol = 5)
  combined_plot
  ggsave(file.path(plot_dir, "16_Simple_linear_regression.pdf"), width=17, height=6)
  

}




CalculateCorrelationMPXeniumSingleCell <- function(seurat_object_path, SampleName_vector, cellid_dir, colors_metaprograms_Xenium) {
  
  message(green('-------------------------------------------------------'))
  message(blue('[1/5] Loading reference single-cell object'))
  
  # load sc/snRNA-seq data
  data_SS2 <- qread(seurat_object_path)
  
  # calculate frequency
  SS2_metaprogram_freq <- data_SS2@meta.data %>% 
    dplyr::group_by(cell_type) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::mutate(freq =(n/sum(n)*100)) 
  SS2_metaprogram_freq
  
  message(green('-------------------------------------------------------'))
  message(blue('[2/5] Loading Xenium object frequency'))
  
  # load Xenium frequency
  
  metaprogram_frequency <- list()
  
  for (i in seq_along(SampleName_vector)) { 
    # read data
    cellID <- suppressWarnings(suppressMessages(read_csv(file.path(cellid_dir, paste0('cell_ID_', SampleName_vector[i], '.csv')))))
    print(paste0("Reading ", SampleName_vector[i]))
    
    # transform in data table with frequency
    metaprogram_frequency[[i]] <- cellID %>%
      group_by(group) %>%
      summarise(n = n()) %>%
      mutate(freq = (n/sum(n)*100))
  }
  
  names(metaprogram_frequency) <- SampleName_vector  
  
  # bind rows together and transform into df
  metaprogram_proportion <- bind_rows(metaprogram_frequency, .id = "Xenium_Region")
  metaprogram_proportion <- as.data.frame(metaprogram_proportion)
  
  metaprogram_proportion <- metaprogram_proportion %>%
    group_by(group) %>%
    summarise(freq = mean(freq))
  metaprogram_proportion
  
  colnames(metaprogram_proportion)[1] <- 'Metaprogram_Xenium'
  colnames(metaprogram_proportion)[2] <- 'freq_Xenium'
  
  SS2_metaprogram_freq <- SS2_metaprogram_freq[, c(1, 3)]
  colnames(SS2_metaprogram_freq)[1] <- 'Metaprogram_SS2'
  colnames(SS2_metaprogram_freq)[2] <- 'freq_SS2'
  
  message(green('-------------------------------------------------------'))
  message(blue('[3/5] Combine dataframes'))
  
  # Combine dataframes SS2 and Xenium
  
  matching_indices <- match(SS2_metaprogram_freq$Metaprogram_SS2, metaprogram_proportion$Metaprogram_Xenium)
  SS2_metaprogram_freq$Metaprogram_Xenium <- metaprogram_proportion$Metaprogram_Xenium[matching_indices]
  SS2_metaprogram_freq$freq_Xenium <- metaprogram_proportion$freq_Xenium[matching_indices]
  
  missing_metaprograms <- metaprogram_proportion[!metaprogram_proportion$Metaprogram_Xenium %in% SS2_metaprogram_freq$Metaprogram_SS2, ]
  
  # Create
  new_rows <- data.frame(
    Metaprogram_SS2 = missing_metaprograms$Metaprogram_Xenium, # Set Metaprogram_SS2 to the missing Xenium ones
    freq_SS2 = 0,                                             # Set freq_SS2 to 0
    Metaprogram_Xenium = missing_metaprograms$Metaprogram_Xenium, # Keep Xenium metaprograms
    freq_Xenium = missing_metaprograms$freq_Xenium            # Retain Xenium frequencies
  )
  
  SS2_metaprogram_freq <- rbind(SS2_metaprogram_freq, new_rows)
  SS2_metaprogram_freq <- SS2_metaprogram_freq %>%
    mutate(Metaprogram_Xenium = ifelse(is.na(Metaprogram_Xenium), Metaprogram_SS2, Metaprogram_Xenium),
           freq_Xenium = ifelse(is.na(freq_Xenium), 0, freq_Xenium))
  
  message(green('-------------------------------------------------------'))
  message(blue('[4/5] Calculate linear regression'))
  
  # calculate linear regression  -----------------------------------    
  model <- lm(freq_Xenium~freq_SS2, data=SS2_metaprogram_freq)
  summary_model <- summary(model)
  summary_model
  
  # Extract R-squared and p-value
  r_squared <- summary_model$r.squared
  p_value <- summary_model$coefficients[2, "Pr(>|t|)"]
  
  
  message(green('-------------------------------------------------------'))
  message(blue('[5/5] Plot linear regression'))
  
  ggplot(SS2_metaprogram_freq, aes(x = freq_SS2, y = freq_Xenium, color = Metaprogram_SS2)) +
    geom_smooth(method='lm', se=T, color='black') +
    geom_point(size = 4) +
    labs(subtitle = paste0('R² = ', round(r_squared, 3), ', p-value: ', signif(p_value, 3)))+ 
    labs(x = "Proportion Smart-seq2, %", y = "Proportion Xenium, %", color = 'Metaprogram') +
    scale_color_manual(values = colors_metaprograms_Xenium) +
    theme_minimal() +
    theme_minimal() + 
    # scale_y_continuous(limits = c(0, NA)) +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
          axis.text.y = element_text(size=12, colour="black"),
          axis.title=element_text(size=12),
          plot.title = element_text(size=10, face="bold"))  
  
}





PlotNicheResults <- function(data_dir, SampleName_vector, colors_metaprogram, colors_niches, niches_k, plot_spatial_maps = TRUE, mid_color = 15, max_col = 30) {
  
  for (i in seq_along(SampleName_vector)) { 
    
    print(paste0('Reading object ', SampleName_vector[i]))
    
    message(green('-------------------------------------------------------'))
    message(blue('[1/5] Loading Xenium data'))
    
    data <- qread(file.path(data_dir, paste0('seurat_obj_', SampleName_vector[i], '.qs')))
    
    if (plot_spatial_maps) {
      message(green('-------------------------------------------------------'))
      message(blue('[2/5] Plotting spatial maps'))
      
      celltype.plot <- ImageDimPlot(data, group.by = 'cell_type', fov = "fov",  cols = colors_metaprogram, border.size = NA, size = 0.5,
                                    dark.background = F) + ggtitle("Cell type")
      niche.plot <- ImageDimPlot(data, group.by = paste0("niche_", niches_k), fov = "fov",  cols = colors_niches, border.size = NA, size = 0.5,
                                 dark.background = F) + ggtitle("Niches")
      celltype.plot | niche.plot
      ggsave(file.path(plot_dir, paste0('1_ImageDimPlot_', SampleName_vector[i], "_", niches_k, '.pdf')), width = 16, height = 5)
    }
    
    message(green('-------------------------------------------------------'))
    message(blue('[3/5] Calculate frequency cell type per niche'))
    
    # Save data table with frequency
    table <- data@meta.data %>% group_by(cell_type, !!rlang::sym(paste0("niche_", niches_k))) %>%
      summarise(n = sum(n()))
    colnames(table)[2] <- 'niches'
    
    # Calculate total number of cells per niche
    table_sum <- table %>% group_by(niches) %>%
      summarise(n = sum(n))
    matching_indices <- match(table$niches, table_sum$niches)
    table$total_cells_per_niche <- table_sum$n[matching_indices]
    
    # Calculate frequency of cell type in niche
    table <- table %>% group_by(cell_type, niches) %>%
      summarise(freq = (n / total_cells_per_niche) * 100)
    
    # Transform into matrix
    cell_frequency <- table %>%
      pivot_wider(names_from = niches, values_from = freq, values_fill = 0)
    cell_frequency2 <- cell_frequency[, !names(cell_frequency) %in% "cell_type"]
    
    cell_frequency_df <- as.matrix(cell_frequency2)
    rownames(cell_frequency_df) <- cell_frequency$cell_type
    write.csv(as.data.frame(cell_frequency_df), row.names = TRUE, file.path(data_dir, paste0('2_niche_frequency_cell_type_', SampleName_vector[i], "_", niches_k, '.csv')))
    
    message(green('-------------------------------------------------------'))
    message(blue('[4/5] Plot heatmap cell frequency per niche'))
    
    # Heatmap annotations
    ann_col <- data.frame(rownames(t(cell_frequency_df)))
    colnames(ann_col) <- c('Niche')
    colours <- colors_niches[1:length(unique(ann_col$Niche))]
    names(colours) <- unique(ann_col$Niche)
    colours <- list('Niche' = colours)
    colAnn <- HeatmapAnnotation(df = ann_col, which = 'col', col = colours, show_annotation_name = FALSE, show_legend = FALSE)
    
    ann_row <- data.frame(colnames(t(cell_frequency_df)))
    colnames(ann_row) <- c('Metaprogram')
    colours <- list('Metaprogram' = colors_metaprogram)
    rowAnn <- HeatmapAnnotation(df = ann_row, which = 'row', col = colours, show_annotation_name = FALSE, show_legend = FALSE)
    
    col_fun = circlize::colorRamp2(c(0, mid_color, max_col), c("steelblue3", "white", "red3"))
    
    pdf(file.path(plot_dir, paste0('2_Heatmap_niche_', SampleName_vector[i], "_", niches_k, '.pdf')),
        width = 5, height = 6)
    draw(Heatmap(cell_frequency_df,
                 show_row_dend = TRUE,
                 show_column_dend = FALSE,
                 cluster_rows = TRUE,
                 cluster_columns = FALSE,
                 col = col_fun,
                 column_names_rot = 45,
                 show_column_names = TRUE,
                 top_annotation = colAnn,
                 left_annotation = rowAnn,
                 row_names_gp = gpar(fontsize = 10), row_names_side = "left", name = '-',
                 heatmap_legend_param = list(legend_direction = "vertical", legend_height = unit(5, "cm")),
                 width = ncol(t(t(cell_frequency_df))) * unit(10, "mm"),
                 height = nrow(t(t(cell_frequency_df))) * unit(10, "mm")))
    dev.off()
    
    message(green('-------------------------------------------------------'))
    message(blue('[5/5] Plot barplot frequency cell type per niche'))
   
    df <-  data@meta.data
    df$cell_type <- factor(df$cell_type, levels = order_metaprograms)
    # Plot niche frequency
    ggplot(df, aes(!!sym(paste0('niche_', niches_k)), fill = factor(cell_type, levels = order_metaprograms))) +
      scale_fill_manual(values = colors_metaprogram) + 
      geom_bar(position = "fill", color = "black") +
      guides(fill = guide_legend(title = "Cell type")) +
      labs(y = 'Proportion, %', x = '') + 
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1, colour = "black"),
            axis.text.y = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12)) 
    ggsave(file.path(plot_dir, paste0('3_BarPlot_niche_', SampleName_vector[i], "_", niches_k, '.pdf')), width = 5, height = 4)
    
    print(paste0('Finished object ', SampleName_vector[i]))
  }
}




HeatmapNicheCorrelation <- function(data_dir, SampleName_vector, colors_metaprogram, colors_niches, niches_k,
                                    nClust, outputname) {
  

  message(green('-------------------------------------------------------'))
  message(blue('[1/3] Read results with frequency of each cell type across each niche'))

  files <- list.files(data_dir, pattern = '2_niche_frequency_cell_type', full.names = T)
  files <- files[sapply(files, function(file) {
    any(sapply(SampleName_vector, function(name) grepl(name, file))) && grepl(paste0(niches_k, "\\.csv$"), file)
  })]
  
  samples <- gsub(paste0('2_niche_frequency_cell_type_|_', niches_k, '.csv'), '', basename(files))
  
  tab <- lapply(files, function(x) {
    read.csv(x, row.names = 1) %>%
      rename_all(~glue(paste0('Niche_{1:', niches_k ,'}'))) %>%
      reshape2::melt()
  })
  
  
  results <- list()
  for (i in seq_along(files)) {
    results[[i]] <- read.csv(files[i], row.names = 1) %>%
      rename_all(~glue(paste0('Niche_{1:', niches_k ,'}'))) %>%
      rownames_to_column('Program') %>%
      reshape2::melt() %>%
      mutate(variable = glue('{variable}-{samples[i]}'))
  }
  results <- do.call(rbind, results)
  
  
  res <- reshape2::acast(results, variable~Program, value.var = "value")
  res[is.na(res)] <- 0
  
  message(green('-------------------------------------------------------'))
  message(blue('[2/3] Plot heatmap'))
  
  
  anno <- data.frame(sample = gsub("\\_.*", "", gsub(paste0('Niche_[1-', niches_k ,']-'), '', rownames(res))),
                    # region = unlist(lapply(strsplit(rownames(res), '-'), `[[`, 2)),
                     row.names = rownames(res))
  
  # add column to annotation with source and location
  anno$sample <- gsub(paste0('Niche_[1-', niches_k ,']-'), '', rownames(anno))
  matching_indices <- match(anno$sample , metadata$SampleName)
  anno$SampleID <- metadata$SampleID[matching_indices]
  anno$Source <- metadata$Source[matching_indices]
  
  
  rowAnno <- HeatmapAnnotation(df = anno, which = 'row', 
                               col = list(sample = structure(clusterExperiment::bigPalette[1:length(unique(anno$sample))], names = unique(anno$sample)),
                                        #  region = structure(clusterExperiment::bigPalette[11:(10+length(unique(anno$region)))], names = unique(anno$region)),
                                          Source = col_sampling,
                                          SampleID = structure(col_patientID[1:length(unique(anno$SampleID))], names = unique(anno$SampleID))
                                          ))
  
  #col_fun = colorRamp2(c(0, 10, 20), c("grey90", "#CE78B3FF", "#573B88FF")) # this works well
  col_fun = colorRamp2(c(0, 7, 15), c("grey90", "#CE78B3FF", "#573B88FF")) # this works well
  
  # reorder columns
  res <- res[, order_metaprograms]
  
  library(dendextend)
  dend = as.dendrogram(hclust(dist(res)), type = "average") # can be used as cluster_rows and columns arg instead of T
  
  ht <- Heatmap(res, right_annotation = rowAnno, cluster_rows = T, 
                col = col_fun,
                column_names_rot = 45, 
                show_row_names = F,
                #row_split = 5, 
                # cluster_rows = T, 
                cluster_columns = F)
  pdf(file.path(plot_dir, paste0('4_Heatmap_frequency_cell_type_', niches_k, "_", outputname,'.pdf')),
      width = 7, height = 5)
  print(ht)
  dev.off()
  

  
  
  ## Hierarchical clustering of nmf factors with correlation among nmf scores
  nmf_factor_hc <- clusterNmfFactors(t(res))
  nmf_factor_cor <- nmf_factor_hc$cor_coef[nmf_factor_hc$hc_obj$order, nmf_factor_hc$hc_obj$order]
  
  ## Heatmap of correlations
  hm_colors <- rev((brewer.pal(n=9, name="RdBu")))
  hm_colors <- colorRampPalette(colors = hm_colors)
  
  Heatmap(nmf_factor_cor, col = hm_colors(100), cluster_rows = F, cluster_columns = F, show_column_dend = F,
          show_row_dend = F, show_row_names = T, show_column_names = F, name = 'cor', row_names_gp = gpar(fontsize = 8),
          layer_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", lwd = 0.5, fill = NA))})
  
  
  
  ## Aggregating programs
  nClust <- nClust
  tmp <- cutree(nmf_factor_hc$hc_obj, nClust)
  
  ## Ordering metaprograms order
  tmp <- tmp[rownames(nmf_factor_cor)]
  tmp <- plyr::mapvalues(tmp, unique(tmp), 1:nClust, warn_missing = FALSE)
  
  ## NMF metaprograms
  nmf_meta_programs <- lapply(sort(unique(tmp)), function(x) names(tmp)[tmp==x])
  names(nmf_meta_programs) <- paste("NMF", sort(unique(tmp)), sep = "_")
  
  anno <- NULL
  for (i in 1:length(nmf_meta_programs)) {
    anno <- rbind(anno, data.frame(nmf = names(nmf_meta_programs)[i], sample = nmf_meta_programs[[i]]))
  }
  rownames(anno) <- anno$sample
  anno$sample <- NULL
  anno$nmf <- factor(anno$nmf, levels = names(nmf_meta_programs))
  
  # add column to annotation with source and location
  anno$sample <- gsub(paste0('Niche_[1-', niches_k ,']-'), '', rownames(anno))
  matching_indices <- match(anno$sample , metadata$SampleName)
  anno$SampleID <- metadata$SampleID[matching_indices]
  anno$Source <- metadata$Source[matching_indices]

  
  colors_to_use <- colors_niches
  annoCol <- list(nmf = structure(colors_to_use[1:length(names(nmf_meta_programs))], names = names(nmf_meta_programs)))
  
  column_ha <- HeatmapAnnotation(df = anno, which = 'column', 
                                 col = list(nmf = annoCol$nmf,
                                            sample = structure(clusterExperiment::bigPalette[1:length(unique(anno$sample))], names = unique(anno$sample)),
                                            #  region = structure(clusterExperiment::bigPalette[11:(10+length(unique(anno$region)))], names = unique(anno$region)),
                                            Source = col_sampling,
                                            SampleID = structure(col_patientID[1:length(unique(anno$SampleID))], names = unique(anno$SampleID)),
                                 show_legend = F, show_annotation_name = F))
  
  row_ha <- HeatmapAnnotation(df = anno[, c('nmf')], which = 'row', 
                               col = list(df = annoCol$nmf),
                              show_legend = F, show_annotation_name = F)
  

  
  ht <- Heatmap(nmf_factor_cor, col = hm_colors(100), cluster_rows = F, cluster_columns = F, show_column_dend = F,
                show_row_dend = F, show_row_names = F, show_column_names = F, name = 'cor',
                top_annotation = column_ha, left_annotation = row_ha,
                row_split = anno$nmf, column_split = anno$nmf, column_title = NULL, row_title = NULL,
                width = ncol(nmf_factor_cor) * unit(1.4, "mm"),
                height = nrow(nmf_factor_cor) * unit(1.4, "mm"),
                layer_fun = function(j, i, x, y, width, height, fill) {
                  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", lwd = 0.5, fill = NA))})
  
  
  pdf(file.path(plot_dir, paste0('5_Heatmap_recurrent_niches_', niches_k, "_", outputname, '.pdf')),
      width = 10, height = 7)
  print(ht)
  dev.off()
  
  message(green('-------------------------------------------------------'))
  message(blue('[3/3] Plot pie chart with cell type distribution'))
  
  
  
  # calculate and plot average frequency of each cell type in each recurrent niche ---------------------------
  # calculate
  niches_programs <- lapply(split(anno, f = anno$nmf), rownames)
  mean_programs <- lapply(niches_programs, function(x) {
    results %>% filter(variable %in% x) %>%
      group_by(Program) %>%
      summarise(mean = mean(value))
  })

  df <- do.call(rbind, lapply(names(mean_programs), function(name) {
    data.frame(value = mean_programs[[name]], group = name, stringsAsFactors = FALSE)
  }))

  write.xlsx(df, file.path(data_dir, paste0('6_recurrent_niche_frequency_niche', niches_k, "_", outputname, '.xlsx')))
  
  
  # Plot 
  piechart_plot <- list()
  
  order <- order_metaprograms
  
  for (i in seq_along(mean_programs)){
    mean_programs[[i]] <- mean_programs[[i]][order(match(mean_programs[[i]]$Program, order)), ]
    
    piechart_plot[[i]] <- ggplot(mean_programs[[i]], aes(x = "", y = mean,  fill = factor(Program, levels = order_metaprograms))) +
      geom_col() + 
      scale_fill_manual(values=colors_metaprogram) + 
      coord_polar(theta = "y") + 
      labs(title = names(mean_programs)[i]) +
      theme_void() + NoLegend() + 
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  plot_grid(plotlist = piechart_plot, ncol = 1) 
  ggsave(file.path(plot_dir, paste0('7_PieChart_', niches_k, "_", outputname, '.pdf')), width=2.5, height=12)
  
  

}


