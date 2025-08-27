#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~ FUNCTIONS MOSTLY RELEVANT TO ASSIGN PROGRAMS STEP~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ReadMarkerGenes function takes in the location of the xenium_panel and reads the information that which gene is a marker for which celltype and 
# returns a list which maps the celltypes to a vector of the genes which are marker for it
# NOTE: Ensure that the annotation column is called Annotation
ReadMarkerGenes <- function(xenium_panel){
  # read the marker genes from the excel
  marker_genes <- read_xlsx(xenium_panel)
  
  # replace the spaces in the Annotation column with underscores
  marker_genes$Annotation<- gsub(" ", "_", marker_genes$Annotation)
  
  # extract the list for all markers as it will also be used
  marker_genes_list <- marker_genes %>%
    filter(!is.na(Annotation)) %>% # remove the NA rows (because these are those genes, which didn't have a label in excel)
    # do some renaming of Annotation values (this is quite general throughout xenium projects hence adding here, in helper functions)
    mutate(Annotation = str_replace(Annotation, 'Neurons.*', 'Neurons')) %>%
    mutate(Annotation = str_replace(Annotation, 'Microglia', 'Myeloid')) %>%
    mutate(Annotation = str_replace(Annotation, 'Endothelial.*', 'Endothelial')) %>%
    mutate(Annotation = str_replace(Annotation, 'Dendritic_cell', 'Myeloid')) %>%
    dplyr::group_by(Annotation) %>%
    dplyr::summarise(Gene_list = list(Gene)) %>%
    deframe()
  
  return (marker_genes_list)
}

# get the annotations of different celltypes using a list of celltype:marker_genes. Returned data object has a column with predictions 
# for each cell (panel_annotations_colname) and an individual column storing score for each celltype for each cell.
# mean_center = TRUE means that for the scores for any celltype, we will mean-center it, so that all the celltypes 
# are at a level playing field before finding which celltype has highest score for a given cell.
GetAnnotationsFromSpatialPanel <- function(data, marker_genes_list, mean_center = FALSE, panel_annotations_colname = 'predicted_label_UCell'){
  
  DefaultAssay(data) <- 'SCT'
  # add module scores with Ucell
  data <- AddModuleScore_UCell(data, features = marker_genes_list)
  U_Cell_signatures <- paste0(names(marker_genes_list), '_UCell')
  # retrieve the UCell scores and print the mean max UCell score for early identification of possible problems
  UCell_scores <- data@meta.data[, U_Cell_signatures]
  max_scores <- apply(UCell_scores, 1, max)
  mean_max_score <- mean(max_scores)
  print(paste0('Mean max UCell score is ', mean_max_score))
  
  if (mean_center == TRUE){
    # perform mean correction on the UCell_scores df, so that the different celltypes' scores are on a level playing field
    UCell_scores <- MeanCenterDF(UCell_scores, center_by = 'col')
    # modify the existing prediction scores to the mean corrected values for each celltype in data@meta.data
    data <- AddMetaData(object = data, metadata = UCell_scores, col.name = colnames(UCell_scores))
  }
  # get the predictions 
  UCell_scores$cell_type_UCell <- apply(UCell_scores, 1, function(row){
    names(UCell_scores)[which.max(row)]
  })
  
  # Remove "_UCell" from each element and add predicted labels to metadata
  UCell_scores$cell_type_UCell <- gsub("_UCell", "", UCell_scores$cell_type_UCell)
  data <- AddMetaData(object = data, metadata = UCell_scores$cell_type_UCell, col.name = panel_annotations_colname)
  return (data)
}

# get annotations from the panel, but just using the seurat's add module score func instead of the UCell variant
GetAnnotationsFromSpatialPanelSimpleAddModule <- function(data, marker_genes_list, mean_center = FALSE, panel_annotations_colname = 'predicted_label_seurat'){
  
  DefaultAssay(data) <- 'SCT'
  # add module scores with AddModuleScore
  data <- AddModuleScore(data, features = marker_genes_list, ctrl = 2) # keeping just 2 ctrl genes, because sometimes we have only very few marker genes for certain celltype and can't have less ctrl genes then
  # now the scores for celltype in index i in names(marker_genes_list) have been added in cols named Clusteri. Hence, loop through each celltype and change the name of the col storing its information
  for (i in seq_along(names(marker_genes_list))){
    # i = 1
    celltype <- names(marker_genes_list)[i]
    # store the scores for current celltype in a new col named celltype_seurat
    data@meta.data[, glue('{celltype}_seurat')] <- data@meta.data[[glue('Cluster{i}')]]
    # remove the col with name Clusteri
    data@meta.data[[glue('Cluster{i}')]] <- NULL
  }
  score_cols <- paste0(names(marker_genes_list), '_seurat')
  # retrieve the scores and print the mean max score for early identification of possible problems
  scores_df <- data@meta.data[, score_cols]
  max_scores <- apply(scores_df, 1, max)
  mean_max_score <- mean(max_scores)
  print(paste0('Mean max annotation score is ', mean_max_score))
  
  if (mean_center == TRUE){
    # perform mean correction on the scores df, so that the different celltypes' scores are on a level playing field
    scores_df <- MeanCenterDF(scores_df, center_by = 'col')
    # modify the existing prediction scores to the mean corrected values for each celltype in data@meta.data
    data <- AddMetaData(object = data, metadata = scores_df, col.name = colnames(scores_df))
  }
  # get the predictions 
  scores_df$cell_type_seurat <- apply(scores_df, 1, function(row){
    names(scores_df)[which.max(row)]
  })
  
  # Remove "_seurat" from each element and add predicted labels to metadata
  scores_df$cell_type_seurat <- gsub("_seurat", "", scores_df$cell_type_seurat)
  data <- AddMetaData(object = data, metadata = scores_df$cell_type_seurat, col.name = panel_annotations_colname)
  return (data)
}

# Updated version of GetAnnotationsFromScMapping, after a thorough investigation into single cell reference objects, 
# batch-correction, and corresponding label transfer. Changes include:
# i) no mean centering variation, as it is likely not the best way to correct for celltype proportion bias.
# ii) by default, recompute_residuals=FALSE and usage of reference.reduction='harmony' in FindTransferAnchors. This allows 
# not requiring a reference SCT model which was a challenge when the reference sc object came post combining from multiple samples 
# This function is only suitable for sc reference objects which had data from different samples and hence required batch correction.
# For the more straightforward sc reference objects which didn't need batch correction, I need to make another function.
GetAnnotationsFromScMappingV2 <- function(data, ref_projection, sc_annotations_colname = 'sc_annot'){
  
  # gracefully read ref_projection (depending on if it is the actual object or the path to it)
  seurat_object <- ReadRefProjection(ref_projection)
  
  # make sure the identity of the sc seurat object is 'cell_type'
  Idents(seurat_object) <- 'cell_type'
  
  # transfer anchors from single cell data
  anchors <- FindTransferAnchors(reference = seurat_object, query = data, normalization.method = "SCT", recompute.residuals = FALSE, reference.reduction = 'harmony')
  # transfer labels
  predictions <- TransferData(anchorset = anchors, refdata = seurat_object$cell_type, weight.reduction = data[["pca"]], dims = 1:20)
  
  # store in a copy, the neatened version (removed the score.max and predicted.id cols, and renamed colnames appropriately) of predictions df. 
  predictions_copy <- NeatenPredictionsDF(predictions)
  
  # get the predicted id according to which col has max value in each row
  predicted_id <- apply(predictions_copy, 1, function(row){
    colnames(predictions_copy)[which.max(row)]
  })
  
  # add the predicted labels to the data and information about each celltypes' prediction scores
  data <- AddMetaData(object = data, metadata = predicted_id, col.name = sc_annotations_colname)
  data <- AddMetaData(object = data, metadata = predictions_copy, col.name = paste0(colnames(predictions_copy), '_sc'))
  return (data)
}

# PlotSpatialMapsCelltypeSplit function takes in the data, colors object, outdir and samplename and makes the spatial maps 
# plot with the celltype split (so that we can see where the different celltypes are located).
PlotSpatialMapsCelltypeSplit <- function(data, colors, out_dir, sample_name, plot_name){
  # plot the spatial plots with celltypes split 
  ImageDimPlot(data, group.by = 'cell_type', border.size = NA,  dark.background = F, cols = colors, split.by = 'cell_type') 
  ggsave(glue('{out_dir}/analysis/2_program_annotation/plots/{plot_name}_{sample_name}.pdf'), width = 12, height = 8)
}

# ExportCellID function takes in the data and the out_dir and sample_name and saves a simple df as a csv file, which stores 
# the cell ids and their celltypes in cols named 'cell_id' and 'group'
ExportCellID <- function(data, out_dir, sample_name){
  # export cell ID
  cell_id <- data.frame(rownames(data@meta.data), data@meta.data$cell_type)
  
  # rename column names
  colnames(cell_id)[colnames(cell_id) == "rownames.data.meta.data."] <- "cell_id"
  colnames(cell_id)[colnames(cell_id) == "data.meta.data.cell_type"] <- "group"
  # save
  write_csv(cell_id, glue('{out_dir}/cell_ID_{sample_name}.csv'))
}

# combine the annotations from the malignant part of seurat object into the whole part (where the malignant cells 
# are annotated as 'Unknown' so far, in the cell_type column). annotations_col_in_data_mal usually takes values out of 
# c('predicted_label_snRNAseq', 'predicted_label_UCell')
# Assumptions: cell_type col stores annotations data. Cancer cells are 'Unknown' so far in data
CombineAnnotationsFromMalAndNormalAndPlot <- function(data, data_mal, annotations_col_in_data_mal, colors){
  # Combine annotations from the malignant part and normal part
  data@meta.data[data@meta.data$cell_type == 'Unknown', 'cell_type'] <- data_mal@meta.data[, annotations_col_in_data_mal] # now that we have the annotations for the cancer cells too, we just have to add the annotation information from data_mal to data
  # Plotting 
  p1 <- DimPlot(data, group.by = 'cell_type', cols = colors)
  p2 <- ImageDimPlot(data, group.by = 'cell_type', cols = colors, dark.background = FALSE)
  plot <- patchwork::wrap_plots(p1, p2, ncol = 2)
  return (plot)
}

# combine the annotations from the malignant part of seurat object into the whole part (where the malignant cells 
# are annotated as 'Unknown' so far, in the cell_type column). annotations_col_in_data_mal usually takes values out of 
# c('predicted_label_snRNAseq', 'predicted_label_UCell')
# Assumptions: cell_type col stores annotations data. Cancer cells are 'Unknown' so far in data
CombineAnnotationsFromMalAndNormal <- function(data, data_mal, annotations_col_in_data_mal, colors){
  warning('Consider using the better named function: CombineAnnotationsFromMalAndNormalAndPlot.')
  # Combine annotations from the malignant part and normal part
  data@meta.data[data@meta.data$cell_type == 'Unknown', 'cell_type'] <- data_mal@meta.data[, annotations_col_in_data_mal] # now that we have the annotations for the cancer cells too, we just have to add the annotation information from data_mal to data
  # Plotting 
  p1 <- DimPlot(data, group.by = 'cell_type', cols = colors)
  p2 <- ImageDimPlot(data, group.by = 'cell_type', cols = colors, dark.background = FALSE)
  plot <- patchwork::wrap_plots(p1, p2, ncol = 2)
  return (plot)
}

# filter the annotations of certain celltypes, which were called in the yaml to be belonging to certain clusters, but 
# likely were only a subset of the cluster. Hence, calling all cells in those clusters with cancer_marker score below a value, to 'Unknown'. 
# marker_genes_list should have the genes for cancer markers in entry named 'cancer_marker'
# since supplying ctrl=5 (control genes for each bin) in AddModuleScore, need atleast 5 genes ig in marker_genes_list['cancer_marker']
FilterAnnotationsOfConfusionCelltypes <- function(data, marker_genes_list, confusion_celltypes = c('Astrocyte'), confident_normal_celltypes = c('Myeloid', 'T-cells', 'Endothelial', 'Oligodendrocytes')){
  
  # first, score the cells for cancer_marker (we use a simple function, seurat's AddModuleScore)
  # data <- AddModuleScore(data, features = marker_genes_list['cancer_marker'], ctrl = 5, name = 'cancer_marker_module_score', nbin = 10)
  data <- AddModuleScore_UCell(data, marker_genes_list['cancer_marker']) # changed to using addmodulescore_UCell because of its ease of use and I think it works fine
  # FeaturePlot(data, features = 'cancer_marker_UCell')
  # now, get the median value of 'cancer_marker_UCell' for 'confident_normal_celltypes' labelled cells
  median_cancer_marker_score_normal <- median(data$cancer_marker_UCell[data$cell_type %in% confident_normal_celltypes])
  # now for the cells belonging to confusion_celltypes and having cancer_marker_UCell > median_cancer_marker_score_normal, we call them Unknown
  data@meta.data[data$cell_type %in% confusion_celltypes & data$cancer_marker_UCell > median_cancer_marker_score_normal, 'cell_type'] <- 'Unknown'
  return (data)
}

# add the normal annotations from the yaml file to data object, and then annotate the remaining cells using the specified 
# algorithm too, which involves  (parameters: 'method' and 'mean_center' determine how the remaining cells will be annotated). 
# method can be: 'sc_ref' or 'panel'. If it is 'sc_ref', value of marker_genes_list is immaterial. If it is 'panel', value 
# of ref_projection is immaterial. So in these cases can supply anything to these variables.
AddNormalAnnotationsAndAnnotateRemainingCells <- function(data, manual_annotation_yaml, malignant_celltypes, ref_projection, marker_genes_list, method, mean_center, recompute_residuals){
  # get the annotations from the yaml file
  data <- AddNormalCellManualAnnotations(data, manual_annotation_yaml, sample_name)
  # at this point, our data may possibly have annotations for celltypes which don't cluster separately too, like Astrocytes. 
  # We just want that only the cells with cancer_marker scores below median of confident normal celltypes are actually called these celltypes
  data <- FilterAnnotationsOfConfusionCelltypes(data, marker_genes_list, confusion_celltypes = c('Astrocyte'), confident_normal_celltypes = c('Myeloid', 'T-cells', 'Endothelial', 'Oligodendrocytes'))
  # subset data to just mal cells
  Idents(data) <- 'cell_type'
  data_mal <- subset(data, idents = 'Unknown')
  # now, according to which method is supplied, we obtain the annotations for the remaining cells
  if (method == 'sc_ref'){
    # read and subset reference to malignant cells
    sc_ref <- qread(ref_projection)
    Idents(sc_ref) <- 'cell_type'
    ref_projection_mal <- subset(sc_ref, idents = intersect(malignant_celltypes, sc_ref$cell_type))
    # get the annotations from the single cell mapping (with mean centering)
    data_mal <- GetAnnotationsFromScMapping(data_mal, ref_projection_mal, mean_center, recompute_residuals = recompute_residuals)
    Idents(data_mal) <- 'predicted_label_snRNAseq'
    # Combine annotations from the malignant part and normal part
    data@meta.data[data@meta.data$cell_type == 'Unknown', 'cell_type'] <- data_mal$predicted_label_snRNAseq # now that we have the annotations for the cancer cells too, we just have to add the annotation information from data_mal to data
  }else if (method == 'panel_ucell'){
    # subset the marker_genes_list to malignant cells
    marker_genes_list_mal <- marker_genes_list[names(marker_genes_list) %in% malignant_celltypes]
    data_mal <- GetAnnotationsFromSpatialPanel(data_mal, marker_genes_list_mal, mean_center)
    Idents(data_mal) <- 'predicted_label_UCell'
    # Combine annotations from the malignant part and normal part
    data@meta.data[data@meta.data$cell_type == 'Unknown', 'cell_type'] <- data_mal$predicted_label_UCell # now that we have the annotations for the cancer cells too, we just have to add the annotation information from data_mal to data
  }
  else if(method == 'panel_simple'){
    # subset the marker_genes_list to malignant cells
    marker_genes_list_mal <- marker_genes_list[names(marker_genes_list) %in% malignant_celltypes]
    # now get annotations in data_mal using seurat's simple AddModuleScore func
    data_mal <- GetAnnotationsFromSpatialPanelSimpleAddModule(data_mal, marker_genes_list_mal, mean_center)
    Idents(data_mal) <- 'predicted_label_seurat' # just calling a different name of the label column and not UCell
    # Combine annotations from the malignant part and normal part
    data@meta.data[data@meta.data$cell_type == 'Unknown', 'cell_type'] <- data_mal$predicted_label_seurat # now that we have the annotations for the cancer cells too, we just have to add the annotation information from data_mal to data
  }else{
    stop("Invalid method. Should be either of 'sc_ref' or 'panel_ucell' or 'panel_simple'")
  }
  Idents(data) <- 'cell_type'
  return (list(data, data_mal))
}

# assign to all cells belonging to same cluster, the celltype which was majority in that cluster till now. 
AssignMajorityCelltypeToClusters <- function(data, celltype_col, cluster_col = 'seurat_clusters'){
  # create two new cols in data, which have the names cluster_col and celltype_col and store the corresponding info. This makes doing dplyr operations easy
  data$celltype_col = data@meta.data[[celltype_col]]
  data$cluster_col = data@meta.data[[cluster_col]]
  
  # get a df which has information about the counts of each celltype for each cluster (in a wide form).
  cluster_celltype_counts_wide = data@meta.data %>% select(cluster_col, celltype_col) %>% 
    group_by(cluster_col, celltype_col) %>% summarise(counts = n(), .groups = 'drop') %>% 
    pivot_wider(id_cols = cluster_col, names_from = celltype_col, values_from = counts, values_fill = 0)
  
  # drop the first col, and keep in a copy so that below apply operations are easier, more readable
  cluster_celltype_counts_wide_without_first_col = cluster_celltype_counts_wide[2:ncol(cluster_celltype_counts_wide)]
  
  # now, use apply func to get the majority celltype for each cluster
  cluster_celltype_counts_wide$majority_celltype <- apply(cluster_celltype_counts_wide_without_first_col, 1, function(row){
    colnames(cluster_celltype_counts_wide_without_first_col)[which.max(row)]
  })
  
  # finally, we get a mapping between cluster and its majority celltype: cluster_celltype_counts_wide %>% select(cluster_col, majority_celltype)
  # apply the mapping to replace the cell_type col in data with updated one, with majority celltypes of each cluster
  mapping_vec <- cluster_celltype_counts_wide %>% select(cluster_col, majority_celltype) %>% deframe()
  data@meta.data[, celltype_col] <- mapping_vec[data$cluster_col]
  
  # remove from data the created cols in this func
  data$celltype_col <- NULL
  data$cluster_col <- NULL
  
  return (data)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~FUNCTIONS MOSTLY RELEVANT TO MANUAL ANNOTATION OF NORMAL CELLTYPES STEP~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Add the annotations for the normal cells, which were manually marked down in a yaml file, commonly named: normal_celltypes_manual_annotation.yaml. 
# The annotations are added in a col named 'cell_type', and the clusters are assumed to be stored in col named 'seurat_clusters'. 
# This function was originally made while working on the Ependymoma project.
AddNormalCellManualAnnotations <- function(data, manual_annotation_yaml, sample_name){
  
  # read the normal_celltypes_manual_annotation.yaml file to and extract the information about which 
  # clusters are normal and which are malignant from that.
  manual_annotation <- yaml.load_file(manual_annotation_yaml)
  data$cell_type <- 'Unknown'
  for (celltype in names(manual_annotation[[sample_name]])){
    clusters_for_current_celltype <- manual_annotation[[sample_name]][[celltype]]
    data$cell_type[data$seurat_clusters %in% clusters_for_current_celltype] <- celltype
  }
  return (data)
}

# version 2 of the original function, in which I make an extra feature plot of just the canonical astrocyte genes 
# (AQP4, GJA1, GFAP) to better assist in finding astrocyte rich clusters. Also, I add a featureplot showing ac-like 
# scores from the single cell data. 
# the specified featureplot_cols or histogram_cols can just be a superset of the cols present in data. It is not 
# necessary for all of them to be actually there in the data
MakePlotsHelpingInManualAnnotationV2 <- function(data, histogram_cols = c('T-cells_UCell', 'Myeloid_UCell', 'Endothelial_UCell', 
        'Neurons_UCell', 'Oligodendrocytes_UCell', 'Astrocyte_UCell', 'canonical_astrocyte_UCell', 'cancer_marker_UCell', 
        'T-cells_sc', 'Myeloid_sc', 'Endothelial_sc', 'Neurons_sc', 'Oligodendrocytes_sc', 'Astrocyte_sc', 'AC-like_sc'), 
      featureplot_cols = NULL, ncol_for_histogram_and_featureplot = 5){ # ncol_for_histogram_and_featureplot indicates ncol when we wrap_plots (histogram and featureplots)
  
  # if no featureplot_cols supplied, then simply use the histogram_cols as featureplot_cols
  if (is.null(featureplot_cols)){
    featureplot_cols = histogram_cols
  }
  
  # first, subset the three types of cols to ones which are available in data
  cols_in_data <- colnames(data@meta.data)
  featureplot_cols <- intersect(featureplot_cols, cols_in_data)
  # top_featureplot_cols <- intersect(top_featureplot_cols, cols_in_data)
  histogram_cols <- intersect(histogram_cols, cols_in_data)
  # right_featureplot_cols <- intersect(right_featureplot_cols, cols_in_data)
  
  # get the df storing the mean value for each histrogram_col for each cluster
  clusterwise_mean_score_hist_cols <- data@meta.data %>% select(all_of(c('seurat_clusters', histogram_cols))) %>% 
    pivot_longer(!seurat_clusters, names_to = 'colname', values_to = 'score') %>% 
    group_by(seurat_clusters, colname) %>% summarise(mean_score = mean(score), .groups = 'drop') %>% 
    pivot_wider(names_from = colname, values_from = mean_score, id_cols = seurat_clusters)
  
  # now make the histogram for each col in histogram_cols
  all_hists_list <- list()
  for (i in seq_along(histogram_cols)){
    col <- histogram_cols[i]
    plot <- CelltypeScoreHistogram(df = clusterwise_mean_score_hist_cols, col)
    all_hists_list[[i]] <- plot
  }
  p1 <- patchwork::wrap_plots(all_hists_list, ncol = ncol_for_histogram_and_featureplot)
  
  # now make another plot by appending, which stores the featureplots we want to make
  all_featureplots_list <- list()
  for (i in seq_along(featureplot_cols)){
    col <- featureplot_cols[i]
    plot <- FeaturePlot(data, features = col) + theme_void()
    all_featureplots_list[[i]] <- plot
  }
  p2 <- patchwork::wrap_plots(all_featureplots_list, ncol = ncol_for_histogram_and_featureplot)
  
  p3 <- DimPlot(data, group.by = 'seurat_clusters', label = TRUE, label.size = 7, cols = 'polychrome') + 
    theme_void() + NoLegend() + ggtitle(NULL) + coord_fixed()
  
  bot_plot <- patchwork::wrap_plots(p2, p3, ncol = 2, widths = c(1.5,1))
  plot <- patchwork::wrap_plots(p1, bot_plot, ncol = 1, heights = c(1.5,1))
  # ggsave('~/plot.pdf', width = 20, height = 15)
  
  return (plot)
}

# function to take in the data object, and other parameters, to perform single cell and spatial annotations on it, and 
# make the dimplot showing seurat clusters, featureplots showing the scores of the annotations from the two methods for 
# the normal celltypes, and histograms showing the frequency of the various mean celltype score for the seurat clusters.
# Added cancer_marker also into the default normal_celltypes vector, because we want the histogram corresponding to it to
# be also made. This is probably gonna be a common theme in most projects, so making an exception for it. It will help a 
# bit in being more confident with the annotations for the normal celltypes
MakePlotsHelpingInManualAnnotation <- function(data, normal_celltypes = c('T-cells', 'Myeloid', 'Endothelial', 'Neurons', 'Oligodendrocytes', 'Astrocyte', 'cancer_marker')){
  
  # print the normal celltypes vector to make the user ensure that they are the intended ones
  message(blue('The following normal celltypes are being used. Ensure that it is a superset of the intended celltypes.'))
  message(blue(paste(normal_celltypes, collapse = ', ')))
  
  # get the column names in data which have annotation scores of the normal celltypes
  panel_scores_colnames <- c(intersect(paste0(normal_celltypes, '_UCell'), colnames(data@meta.data)), 'seurat_clusters')
  sc_scores_colnames <- c(intersect(paste0(normal_celltypes, '_sc'), colnames(data@meta.data)), 'seurat_clusters')
  
  # get the df storing the mean score for each celltype for each cluster
  clusterwise_mean_score_normal_celltypes <- data@meta.data %>% select(all_of(panel_scores_colnames)) %>%
    tidyr::gather(key = 'cell_type', value = 'score', -seurat_clusters) %>% 
    group_by(seurat_clusters, cell_type) %>% summarise(mean = mean(score)) %>% 
    tidyr::spread(key = 'cell_type', value = 'mean')
  
  # now make the histogram for each celltype (for panel scorings)
  cols_for_histogram <- intersect(paste0(normal_celltypes, '_UCell'), colnames(data@meta.data))
  for (i in seq_along(cols_for_histogram)){
    col <- cols_for_histogram[i]
    if (i == 1){
      all_hists <- CelltypeScoreHistogram(df = clusterwise_mean_score_normal_celltypes, col)
    }else{
      plot <- CelltypeScoreHistogram(df = clusterwise_mean_score_normal_celltypes, col)
      all_hists <- all_hists + plot
    }
  }
  # depending on the length of all_hists, we will add more featureplots. Things to possibly include: 
  # cluster_homogeneity scores if available, featureplots of Neurons' scores, Oligos scores, Astrocyte scores.
  # Target shape of the first section of plot (which we are calling all_hists) is 3x3
  
  # add plots according to preference order
  if ('channel_1_mean' %in% colnames(data@meta.data)){
    all_hists <- all_hists + FeaturePlot(data, features = 'channel_1_mean') # IF data telling which cells have h3k27m (driver mutation for dmg)
  }
  if ('cancer_marker_UCell' %in% colnames(data@meta.data)){ # important to tell normal celltypes from cancerous
    all_hists <- all_hists + FeaturePlot(data, features = 'cancer_marker_UCell')
  }
  if ('Cycling_UCell' %in% colnames(data@meta.data)){ # important to tell normal celltypes from cancerous
    all_hists <- all_hists + FeaturePlot(data, features = 'Cycling_UCell')
  }
  if ('cluster_homogeneity' %in% colnames(data@meta.data)){ # important because in general helpful in telling normal cells from cancer
    all_hists <- all_hists + FeaturePlot(data, features = 'cluster_homogeneity')
  }
  if ('Astrocyte_UCell' %in% colnames(data@meta.data)){ # important because generally, there is confusion between Astrocyte and AC-like
    all_hists <- all_hists + FeaturePlot(data, features = 'Astrocyte_UCell')
  }
  if ('AC-like_UCell' %in% colnames(data@meta.data)){ # important because sometimes just looking at Astrocytes distribution may not give an idea of which clusters have Astrocytes/ac-like (contributed to by the fact that the two important genes for them, aqp4 and gja1 are ac-like markers and not astrocyte markers in many panels)
    all_hists <- all_hists + FeaturePlot(data, features = 'AC-like_UCell')
  }
  if ('Neurons_UCell' %in% colnames(data@meta.data)){ # important because generally, we don't have informatoin about Neuronal annotations in single cell data
    all_hists <- all_hists + FeaturePlot(data, features = 'Neurons_UCell')
  }
  if ('Oligodendrocytes_UCell' %in% colnames(data@meta.data)){ # important because sometimes sc data isn't able to tell oligos, and panel information can be helpful in these cases
    all_hists <- all_hists + FeaturePlot(data, features = 'Oligodendrocytes_UCell')
  }
  
  # make the relevant plots which will help in doing the manual annotation
  p1 <- DimPlot(data, group.by = 'seurat_clusters', label = TRUE, label.size = 7, cols = 'polychrome') + NoLegend()
  p2 <- FeaturePlot(data, features = intersect(paste0(normal_celltypes, '_sc'), colnames(data@meta.data)))
  top <- all_hists
  bottom <- p1 + p2
  plot <- top/bottom + plot_layout(heights = c(3,2))
  return (plot)
}

# function to make the histogram for a particular column in clusterwise_mean_score_normal_celltypes df
CelltypeScoreHistogram <- function(df, col, clustering_col = 'seurat_clusters'){
  
  # just naming the col of interest as temp so that below code is more readable
  df[,'temp'] <- df[,col]
  
  # just subset to relevant cols for cleanliness
  df <- df %>%
    select(c(all_of(clustering_col), temp))
  
  # Calculating the Sturges bins
  breaks <- pretty(range(df$temp), 
                   n = nclass.Sturges(df$temp),
                   min.n = 1)
  
  # now need to add information about x and y coordinates of the text (seurat_clusters names) in the df
  # get midpoints of breaks
  midpoints <- (breaks[-length(breaks)] + breaks[-1])/2
  
  # for each entry in df, map the the score to the correct bin midpoint (it will just be the nearest point)
  df$xmid <- cut(df$temp, breaks = breaks, labels = midpoints)
  
  # then create col storing information that which entry that is, for its midpoint
  df$text_y <- ave(seq_along(df$xmid), df$xmid, FUN = seq_along)
  
  # then fill y values using that information
  df$text_y <- df$text_y + 0.5
  
  # also make the xmid col numeric because otherwise plot gives error
  df$xmid <- as.numeric(as.character(df$xmid))
  
  # now make the plot with the appropriate bin breaks
  plot <- ggplot(df) + 
    geom_histogram(color = 1, aes(x = temp), breaks = breaks, fill = 'grey') + 
    xlab(col) + 
    geom_text(aes(x = xmid, y = text_y, label = !!sym(clustering_col)), size = 5) + 
    theme_minimal() + 
    ylab(NULL)
  
  return (plot)
}

# function which will take a seurat object, and two vectors, which hold values from the annotations col 
# of the object, and returns the middle point of the cluster_homogeneity means of the two groups
FindBestClusterHomogeneitySplit <- function(data, group1, group2, annotations_col = 'cell_type'){
  
  # get the first group's cluster homogeneity mean
  group1_df <- data@meta.data %>%
    filter(!!sym(annotations_col) %in% group1)
  mean1 <- mean(group1_df[,'cluster_homogeneity']) # assuming that the col storing cluster homogeneity values will be called cluster homogeneity only
  
  # get the second group's cluster homogeneity mean  
  group2_df <- data@meta.data %>%
    filter(!!sym(annotations_col) %in% group2)
  mean2 <- mean(group2_df[,'cluster_homogeneity'])
  
  # return the middle point of the two means
  return ((mean1 + mean2)/2)
}



