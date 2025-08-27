# NicheAnalysis function takes in the location of the annotated seurat object, the minimum and maximum number of niches 
# and builds niche assay for all the values of niche numbers between min and max, and the given number of neighbors for 
# the same, and saves the seurat object appropriately.
NicheAnalysis <- function(data, niches_min, niches_max, n_neighbors) {
  
  # performing niche analysis
  for (i in seq(niches_min, niches_max, 1)){
    print(glue('Calculating niche {i}')) 
    data <- BuildNicheAssay(object = data, fov = "fov", group.by = "cell_type", niches.k = i, 
                            neighbors.k = n_neighbors, cluster.name = glue('niche_{i}'))
    print(glue('Finished calculating niche {i}')) 
  }

  return (data)
}

# make the celltype and niche grouped imagedimplot of a sample
PlotNicheSpatialMaps <- function(data, niche_resolution, colors_metaprogram, colors_niches){
  
  # make the imagedimplot with celltypes as grouping
  celltype.plot <- ImageDimPlot(data, group.by = 'cell_type', fov = "fov",  cols = colors_metaprogram, border.size = NA, size = 0.5,
                                dark.background = F) + ggtitle("Cell type")
  # make the imagedimplot with niches as grouping
  niche.plot <- ImageDimPlot(data, group.by = paste0("niche_", niche_resolution), fov = "fov",  cols = colors_niches, border.size = NA, size = 0.5,
                             dark.background = F) + ggtitle("Niches")
  # return the plot
  return (celltype.plot + niche.plot)
}

# niche_resolution is the niche resolution which we are interested in (generally one of 4,5,6)
# returns a df with cols: cell_id, group, niche
GetNicheInformation <- function(sample_name, niche_dir, niche_resolution){
  
  # read the csv (cols: cell_id, group, niche_4, niche_5, niche_6, etc...)
  niche_info = read.csv(glue('{niche_dir}/{sample_name}.csv'))
  
  # limit to a particular niche resolution
  niche_info <- niche_info %>% select(c(cell_id, group, glue('niche_{niche_resolution}')))
  
  # rename the column with niche information by 'niche'
  niche_info[, 'niche'] <- niche_info[, glue('niche_{niche_resolution}')]
  niche_info[, glue('niche_{niche_resolution}')] <- NULL
  
  # return df in its current, long form so that user of this function has more flexibility in modifying it in whichever way they like
  return (niche_info)
}

# get the normalized (proportion of each celltype in each niche) niche information in wide form (colnames: group,1,2,3...)
# form is either wide or long
GetNicheInformationNormalized <- function(sample_name, niche_dir, niche_resolution, form = 'wide'){
  
  #### Obtaining Niche information for sample 
  # get the niche information (cols: cell_id, group, niche)
  niche_info <- GetNicheInformation(sample_name, niche_dir, niche_resolution)
  
  # get count of each celltype in each niche
  niche_counts <- niche_info %>% group_by(group, niche) %>% summarise(counts = n(), .groups = 'drop')
  # get total number of cells in each niche (cols: niche, counts)
  niche_total_counts <- niche_info %>% group_by(niche) %>% summarise(niche_counts = n())
  
  # left join niche_counts with niche_total_counts
  niche_counts <- niche_counts %>% left_join(niche_total_counts, by = 'niche')
  
  # Calculate frequency (proportion) of cell type in niche. After this, table has cols: cell_type, niches, freq
  niche_counts <- niche_counts %>% mutate(proportion = (counts/niche_counts) * 100)
  
  # return the dataframe according to supplied argument (wide form or long form)
  if (form == 'wide'){
    # transform into wide form, such that the rownames are celltypes and colnames are the niche numbers. Each 
    # entry tells the proportion of a particular celltype in a particular niche.
    cell_frequency <- niche_counts %>% pivot_wider(names_from = niche, values_from = proportion, values_fill = 0, id_cols = group)
  }else if (form == 'long') {
    # if long form, then just select the appropriate columns and return
    cell_frequency <- niche_counts %>% select(group, niche, counts, proportion) # counts is useful for correct averaging when we do niche clustering
  }else{
    stop ("Incorrect value of form. It is one of 'wide'/'long'")
  }
  
  return (cell_frequency)
}

# making the plots related to niche- heatmap for each sample, showing proportion of celltypes in each niche
# we save the plot from within the function because then specifying the dimensions of the plot is easier
PlotNicheResults <- function(sample_name, niche_dir, niche_resolution, colors_metaprogram, colors_niches, out_dir) {
  
  # get the normalized (proportion of each celltype in each niche) niche information in wide form (colnames: group,1,2,3...)
  cell_frequency <- GetNicheInformationNormalized(sample_name, niche_dir, niche_resolution)
  
  # obtain the matrix of same information as cell_frequency df
  mat <- as.matrix(cell_frequency[, 2:ncol(cell_frequency)])
  rownames(mat) <- cell_frequency$group

  #### Plot heatmap cell frequency per niche
  
  # initialize the colors for the heatmap
  heatmap_color = circlize::colorRamp2(c(0, max(mat)), c("white", "gold3"))
  annotations_color = list(niche = colors_niches, celltype = colors_metaprogram)
  
  # make the HeatmapAnnotation objects
  colAnn <- HeatmapAnnotation(niche = colnames(mat), which = 'col', col = annotations_color, show_legend = TRUE, show_annotation_name = FALSE)
  rowAnn <- HeatmapAnnotation(celltype = rownames(mat), which = 'row', col = annotations_color, show_legend = FALSE, show_annotation_name = FALSE)
  
  # make the heatmap
  # construct the heatmap-class object which stores the information about the plot
  plot <- Heatmap(mat,
                  cluster_rows = TRUE,
                  show_row_dend = TRUE,
                  cluster_columns = FALSE,
                  show_column_dend = FALSE,
                  col = heatmap_color,
                  column_names_rot = 45,
                  top_annotation = colAnn,
                  left_annotation = rowAnn,
                  row_names_gp = gpar(fontsize = 10),
                  heatmap_legend_param = list(legend_direction = "vertical", legend_height = unit(5, "cm"), title = 'Frequency (%)'),
                  width = ncol(mat) * unit(10, "mm"),
                  height = nrow(mat) * unit(10, "mm")
                  )
  SaveHeatmap(plot, out_dir, width = ncol(mat)*unit(2, "mm"), height = nrow(mat)*unit(1, "mm"))
  # draw(plot)
  return (plot) # also returning just in case we need to check the plot from outside the function
}

# make the big heatmap showing the correlations between niches from the different samples, and the pi charts showing 
# average proportions of the different celltypes in each niche.
# return a list of heatmap object and pi chart object
# plot_suffix stores a string to append at the end in the name of plots (was used when we made plots for just Primary and just Reccurence samples)
HeatmapNicheCorrelation <- function(metadata, niche_dir, niche_resolution, n_clusters, plot_dir, plot_suffix = ''){
  library(dendextend)
  library(progress)
  library(scales)
  
  # get the sample_names from metadata
  sample_names <- metadata$SampleName
    
  # get the niche info from each sample (has cols: group, niche, proportion, sample_name (which we add in below loop))
  niche_info <- data.frame()
  pb <- progress_bar$new(format = 'Loading niche results [:bar] :percent', total = length(sample_names), clear = FALSE)
  for (sample_name in sample_names){
    cell_frequency <- GetNicheInformationNormalized(sample_name, niche_dir, niche_resolution, form = 'long')
    cell_frequency$sample_name <- sample_name
    niche_info <- rbind(niche_info, cell_frequency)
    pb$tick()
  }
  
  # create a new column in niche_info, which stores the combination of sample_name and niche. We will call this as the 'Identifier' column
  niche_info <- niche_info %>% mutate(Identifier = glue('Niche_{niche}-{sample_name}'))
  
  # convert the niche_info df into wide form so that we have proportion for all celltypes for all identifiers (niche-sample pair)
  niche_info_wider <- niche_info %>% pivot_wider(values_from = proportion, values_fill = 0, names_from = group, id_cols = Identifier)
  
  # now, get a matrix from niche_info_wider, storing just the numerical information from it, and rename the rownames and colnames accordingly
  niche_info_wider_mat <- as.matrix(niche_info_wider[2:ncol(niche_info_wider)])
  rownames(niche_info_wider_mat) <- niche_info_wider$Identifier
  
  # now perform hierarchical clustering of nmf factors with correlation among nmf scores
  nmf_factor_hc <- clusterNmfFactors(t(niche_info_wider_mat))
  # reorder the columns and rows for visualization
  nmf_factor_cor <- nmf_factor_hc$cor_coef[nmf_factor_hc$hc_obj$order, nmf_factor_hc$hc_obj$order]
  
  # split the nmf factors into appropriate groupings
  identifier_cluster_vec <- cutree(nmf_factor_hc$hc_obj, n_clusters)
  # save metaniche assignment csv
  # write_csv(data.frame(identifier = names(identifier_cluster_vec), metaniche = identifier_cluster_vec), '~/temp.csv')
  
  # at this point, identifier_cluster_vec has names of identifiers in the order of niche_info_wider and not in order of 
  # row/colnames in nmf_factor_cor, which is what we want. Hence, reorder the entries in identifier_cluster_vec accordingly
  identifier_cluster_vec_reordered <- identifier_cluster_vec[rownames(nmf_factor_cor)]
  # we also want that class 1 appears first time before 2, 2 appears first time before 3 and so on, so that in heatmap splitting, the sequence of rows and cols doesn't change when we call split (split seemingly orders the rows and cols according to their sorted values, idk why)
  identifier_cluster_vec_reordered <- plyr::mapvalues(identifier_cluster_vec_reordered, from = unique(identifier_cluster_vec_reordered), to = sort(unique(identifier_cluster_vec_reordered)))
  
  # get a vector of sample_names in appropriate order as the heatmap's rows, so that we can get metadata's entries accordingly
  annotations_df <- data.frame(Identifier = names(identifier_cluster_vec_reordered))
  annotations_df <- annotations_df %>% mutate(SampleName = gsub('Niche_.?-', '', Identifier))
  annotations_df <- annotations_df %>% left_join(metadata, by = 'SampleName') %>% select(Identifier, SampleName, SampleID, Source) # left join from metadata to get the relevant information about each sample. Selecting few columns for ease of handling.
  # also get information about total size of each niche (identified by identifier)
  niche_info_sizes <- niche_info %>% group_by(Identifier) %>% summarize(total_counts = sum(counts)) # df with cols: Identifier, total_counts (total cells belonging to that identifier (particular niche in a sample))
  annotations_df <- annotations_df %>% left_join(niche_info_sizes, by = 'Identifier')
  
  # set up colors for the heatmap
  # hm_colors <- circlize::colorRamp2(c(min(nmf_factor_cor), 0, max(nmf_factor_cor)), c("steelblue3", "white", "red3"))
  hm_colors <- rev((brewer.pal(n=9, name="RdBu")))
  hm_colors <- colorRampPalette(colors = hm_colors)
  cluster_colors <- hue_pal()(n_clusters) # making a named vector storing colors we will allocate to clusters
  names(cluster_colors) <- seq(from = 1, to = n_clusters)
  col_annot_colors <- list(Source = col_sampling, cluster = cluster_colors, sample = col_sampleid_all)
  
  # Annotations
  # get a vector storing the size of 
  row_annot <- HeatmapAnnotation(cluster = as.factor(identifier_cluster_vec_reordered), which = 'row', show_legend = F, show_annotation_name = F, col = col_annot_colors)
  col_annot <- HeatmapAnnotation(
                                 # niche_size = annotations_df$total_counts,
                                 cluster = as.factor(identifier_cluster_vec_reordered),
                                 # Source = annotations_df$Source, 
                                 sample = annotations_df$SampleID,
                                 which = 'col', show_legend = TRUE, show_annotation_name = FALSE, col = col_annot_colors)
  
  # finally, make the heatmap
  ht <- Heatmap(
    nmf_factor_cor,
    col = hm_colors(100),
    top_annotation = col_annot,
    left_annotation = row_annot,
    # bottom_annotation = col_annot,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    name = 'correlation', # name of the heatmap which is used as title of heatmap legend
    row_split = identifier_cluster_vec_reordered,
    column_split = identifier_cluster_vec_reordered, 
    row_title = NULL,
    column_title = NULL,
    # width = ncol(nmf_factor_cor) * unit(0.8, "mm"), 
    # height = nrow(nmf_factor_cor) * unit(0.8, "mm"), # height and width we will supply during saving the plot
    # layer_fun = function(j, i, x, y, width, height, fill) {
    #   grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = 'grey', lwd = 0.5, fill = NA))
    # }
  )
  # draw(ht)
  
  message(blue('Making pie chart with cell type distribution'))
  
  # need to obtain the average frequency of each celltype in each niche cluster. For this, will add all information to annotations_df
  identifier_cluster_vec_reordered_df <- data.frame(Identifier = names(identifier_cluster_vec_reordered), cluster = identifier_cluster_vec_reordered)
  identifier_cluster_vec_reordered_df <- niche_info %>% left_join(identifier_cluster_vec_reordered_df, by = 'Identifier') # has each row storing info about a particular celltype in a particular niche in a particular sample
  pie_chart_df <- identifier_cluster_vec_reordered_df %>% group_by(cluster, group) %>% summarise(celltype_counts_in_cluster = sum(counts), .groups = 'drop') # has cols: cluster, group (celltype), celltype_counts_in_cluster. Here, cluster refers to niche cluster
  pie_chart_df <- GroupAndAddProportionsCol(df = pie_chart_df, grouping_col = 'cluster', class_col = 'group', counts_col = 'celltype_counts_in_cluster') # group by cluster, and then compute proportions of each celltype and accordingly add new col named 'proportions'
  
  # now we can simply make the pie_charts
  pie_plot <- ggplot(pie_chart_df, aes(x = "", y = proportions, fill = group)) + 
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = colors_metaprograms_Xenium) + 
    coord_polar(theta = 'y') + 
    facet_wrap(~cluster) + 
    theme_void()
  
  # save the plots
  ggsave(glue('{plot_dir}/Proportion_pie_charts_res_{niche_resolution}_nclust_{n_clusters}{plot_suffix}.pdf'), width = 7, height = 5)
  # SaveHeatmap(ht, glue('~/plot.pdf'), width = ncol(nmf_factor_cor) * unit(0.05, "mm"), height = nrow(nmf_factor_cor) * unit(0.04, "mm"))
  SaveHeatmap(ht, glue('{plot_dir}/Correlation_heatmap_res_{niche_resolution}_nclust_{n_clusters}{plot_suffix}.pdf'), 
              width = ncol(nmf_factor_cor) * unit(0.1, "mm"), height = nrow(nmf_factor_cor) * unit(0.08, "mm"))
  
  # return (pie_chart_df) # this is used for downstream analysis, to find out which niches have interesting proportions of celltypes
  # return (identifier_cluster_vec_reordered_df) # this is used to relabel the niches in all the samples according to the metaniches
}


