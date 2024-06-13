library(infercnv)

addNormalControl <- function(ext_ctrl, inferCNV_analysis_folder, cm_raw, orig_samples, type){
  ## Read ctrl cm
  message("Loading normal control count matrix...")
  if (type == 'frozen'){
    suffix <- "_nuc_premrna_counts.rds"
  } else if (type == 'fresh') {
    suffix <- "_fresh_counts.rds"
  }

  mg_cm <- readRDS(file.path(ext_ctrl, paste0("mg", suffix)))
  od_cm <- readRDS(file.path(ext_ctrl, paste0("od", suffix)))
  ext_ctrl_cm <- cbind(mg_cm, od_cm)

  ## Generate control sample list
  message("Loading normal control annotations...")
  normal_samples <- structure(c(rep("mg", ncol(mg_cm)), rep("od", ncol(od_cm))), names = c(colnames(mg_cm), colnames(od_cm)))

  ## Reorder the ext_ctrl cm and subset genes
  message("Reordering cell names for normal control cm and annotations")
  ext_ctrl_cm <- ext_ctrl_cm[,names(normal_samples)]
  ext_ctrl_cm <- ext_ctrl_cm[rownames(cm_raw),]

  ## Concatenate samples and normal controls
  message("Concatenate cm and annotations of normal controls to patient data...")
  cm <- cbind(cm_raw, ext_ctrl_cm)
  samples <- c(orig_samples, normal_samples)
  message("Saving data...")
  saveRDS(cm, file.path(inferCNV_analysis_folder, "cm_exp_ctrl.rds"))
  saveRDS(samples, file.path(inferCNV_analysis_folder, "samples_exp_ctrl.rds"))
}

#addNormalControl <- function(ext_ctrl, inferCNV_analysis_folder,
#                             cm_raw, orig_samples,
#                             normal_cm_name="K27M_20190124normal.expr.txt"){
#  message("Loading normal control count matrix...")
#  ext_ctrl_cm = read.table(paste0(ext_ctrl, normal_cm_name), sep = "\t", header=T, stringsAsFactors = F, row.names = 1)
#  ext_ctrl_cm = subset(ext_ctrl_cm, select=-c(id, coding))
#  colnames(ext_ctrl_cm) = cleanCellName(colnames(ext_ctrl_cm))
#
#  ## Read and process id of the ctrl data
#  message("Loading normal control annotations...")
#  microglia_cell = readLines(paste0(ext_ctrl, "Microglia.txt"))
#  oc_cell = readLines(paste0(ext_ctrl, "Oligodendrocytes.txt"))
#  microglia_cell = cleanCellName(microglia_cell)
#  oc_cell = cleanCellName(oc_cell)
#  normal_samples = c(rep("mg", length(microglia_cell)), rep("od", length(oc_cell)))
#  names(normal_samples) = c(microglia_cell, oc_cell)
#
#  ## Reorder the ext_ctrl cm and subset genes
#  message("Reordering cell names for normal control cm and annotations")
#  ext_ctrl_cm = ext_ctrl_cm[,names(normal_samples)]
#  ext_ctrl_cm = ext_ctrl_cm[rownames(cm_raw),]
#
#  ## Concatenate samples and normal controls
#  message("Concatenate cm and annotations of normal controls to patient data...")
#  cm = cbind(cm_raw, ext_ctrl_cm)
#  samples = c(orig_samples, normal_samples)
#  message("Saving data...")
#  saveRDS(cm, file.path(inferCNV_analysis_folder, "cm_exp_ctrl.rds"))
#  saveRDS(samples, file.path(inferCNV_analysis_folder, "samples_exp_ctrl.rds"))
#}

#addNormalControlNuc <- function(ext_ctrl, inferCNV_analysis_folder, cm_raw,
#                                orig_samples, premrna = F){
#  ## Read ctrl cm
#  message("Loading normal control count matrix...")
#  if (premrna){
#    suffix = "_premrna_counts.rds"
#  }else{
#    suffix = ".rds"
#  }
#  mg_cm = readRDS(paste0(ext_ctrl, "mg_nuc", suffix))
#  ##tc_cm = readRDS(file.path(ext_ctrl, "tc_nuc", suffix))
#  od_cm = readRDS(paste0(ext_ctrl, "od_nuc", suffix))
#  ext_ctrl_cm = cbind(mg_cm, od_cm)
#
#  ## Generate control sample list
#  message("Loading normal control annotations...")
#  normal_samples = c(rep("mg", dim(mg_cm)[2]), rep("od", dim(od_cm)[2]))
#  names(normal_samples) = c(colnames(mg_cm), colnames(od_cm))
#
#  ## Reorder the ext_ctrl cm and subset genes
#  message("Reordering cell names for normal control cm and annotations")
#  ext_ctrl_cm = ext_ctrl_cm[,names(normal_samples)]
#  ext_ctrl_cm = ext_ctrl_cm[rownames(cm_raw),]
#
#  ## Concatenate samples and normal controls
#  message("Concatenate cm and annotations of normal controls to patient data...")
#  cm = cbind(cm_raw, ext_ctrl_cm)
#  samples = c(orig_samples, normal_samples)
#  message("Saving data...")
#  saveRDS(cm, file.path(inferCNV_analysis_folder, "cm_exp_ctrl.rds"))
#  saveRDS(samples, file.path(inferCNV_analysis_folder, "samples_exp_ctrl.rds"))
#}

## Compute score for oligodendrocytes or immune cells
## @param cm_norm count matrix with log transformed but uncentered data
## @param cm_center count matrix with log transformed and centered data
## @param macrophage_marker marker genes for microglias
## @param t_cell_marker marker genes for T cells
## @param oc_marker markger genes for oligodendrocytes
## @param simple whether use simple average, default = FALSE
## @return a df with scores for these cell types
## How to deal with marker genes not found in this dataset?
computeNormalCellScore <- function(cm_norm, cm_center, gene_mean, scoreName=c("microglia","t_cell","oc"),
                                   microglia_marker=c("CD14", "AIF1", "FCER1G", "FCGR3A", "TROBP", "CSF1R"),
                                   t_cell_marker=c("CD2", "CD3D", "CD3E", "CD3G"),
                                   oc_marker=c("MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11"),
                                   simple=FALSE,
                                   verbose=TRUE){
  ## Check if raw and normalized data have the same cell names
  try (if(any(rownames(cm_norm) != rownames(cm_center))) stop("normalized and centered cm should have identical cell names and orders"))

  ## Compute mean gene expression
  ##message("Compute mean gene expressions...")
  ##gene_mean = rowMeans(cm_norm)

  ## Compute each cell score
  message("Compute single cell scores...")
  microglia_score = scoreSignature(cm_center, gene_mean, microglia_marker, n=100, simple=simple, verbose=verbose)
  t_cell_score = scoreSignature(cm_center, gene_mean, t_cell_marker, n=100, simple=simple, verbose=verbose)
  oc_score = scoreSignature(cm_center, gene_mean, oc_marker, n=100, simple=simple, verbose=verbose)

  ## cbind all scores into a df and result
  message("Concatenate scores into a df...")
  result = data.frame(cbind(microglia_score, t_cell_score, oc_score))
  colnames(result) = c("mg", "tc", "od")
  rownames(result) = colnames(cm_norm)
  return(result)
}

## Add sample annotation
## @param sc_scores single cell scores for microglia, t-cells and oligodendrocytes
## @param sample_names the vector of all patient/sample names
## @param threshold single cell threshold to call a cell any cell type
## @returns a df with rownames = cell_name and one column = malignancy status
generateMalignancyAnnotation <- function(sc_scores, sample_names, useCellType = F, threshold = 4){
  ## Check cell names and orders match
  try (if(any(rownames(sc_scores) != names(sample_names))) stop("sc scores and cell names should have identical cell names and orders"))
  ## Assign maligancy status based on scores and threshold
  if (useCellType){
    sc_scores$malig_status = ifelse((apply(sc_scores, 1, max) < threshold),
                                    paste("malignant", sample_names, sep="_"),
                                    apply(sc_scores, 1, function(x) names(x)[x==max(x)])
    )
  }
  else{
    sc_scores$malig_status = ifelse((apply(sc_scores, 1, max) >= threshold),
                                    "non-malignant",
                                    paste("malignant", sample_names, sep="_"))
  }

  result = sc_scores[,c("malig_status")]
  names(result) = rownames(sc_scores)
  return(data.frame(result))
}

computeCnvScore <- function(cnv_values, squared=T){
  if (squared){
    return (colMeans(cnv_values^2))
  }
  else{
    return (colMeans(abs(cnv_values)))
  }
}

## Compute sample average CNV values
## @param cnv_values cnv_values from outputs of inferCNV
## @param samples a vector of sample names
## @param malig_status malignancy status of each cell (reuse inputs for inferCNV)
## @param min_tumor_cell minium number of tumor cells required to calculate sample mean CNV value
calculateAverageCnvValue <- function(cnv_values, samples, malig_status, min_tumor_cell=3){
  try (if(any(colnames(cnv_values) != names(samples) | any(colnames(cnv_values) != rownames(malig_status))))
    stop("CNV values, sample names, and malignancy status should have identical cell names and orders"))

  result = list()
  unique_sample_names = names(table(samples))
  message("Subsetting cnv_values and samples with only expressionally defined tumor cells...")
  ##cnv_values_tumor = cnv_values[,which(malig_status$result != "non-malignant")]
  putative_malig = sapply(malig_status, function(x) grepl("^malignant", x))
  cnv_values_tumor = cnv_values[,putative_malig]
  ##samples = samples[malig_status$result != "non-malignant"]
  samples = samples[putative_malig]

  try(if(dim(cnv_values_tumor)[2] != length(samples)) stop("Subsetted cnv values and sample names should have equal length"))

  message("Computing average CNV values for each sample...")
  for (sample_name in unique_sample_names){
    tmp_cnv_values = cnv_values_tumor[,which(samples == sample_name)]
    message("The number of putative tumor cells in sample ", sample_name, " is ", dim(tmp_cnv_values)[2])
    if (dim(tmp_cnv_values)[2] >= min_tumor_cell){
      result[[sample_name]] = rowMeans(tmp_cnv_values)
    }
    else{
      result[[sample_name]] = NA
    }
  }
  return(result)
}

## Compute correlation between CNV values of each cell and average CNV values of all putative tumor cells from each sample
## @param cnv_values cnv_values from outputs of inferCNV
## @param mean_cnv_values average CNV values of all putative tumor cells from a sample
## @param samples a vector of sample names
calculateCnvCor <- function(cnv_values, mean_cnv_values, samples){
  try (if(any(colnames(cnv_values) != names(samples)))
    stop("CNV values and sample names should have identical cell names and orders"))

  cor_coef = NULL
  for (i in seq(dim(cnv_values)[2])){
    sample = as.character(samples[i])
    #print(sample)
    #print(head(mean_cnv_values$sample))
    if (sum(is.na(mean_cnv_values[[sample]]))>0){
      cor_coef = c(cor_coef, 0)
    }else{
      ##cor_coef = c(cor_coef, cor(cnv_values[,i], mean_cnv_values[[sample]], method = "spearman"))
      cor_obj = cor.test(cnv_values[,i], mean_cnv_values[[sample]],
                         alternative = "two.sided", method = "spearman", exact = FALSE)
      if (cor_obj$p.value > 0.05){
        cor_coef = c(cor_coef, 0)
      }
      else{
        cor_coef = c(cor_coef, cor_obj$estimate)
      }
    }
  }

  names(cor_coef) = names(samples)
  return(cor_coef)
}

## Distinguish malignant and non-malignant cells
## @param cnv_score a vector of cnv scores (output from calculateAverageCnvValue)
## @param cor_coef a vector of correlation coefficient (output from calculateCnvCor)
callTumorCell <- function(cnv_values, samples, malig_status, cnv_score_cutoff=0.03, cnv_coef_cutoff=0.2){
  ## Calculate CNV score
  message("Calculating CNV scores...")
  cnv_score = colMeans(cnv_values^2)
  ## Calculate CNV correlation
  message("Calculating CNV correlations...")
  mean_cnv_value = calculateAverageCnvValue(cnv_values, samples, malig_status)
  cor_coef = calculateCnvCor(cnv_values, mean_cnv_value, samples)
  ## Call malignant cells
  message("Call malignant cells...")
  result = ifelse((cnv_score < cnv_score_cutoff & cor_coef < cnv_coef_cutoff), "non-malignant", "malignant")
  names(result) = names(samples)
  return(result)
}

plot_subcluster_10X <- function(folder, sample, cm, samples, normal_type, cutree_res,
                                gene_order_file, cutoff = 0.1, window_length = 201, HMM=FALSE){
  ## Create dir to store results
  if (!dir.exists(folder)){
    dir.create(folder)
  }

  ## store cell names in each cluster to a list
  group_num = unique(cutree_res[[sample]])
  group_name = lapply(group_num, function(x){
    tmp = names(cutree_res[[sample]][cutree_res[[sample]] == x])
    tmp = tmp[tmp %in% names(samples)[samples == sample]]})
  names(group_name) = group_num
  group_name = group_name[sapply(group_name, length)>1]


  cm_sub = cm[,c(unlist(group_name), names(samples)[samples %in% normal_type])]
  label = c(unlist(sapply(names(group_name), function(x) rep(paste0("malignant_", x), length(group_name[[x]])))),
            unlist(sapply(normal_type, function(x) rep(x, sum(samples==x)))))
  names(label) = colnames(cm_sub)

  write.table(as.matrix(cm_sub), file = file.path(folder, "counts.matrix"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(label, file = file.path(folder, "cellAnnotations.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=file.path(folder, "counts.matrix"),
                                      annotations_file=file.path(folder, "cellAnnotations.txt"),
                                      delim="\t",
                                      gene_order_file=gene_order_file,
                                      ref_group_names=normal_type)

  # perform infercnv operations to reveal cnv signal
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=cutoff,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir=file.path(folder, "out_dir"),  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               hclust_method="ward.D2",
                               window_length = window_length,
                               plot_steps = F,
                               denoise=T,
                               HMM=HMM
  )

}

call_cnv_hc <- function(cnv_cutree_table, cnv_cutree_each_sample, samples, obs_ref, inferCNV_analysis_folder ,norm_cutoff=3){
  ## Identify clusters of normal cells in each cluster
  ## Defined as group of samples that clusters with >= norm_cutoff normal control cells
  normal_clusters = lapply(cnv_cutree_table, function(x) which((x["mg",] + x["od",])>=norm_cutoff))
  try (if(any(names(cnv_cutree_each_sample) != names(normal_clusters))) stop("cutree results list and normal clusters list should have identical names for elements"))

  ## Assign CNV to cells based on whether cells in normal clusters or not
  call_cnv_res = NULL
  for (i in seq(length(normal_clusters))){
    ## Identify current sample and subset normal clusters and cutree results for that sample
    sample = names(normal_clusters)[i]
    normal_clust = normal_clusters[[i]]
    hc_clust_res = cnv_cutree_each_sample[[i]]
    ## Note: this way of extracting sample names exclude normal controls from the final result
    sample_names = sapply(names(hc_clust_res), function(x) unlist(strsplit(x, split = ".", fixed = T))[1])
    ## For each sample, only keep cells in that sample, not normal controls
    hc_clust_res = hc_clust_res[sample_names == sample]
    tmp = ifelse(hc_clust_res %in% normal_clust, F, T)
    names(tmp) = names(hc_clust_res)[sample_names == sample]
    call_cnv_res = c(call_cnv_res, tmp)
  }

  ## Make a vector of mg and od; they are both normal cells
  mg_names = names(samples)[samples == "mg"]
  mg_res = rep(F, length(mg_names))
  names(mg_res) = mg_names
  od_names = names(samples)[samples == "od"]
  od_res = rep(F, length(od_names))
  names(od_res) = od_names

  ## Concatenate samples and normal controls into a single vector
  call_cnv_res_complete = c(call_cnv_res, mg_res, od_res)
  call_cnv_res_complete = call_cnv_res_complete[names(samples)]
  call_cnv_res = call_cnv_res_complete[names(samples)[obs_ref == "obs"]]

  ## Save CNV results (w/o normal controls)
  saveRDS(call_cnv_res, paste0(inferCNV_analysis_folder, "call_cnv.rds"))
  saveRDS(call_cnv_res_complete, paste0(inferCNV_analysis_folder, "call_cnv_w_ctrl.rds"))
}

# plot_subcluster <- function(folder, sample, group_1, group_2, group_1_num, group_2_num,
#                             cm, samples, gene_order_file, cnv_hc_each_sample, window_length = 201){
#     if (!dir.exists(folder)){
#         dir.create(folder)
#     }
#
#     group_1 = cnv_hc_each_sample[[sample]]$labels[cnv_cutree_each_sample[[sample]] == group_1][1:group_1_num]
#     group_2 = cnv_hc_each_sample[[sample]]$labels[cnv_cutree_each_sample[[sample]] == group_2][1:group_2_num]
#
#     cm = cm[,c(group_1, group_2, names(samples)[samples=="od" | samples == "mg"])]
#     label = c(rep("malignant_1", length(group_1)), rep("malignant_2", length(group_2)), rep("mg", 96), rep("od", 94))
#     names(label) = colnames(cm)
#
#     write.table(cm, file = paste0(folder, "counts.matrix"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
#     write.table(label, file = paste0(folder, "cellAnnotations.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
#
#     infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(folder, "counts.matrix"),
#                                         annotations_file=paste0(folder, "cellAnnotations.txt"),
#                                         delim="\t",
#                                         gene_order_file=gene_order_file,
#                                         ref_group_names=c("mg", "od"))
#
#     # perform infercnv operations to reveal cnv signal
#     infercnv_obj = infercnv::run(infercnv_obj,
#                                  cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
#                                  out_dir=paste0(folder, "out_dir"),  # dir is auto-created for storing outputs
#                                  cluster_by_groups=T,   # cluster
#                                  hclust_method="ward.D2",
#                                  window_length = window_length,
#                                  plot_steps = F,
#                                  include.spike=T
#     )
# }

# plot_subcluster_nuc <- function(folder, sample, group_1, group_2, group_1_num, group_2_num,
#                                 cm, samples, gene_order_file, cnv_hc_each_sample, window_length = 201){
#     if (!dir.exists(folder)){
#         dir.create(folder)
#     }
#
#     group_1 = cnv_hc_each_sample[[sample]]$labels[cnv_cutree_each_sample[[sample]] == group_1][1:group_1_num]
#     group_2 = cnv_hc_each_sample[[sample]]$labels[cnv_cutree_each_sample[[sample]] == group_2][1:group_2_num]
#
#     cm = cm[,c(group_1, group_2, names(samples)[samples=="od" | samples == "tc" | samples == "mg"])]
#     label = c(rep("malignant_1", length(group_1)), rep("malignant_2", length(group_2)),
#               rep("mg", 243), rep("tc", 22), rep("od", 38))
#     names(label) = colnames(cm)
#
#     write.table(cm, file = paste0(folder, "counts.matrix"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
#     write.table(label, file = paste0(folder, "cellAnnotations.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
#
#     infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(folder, "counts.matrix"),
#                                         annotations_file=paste0(folder, "cellAnnotations.txt"),
#                                         delim="\t",
#                                         gene_order_file=gene_order_file,
#                                         ref_group_names=c("mg", "tc", "od"))
#
#     # perform infercnv operations to reveal cnv signal
#     infercnv_obj = infercnv::run(infercnv_obj,
#                                  cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
#                                  out_dir=paste0(folder, "out_dir"),  # dir is auto-created for storing outputs
#                                  cluster_by_groups=T,   # cluster
#                                  hclust_method="ward.D2",
#                                  window_length = window_length,
#                                  plot_steps = F,
#                                  include.spike=T
#     )
#
#     infercnv::plot_cnv(infercnv_obj,
#                        out_dir=paste0(folder, "out_dir"),
#                        cluster_by_groups=T,
#                        color_safe_pal=FALSE,
#                        x.center=1,
#                        x.range=c(0.5,1.5),
#                        title="inferCNV",
#                        obs_title="Observations (Cells)",
#                        ref_title="References (Cells)",
#                        output_filename="heatmap",
#                        hclust_method = "ward.D2",
#                        write_expr_matrix = F
#     )
# }


save_cnv_files <- function(infercnv_obj, k_obs_groups = 1, cluster_by_groups = T, cluster_references = T,
                           out_dir, x.center = 1, x.range = 'auto', output_filename = "infercnv") {
  plot_data <- infercnv_obj@expr.data

  expr_dat_file <- paste(out_dir, paste("expr.", output_filename, ".dat.gz", sep = ""), sep = "/")
  fwrite(as.data.frame(plot_data), file = expr_dat_file, quote = FALSE, sep = ',', row.names = T)

  quantiles <- quantile(plot_data[plot_data != x.center], c(0.01, 0.99))
  delta <- max(abs(c(x.center - quantiles[1], quantiles[2] - x.center)))
  low_threshold <- x.center - delta
  high_threshold <- x.center + delta
  x.range <- c(low_threshold, high_threshold)
  plot_data[plot_data < low_threshold] <- low_threshold
  plot_data[plot_data > high_threshold] <- high_threshold
  infercnv_obj@expr.data <- plot_data

  contigs <- infercnv_obj@gene_order[['chr']]
  unique_contigs <- unique(contigs)
  n_contig <- length(unique_contigs)
  ct.colors <- colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(n_contig)
  names(ct.colors) <- unique_contigs

  custom_pal <- color.palette(c("darkblue", "white", "darkred"), c(2, 2))

  ref_idx <- NULL
  ref_idx <- unlist(infercnv_obj@reference_grouped_cell_indices)
  ref_idx <- ref_idx[order(ref_idx)]

  contig_tbl <- table(contigs)[unique_contigs]
  col_sep <- cumsum(contig_tbl)
  col_sep <- col_sep[-1 * length(col_sep)]
  contig_labels = names(contig_tbl)
  contig_names = unlist(lapply(contig_labels, function(contig_name) {
    rep(contig_name, contig_tbl[contig_name])
  }))
  grouping_key_coln <- c()
  obs_annotations_names <- names(infercnv_obj@observation_grouped_cell_indices)
  obs_annotations_groups = rep(-1, length(colnames(infercnv_obj@expr.data)))
  names(obs_annotations_groups) = colnames(infercnv_obj@expr.data)
  obs_index_groupings = infercnv_obj@observation_grouped_cell_indices
  counter <- 1
  for (obs_index_group in obs_index_groupings) {
    obs_annotations_groups[obs_index_group] <- counter
    counter <- counter + 1
  }
  obs_annotations_groups <- obs_annotations_groups[-ref_idx]

  nobs <- length(unlist(infercnv_obj@observation_grouped_cell_indices))
  grouping_key_coln[1] <- floor(123/(max(nchar(obs_annotations_names)) + 6))
  name_ref_groups <- names(infercnv_obj@reference_grouped_cell_indices)
  grouping_key_coln[2] <- floor(123/(max(nchar(name_ref_groups)) + 6))

  obs_data <- infercnv_obj@expr.data
  obs_data <- plot_data[, -ref_idx, drop = FALSE]
  obs_data <- t(obs_data)

  ref_data_t <- NULL
  updated_ref_groups <- list()
  current_ref_count <- 1
  current_grp_idx <- 1
  plot_data <- infercnv_obj@expr.data
  ref_groups = infercnv_obj@reference_grouped_cell_indices
  for (ref_grp in ref_groups) {
    ref_data_t <- cbind(ref_data_t, plot_data[, ref_grp,
                                              drop = FALSE])
    updated_ref_groups[[current_grp_idx]] = seq(current_ref_count,
                                                current_ref_count + length(ref_grp) - 1)
    current_ref_count <- current_ref_count + length(ref_grp)
    current_grp_idx <- current_grp_idx + 1
  }
  ref_groups <- updated_ref_groups
  nb_breaks <- 16
  breaksList_t <- seq(x.range[1], x.range[2], length.out = nb_breaks)
  gene_position_breaks = NULL

  ## saving observations file
  save_cnv_observation(file_base_name = out_dir, output_filename_prefix = output_filename, obs_data, infercnv_obj, obs_annotations_names, breaksList = breaksList_t)
  obs_data <- NULL

  ## saving references file
  save_cnv_ref(infercnv_obj = infercnv_obj, ref_data = ref_data_t, file_base_name = out_dir, output_filename_prefix = output_filename, ref_groups = ref_groups, cluster_references = cluster_references, name_ref_groups = name_ref_groups)
}


save_cnv_observation <- function(file_base_name, output_filename_prefix, obs_data, infercnv_obj, obs_annotations_names, breaksList) {
  observation_file_base <- paste(file_base_name, sprintf("%s.observations.txt.gz", output_filename_prefix), sep=.Platform$file.sep)

  # Output dendrogram representation as Newick
  # Need to precompute the dendrogram so we can manipulate
  # it before the heatmap plot
  ## Optionally cluster by a specific contig
  hcl_desc <- "General"
  hcl_group_indices <- seq_len(ncol(obs_data))

  obs_dendrogram <- list()
  ordered_names <- NULL
  isfirst <- TRUE
  hcl_obs_annotations_groups <- vector()
  obs_seps <- c()
  sub_obs_seps <- c()  # never use at this time? available if we want to add splits in the heatmap for subclusters

  for (i in seq_along(obs_annotations_names)) {
    if (!is.null(infercnv_obj@tumor_subclusters$hc[[ obs_annotations_names[i] ]])) {

      obs_dendrogram[[i]] = as.dendrogram(infercnv_obj@tumor_subclusters$hc[[ obs_annotations_names[i] ]])
      ordered_names <- c(ordered_names, infercnv_obj@tumor_subclusters$hc[[ obs_annotations_names[i] ]]$labels[infercnv_obj@tumor_subclusters$hc[[ obs_annotations_names[i] ]]$order])
      obs_seps <- c(obs_seps, length(ordered_names))
      hcl_obs_annotations_groups <- c(hcl_obs_annotations_groups, rep(i, length(infercnv_obj@tumor_subclusters$hc[[ obs_annotations_names[i] ]]$order)))

      if (isfirst) {
        write.tree(as.phylo(infercnv_obj@tumor_subclusters$hc[[ obs_annotations_names[i] ]]),
                   file=paste(file_base_name, sprintf("%s.observations_dendrogram.txt", output_filename_prefix), sep=.Platform$file.sep))
        isfirst <- FALSE
      }
      else {
        write.tree(as.phylo(infercnv_obj@tumor_subclusters$hc[[ obs_annotations_names[i] ]]),
                   file=paste(file_base_name, sprintf("%s.observations_dendrogram.txt", output_filename_prefix), sep=.Platform$file.sep), append=TRUE)
      }
    }
    else { ## should only happen if there is only 1 cell in the group so the clustering method was not able to generate a hclust
      #### actually happens with 2 cells only too, because can't cluster 2
      if ((length(unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]])) == 2) || (length(unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]])) == 1)) {
        if (length(unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]])) == 2) {
          obs_dendrogram[[i]] <- .pairwise_dendrogram(colnames(infercnv_obj@expr.data[, unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]), drop=FALSE]))
        }
        else { ## == 1
          obs_dendrogram[[i]] <- .single_element_dendrogram(colnames(infercnv_obj@expr.data[, unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]), drop=FALSE]))
        }
        ordered_names <- c(ordered_names, colnames(infercnv_obj@expr.data[, unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]), drop=FALSE]))
        obs_seps <- c(obs_seps, length(ordered_names))
        hcl_obs_annotations_groups <- c(hcl_obs_annotations_groups, rep(i, length(unlist(infercnv_obj@tumor_subclusters$subclusters[[obs_annotations_names[i]]]))))
      }
      else {
        flog.error("Unexpected error, should not happen.")
        stop("Error")
      }
    }
  }

  obs_dendrogram <- do.call(merge, obs_dendrogram)
  split_groups <- rep(1, dim(obs_data)[1])
  names(split_groups) <- ordered_names

  for(subtumor in infercnv_obj@tumor_subclusters$subclusters[[ obs_annotations_names[i] ]]) {
    sub_obs_seps <- c(sub_obs_seps, (sub_obs_seps[length(sub_obs_seps)] + length(subtumor)))
  }

  obs_seps <- obs_seps[length(obs_seps)] - obs_seps[(length(obs_seps) - 1):1]

  # Output HCL group membership.
  # Record locations of seperations

  # Make colors based on groupings
  row_groupings <- colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(length(table(split_groups)))[split_groups]
  row_groupings <- cbind(row_groupings, colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(length(table(hcl_obs_annotations_groups)))[hcl_obs_annotations_groups])

  # Make a file of coloring and groupings
  groups_file_name <- file.path(file_base_name, sprintf("%s.observation_groupings.txt", output_filename_prefix))
  # file_groups <- rbind(split_groups,row_groupings)
  file_groups <- cbind(split_groups, row_groupings[,1], hcl_obs_annotations_groups, row_groupings[,2])
  #row.names(file_groups) <- c("Group","Color")
  colnames(file_groups) <- c("Dendrogram Group", "Dendrogram Color", "Annotation Group", "Annotation Color")
  # write.table(t(file_groups), groups_file_name)
  write.table(file_groups, groups_file_name)

  # reorder expression matrix based on dendrogram ordering
  obs_data <- obs_data[ordered_names, ]

  # Remove row/col labels, too cluttered
  # and print.
  orig_row_names <- row.names(obs_data)
  row.names(obs_data) <- rep("", nrow(obs_data))

  heatmap_thresholds_file_name <- file.path(file_base_name, sprintf("%s.heatmap_thresholds.txt", output_filename_prefix))
  write.table(breaksList, heatmap_thresholds_file_name, row.names=FALSE, col.names=FALSE)

  # Write data to file.
  row.names(obs_data) <- orig_row_names
  # Rowv inherits dendrogram, Colv is FALSE
  # rowInd = seq_len(nrow(ref_data)) == everything in normal order
  # colInd = seq_len(ncol(ref_data)) == everything in normal order
  #write.table(as.matrix(t(obs_data[labels(obs_dendrogram),])), file = observation_file_base)
  fwrite(as.data.frame(t(obs_data[labels(obs_dendrogram),])), file = observation_file_base, quote = FALSE, sep = ',', row.names = T)
}


save_cnv_ref <- function(infercnv_obj, ref_data, file_base_name, output_filename_prefix, ref_groups, cluster_references, name_ref_groups) {
  number_references <- ncol(ref_data)
  reference_ylab <- NA
  reference_data_file <- paste(file_base_name, sprintf("%s.references.txt.gz", output_filename_prefix), sep=.Platform$file.sep)

  ref_seps <- c()

  # Handle reference groups
  # If there is more than one reference group, visually break
  # up the groups with a row seperator. Also plot the rows in
  # order so the current groups are shown and seperated.

  ordered_names <- c()

  if (cluster_references) {
    split_groups <- c()
    for (i in seq_along(name_ref_groups)) {
      if (!is.null(infercnv_obj@tumor_subclusters$hc[[ name_ref_groups[i] ]])) {
        ordered_names <- c(ordered_names, infercnv_obj@tumor_subclusters$hc[[ name_ref_groups[i] ]]$labels[infercnv_obj@tumor_subclusters$hc[[ name_ref_groups[i] ]]$order])
        ref_seps <- c(ref_seps, length(ordered_names))
        split_groups <- c(split_groups, rep(i, length(infercnv_obj@tumor_subclusters$hc[[ name_ref_groups[i] ]]$order)))
      }
      else {  ## should only happen if there is only 1 cell in the group so the clustering method was not able to generate a hclust
        #### actually happens with 2 cells only too, because can't cluster 2
        if ((length(unlist(infercnv_obj@tumor_subclusters$subclusters[[name_ref_groups[i]]])) == 2) || (length(unlist(infercnv_obj@tumor_subclusters$subclusters[[name_ref_groups[i]]])) == 1)) {
          ordered_names <- c(ordered_names, colnames(infercnv_obj@expr.data[, unlist(infercnv_obj@tumor_subclusters$subclusters[[ name_ref_groups[i] ]]), drop=FALSE]))
          ref_seps <- c(ref_seps, length(ordered_names))
          split_groups <- c(split_groups, rep(i, length(length(unlist(infercnv_obj@tumor_subclusters$subclusters[[name_ref_groups[i]]])))))
        }
        else {
          stop("Error")
        }
      }
    }
  }
  ref_data <- ref_data[, ordered_names, drop=FALSE]

  # Transpose data.
  ref_data <- t(ref_data)

  # Remove labels if too many.
  ref_orig_names <- row.names(ref_data)
  row.names(ref_data) <- rep("", number_references)

  row_groupings <- as.matrix(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(length(table(split_groups)))[split_groups])
  annotations_legend <- cbind(name_ref_groups, colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(length(name_ref_groups)))

  # Write data to file
  row.names(ref_data) <- ref_orig_names
  ## Rowv is FALSE, Colv is FALSE
  # colInd = seq_len(ncol(ref_data)) == everything in normal order
  fwrite(as.data.frame(t(ref_data)), file = reference_data_file, quote = FALSE, sep = ',', row.names = T)
}
