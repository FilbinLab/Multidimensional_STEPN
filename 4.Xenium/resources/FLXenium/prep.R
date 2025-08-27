# Function to preprocess the xenium sequencing output (creates seurat object, find markers for 
# different resolutions and score the signatures)
# spatially_subset_coords is a vector storing 4 numerics (xmin, xmax, ymin, ymax) indicating how to spatially subset the sample
# This function both saves the processed file to the specified location and returns the object for the processed file
RunFullXenium <- function(smp, rawDir, outFile, spatially_subset_coords = NULL) {
  
  ## Preprocessing seurat objects
  message(green('-------------------------------------------------------'))
  message(blue(glue('Starting pre-processing for sample: {yellow$underline$bold(smp)}!!')))
  message(blue('[1/5] Loading Xenium raw files...'))
  
  suppressWarnings(suppressMessages(xenium.obj <- LoadXenium(rawDir, fov = "fov", segmentations = 'cell')))
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 5) # atleast 6 transcript should be there in the cell, else we capture badqual or non-cells too
  
  # spatially subset the object if spatially_subset_coords != NULL
  if (!is.null(spatially_subset_coords)){
    xenium.obj <- SubsetToRect(xenium.obj, 
                               xmin = spatially_subset_coords[1], 
                               xmax = spatially_subset_coords[2], 
                               ymin = spatially_subset_coords[3], 
                               ymax = spatially_subset_coords[4])
  }
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
  
  qsave(xenium.obj, outFile)
  
  # also return the data object, so that we can make some plots visualizing the quality of data
  return (xenium.obj)
}

# make the clustering columns (like SCT_snn_res.0.1) in seurat_obj into factor type. Hard-coded to handle 
# resolutions from 0.1-1
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

# get the best number of principle components. Originally written by Carlos.
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


