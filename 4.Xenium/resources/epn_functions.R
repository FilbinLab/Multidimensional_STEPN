# In this file, I store ependymoma project specific functions (used for analysis of xenium data).

# set up the global vars of the epn project, without any sample_specific variables
# home_dir: where the code, and light weight files, like metadata etc are stored
# data_dir: where the heavy files, like raw and processed data are stored
SetUpEpendymomaGlobalVarsGeneral <- function(home_dir = '/home/shk490/Multidimensional_STEPN/4.Xenium', 
                                             data_dir = '/n/data1/dfci/pedonc/filbin/lab/users/shk490/ependymoma'){
  ### get metadata
  metadata <- read_xlsx(glue('{home_dir}/resources/SampleIdentifier.xlsx'))
  
  ### read the panels for zfta-rela and non-canonical subtypes
  marker_genes_zfta_rela <- read_xlsx(glue('{home_dir}/resources/Xenium_panel.xlsx'))
  marker_genes_non_canonical <- read_xlsx(glue('{home_dir}/resources/Xenium_panel_Clusters.xlsx'))
  
  ### get some directories
  preparation_dir <- glue('{data_dir}/analysis/1_preparation')
  preprocessing_dir <- glue('{data_dir}/analysis/2_preprocessing')
  manual_annotation_plots_dir <- glue('{data_dir}/analysis/3_manual_annotation_plots')
  annotation_dir <- glue('{data_dir}/analysis/3_program_annotation')
  cellid_dir <- glue('{annotation_dir}/data')
  niche_dir <- glue('{data_dir}/analysis/4_niche')
  coherence_dir <- glue('{data_dir}/analysis/6_coherence')
  nhood_dir <- glue('{data_dir}/analysis/8_neighborhood')
  manual_annotation_yaml <- glue('{home_dir}/resources/normal_celltypes_manual_annotation.yaml')
  raw_data_dir <- glue('{data_dir}/raw_data/xenium_folders')
  sc_mal_path <- glue('{preparation_dir}/data/seurat_obj_malignant_annotated2.qs') # mal only single cell object path
  sc_full_path <- glue('{preparation_dir}/data/seurat_obj_ST_normal_malig_annotated.qs') # all cells (mal and normal and from all subtypes) object
  sc_reference <- glue('{preparation_dir}/data/Nowakowski_Eze.qs') # single cell reference object path (Nowakowski + Eze)
  zfta_reference <- glue('{data_dir}/analysis/1_preparation/data/ZR_Xenium_projection.qs')
  nc_reference <- glue('{data_dir}/analysis/1_preparation/data/NC_Xenium_projection.qs')
  yap1_reference <- glue('{data_dir}/analysis/1_preparation/data/YAP1_Xenium_projection.qs')
  
  # the malignant celltypes we expect in our data 
  malignant_celltypes <- c("Neuroepithelial-like", "Radial-glia-like", "Embryonic-neuronal-like", "Neuronal-like" ,"Ependymal-like", "MES-like", "Embryonic-like")
  non_malignant_celltypes <- c("T-cells", "Myeloid",  "Endothelial", "Oligodendrocytes", 'Neurons')
  # this variable is required in many of the plotting functions, because we want the celltypes ordered in a particular order
  order_metaprograms <- c("Neuroepithelial-like", "Embryonic-like", "Radial-glia-like", "Embryonic-neuronal-like", 
                          "Neuronal-like" ,"Ependymal-like", "MES-like", "T-cells", "Myeloid",  "Endothelial",  
                          "Oligodendrocytes", 'Neurons')
  
  ### source the color palette and get the colors object
  source(glue('{home_dir}/resources/color_palette.R'))
  colors <- colors_metaprograms_Xenium
  
  # make the list in which to store all the attributes and return the list
  all_vars <- list()
  all_vars$home_dir <- home_dir
  all_vars$data_dir <- data_dir
  all_vars$metadata <- metadata
  all_vars$marker_genes_zfta_rela <- marker_genes_zfta_rela
  all_vars$marker_genes_non_canonical <- marker_genes_non_canonical
  all_vars$preparation_dir <- preparation_dir
  all_vars$preprocessing_dir <- preprocessing_dir
  all_vars$manual_annotation_plots_dir <- manual_annotation_plots_dir
  all_vars$annotation_dir <- annotation_dir
  all_vars$cellid_dir <- cellid_dir
  all_vars$niche_dir <- niche_dir
  all_vars$coherence_dir <- coherence_dir
  all_vars$nhood_dir <- nhood_dir
  all_vars$manual_annotation_yaml <- manual_annotation_yaml
  all_vars$raw_data_dir <- raw_data_dir
  all_vars$sc_mal_path <- sc_mal_path
  all_vars$sc_full_path <- sc_full_path
  all_vars$sc_reference <- sc_reference
  all_vars$zfta_reference <- zfta_reference
  all_vars$nc_reference <- nc_reference
  all_vars$yap1_reference <- yap1_reference
  all_vars$malignant_celltypes <- malignant_celltypes
  all_vars$non_malignant_celltypes <- non_malignant_celltypes
  all_vars$order_metaprograms <- order_metaprograms
  all_vars$colors <- colors
  return (all_vars)
}

# set up the sample-specific (and general) global vars for the ependymoma project
SetUpEpendymomaGlobalVars <- function(sample_name, home_dir = '/home/shk490/Multidimensional_STEPN/4.Xenium', 
                                      data_dir = '/n/data1/dfci/pedonc/filbin/lab/users/shk490/ependymoma'){
  
  # sample_name <- 'STEPN06_Region_1' # for testing code
  # first get the general variables using the SetUpEpendymomaGlobalVarsGeneral function
  all_vars <- SetUpEpendymomaGlobalVarsGeneral(home_dir, data_dir)
  
  # get the subtype of given sample
  metadata <- all_vars$metadata
  metadata_sample <- metadata[metadata[['SampleName']] == sample_name, ] # the subset of metadata specific to the current sample
  diagnosis <- metadata_sample[['Subtype']]
  
  # now according to diagnosis, get the sample specific variables like ref_projection, annotation_file, 
  # preprocessed_file, etc.
  if (diagnosis == "ZFTA-RELA"){ref_projection <- all_vars$zfta_reference}
  if (diagnosis %in% c("ZFTA-Cluster 1", "ZFTA-Cluster 2", "ZFTA-Cluster 3", "ZFTA-Cluster 4")){ref_projection <- all_vars$nc_reference}
  if (diagnosis == "ST-YAP1"){ref_projection <- all_vars$yap1_reference}
  preprocessing_file <- glue('{all_vars$preprocessing_dir}/{sample_name}/seurat.qs')
  annotation_file <- glue('{all_vars$annotation_dir}/data/{sample_name}.qs')
  niche_file <- glue('{all_vars$niche_dir}/data/{sample_name}.qs')
  
  # get the raw data directory for this sample
  xenium_folder <- metadata_sample[['RawDataPath']]
  sample_raw_data_dir <- glue('{all_vars$raw_data_dir}/{xenium_folder}')
  h5_file_cell_feature_mat <- glue('{sample_raw_data_dir}/cell_feature_matrix.h5')
  
  # get the marker genes list, filtered and processed appropriately
  if (diagnosis == 'ZFTA-RELA'){
    marker_genes_list <- ReadMarkerGenesEPN(all_vars$marker_genes_zfta_rela, metadata_sample)
  } else{
    marker_genes_list <- ReadMarkerGenesEPN(all_vars$marker_genes_non_canonical, metadata_sample)
  }
  
  # add the sample-specific information obtained in this function to all_vars and return it
  all_vars$metadata_sample <- metadata_sample
  all_vars$diagnosis <- diagnosis
  all_vars$ref_projection <- ref_projection
  all_vars$preprocessing_file <- preprocessing_file
  all_vars$annotation_file <- annotation_file
  all_vars$niche_file <- niche_file
  all_vars$sample_raw_data_dir <- sample_raw_data_dir
  all_vars$h5_file_cell_feature_mat <- h5_file_cell_feature_mat
  all_vars$marker_genes_list <- marker_genes_list
  return (all_vars)
}

# ReadMarkerGenesEPN is the ReadMarkerGenes function specific to ependymoma project. It returns a list which maps 
# the celltypes to a vector of the genes which are marker for it
# NOTE: Here, we are assuming that the annotation column is called Annotation_Sara
ReadMarkerGenesEPN <- function(marker_genes, metadata){
  
  # replace the spaces in the Annotation column with underscores
  marker_genes$Annotation_Sara <- gsub(" ", "_", marker_genes$Annotation_Sara)
  
  # just limit to the columns we are interested in
  marker_genes <- marker_genes[, c('Gene', 'Annotation_Sara')]
  
  # if the sample doesn't have zfta-rela fusion genes, then remove them from marker_genes 
  if (metadata$ZFTAFusionPositive == 0){
    marker_genes <- marker_genes %>%
      filter(!(Gene %in% c('ZFTA_RELA_Fusion1', 'ZFTA_RELA_Fusion2', 'ZFTA_RELA_Fusion3')))
  }
  
  # extract the list for all markers as it will also be used
  marker_genes_list <- marker_genes %>%
    filter(!is.na(Annotation_Sara)) %>% # remove the NA rows (because these are those genes, which didn't have a label in excel)
    # do some renaming of Annotation values
    mutate(Annotation_Sara = str_replace(Annotation_Sara, 'Neurons.*', 'Neurons')) %>%
    mutate(Annotation_Sara = str_replace(Annotation_Sara, 'Microglia', 'Myeloid')) %>%
    mutate(Annotation_Sara = str_replace(Annotation_Sara, 'Endothelial.*', 'Endothelial')) %>%
    mutate(Annotation_Sara = str_replace(Annotation_Sara, 'Dendritic_cell', 'Myeloid')) %>%
    dplyr::group_by(Annotation_Sara) %>%
    dplyr::summarise(Gene_list = list(Gene)) %>%
    deframe()
  
  # remove astrocyte marker genes from the list because in ependymoma, we don't expect much astrocytes and 
  # with these genes, we were observing many astrocytes where we expect cancer cells
  marker_genes_list$Astrocyte = NULL
  
  return (marker_genes_list)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~ FUNCTIONS WHILE PREPARING SINGLE CELL OBJECT FOR SPATIAL ANALYSIS ~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to process a single cell object by performing the standard steps like pca, clustering etc.
RunFullSeurat_v5 <- function(cm, metadata, do.scale = FALSE, do.center = TRUE, doBatch = F, var2batch = NULL, 
                             batchMethod = 'harmony', resolution = 1, dims = 20, pca_dims = 100, n.neighbors = 30, 
                             tsne_perplexity = 30, project, norm.type = 'RNA', verbose = FALSE){
  #seurat_obj <- preprocessSeuratObject_v5(cm, project = project, min.cells = 0, min.genes = 0, scale.factor = 1E5, do.scale = do.scale, do.center = do.center, norm.type = norm.type, mt.pattern = mt.pattern, verbose = verbose)
  seurat_obj <- preprocessSeuratObject_v5(cm, project = project, min.cells = 0, min.genes = 0, scale.factor = 1E5, do.scale = do.scale, do.center = do.center, norm.type = norm.type, verbose = verbose)
  
  if (is.data.frame(metadata) | is.matrix(metadata)) {
    seurat_obj <- AddMetaData(seurat_obj, metadata)
  } else if (is.vector(metadata)) {
    seurat_obj <- AddMetaData(seurat_obj, metadata, col.name = 'sample')
  }
  
  if (norm.type == 'RNA') {
    seurat_obj <- RunPCA(object = seurat_obj, features = VariableFeatures(seurat_obj), npcs = pca_dims, ndims.print  = 1:5, nfeatures.print = 5, verbose = verbose)
    if (dims <= 1) {
      dims <- optimizePCA(seurat_obj, dims)
      #message("Using ", dims, " PCs as optimal number!")
    }
  }
  
  if (doBatch){
    seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj@meta.data %>% pull(var2batch))
    
    if (norm.type == 'RNA') {
      seurat_obj <- suppressWarnings(NormalizeData(seurat_obj, verbose = verbose)) %>%
        FindVariableFeatures(verbose = verbose) %>%
        ScaleData(verbose = verbose) %>%
        RunPCA(npcs = pca_dims, verbose = verbose) %>%
        FindNeighbors(dims = 1:dims, reduction = "pca", verbose = verbose) %>%
        FindClusters(resolution = resolution, cluster.name = "unintegrated_clusters", verbose = verbose) %>%
        RunUMAP(dims = 1:dims, reduction = "pca", reduction.name = "umap.unintegrated", min.dist = 0.8, verbose = verbose, return.model = TRUE) %>%
        RunTSNE(dims = 1:dims, reduction = "pca", reduction.name = "tsne.unintegrated", num_threads = 8, perplexity = tsne_perplexity, check_duplicates = F, verbose = verbose)
    } else if (norm.type == 'SCT') {
      seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2", method = "glmGamPoi", verbose = verbose)
      seurat_obj <- RunPCA(seurat_obj, npcs = pca_dims, verbose = verbose)
      if (dims <= 1) {dims <- optimizePCA(seurat_obj, dims)}
      seurat_obj <- seurat_obj %>%
        FindNeighbors(reduction = "pca", dims = 1:dims, verbose = verbose) %>%
        FindClusters(resolution = resolution, cluster.name = "unintegrated_clusters", verbose = verbose) %>%
        RunUMAP(reduction = "pca", dims = 1:dims, reduction.name = "umap.unintegrated", min.dist = 0.8, return.model = TRUE, verbose = verbose) %>%
        RunTSNE(dims = 1:dims, reduction = "pca", reduction.name = "tsne.unintegrated", num_threads = 8, perplexity = tsne_perplexity, check_duplicates = F, verbose = verbose)
    }
    
    if (batchMethod != 'all') {
      message(glue("Performing batch correction using {batchMethod}"))
      met <- switch(batchMethod, cca = CCAIntegration, rpca = RPCAIntegration, harmony = HarmonyIntegration, mnn = FastMNNIntegration, NA)
      reduc <- switch(batchMethod, cca = 'integrated.cca', rpca = 'integrated.rpca', harmony = 'harmony', mnn = 'integrated.mnn', NA)
      
      #if (sum(seurat_obj@meta.data %>% group_by_at(var2batch) %>% summarise(n = n()) %>% arrange(n) %>% pull(n) < 200) > 0 & batchMethod %in% c('cca', 'rpca')) {
      #  seurat_obj <- IntegrateLayers(object = seurat_obj, method = met, orig.reduction = "pca", new.reduction = reduc, verbose = verbose, k.anchor = 10)
      #} else {
      #  seurat_obj <- IntegrateLayers(object = seurat_obj, method = met, orig.reduction = "pca", new.reduction = reduc, verbose = verbose)
      #}
      
      seurat_obj <- IntegrateLayers(object = seurat_obj, method = met, orig.reduction = "pca", new.reduction = reduc, verbose = verbose)
      seurat_obj <- seurat_obj %>%
        FindNeighbors(reduction = reduc, dims = 1:dims, verbose = verbose) %>%
        FindClusters(resolution = resolution, cluster.name = glue("{batchMethod}_clusters"), verbose = verbose) %>%
        RunUMAP(reduction = reduc, dims = 1:dims, reduction.name = glue("umap.{batchMethod}"), min.dist = 0.8, verbose = verbose, return.model = TRUE) %>%
        RunTSNE(reduction = reduc, dims = 1:dims, reduction.name = glue("tsne.{batchMethod}"), num_threads = 8, perplexity = tsne_perplexity, check_duplicates = F, verbose = verbose)
    } else if (batchMethod == 'all') {
      integration_methods <- c('cca', 'rpca', 'harmony', 'mnn')
      for (intMet in integration_methods) {
        message(glue("Performing batch correction using {intMet}"))
        met <- switch(intMet, cca = CCAIntegration, rpca = RPCAIntegration, harmony = HarmonyIntegration, mnn = FastMNNIntegration, NA)
        reduc <- switch(intMet, cca = 'integrated.cca', rpca = 'integrated.rpca', harmony = 'harmony', mnn = 'integrated.mnn', NA)
        
        #if (sum(seurat_obj@meta.data %>% group_by_at(var2batch) %>% summarise(n = n()) %>% arrange(n) %>% pull(n) < 200) > 0 & intMet %in% c('cca', 'rpca')) {
        #  seurat_obj <- IntegrateLayers(object = seurat_obj, method = met, orig.reduction = "pca", new.reduction = reduc, verbose = verbose, k.anchor = 10)
        #} else {
        #  seurat_obj <- IntegrateLayers(object = seurat_obj, method = met, orig.reduction = "pca", new.reduction = reduc, verbose = verbose)
        #}
        
        seurat_obj <- IntegrateLayers(object = seurat_obj, method = met, orig.reduction = "pca", new.reduction = reduc, verbose = verbose)
        seurat_obj <- seurat_obj %>%
          FindNeighbors(reduction = reduc, dims = 1:dims, verbose = verbose) %>%
          FindClusters(resolution = resolution, cluster.name = glue("{intMet}_clusters"), verbose = verbose) %>%
          RunUMAP(reduction = reduc, dims = 1:dims, reduction.name = glue("umap.{intMet}"), min.dist = 0.8, verbose = verbose, return.model = TRUE) %>%
          RunTSNE(reduction = reduc, dims = 1:dims, reduction.name = glue("tsne.{intMet}"), num_threads = 8, perplexity = tsne_perplexity, check_duplicates = F, verbose = verbose)
      }
    }
    if (norm.type == 'RNA') {
      seurat_obj <- JoinLayers(seurat_obj)
    } else if (norm.type == 'SCT') {
      seurat_obj <- PrepSCTFindMarkers(seurat_obj, verbose = verbose)
    }
  } else {
    if (norm.type == 'RNA') {
      seurat_obj <- seurat_obj %>%
        FindNeighbors(reduction = 'pca', dims = 1:dims, verbose = verbose) %>%
        FindClusters(resolution = resolution, verbose = verbose) %>%
        RunUMAP(reduction = 'pca', dims = 1:dims, n.neighbors = n.neighbors, min.dist = 0.8, verbose = verbose, return.model = TRUE) %>%
        RunTSNE(reduction = 'pca', dims = 1:dims, num_threads = 8, perplexity = tsne_perplexity, check_duplicates = F, verbose = verbose)
    } else if (norm.type == 'SCT') {
      seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2", method = "glmGamPoi", verbose = verbose)
      seurat_obj <- RunPCA(seurat_obj, npcs = pca_dims, verbose = verbose)
      if (dims <= 1) {dims <- optimizePCA(seurat_obj, dims)}
      seurat_obj <- seurat_obj %>%
        FindNeighbors(reduction = "pca", dims = 1:dims, verbose = verbose) %>%
        FindClusters(resolution = resolution, verbose = verbose) %>%
        RunUMAP(reduction = "pca", dims = 1:dims, min.dist = 0.8, return.model = TRUE, verbose = verbose) %>%
        RunTSNE(reduction = 'pca', dims = 1:dims, num_threads = 8, perplexity = tsne_perplexity, check_duplicates = F, verbose = verbose)
    }
  }
  return(seurat_obj)
}

# preprocessing helper function used in RunFullSeurat_v5 function
preprocessSeuratObject_v5 <- function(cm, project, min.cells = 0, min.genes = 0, scale.factor = 1E5, do.scale = F, do.center = T, norm.type = 'RNA', verbose){
  if (norm.type == 'RNA') {
    ## Create Seurat Obj, Normalize, variable genes, scaling data,
    message("Creating Seurat object, normalizing, variable genes and scaling data...")
    suppressWarnings(seurat_obj <- CreateSeuratObject(Matrix::Matrix(as.matrix(cm),sparse = T), min.cells = min.cells, min.features = min.genes, project = project) %>%
                       NormalizeData(normalization.method = "LogNormalize", scale.factor = scale.factor, verbose = verbose) %>%
                       FindVariableFeatures(mean.function = ExpMean, dispersion.function = LogVMR, verbose = verbose) %>%
                       ScaleData(do.scale = do.scale, do.center = do.center, verbose = verbose))
  } else if (norm.type == 'SCT') {
    message("Creating Seurat object...")
    suppressWarnings(seurat_obj <- CreateSeuratObject(Matrix::Matrix(as.matrix(cm),sparse = T), min.cells = min.cells, min.features = min.genes, project = project) %>%
                       NormalizeData(normalization.method = "LogNormalize", scale.factor = scale.factor, verbose = verbose))
  }
  return(seurat_obj)
}

# create the reference from base single cell object. Specifically used in epn project.
CreateReferenceFromBaseSC_EPN <- function(seurat_object, seurat_object_mal, subset){
  
  # for both objects, set the identity as subtype, and subset according to which subtype's reference we are creating
  Idents(seurat_object) <- 'Subtype'
  Idents(seurat_object_mal) <- 'Subtype'
  if (subset == "ZR"){
    seurat_object_subset <- subset(seurat_object, idents = 'ZFTA-RELA')
    seurat_object_subset_mal <- subset(seurat_object_mal, idents = 'ZFTA-RELA')
  }else if (subset == "NC"){
    seurat_object_subset <- subset(seurat_object, idents = c("ZFTA-Cluster 1", "ZFTA-Cluster 2", "ZFTA-Cluster 3", "ZFTA-Cluster 4"))
    seurat_object_subset_mal <- subset(seurat_object_mal, idents = c("ZFTA-Cluster 1", "ZFTA-Cluster 2", "ZFTA-Cluster 3", "ZFTA-Cluster 4"))
  }else if (subset == "YAP1"){
    seurat_object_subset <- subset(seurat_object, idents = 'ST-YAP1')
    seurat_object_subset_mal <- subset(seurat_object_mal, idents = 'ST-YAP1')
  }else{
    stop('Wrong value of argument "subset" supplied. It is should be one of "ZR", "NC", "YAP1"')
  }
  
  # subset seurat_object_subset to normal cells only
  Idents(seurat_object_subset) <- 'malignant'
  seurat_object_subset_normal <- subset(seurat_object_subset, idents = 'Non-malignant')
  sort(unique(seurat_object_subset_normal$Sample_deID))
  
  # concatenate the normal and mal subparts
  seurat_object <- merge(seurat_object_subset_normal, seurat_object_subset_mal) # after this point, the expression 
  # information is still split into different layers (which is useful if doing integration). Since we 
  # are not doing that, we rejoin the layers by calling the JoinLayers function
  table(seurat_object$cell_type)
  seurat_object <- JoinLayers(seurat_object)
  
  # remove unclassified
  Idents(seurat_object) <- 'cell_type'
  seurat_object <- subset(seurat_object, idents = 'Unclassified', invert = T)
  
  # reprocess normal+malignant object with RunFullseurat
  metadata <- seurat_object@meta.data
  cm <- seurat_object[["RNA"]]$counts
  cm_norm <- as.matrix(log2(cm/10+1))
  cm_mean <- log2(Matrix::rowMeans(cm)+1)
  cm_center <- cm_norm - rowMeans(cm_norm)
  seurat_object <- RunFullSeurat_v5(cm = cm, metadata = metadata,  doBatch = F,  project = 'EPN')
  
  # make plots and just display
  p1 <- DimPlot(seurat_object, reduction = 'umap', group.by = 'cell_type', cols = colors, label.size = 4, label = T)
  p2 <- DimPlot(seurat_object, reduction = 'umap', group.by = 'Sample_deID', label.size = 4, label = T)
  show(patchwork::wrap_plots(p1, p2, ncol = 2))
  
  # run SCT 
  seurat_object <- SCTransform(seurat_object,  verbose = T)
  
  # return object
  return (seurat_object)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS FOR PREPROCESSING XENIUM DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to preprocess the xenium sequencing output. Steps include: sctransform, pca, umap, clustering
RunFullXeniumEPN <- function(smp, rawDir, outFile) {
  
  ## Preprocessing seurat objects
  message(green('-------------------------------------------------------'))
  message(blue(glue('Starting pre-processing for sample: {yellow$underline$bold(smp)}!!')))
  message(blue('[1/5] Loading Xenium raw files...'))
  
  suppressWarnings(suppressMessages(xenium.obj <- LoadXenium(rawDir, fov = "fov", segmentations = 'cell')))
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
  qsave(xenium.obj, outFile)
  
  # also return the data object, so that we can make some plots visualizing the quality of data
  return (xenium.obj)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS FOR ANNOTATING XENIUM DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Perform annotation on the spatial sample. Some important arguments of it are: i) a yaml file with information about the 
# manual annotation of the normal celltypes ii) path to the single cell reference iii) marker genes list, etc. This is one 
# of the most important functions.
# data_dir is where we save the output data and csv file with cell annotations
# plot_dir is where we save the plots of annotation (dimplots and imagedimplots)
XeniumCellAssignment <- function(data, ref_projection, sample_name, marker_genes_list, malignant_celltypes, 
                                 manual_annotation_yaml, colors, data_dir, plot_dir) {
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[1/6] Add the manual, normal celltype annotations to the seurat obj')))
  data <- AddNormalCellManualAnnotations(data, manual_annotation_yaml, sample_name)
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[2/6] Get annotations from scMapping on the malignant subset')))
  # subset data to just mal cells
  Idents(data) <- 'cell_type'
  data_mal <- subset(data, idents = 'Unknown')
  # read and subset reference to malignant cells
  sc_ref <- qread(ref_projection)
  Idents(sc_ref) <- 'cell_type'
  ref_projection_mal <- subset(sc_ref, idents = intersect(malignant_celltypes, sc_ref$cell_type))
  # get the annotations from the single cell mapping (with mean centering)
  data_mal <- GetAnnotationsFromScMappingEPN(data_mal, ref_projection_mal, mean_center = TRUE)
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[3/6] Also get annotations from panel on the malignant subset, just for plotting')))
  marker_genes_list_mal <- marker_genes_list[intersect(malignant_celltypes, names(marker_genes_list))]
  data_mal <- GetAnnotationsFromSpatialPanelEPN(data_mal, marker_genes_list_mal, mean_center = TRUE)
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[4/6] Combine annotations from the malignant part and normal part')))
  data@meta.data[data@meta.data$cell_type == 'Unknown', 'cell_type'] <- data_mal$predicted_label_snRNAseq # now that we have the annotations for the cancer cells too, we just have to add the annotation information from data_mal to data
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[5/6] Plot  cell assignment')))
  # UMAPs
  p1 <- DimPlot(data_mal, reduction = 'umap', group.by = 'predicted_label_snRNAseq',  cols = colors) + labs(title = 'scMapping labels on the cancer subset')
  p2 <- DimPlot(data_mal, reduction = 'umap', group.by = 'predicted_label_UCell',  cols = colors) + labs(title = 'Panel labels on the cancer subset')
  p3 <- DimPlot(data, reduction = 'umap', group.by = 'cell_type',  cols = colors) + labs(title = 'Final cell type labels')
  p1+p2+p3
  ggsave(file.path(plot_dir, paste0('2_UMAP_classification_', sample_name, '.pdf')), width = 20, height = 6)
  # Spatial Plots
  dot_size <- 0.4 # dot_size suitable for visualizing samples with less cells and those with more cells
  ImageDimPlot(data, group.by = 'cell_type', size = dot_size, border.size = NA, dark.background = F, cols = colors)
  ggsave(file.path(plot_dir, paste0('3_SpatialMaps_', sample_name, '.pdf')), width = 7, height = 6)
  
  message(green('-------------------------------------------------------'))
  message(blue(paste0('[6/6] Save annotated object')))
  # save qs file
  qsave(data, file.path(data_dir, paste0('seurat_obj_', sample_name, '.qs')))
  # export cell ID
  cell_id <- data.frame(rownames(data@meta.data), data@meta.data$cell_type)
  colnames(cell_id)[colnames(cell_id) == "rownames.data.meta.data."] <- "cell_id"
  colnames(cell_id)[colnames(cell_id) == "data.meta.data.cell_type"] <- "group"
  write_csv(cell_id, file.path(data_dir, paste0('cell_ID_', sample_name, ".csv")))
}

# GetAnnotationsFromScMapping function takes in the data, the location of the single cell reference 
# object and using the anchor finding and transfer anchors functions, determines the celltypes for each cell
# mean_center = TRUE means that after we get scores for each cell for each celltype, we will mean center the 
# scores values for a given celltype. 
# NOTE: ENSURE THE COLUMN IN THE REFERENCE SEURAT OBJECT WITH ANNOTATIONS IS CALLED cell_type
GetAnnotationsFromScMappingEPN <- function(data, ref_projection, mean_center = FALSE, sc_annotations_colname = 'predicted_label_snRNAseq'){
  
  # gracefully read ref_projection (depending on if it is the actual object or the path to it)
  seurat_object <- ReadRefProjection(ref_projection)
  
  # make sure the identity of the sc seurat object is 'cell_type'
  Idents(seurat_object) <- 'cell_type'
  
  # transfer anchors from single cell data
  anchors <- FindTransferAnchors(reference = seurat_object, query = data, normalization.method = "SCT", npcs = 50)
  # transfer labels
  predictions <- TransferData(anchorset = anchors, refdata = seurat_object$cell_type, weight.reduction = data[["pca"]], dims = 1:20)
  
  # store in a copy, the neatened version (removed the score.max and predicted.id cols, and renamed colnames appropriately) of predictions df. 
  predictions_copy <- NeatenPredictionsDF(predictions)
  
  if (mean_center == TRUE){
    # mean-center the prediction scores so that celltypes which have high score in general get removed (like neuroepithelial)
    predictions_copy <- MeanCenterDF(predictions_copy, center_by = 'col')
  }
  # get the predicted id according to which col has max value in each row
  predicted_id <- apply(predictions_copy, 1, function(row){
    colnames(predictions_copy)[which.max(row)]
  })
  
  # add the predicted labels to the data and information about each celltypes' prediction scores
  data <- AddMetaData(object = data, metadata = predicted_id, col.name = sc_annotations_colname)
  data <- AddMetaData(object = data, metadata = predictions_copy, col.name = paste0(colnames(predictions_copy), '_sc'))
  return (data)
}

# get the annotations of different celltypes using a list of celltype:marker_genes. Returned data object has a column with predictions 
# for each cell (panel_annotations_colname) and an individual column storing score for each celltype for each cell.
# mean_center = TRUE means that for the scores for any celltype, we will mean-center it, so that all the celltypes 
# are at a level playing field before finding which celltype has highest score for a given cell.
GetAnnotationsFromSpatialPanelEPN <- function(data, marker_genes_list, mean_center = FALSE, panel_annotations_colname = 'predicted_label_UCell'){
  
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

# version of AddNormalAnnotationsAndAnnotateRemainingCells specific to EPN project, which doesnt have step of annotating astrocytes
AddNormalAnnotationsAndAnnotateRemainingCellsEPN <- function(data, manual_annotation_yaml, malignant_celltypes, ref_projection, marker_genes_list, method, mean_center){
  # get the annotations from the yaml file
  data <- AddNormalCellManualAnnotations(data, manual_annotation_yaml, sample_name)
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
    data_mal <- GetAnnotationsFromScMapping(data_mal, ref_projection_mal, mean_center)
    # Combine annotations from the malignant part and normal part
    data@meta.data[data@meta.data$cell_type == 'Unknown', 'cell_type'] <- data_mal$predicted_label_snRNAseq # now that we have the annotations for the cancer cells too, we just have to add the annotation information from data_mal to data
  }else if (method == 'panel'){
    # subset the marker_genes_list to malignant cells
    marker_genes_list_mal <- marker_genes_list[names(marker_genes_list) %in% malignant_celltypes]
    data_mal <- GetAnnotationsFromSpatialPanel(data_mal, marker_genes_list_mal, mean_center)
    # Combine annotations from the malignant part and normal part
    data@meta.data[data@meta.data$cell_type == 'Unknown', 'cell_type'] <- data_mal$predicted_label_UCell # now that we have the annotations for the cancer cells too, we just have to add the annotation information from data_mal to data
  }else{
    stop("Invalid method. Should be either of 'sc_ref' or 'panel'")
  }
  return (data)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS FOR COHERENCE ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Prepare the data to perform Spatial Coherence. It basically modifies the cells_info df slightly by renaming cols etc, 
# and filters some cells and genes and returns the updated cells_info df, and cell_feature_matrix matrix. 
# cell_feature_matrix is the Matrix of raw counts (genes x cells). cells_info is a df with 3 cols: cell_id, 
# x_centroid, y_centroid
prepareData <- function(cell_feature_matrix, cells_info, scale = 1, filter_genes = 10, filter_cells = 10) {
  
  # just modify the cell_feature_matrix slightly - rename the columns and make the cell_ids as rownames
  spatial <- cells_info %>% dplyr::select(c("cell_id", "x_centroid", "y_centroid")) %>%
    column_to_rownames("cell_id") %>%
    rename_all(~c("imagecol", "imagerow"))
  
  # multiply by the scale value, the coordinates of each cell
  spatial[["imagecol"]] = spatial[["imagecol"]] * scale
  spatial[["imagerow"]] = spatial[["imagerow"]] * scale
  
  # subset cell_feature_matrix to just those genes which are expressed in more than filter_genes cells.
  cell_feature_matrix <- cell_feature_matrix[rowSums(cell_feature_matrix) >= filter_genes, ]
  
  # subset cell_feature_matrix to just those cells which have more than filter_cells transcripts. 
  cell_feature_matrix <- cell_feature_matrix[, colSums(cell_feature_matrix) >= filter_cells]
  
  # subset the spatial dataframe to just those cells which satisfied the above criteria
  spatial <- spatial[colnames(cell_feature_matrix), ]
  
  # return the list of spatial df, and the cell_feature_matrix, both filtered to have good cells and genes
  return(list(spatial = spatial, cell_feature_matrix = cell_feature_matrix))
}

# make a grid over the spatial data and get the information about x and y coordinates of each valid square (square in the 
# grid with atleast one cell beneath), and the celltype assigned to that square (the most-common celltype beneath that square)
# norm_data: expression matrix of data (genes x cells), spatial: df storing x and y coordinates of each cell
# variable: vector of celltypes for each cell,          nbins: number of horizontal and vertical bins in the grid
grid_spatial <- function(norm_data, spatial, variable, nbinsrow, nbinscol) {
  library(gplots)
  message(glue('{nbinsrow} by {nbinscol} has this many spots: {nbinsrow*nbinscol}'))
  
  n_squares = nbinsrow * nbinscol
  cell_bcs = colnames(norm_data)
  xs <- as.integer(unname(spatial$imagecol))
  ys <- as.integer(unname(spatial$imagerow))
  
  # make a 2d histogram which takes in two vectors of values in x and y dimensions, and the 
  # number of bins in each dimension, and gives a histogram with the counts of cells (pair 
  # of x and y coordinates) lying in each square of the histogram. We do this to get the x 
  # breaks and y breaks which help in gridding.
  h2 <- hist2d(xs, ys, nbins = c(nbinscol, nbinsrow), show = F)
  grid_counts <- h2$counts
  xedges <- h2$x.breaks
  yedges <- h2$y.breaks
  
  grid_expr <- matrix(data = 0, nrow = n_squares, ncol = nrow(norm_data))
  grid_coords <- matrix(data = 0, nrow = n_squares, ncol = 2)
  grid_cell_counts <- rep(0, n_squares) # how many cells are there in each square
  
  cell_labels <- as.character(variable)
  cell_set = sort(unique(cell_labels))
  cell_info <- matrix(data = 0, nrow = n_squares, ncol = length(cell_set))
  
  pb <- txtProgressBar(min = 0, max = nbinscol, style = 3, width = 50, char = "=")
  for (x_coor in 1:nbinscol) {
    # x_coor <- 60
    # y_coor <- 50
    x_left <- xedges[x_coor]
    x_right <- xedges[x_coor + 1]
    for (y_coor in 1:nbinsrow) {
      # which square position we are in, in the flattened vector of bins
      n <- ((y_coor-1) * nbinscol) + x_coor
      y_down <- yedges[y_coor] # just the value of the current square's bottom edge
      y_up <- yedges[y_coor + 1] # just the value of the current square's top edge
      # get the center value of the current square and add to grid_coords matrix
      grid_coords[n, 1] = (x_right + x_left) / 2
      grid_coords[n, 2] = (y_up + y_down) / 2
      
      # Now determining the cells within the gridded area 
      #if ((i != nbins-1) & (j == nbins-1)) { # top left corner
      if ((x_coor != nbinscol-1) & (y_coor == nbinsrow)) { # top left corner
        x_true = (xs >= x_left) & (xs < x_right)
        y_true = (ys <= y_up) & (ys > y_down)
        #} else if ((i == nbins - 1) & (j != nbins)) { # bottom right corner
      } else if ((x_coor == nbinscol - 1) & (y_coor != nbinsrow - 1)) { # bottom right corner
        x_true = (xs > x_left) & (xs <= x_right)
        y_true = (ys < y_up) & (ys >= y_down)
      } else {   # average case
        x_true = (xs >= x_left) & (xs < x_right) # boolean of length num_cells, having true for cells with x coor within current square's x breaks
        y_true = (ys < y_up) & (ys >= y_down) # boolean of length num_cells, having true for cells with y coor within current square's y breaks
      }
      
      cell_bool = x_true & y_true # boolean of length num_cells, having true for cells with x and y coor within current square's x and y breaks
      grid_cells = cell_bcs[cell_bool] # barcodes of cells which are within current square
      
      grid_cell_counts[n] = length(grid_cells) # how many cells were there in the current square
      
      # Summing the expression across these cells to get the grid expression
      if (length(grid_cells) > 0) {
        if (length(grid_cells) != 1) {
          grid_expr[n,] = rowSums(norm_data[, cell_bool])
        } else {
          grid_expr[n,] = norm_data[, cell_bool]
        }
      }
      
      # If we have cell type information, will record
      if (length(grid_cells) > 0) {
        grid_cell_types <- cell_labels[cell_bool] # vector of celltypes of the cells in the current square
        tmp <- c() # vector of length=cell_set which will store the proportions of each celltype in current square
        for (ct in cell_set) {
          tmp <- c(tmp, length(which(grid_cell_types == ct)) / length(grid_cell_types))
        }
        cell_info[n, ] <- tmp
      }
    }
    setTxtProgressBar(pb, x_coor)
  }
  close(pb)
  # at the end of this loop, we have 
  # grid_coords matrix filled (coordinate of each square), 
  # grid_expr matrix filled (expression of each gene in each square), 
  # grid_cell_counts vector filled (num_cells in each square), 
  # cell_info matrix filled (proportion of each celltype in each square)
  
  # make grid_expr into a df from a matrix with appropriate rownames and colnames
  grid_expr <- grid_expr %>% as.data.frame()
  rownames(grid_expr) <- glue('grid_{1:n_squares}')
  colnames(grid_expr) <- rownames(norm_data)
  
  # in a list, combine all the data and store
  grid_data <- list(expr = grid_expr,
                    image_col = grid_coords[,1],
                    image_row = grid_coords[,2],
                    n_cells = grid_cell_counts,
                    grid_coords = grid_coords)
  
  # make cell_info into a df from a matrix with appropriate rownames and colnames
  cell_info <- cell_info %>% as.data.frame()
  rownames(cell_info) <- rownames(grid_expr)
  colnames(cell_info) <- cell_set
  
  # get the celltype associated to each square
  max_indices <- apply(cell_info, 1, function(x) names(which.max(x)))
  max_indices <- max_indices[grid_data$n_cells > 0] # limit to just those squares which have atleast one cell inside
  max_indices <- as.character(max_indices) # convert to a vector
  
  # subset the entries in grid data to just those squares which have atleast on cell inside them
  grid_data$expr <- grid_data$expr[grid_data$n_cells > 0, ]
  grid_data$image_col <- grid_data$image_col[grid_data$n_cells > 0]
  grid_data$image_row <- grid_data$image_row[grid_data$n_cells > 0]
  
  # make a dataframe storing the squares' x and y coordinates, and the most prominent celltypes in 
  # them (only for the squares having atleast one cell inside)
  dt <- data.frame(x = grid_data$image_col, y = grid_data$image_row, Metaprogram = max_indices)
  
  return(dt)
}

# modification of function coherence_score. It has option of giving to pixels, a +1 even if their neighbor is empty pixel (equivalent 
# to taking percentage of neighboring cells same as itself)
coherence_score <- function(grid_df, variable, empty_neighbors_increase_coherence = FALSE) {
  # get a dataframe, which will have the squares' y indices as rownames, their x indices as colnames, 
  # and their associated celltype as entry. (There will be some entries with NA too, as not all squares were valid) 
  mat <- grid_df %>% mutate(x = as.character(x)) %>% mutate(y = as.character(y))
  mat <- reshape2::acast(mat, y~x, value.var = variable) %>% as.data.frame()
  # IMPORTANT: here, we should reorder the colnames and rownames to be 'numeric' ascending
  mat <- mat[,as.character(sort(as.numeric(colnames(mat))))]
  mat <- mat[as.character(sort(as.numeric(rownames(mat)))),]
  rownames(mat) <- 1:nrow(mat)
  colnames(mat) <- 1:ncol(mat)
  
  programs <- unique(grid_df[[variable]]) # vector of unique cell_types
  results <- list()
  
  # loop through each program and for each, loop through grid to find the pixels for it and find their coherence
  for (program in programs) {
    tmp <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = list(rownames(mat), colnames(mat)))
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        # if the current square has a celltype associated to it (and not NA)
        if (!is.na(mat[i,j])) {
          if (mat[i,j] == program) {
            # i = 50
            # j = 87
            # mat[i,j]
            # mat[49:51, 86:88]
            # tmp[49:51, 86:88]
            # perform the check for all neighbors (except itself, meaning i,j) and accordingly update tmp
            tmp[i,j] <- tmp[i,j] + CheckIfPixelIsACelltype(i-1, j-1, mat, program, empty_neighbors_increase_coherence)
            tmp[i,j] <- tmp[i,j] + CheckIfPixelIsACelltype(i-1, j, mat, program, empty_neighbors_increase_coherence)
            tmp[i,j] <- tmp[i,j] + CheckIfPixelIsACelltype(i-1, j+1, mat, program, empty_neighbors_increase_coherence)
            
            tmp[i,j] <- tmp[i,j] + CheckIfPixelIsACelltype(i, j-1, mat, program, empty_neighbors_increase_coherence)
            tmp[i,j] <- tmp[i,j] + CheckIfPixelIsACelltype(i, j+1, mat, program, empty_neighbors_increase_coherence)
            
            tmp[i,j] <- tmp[i,j] + CheckIfPixelIsACelltype(i+1, j-1, mat, program, empty_neighbors_increase_coherence)
            tmp[i,j] <- tmp[i,j] + CheckIfPixelIsACelltype(i+1, j, mat, program, empty_neighbors_increase_coherence)
            tmp[i,j] <- tmp[i,j] + CheckIfPixelIsACelltype(i+1, j+1, mat, program, empty_neighbors_increase_coherence)
          }
        }
      }
    }
    results[[program]] <- tmp
  }
  
  counts <- grid_df %>% group_by(Metaprogram) %>% summarise(n = n())
  
  # get the mean coherence for each celltype
  tmp <- data.frame()
  for (metaprogram in counts[[variable]]) {
    if (sum(results[[metaprogram]]) == 0) {
      tmp <- rbind(tmp, data.frame(metaprogram = metaprogram, average = 0))
    } else {
      tmp <- rbind(tmp, data.frame(metaprogram = metaprogram, average = (sum(results[[metaprogram]]))/counts$n[counts[[variable]] == metaprogram]))
    }
  }
  
  # get the mean coherence for whole sample (Earlier, it was computed by simply taking average of mp_avg_coherences, but that is inaccurate)
  total_pixels_with_cells <- sum(counts$n)
  total_coherence_celltype_list <- lapply(results, sum) # list storing celltype:total_coherence
  total_coherence <- sum(as.numeric(total_coherence_celltype_list))
  overall_mean_coherence <- total_coherence/total_pixels_with_cells
  return(list(coherence_score = overall_mean_coherence, coherence_score_program = tmp, results_df = results, counts = counts))
}

# check a specific neighbor of a pixel and perform relevant checks (like being outside of grid), and return 1 if 
# that neighbor is causing coherence of the center pixel to increase (generally meaning its same celltype as center 
# pixel) and 0 otherwise. specific to usage within coherence_score function. Used when we check the neighbors of a 
# pixel in the grid, and assign it a coherence score.
# empty_neighbors_increase_coherence = TRUE means we increase coherence if there is empty neighbor (return 1), else we dont.
# tmp is the matrix in which we are storing the coherence scores for current program
CheckIfPixelIsACelltype <- function(coor1, coor2, mat, program, empty_neighbors_increase_coherence = FALSE){
  # coor1 = 49
  # coor2 = 87
  # ensure its not outside grid. If its outside, simply return tmp
  if (coor1 < 1 | coor1 > nrow(mat)){return (0)}
  if (coor2 < 1 | coor2 > ncol(mat)){return (0)}
  
  # check if mat[coor1, coor2] is empty or non-empty (meaning there was some cell beneath the pixel or not)
  if (!is.na(mat[coor1,coor2])) {
    # if the pixel is not NA, then check its program, if same as given
    if (mat[coor1,coor2] == program){ 
      return (1)
    }
  }else{
    # if it was NA (empty pixel), then we change tmp depending on the argument empty_neighbors_increase_coherence
    if (empty_neighbors_increase_coherence == TRUE){
      return (1)
    }
  }
  return (0)
}

# plot the gridding of the coherence and save at the appropriate location. Can benefit from some reorganization.
gridding_coherence_density <- function(gridding, coherence, sample_name, color_coherence_density, colors_metaprograms_Xenium) {
  
  # make the gridding plot
  p1 <- ggplot(gridding, aes(x = x, y = y, color = Metaprogram)) +
    geom_point(size = 0.2) +
    scale_color_manual(values = colors_metaprograms_Xenium) +
    theme_void() +
    guides(colour = guide_legend(override.aes = list(size = 4)))
  
  mat <- gridding %>% mutate(x = as.character(x)) %>% mutate(y = as.character(y))
  mat <- reshape2::acast(mat, y~x, value.var = 'Metaprogram') %>% as.data.frame()
  # IMPORTANT: reorder the columns and rows numerically and not alphabetically
  mat <- mat[,as.character(sort(as.numeric(colnames(mat))))]
  mat <- mat[as.character(sort(as.numeric(rownames(mat)))),]
  
  # Sum coherence score across all MPs for each position (regions with high coherence will score high)
  summed_matrix <- Reduce(`+`, coherence$results_df)
  
  #change x/y values with actual positions
  rownames(summed_matrix) <- rownames(mat) 
  colnames(summed_matrix) <- colnames(mat) 
  summed_matrix <- reshape2::melt(summed_matrix)
  
  df <- summed_matrix 
  # y comes before because when we reshaped summed_matrix above, the rownames (which had y-values) became the first col after reshaping
  colnames(df) <- c("y", "x", "value")  # Rename columns
  
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
  
  # plot combined
  plot <- patchwork::wrap_plots(p1, p2, ncol = 2)
  return (plot)
}

# get the combined df of the average coherence values of all the samples in sample_names. coherence dir stores the coherence 
# qs files (lists). Returns a df with cols Identifier, coherence, and coherence_celltype corresponding to each celltype in the sample. 
# Identifier col stores SampleNames or SampleIDs depending on average_by_sample_id.
# average_by_sample_id indicates if we want to return the average coherences averaged by replicates. This was required 
# in EPN project. sample_ids is a vector of length same as sample_names which holds the ids for each sample.
# NOTE: this function is highly specific to EPN project.
GetAverageSampleCoherencesGridVersion <- function(sample_names, coherence_dir, average_by_sample_id = FALSE, sample_ids = NULL, fill_na = 0){
  
  # ensure that if average_by_sample_id is TRUE, sample_ids is supplied 
  if (average_by_sample_id == TRUE & is.null(sample_ids)) stop('If average_by_sample_id is TRUE, sample_ids should be supplied!')
  
  # read the coherence list of each sample and append into df with identifier col called Identifier
  average_coherence <- data.frame()
  metaprogram_wise_coherence <- data.frame() # this will be long form dataframe, with cols: sample_name, metaprogram, average
  metaprogram_grid_proportion_and_counts <- data.frame() # this will be a dataframe holding the information for proportion of squares in the grid of a particular celltype
  for (sample_name in sample_names){
    # sample_name = 'STEPN12_Region_3'
    results <- qread(glue('{coherence_dir}/results_{sample_name}.qs'))
    # get the sample_coherence by averaging the MP's coherences. This was found to more robustly show the correlation between hypoxia and coherence.
    sample_coherence <- mean(results$coherence_score_program$average)
    average_coherence <- rbind(average_coherence, data.frame(Identifier = sample_name, coherence = sample_coherence))
    # now append the metaprogram wise coherence scores to metaprogram_wise_coherence df
    results$coherence_score_program$Identifier <- sample_name # add the SampleName col so that can make the concatenated df into wide form later (else there will be multiple rows with same celltype without any way to identify them)
    metaprogram_wise_coherence <- rbind(metaprogram_wise_coherence, results$coherence_score_program)
    # now append the metaprogram proportion from the grid squares into metaprogram_grid_proportion_and_counts
    metaprogram_grid_proportion_and_counts <- rbind(metaprogram_grid_proportion_and_counts, results$counts %>% mutate(proportions = (100*n)/sum(n)) %>% mutate(Identifier = sample_name))
  }
  # make metaprogram_wise_coherence into wide form, and join to average_coherence
  metaprogram_wise_coherence <- metaprogram_wise_coherence %>% pivot_wider(names_from = 'metaprogram', values_from ='average', names_prefix = 'coherence_', values_fill = fill_na) 
  # perform left join with average_coherence to get one df with all information (sample level coherence avg, and metaprogram level) 
  average_coherence <- average_coherence %>% left_join(metaprogram_wise_coherence, by = 'Identifier')
  
  # make metaprogram_grid_proportion from metaprogram_grid_proportion_and_counts: remove counts information and make into wide form 
  metaprogram_grid_proportion <- metaprogram_grid_proportion_and_counts %>% pivot_wider(names_from = Metaprogram, values_from = proportions, names_prefix = 'pixel_proportion_', id_cols = Identifier, values_fill = fill_na) # specifying id_cols automatically removes the col n, which stores the actual counts of the squares
  # perform left join with average_coherence to get one df with all information 
  average_coherence <- average_coherence %>% left_join(metaprogram_grid_proportion, by = 'Identifier')
  
  # make metaprogram_grid_counts from metaprogram_grid_proportion_and_counts: remove proportions information and make into wide form 
  metaprogram_grid_counts <- metaprogram_grid_proportion_and_counts %>% pivot_wider(names_from = Metaprogram, values_from = n, names_prefix = 'pixel_counts_', id_cols = Identifier, values_fill = fill_na) # specifying id_cols automatically removes the col proportions, which stores the proportions of the squares of a metaprogram
  # perform left join with average_coherence to get one df with all information 
  average_coherence <- average_coherence %>% left_join(metaprogram_grid_counts, by = 'Identifier')
  
  # if we want to average the coherence values by sample id (average over the replicates), then add a new column 
  # for sample_ids and groupby. In this case, the resultant df will have cols SampleID, coherence, instead of SampleName, coherence
  if (average_by_sample_id == TRUE){
    average_coherence$Identifier <- sample_ids
    average_coherence <- average_coherence %>% group_by(Identifier) %>% summarise(across(where(is.numeric), mean)) # where(is.numeric) like functions are selection_helpers and are designed to be used only with functions like across(), select(), rename()
  }
  return (average_coherence)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS FOR PLOTTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make big barplot ordered by coherence scores of samples. Made originally for EPN. Requires some work.
# NOTE: this only processes only samples which are rows in metadata, and plots them in same order
# if we do order_by_coherence = TRUE, then the bars are ordered from left to right in ascending order of coherence, 
# else they are ordered in same order as in metadata, be it average_by_sample_id = TRUE or FALSE (Identifier is sample_name or sample_id)
PlotXeniumMetaprogramsCoherence <- function(metadata, cellid_dir, coherence_dir, order_metaprograms, colors, 
                                            malignant_vec, average_by_sample_id, order_by_coherence = TRUE){
  
  # since we don't want discrepancies between metadata and sample_names, we are extracting sample_names from metadata only. 
  # Hence, one needs to ensure that metadata has only those entries which need to be worked on
  sample_names <- metadata$SampleName
  
  # if we are averaging by sample_id, get a average the metadata, so that we have just one row for each sample_id
  # first limit to just those columns which can be averaged with sample_id
  if (average_by_sample_id == TRUE){
    metadata_processed <- metadata %>% select(SampleID, Subtype, Source)
    # rename to Identifier
    metadata_processed$Identifier <- metadata_processed$SampleID
    metadata_processed$SampleID <- NULL
    # get unique of metadata
    metadata_processed <- unique(metadata_processed)
  } else{
    metadata_processed <- metadata %>% select(SampleName, Subtype, Source)
    # rename to Identifier
    metadata_processed$Identifier <- metadata_processed$SampleName
    metadata_processed$SampleName <- NULL
  }
  
  # get coherence information. Here, loading the grid version coherences, because in ependymoma, it better explained the relationship between coherence sample and MES-like proportion
  coherences <- GetAverageSampleCoherencesGridVersion(metadata$SampleName, coherence_dir, average_by_sample_id, metadata$SampleID)
  metadata_processed <- metadata_processed %>% left_join(coherences, by = 'Identifier')
  # order metadata according to coherence scores
  if (order_by_coherence == TRUE) metadata_processed <- metadata_processed %>% arrange(coherence)
  # get the order in which Identifier col of metadata_processed and metaprogram_proportion will be ordered
  datapoints_order <- unique(metadata_processed$Identifier)
  # now metadata_processed$Identifier will be a factor with levels same as datapoints_order
  metadata_processed$Identifier <- factor(metadata_processed$Identifier, levels = datapoints_order)
  
  # get the long df which has proportions of the different metaprograms in each sample_name/sample_id. (cols = Identifier, group, counts, proportions)
  metaprogram_proportion <- GetMetaprogramProportions(sample_names = sample_names, cellid_dir, average_by_sample_id, sample_ids = metadata$SampleID)
  metaprogram_proportion <- left_join(metaprogram_proportion, coherences, by = 'Identifier')
  # make metaprogram_proportion$Identifier as be a factor with levels same as datapoints_order
  metaprogram_proportion$Identifier <- factor(metaprogram_proportion$Identifier, levels = datapoints_order)
  
  # add metadata info in metaprogram_proportion by left joining by the Identifier column 
  metaprogram_proportion <- left_join(metaprogram_proportion, metadata_processed, by = 'Identifier')
  
  # reorder metaprograms so that it is the same as in the paper
  metaprogram_proportion$group <- factor(metaprogram_proportion$group, levels = order_metaprograms)
  # add information about malignant/malignant celltype
  metaprogram_proportion$malignant <- ifelse(metaprogram_proportion$group %in% malignant_vec, "Malignant", "Non-malignant")
  # trick for plotting negative axis
  metaprogram_proportion <- metaprogram_proportion %>% 
    mutate(proportions = if_else(malignant == "Non-malignant", -proportions, proportions)) 
  
  # Plot a stacked bar plot. Note that the order in which the bars are is according to the order of levels in Identifier col of df, which should be of factor type now
  p1 <- ggplot(metaprogram_proportion, aes(x = Identifier, y = proportions, fill = group)) +
    scale_fill_manual(values = colors, name = 'cell_type') +
    geom_col() + 
    theme(panel.spacing.x = unit(0, "mm")) +
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +
    scale_y_continuous(breaks = seq(-100, 100, 20), 
                       labels = abs(seq(-100, 100, 20))) +
    labs(y = 'Proportion, %', x = '') + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          # legend.position = 'bottom'
    )
  
  # plot patient info
  p2 <- metadata_processed %>%
    ggplot(aes(x = Identifier, y = 1)) +
    geom_tile(aes(fill = Subtype), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = col_subtype , na.value = "grey95", guide = guide_legend(ncol = 2))  +
    theme_void() +
    theme(axis.title.y = element_text(),
          # legend.position = 'bottom'
    ) +
    labs(y = 'Subtype')
  
  p3 <- metadata_processed %>%
    ggplot(aes(x = Identifier, y = 1)) +
    geom_tile(aes(fill = Source), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = col_sampling, na.value = "grey95", guide = guide_legend(ncol = 2))  +
    theme_void() +
    theme(axis.title.y = element_text(), 
          # legend.position = 'bottom'
    ) +
    labs(y = 'Source')
  
  p4 <- metadata_processed %>%
    ggplot(aes(x = Identifier, y = 1)) +
    geom_tile(aes(fill = coherence), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_gradientn(colours = color_coherence_density_without_white) + 
    theme_void() +
    theme(axis.title.y = element_text(), 
          axis.text.x= element_text(color = 'black', hjust=1, angle=90)) +
    labs(y = 'Coherence')
  
  # plot combined plot. Editing the line below allows us to either plot the subtype information about each sample, or the source information
  plot <- patchwork::wrap_plots(list(p1, p3, p4), ncol = 1) + plot_layout(heights = c(1, 0.1, 0.1), guides = 'collect')
  return (plot)
}

# read in the coherence information for the samples and proportions of different celltypes in them, and 
# compute linear regression and make the corresponding plots.
# NOTE: Ensure that metadata has only those entries which need to be worked on
# On 10-21, added new argument transform_proportions, which specifies if we want to transform the proportion values 
# before the LR and making the plot. This was required for the second review when we performed a multiple linear 
# regression using all celltypes
OverallCoherenceCelltypeProportionLinearRegression <- function(metadata, cellid_dir, coherence_dir, average_by_sample_id, colors, transform_proportions = F) {
  library(reshape)
  # since we don't want discrepancies between metadata and sample_names, we are extracting sample_names from metadata only. 
  # Hence, one needs to ensure that metadata has only those entries which need to be worked on
  sample_names <- metadata$SampleName
  
  # read proportion of each cell type in each tumor 
  message(green('-------------------------------------------------------'))
  message(blue('[1/3] Reading metaprogram proportions'))
  
  # get the long df which has proportions of the different metaprograms in each sample_name/sample_id. (cols = Identifier, group, counts, proportions)
  metaprogram_proportion <- GetMetaprogramProportions(sample_names = sample_names, cellid_dir, average_by_sample_id, sample_ids = metadata$SampleID)
  
  # obtain a vector storing the metaprograms we have in our data
  cell_types <- sort(unique(metaprogram_proportion$group))
  # limiting to just the relevant columns
  metaprogram_proportion <- metaprogram_proportion[ , c("Identifier", "group", "proportions")]
  # ensure that every combination of Identifier and group is present in the data, and if they were not present till now, their proportions will be 0
  metaprogram_proportion_complete <- metaprogram_proportion %>% 
    complete(Identifier, group, fill = list(proportions = 0))
  # here, according to transform_proportions argument value, we transform the proportion values in metaprogram_proportion_complete using logit transformation
  if (transform_proportions == T){
    library(car)
    metaprogram_proportion_complete$proportions <- logit(metaprogram_proportion_complete$proportions/100, percents = F, adjust = 0.025)
  }
  # cast the just created df into wide form so that we have Identifier x group shape of df.
  metaprogram_proportion_wide <- cast(metaprogram_proportion_complete, Identifier~group, value = 'proportions')
  
  message(green('-------------------------------------------------------'))
  message(blue('[2/3] Get spatial coherence scores'))
  
  # get the df storing coherence values for each SampleName or SampleID. colnames: Identifier, coherence
  average_spatial_coherence_df <- GetAverageSampleCoherencesGridVersion(sample_names, coherence_dir, average_by_sample_id, sample_ids = metadata$SampleID)
  # just keep the whole-sample coherence information as only that is required in this function
  average_spatial_coherence_df <- average_spatial_coherence_df %>% select(Identifier, coherence)
  # calculate scaled score (get values to 0-1 range and add as a new col)
  average_spatial_coherence_df <- average_spatial_coherence_df %>% mutate(scaled_spatial_coherence = MinMaxScaleVector(coherence))
  
  # write.csv(average_spatial_coherence_df, '~/temp.csv')
  
  message(green('-------------------------------------------------------'))
  message(blue('[3/3] Perform linear regression between coherence and each MP'))
  
  # combine the information about celltypes proportions into average_spatial_coherence_df
  linear_regression_df <- average_spatial_coherence_df %>% 
    left_join(metaprogram_proportion_wide, by = 'Identifier')
  # make the cols in linear_regression_df to have dots instead of spaces because as.formula function which we 
  # use below doesn't handle columns with hyphens well
  cell_types <- gsub('-', '.', cell_types)
  colnames(linear_regression_df) <- gsub('-', '.', colnames(linear_regression_df))
  # write.csv(as.data.frame(linear_regression_df), file.path(plot_dir, '13_Linear_regression_input.csv'))
  
  # loop through each predictor (celltype which drives coherence value), perform regression on scaled_spatial_coherence
  # with just that, and save the corresponding statistics in a df, and also save the corresponding linear regression plots.
  
  # Initialize a dataframe to store regression results
  results_df <- data.frame()
  results_df <- data.frame(Predictor = character(), Estimate = numeric(), StdError = numeric(), 
                           tValue = numeric(), pValue = numeric(), stringsAsFactors = FALSE)
  plot_list <- list()
  # Loop over each predictor and perform regression using just that
  for (cell_type in cell_types){
    # cell_type <- 'T.cells'
    # cell_type = 'Embryonic.like'
    # Dynamically create the formula. Make quoted because else as.formula() function splits string open at the hyphen 
    formula <- as.formula(paste("scaled_spatial_coherence ~", cell_type))
    # Fit the linear model
    model <- lm(formula, data = linear_regression_df)
    # Extract coefficients summary
    summary_model <- summary(model)
    coefficients <- summary_model$coefficients
    # printing equation
    slope = coefficients[2, 'Estimate']
    intercept = coefficients[1, 'Estimate']
    if (intercept > 0) print(glue("{cell_type}: Y = {signif(slope,3)}*X + {signif(intercept,3)}")) else print(glue("{cell_type}: Y = {signif(slope,3)}*X - {-signif(intercept,3)}"))
    # Extract R-squared and p-value
    r_squared <- summary_model$r.squared
    p_value <- coefficients[2, "Pr(>|t|)"]  # p-value of the predictor term
    # also obtain the pearson correlation because the reviewer asked for it
    correlation <- cor(x = linear_regression_df[, cell_type], y = linear_regression_df$scaled_spatial_coherence)
    # Append results to dataframe
    results_df <- rbind(results_df, data.frame(Predictor = cell_type, Estimate = coefficients[2, "Estimate"], StdError = coefficients[2, "Std. Error"], 
                                               tValue = coefficients[2, "t value"], pValue = p_value, correlation = correlation))
    # Create and save the plot
    plot_list[[cell_type]] <- ggplot(data.frame(linear_regression_df, celltype = cell_type), aes(x = scaled_spatial_coherence, y = .data[[cell_type]], color = celltype)) +
      # labs(x = "Scaled spatial coherence score", y = cell_type) +
      geom_smooth(method = 'lm', se = TRUE, color = 'black') +
      geom_point(size = 4) +
      scale_color_manual(values = colors) + 
      theme_classic() + 
      theme(legend.position = 'none') +
      labs(title = paste0(cell_type, '\n', 'R = ', round(correlation, 2), ', R² = ', round(r_squared, 2), '\nAdj p-value: ', signif(min(p_value*length(cell_types), 1), 3), ', n = ', nrow(linear_regression_df))) +
      ylab('Transformed Proportion') + xlab('Spatial coherence score') + 
      theme(plot.title = element_text(hjust = 0.5), 
            axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
            axis.text = element_text(size = 8.5),
            axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))
  }
  # Export results to a CSV file
  # write.csv(results_df, file.path(plot_dir, "15_Simple_linear_regression_results.csv"), row.names = FALSE)
  # plot 
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = 4)
  return(combined_plot)
}

# make the set of boxplots showing the distribution of coherence values for each metaprogram in a set of samples
MetaprogramCoherenceBoxplots <- function(metadata, coherence_dir, average_by_sample_id, colors){
  # get the sample names from metadata
  sample_names <- metadata$SampleName
  
  # get the df storing coherence values for each SampleName or SampleID. colnames: Identifier, coherence, 
  # coherence_celltype1, coherence_celltype2, etc. Using the grid version of coherences for the ependymoma 
  # project, as it seemed to work better for it.
  average_spatial_coherence_df <- GetAverageSampleCoherencesGridVersion(sample_names, coherence_dir, average_by_sample_id, sample_ids = metadata$SampleID)
  # now, we make it into a long df, with cols: celltype, coherence.
  # now, we can exclude the sample level coherence as that is not required for making metaprogram-coherence plot
  df <- average_spatial_coherence_df %>% 
    select(!c(Identifier, coherence)) %>% # these cols are not required. Just the cols storing coherence for celltype are required
    pivot_longer(everything(), names_to = 'celltype', values_to = 'coherence')
  
  # rename the values in celltype col
  df$celltype <- gsub('coherence_', '', df$celltype)
  # remove the rows which have pixel_proportion values. They store that out of all the pixels, how many pixels were of a particular celltype.
  df <- df %>% filter(!grepl('pixel', celltype))
  
  # rename the MES-like entries to MES/Hypoxia
  df$celltype <- gsub('MES-like', 'MES/Hypoxia', df$celltype)
  
  # now make the plot
  plot <- ggplot(df, aes(x = fct_reorder(celltype, coherence, .fun = median, .desc = TRUE), y = coherence, fill = celltype)) + 
    geom_boxplot(outliers = F) + 
    geom_jitter(size = 0.5) +
    scale_fill_manual(values = colors) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    xlab('Metaprogram')
  
  return (plot)
}

# get the plots showing correlation of a celltype's proportion with coherence of everything apart from it. 
# Originally made with respect to MES-like celltype
OverallCoherenceCelltypeProportionLinearRegressionAfterRemoval <- function(celltype_to_remove, metadata, cellid_dir, coherence_dir, average_by_sample_id){
  library(reshape)
  # since we don't want discrepancies between metadata and sample_names, we are extracting sample_names from metadata only. 
  # Hence, one needs to ensure that metadata has only those entries which need to be worked on
  sample_names <- metadata$SampleName
  
  # read proportion of each cell type in each tumor 
  message(green('-------------------------------------------------------'))
  message(blue('[1/3] Reading metaprogram proportions'))
  
  # get the long df which has proportions of the different metaprograms in each sample_name/sample_id. (cols = Identifier, group, counts, proportions)
  metaprogram_proportion <- GetMetaprogramProportions(sample_names = sample_names, cellid_dir, average_by_sample_id, sample_ids = metadata$SampleID)
  # limiting to just the relevant columns
  metaprogram_proportion <- metaprogram_proportion[ , c("Identifier", "group", "proportions")]
  # ensure that every combination of Identifier and group is present in the data, and if they were not present till now, their proportions will be 0
  metaprogram_proportion_complete <- metaprogram_proportion %>% 
    complete(Identifier, group, fill = list(proportions = 0))
  # cast the just created df into wide form so that we have Identifier x group shape of df.
  metaprogram_proportion_wide <- cast(metaprogram_proportion_complete, Identifier~group, value = 'proportions')
  
  message(green('-------------------------------------------------------'))
  message(blue('[2/3] Get spatial coherence scores (after removal of celltype to remove)'))
  
  # get the df storing coherence values for each SampleName or SampleID. colnames: Identifier, coherence
  average_spatial_coherence_df <- GetAverageSampleCoherencesAfterRemovingCelltype(celltype_to_remove, sample_names, coherence_dir, average_by_sample_id, sample_ids = metadata$SampleID)
  # calculate scaled score (get values to 0-1 range and add as a new col)
  average_spatial_coherence_df <- average_spatial_coherence_df %>% mutate(coherence = MinMaxScaleVector(coherence)) # calling the new col also coherence, because below, that general name is used of that column
  
  message(green('-------------------------------------------------------'))
  message(blue('[3/3] Organize data for performing Linear Regression'))
  
  # combine the information about celltypes proportions into average_spatial_coherence_df
  linear_regression_df <- average_spatial_coherence_df %>% 
    left_join(metaprogram_proportion_wide, by = 'Identifier')
  # make the cols in linear_regression_df to have dots instead of spaces because as.formula function which we 
  # use below doesn't handle columns with hyphens well
  celltype_to_remove <- gsub('-', '.', celltype_to_remove)
  colnames(linear_regression_df) <- gsub('-', '.', colnames(linear_regression_df))
  # just perform linear regression between the given celltype's proportion and coherence
  result <- PerformLinearRegression(linear_regression_df, col1 = 'coherence', col2 = celltype_to_remove)
  plot <- PlotLinearRegression(linear_regression_df, result, x_label = 'Spatial coherence score \nof everything remaining', y_label = glue('Proportion of {celltype_to_remove}'))
  return (plot)
}

# plot the linear regression plot showing comparison of celltypes proportion in sc and spatial data. 
# metadata should have only those samples which are to be included in the analysis (in EPN project, we had 
# multiple plots for the different cancertypes, hence this point).
# sc_data should have the single cell data (which we need the proportions from only) with annotations in the cell_type col.
PlotCompositionComparisonSpatialSC <- function(sc_data, metadata, cellid_dir, transform_proportions = F){
  # get the sample_names, sample_ids from metadata
  sample_names <- metadata$SampleName
  sample_ids <- metadata$SampleID
  
  ### we just need the overall avg of each celltype in both datasets. Hence, no need to mapping sample to sample
  # get the proportions of the different celltypes from single cell
  sc_props <- sc_data@meta.data %>% group_by(cell_type) %>% summarise(counts_sc = n()) %>% mutate(proportions_sc = 100*counts_sc/sum(counts_sc))
  # get the proportions of the celltypes in the different xenium samples
  metaprogram_proportions <- GetMetaprogramProportions(sample_names, cellid_dir, average_by_sample_id = TRUE, sample_ids)
  spatial_props <- metaprogram_proportions %>% group_by(group) %>% summarise(counts_spatial = sum(counts)) %>% mutate(proportions_spatial = 100*counts_spatial/sum(counts_spatial))
  colnames(spatial_props)[1] <- 'cell_type' # rename so that we can do the join in next step
  
  # perform left join and combine into one df
  linear_regression_df <- spatial_props %>% full_join(sc_props, by = 'cell_type')
  linear_regression_df[is.na(linear_regression_df)] <- 0 # fill the na values by 0
  
  # at this point, according to the value of the parameter transform_proportions, we perform a logit transformation
  # on the proportions_spatial and proportions_sc in proportions in linear_regression_df
  if (transform_proportions == T){
    library(car)
    linear_regression_df$proportions_sc <- logit(linear_regression_df$proportions_sc/100, percent = F, adjust = 0.025)
    linear_regression_df$proportions_spatial <- logit(linear_regression_df$proportions_spatial/100, percent = F, adjust = 0.025)
  }
  
  # now do the regression
  linear_reg_result <- PerformLinearRegression(linear_regression_df, 'proportions_sc', 'proportions_spatial')
  plot <- PlotLinearRegression(linear_regression_df, linear_reg_result, x_label = 'Proportion Smart-seq2', y_label = 'Proportion Xenium', color = colors_metaprograms_Xenium, color_col = 'cell_type')
  return (plot)
}









