# Load libraries

## Libraries for basic preprocessing
library(reshape2)
library(dplyr)

## Single cell libraries
library(Seurat)
#library(pagoda2)
##library(conos)

## Libraries for plotting
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)

## Annotation and pathway analysis library
##library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

## Others
library(Rtsne)
library(harmony)
library(parallel)
library(glue)
library(SeuratWrappers)

# Define functions
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


## Functions for preprocessing count matrix

## Function to make gene_id as rownames for a cm
## @para: df (a raw cm still with a gene_id column)
## @returns: df with gene_id removed, and rownames as gene_id
addRownames <- function(df){
  ensembl_id = sapply(df$gene_id, function(x) unlist(strsplit(x, "\\."))[1])
  df$ensembl_id = ensembl_id
  df = distinct(df, ensembl_id, .keep_all = TRUE)
  rownames(df) = df$ensembl_id
  df = subset(df, select=-c(gene_id, ensembl_id))
}

## Function to convert gene_id into gene symbols using biomart
## @para: df, a cm still with gene_id as rownames
## @para: dataset, biomart dataset ("hsapiens_gene_ensembl" or "mmusculus_gene_ensembl")
## @returns: df with gene symbol as rownames
ensembl_to_symbol_biomart <- function(df, dataset = c("hsapiens_gene_ensembl", "mmusculus_gene_ensembl")){
  ## Initiate ensembl db
  ensembl = useMart("ensembl")
  mart <- useDataset(dataset, useMart("ensembl"))
  ## Look up HUGO from ensembl gene id (with version) and remove duplicates
  ids <- rownames(df)
  genes <- getBM(filters="ensembl_gene_id_version", attributes=c("ensembl_gene_id_version", "external_gene_name"), values=ids, mart=mart)
  genes <- genes[which(duplicated(genes$external_gene_name) == FALSE), ]
  ##merged_gene_ids <- merge(x=df, y=genes, by.x="row.names", by.y="ensembl_gene_id")
  ##rownames(merged_gene_ids) = merged_gene_ids$external_gene_name
  ##merged_gene_ids = subset(merged_gene_ids, select=-c(Row.names, external_gene_name))
  ## Make a vector of HUGO, named by ensembl id
  gene_symbols = genes$external_gene_name
  names(gene_symbols) = genes$ensembl_gene_id_version
  ## Remove na and duplicates
  gene_symbols = gene_symbols[!is.na(gene_symbols)]
  gene_symbols = gene_symbols[!duplicated(gene_symbols)]
  ## Convert rownames of df (cm) from ensembl id to HUGO
  df = df[names(gene_symbols),]
  rownames(df) = gene_symbols
  return(df)
}

## Function to convert gene_id into gene symbols using AnnotationDB
## @para: cm, a cm still with gene_id as rownames
ensembl_to_symbol_annotationdb <- function(cm){
  ## Remove duplicated gene ids if version number is removed
  gene_names = sapply(rownames(cm), function(x) unlist(strsplit(x, split = ".", fixed=T))[1])
  names(gene_names) = rownames(cm)
  gene_names = gene_names[!duplicated(gene_names)]
  ## Convert gene_names from gene ids into gene symbols
  cm = cm[names(gene_names),]
  rownames(cm) = gene_names
  gene_symbols = mapIds(org.Hs.eg.db, keys = rownames(cm), column = "SYMBOL", multiVals = "first", keytype = "ENSEMBL")
  ## Remove NA and duplicates from gene_symbols and make them rownames of cm
  gene_symbols = gene_symbols[!is.na(gene_symbols)]
  gene_symbols = gene_symbols[!duplicated(gene_symbols)]
  cm = cm[names(gene_symbols),]
  rownames(cm) = gene_symbols
  return(cm)
}

## Remove "_" and "." from cell names
## @para: cell_name, a vector of cell names
cleanCellName <- function(cell_name){
  cell_name = gsub("_", "", cell_name, fixed = T)
  cell_name = gsub("-", ".", cell_name, fixed = T)
  return(cell_name)
}

## Construct a seurat object
## @para: count_matrix, cm
## @para: mito_symbol, mitocondrial genes in gene name/symbol format
## @para: min.cells, minimal number of cells that express a certain gene for that gene to be considered as expressed (default = 1)
## @para: is.expr, gene expression detection limit; default = 0 (for UMI); reads should use 5
## @para: project, project name for this seurat object
seuratObj <- function(count_matrix, mito_symbol, min.cells = 1, is.expr = 0, project){
  seurat_obj <- CreateSeuratObject(round(count_matrix), min.cells = min.cells, is.expr = is.expr, project = project)
  mito_symbol = mito_symbol[mito_symbol %in% rownames(seurat_obj@raw.data)]
  percent.mito <- Matrix::colSums(seurat_obj@raw.data[mito_symbol, ]) / Matrix::colSums(seurat_obj@raw.data)
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.mito, col.name = "percent.mito")
}

## normalize and filtering tpm dataset
## @param cm count matrix
## @param filter_cell whether to filter cells, default = TRUE
## @param scale_factor scaling factor to shrink tpm, default = 10
## @param gene_detect gene detection limit, default = 1
## @param min_genes minimal number of genes for calling a cell, default = 2500
## @param min_mean_hg minimal mean expression of house keeping genes for calling a cell, default = 2.5
## @param hg_list the list of house keeping genes, default =
## @param gene_cutoff average expression cutoff to keep a gene
## @param centering whether to center filtered gene expressions, default = TRUE
## @param verbose whether to print summary statistics, default = TRUE
normalizeTPM <- function(cm, filter_cell = TRUE, log_base=2, scale_factor=10,
                         gene_detect=1, min_genes=2500,
                         min_mean_hg=2.5, hg_list=NULL,
                         filter_gene_relaxed=F, gene_cutoff=4, gene_cell_cutoff=10, verbose = TRUE){
  cm = as.matrix(cm)
  message("The number of cells in the unfilered cm:", dim(cm)[2])
  if (filter_cell){
    ## Filter cells with low total number of genes detected
    message("Filtering based on total number of detected genes")
    cm = cm[,which(colSums(cm >= gene_detect) >= min_genes)]
    message("Number of cells left after total number of gene filtering: ", dim(cm)[2])
  }
  ## Log2 transform
  message("Scaling and log2 transformation...")
  cm_norm = log(cm/scale_factor+1, base=log_base)
  if (filter_cell){
    ## Filter cells with low house keeping gene expressions
    message("Filtering based on expressions of house keeping genes...")
    hg_list = hg_list[hg_list %in% rownames(cm)]
    hg_filtering = colMeans(cm_norm[hg_list,]) >= min_mean_hg
    cm_norm = cm_norm[,hg_filtering]
    cm = cm[,hg_filtering]
    message("Number of cells left after house-keeping gene filtering: ", dim(cm)[2])
  }
  ## Gene filtering
  message("Filtering gene with low expressions...")
  if (filter_gene_relaxed){
    genes_to_keep = rowSums(cm >= 2^gene_cutoff) >= gene_cell_cutoff
  }else{
    genes_to_keep = log2(rowMeans(cm)+1) >= gene_cutoff
  }
  cm = cm[genes_to_keep,]
  cm_norm = cm_norm[genes_to_keep,]
  message("The number of genes left after gene filtering: ", dim(cm)[1])
  ## Centering
  message("Centering gene expressions...")
  cm_norm_center = cm_norm-rowMeans(cm_norm)
  ## Check # of genes detected per cell
  if (verbose){
    message(paste0("Cell: ", dim(cm_norm)[2], "; ", "Gene: ", dim(cm_norm)[1]))
    message("Printing the number of genes expressed in all cells...")
    cat(summary(colSums(cm >= gene_detect)))
    cat("\n")
  }
  result = list()
  result[["raw_data"]] = cm
  result[["norm_data"]] = cm_norm
  result[["center_data"]] = cm_norm_center
  return (result)
}

## QC
## @param cm count matrix (raw or log transformed)
## @param hg_list list of house-keeping genes
## @param scale_factor factor to scale tpm, default = 10
## @param log_base base of log, default = 2
## @param gene_detect detection limit for gene, default = 1
calculateQcMetrics <- function(cm, hg_list, scale_factor=10, log_base=2, gene_detect=1){
  ## Info for sample, plate and well
  sample = sapply(colnames(cm), function(x) unlist(strsplit(x, split=".", fixed=T))[1])
  plate = sapply(colnames(cm), function(x) unlist(strsplit(x, split=".", fixed=T))[2])
  well = sapply(colnames(cm), function(x) unlist(strsplit(x, split=".", fixed=T))[3])
  ## Total number of genes
  gene = colSums(cm >= gene_detect)
  cm_norm = log2(cm/scale_factor+1)
  ## Average expression of house keeping genes
  hg_list = hg_list[hg_list %in% rownames(cm)]
  hk = colMeans(cm_norm[hg_list,])
  qc_df = cbind.data.frame(sample, plate, well, gene, hk)
  return(qc_df)
}

## Calculate pairwise cell correlation and tSNE embedding
## @param cm_norm normalized count matrix
tSNEbyCor <- function(cm_norm, hc_method = "complete"){
  ## Compute pairwise pearson correlation and distance (1-r)
  message("Computing pairwise correlation and distance...")
  pairwise_cor = cor(cm_norm, method = "pearson")
  pairwise_dist = 1-pairwise_cor
  ## hirarchical clustering
  message("Computing hierachical clustering...")
  hc = hclust(as.dist(pairwise_dist), method = hc_method)
  ## tSNE
  message("Computing tSNE embedding...")
  tsne_cor = Rtsne(as.dist(pairwise_dist), pca = F, is_distance=T)
  result = list()
  result[["pairwise_cor"]] = pairwise_cor
  result[["pairwise_dist"]] = pairwise_dist
  result[["hc_obj"]] = hc
  result[["tsne_obj"]] = tsne_cor
  return(result)
}

## Function to plot and store tSNE plots
## @param tsne_obj tSNE object
## @param cm_norm log transformed cm
## @param color color scheme for tSNE
## @param title figure title
## @param filename file name of the figure
## @param path path to store the figure
## @param dot_size size of dot on tSNE
## @param axis_size font size of axis labels
## @param legend_size font size of legend
## @param save_obj whether to save the ggplot2 obj
## @param save_plot whether to save the plot
plotTsne <- function(tsne_obj, cm_norm, color, title="tSNE plot", filename="tSNE_cor.jpg", path="figures/",
                     dot_size=4, title_size=24, axis_size=18, legend_size=18, save_obj=FALSE, save_plot=TRUE){
  ## Check if length of tsne_obj and cm_norm matches
  try (if(dim(tsne_obj$Y)[1] != dim(cm_norm)[2]) stop("tsne_obj and cm_norm should have equal lengths..."))
  
  ## Generate a df for tSNE and color
  tsne_mat = data.frame(tsne_obj$Y)
  rownames(tsne_mat) = colnames(cm_norm)
  colnames(tsne_mat) = c("tSNE1", "tSNE2")
  tsne_mat$color = color
  ## Generate ggplot2 obj
  g = ggplot(tsne_mat, aes_string(x= "tSNE1", y = "tSNE2", colour="color")) +
    geom_point(size = dot_size) + theme_bw() + ggtitle(title) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = title_size),
          axis.title = element_text(size = axis_size), legend.text=element_text(size=legend_size))
  ## Display tSNE plot
  print(g)
  ## Save ggplot2 obj if save_obj = TRUE
  if (save_obj){
    return(g)
  }
  ## Save plot if save_plot = TRUE
  ggsave(filename, path = path, height = 8, width = 12)
}

## Identify highly variable genes by pagoda2
## @param cm count matrix
## @param gam.k parameter for pagoda's adjustVariance, default = 10
## @param plot_var plot variance graph or not, default = TRUE
## @param n.cores # of cores to calculate variance, default = 6
hvgPagoda <- function(cm, n.OdGenes=3000, gam.k=10, plot_var=TRUE, n.cores=6, verbose=TRUE){
  pagoda_obj <- Pagoda2$new(x = Matrix(as.matrix(cm), sparse=TRUE), n.cores=n.cores)
  ## Adjust the variance
  pagoda_obj$adjustVariance(plot=plot_var, gam.k=gam.k)
  ##pagoda_hvg = pagoda_obj$getOdGenes(n.odgenes=n.OdGenes)
  pagoda_hvg = pagoda_obj$misc$varinfo
  pagoda_hvg = pagoda_hvg[order(pagoda_hvg$lpa), ]
  if (verbose){
    cat(paste0("The number of highly variable genes identified by Pagoda2: ", length(pagoda_obj$misc$odgenes), '\n'))
  }
  return(pagoda_hvg)
}

## Identify highly variable genes by seurat
## @param cm count matrix
## @param scale.factor default = 100000
## @param x.low.cutoff lower bound for expressions of highly variable genes, default = 0.0125
## @param x.high.cutoff upper bound for expressions of highly variable genes, default = 8
## @param y.cutoff lower bound for variance of highly variable genes, default = 1
hvgSeurat <- function(cm, scale.factor=1E5, x.low.cutoff=0.0125, x.high.cutoff=8, y.cutoff=1, verbose=TRUE, variableFeatures=2000){
  seurat_obj = CreateSeuratObject(cm)
  ## Log normalize
  seurat_obj <- NormalizeData(
    object = seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = scale.factor
  )
  ## Detection of variable genes across the single cells
  seurat_obj <- FindVariableFeatures(
    object = seurat_obj,
    mean.function = ExpMean,
    dispersion.function = LogVMR,
    x.low.cutoff = 0.0125, ## Default is 0.0125
    x.high.cutoff = 8, ## Default is 8
    y.cutoff = 1, ## Default is 1
    nfeatures = variableFeatures
  )
  seurat_hvg = seurat_obj@assays$RNA@var.features
  cat(paste0("The number of highly variable genes identified by seurat: ", length(seurat_hvg)))
  return(seurat_hvg)
}

## Compute plain average expression or control corrected signature score
## @param X.center centered relative expression
## @param X.mean average of relative expression of each gene (log2 transformed)
## @param n number of genes with closest average expression for control genesets, default = 100
## @param simple whether use average, default  = FALSE
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

gene_symbol_to_ensembl_id <- function(gene, dataset){
  ensembl = useMart("ensembl")
  mart <- useDataset(dataset, useMart("ensembl"))
  ids <- getBM(filters="external_gene_name", attributes=c("ensembl_gene_id", "external_gene_name"), values=gene, mart=mart)
  non_duplicates <- which(duplicated(ids$external_gene_name) == FALSE)
  ids <- ids[non_duplicates, ]
}

go_analysis <- function(sigOE_genes, allOE_genes, keyType="ENSEMBL",
                        OrgDb="org.Hs.eg.db", ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05){
  ego <- enrichGO(gene = sigOE_genes,
                  universe = allOE_genes,
                  keyType = keyType,
                  OrgDb = OrgDb,
                  ont = ont,
                  pAdjustMethod = pAdjustMethod,
                  qvalueCutoff = qvalueCutoff,
                  readable = TRUE)
  
  cluster_summary <- subset(data.frame(ego), select = -c(geneID))
  cluster_summary <- data.frame(ego)
  return(list(ego=ego, cluster_summary=cluster_summary))
}

## Create seurat obj, normalize, find variable genes, and scale data
## @param cm, count matrix
## @param project project name
## @param min.cells minimal number of cells required to express a gene
## @param min.genes minimal number of genes required to be expressed in a cell
## @param is.expr detection limit for a gene
## @param scale.factor scale factor normalization, default=1E5
## @param do.scale whether to scale the data (divided by sd), default=F
## @param do.center whether to center the data (substract mean), default=T
preprocessSeuratObject <- function(cm, project, min.cells = 0, min.genes = 0,
                                   scale.factor = 1E5, do.scale = F, do.center = T, version = 'v4'){
  ## Create Seurat Obj, Normalize, variable genes, scaling data,
  if (version == 'v4') {
    seurat_obj <- CreateSeuratObject(cm, min.cells = min.cells, min.features = min.genes, project = project)
  } else if (version == 'v5') {
    seurat_obj <- CreateSeuratObject(Matrix::Matrix(as.matrix(cm),sparse = T), min.cells = min.cells, min.features = min.genes, project = project)
  }
  seurat_obj <- seurat_obj %>%
    NormalizeData(normalization.method = "LogNormalize", scale.factor = scale.factor, verbose = F) %>%
    FindVariableFeatures(mean.function = ExpMean, dispersion.function = LogVMR, verbose = F) %>%
    ScaleData(do.scale = do.scale, do.center = do.center, verbose = F)
  
  return(seurat_obj)
}


#preprocessSeuratObject_v5 <- function(cm, project, min.cells = 0, min.genes = 0, scale.factor = 1E5, do.scale = F, do.center = T, verbose){
#  ## Create Seurat Obj, Normalize, variable genes, scaling data,
#  message("Creating Seurat object, normalizing, variable genes and scaling data...")
#  suppressWarnings(seurat_obj <- CreateSeuratObject(Matrix::Matrix(as.matrix(cm),sparse = T), min.cells = min.cells, min.features = min.genes, project = project) %>%
#                     NormalizeData(normalization.method = "LogNormalize", scale.factor = scale.factor, verbose = verbose) %>%
#                     FindVariableFeatures(mean.function = ExpMean, dispersion.function = LogVMR, verbose = verbose) %>%
#                     ScaleData(do.scale = do.scale, do.center = do.center, verbose = verbose))
#
#  return(seurat_obj)
#}

#preprocessSeuratObject_v5 <- function(cm, project, min.cells = 0, min.genes = 0, scale.factor = 1E5, do.scale = F, do.center = T, norm.type = 'RNA', mt.pattern = '^MT-', verbose){
#  if (norm.type == 'RNA') {
#    ## Create Seurat Obj, Normalize, variable genes, scaling data,
#    message("Creating Seurat object, normalizing, variable genes and scaling data...")
#    suppressWarnings(seurat_obj <- CreateSeuratObject(Matrix::Matrix(as.matrix(cm),sparse = T), min.cells = min.cells, min.features = min.genes, project = project))
#    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt.pattern)
#    suppressWarnings(seurat_obj <- seurat_obj %>%
#                       NormalizeData(normalization.method = "LogNormalize", scale.factor = scale.factor, verbose = verbose) %>%
#                       FindVariableFeatures(mean.function = ExpMean, dispersion.function = LogVMR, verbose = verbose) %>%
#                       ScaleData(do.scale = do.scale, do.center = do.center, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"), verbose = verbose) %>% 
#                       CellCycleScoring(s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE, search = T) %>% 
#                       NormalizeData(normalization.method = "LogNormalize", scale.factor = scale.factor, verbose = verbose) %>%
#                       FindVariableFeatures(mean.function = ExpMean, dispersion.function = LogVMR, verbose = verbose) %>%
#                       ScaleData(do.scale = do.scale, do.center = do.center, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score"), verbose = verbose))
#  } else if (norm.type == 'SCT') {
#    message("Creating Seurat object...")
#    suppressWarnings(seurat_obj <- CreateSeuratObject(Matrix::Matrix(as.matrix(cm),sparse = T), min.cells = min.cells, min.features = min.genes, project = project) %>%
#                       NormalizeData(normalization.method = "LogNormalize", scale.factor = scale.factor, verbose = verbose))
#  }
#  return(seurat_obj)
#}

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

## Preprocess ss2 data for pagoda analysis
## @param cm count matrix
## @param total_gene_cutoff minimum total number of genes in a cell
## @param gene_expr gene detection limit for each cell
## @param gene_per_cell minimum number of cells with expressions of a gene
## @return filtered count matrix
preprocessing_ss2_data <- function(cm, total_gene_cutoff = 500, gene_expr = 5, gene_per_cell = 2){
  ## Filter out cells with total number of genes < 500
  cm <- cm[,which(apply(cm, 2, function(x) sum(x > gene_expr) > total_gene_cutoff))]
  ## Filter out genes that are expressed in 0 or 1 cell
  filter_genes <- apply(cm, 1, function(x) length(x[x > gene_expr]) >= gene_per_cell)
  cm <- cm[filter_genes,]
  ## Filter out duplicated gene names (usually do not exist after converting ensembl-ID to gene symbols)
  rownames(cm) <- make.unique(rownames(cm))
  ## Convert df to matrix and return
  cm <- as.matrix(cm)
  return(cm)
}

Head<-function(x){x[1:3,1:3]}

## Normalize and center count matrix, for input into ScoreSignature
NormCenter<-function(cm, scale_factor=10, log_base=2){
  cm = as.matrix(cm)
  cm_norm = log(cm/scale_factor+1, base=log_base)
  cm_norm_center = cm_norm-rowMeans(cm_norm)
  result = list()
  result[["raw_data"]] = cm
  result[["norm_data"]] = cm_norm
  result[["center_data"]] = cm_norm_center
  return (result)
}

## Convert v3 seurat object to v4 seurat object
## from raw data --> UMAP
ConvertV3ToV4<- function(seurat_v3, project, min.cells=0, min.genes=0,
                         scale.factor=1E5, do.scale=F, do.center=T,
                         harmonyTheta=2, resolution=0.8, dims=20){
  ## remove any loaded seurat package, reload v4 seurat
  detach("package:Seurat", unload=TRUE)
  library(Seurat, lib.loc = "C:/Users/jenna/OneDrive/Documents/R/win-library/4.0/Seurat/V4")
  
  ## extract counts and meta from v3 seurat object
  cm<- seurat_v3@assays$RNA@counts
  cm<-as.matrix(cm)
  meta<-seurat_v3@meta.data
  meta<- meta[,!(colnames(meta) %in% c("orig.ident", "nCount_RNA", "nFeature_RNA"))]
  
  ## Create new seurat object
  seurat_v4<- CreateSeuratObject(counts=cm, min.features = 0, min.cells = 0,
                                 meta.data = meta)
  
  ## Normalize data
  seurat_v4 <- NormalizeData(
    object = seurat_v4,
    normalization.method = "LogNormalize",
    scale.factor = scale.factor
  )
  
  ## Detection of variable genes across the single cells
  seurat_v4 <- FindVariableFeatures(
    object = seurat_v4,
    mean.function = ExpMean,
    dispersion.function = LogVMR
  )
  ##length(x = seurat_obj@var.genes)
  
  ## Scaling the data and removing unwanted sources of variation
  seurat_v4 <- ScaleData(
    object = seurat_v4,
    do.scale = do.scale,
    do.center = do.center
    #vars.to.regress = c("nUMI")
  )
  
  ## PCA
  seurat_v4 <- RunPCA(
    object = seurat_v4,
    features = seurat_v4@assays$RNA@var.features,
    npcs = 100,
    verbose = TRUE,
    ndims.print  = 1:5,
    nfeatures.print = 5
  )
  
  ## Harmony
  seurat_v4 = RunHarmony(seurat_v4, "sample", theta = harmonyTheta,
                         max.iter.harmony = 50, plot_convergence = TRUE)
  
  ## Run UMAP
  seurat_v4 <- RunUMAP(seurat_v4, reduction = "harmony", dims = 1:dims)
  
  
  seurat_v4 <- FindNeighbors(seurat_v4,
                             reduction = "harmony",
                             dims = 1:dims,
                             force.recalc = TRUE)   %>%
    FindClusters(resolution = resolution)
  
  detach("package:Seurat", unload=TRUE)
  return(seurat_v4)
}

## Full seurat v3 analysis
RunFullSeurat <- function(cm, do.scale = FALSE, do.center = TRUE, RunHarmony, samples,
                          resolution = .8, dims = 20, pca_dims = 100, n.neighbors = 30, project, version = 'v4'){
  seurat_obj <- preprocessSeuratObject(cm, project = project, min.cells = 0, min.genes = 0,
                                       scale.factor = 1E5, do.scale = do.scale, do.center = do.center, version = version)
  seurat_obj$sample <- samples
  
  seurat_obj <- RunPCA(
    object = seurat_obj,
    features = VariableFeatures(seurat_obj),
    npcs = pca_dims,
    verbose = F,
    ndims.print  = 1:5,
    nfeatures.print = 5)
  
  if (dims <= 1) {
    dims <- optimizePCA(seurat_obj, dims)
    #message("Using ", dims, " PCs as optimal number!")
  }
  
  if (RunHarmony){seurat_obj <- RunHarmony(seurat_obj, "sample", theta = 2, max.iter.harmony = 50, plot_convergence = TRUE)}
  reduction_use <- ifelse(RunHarmony, 'harmony', 'pca')
  
  seurat_obj <- seurat_obj %>%
    RunUMAP(reduction = reduction_use, dims = 1:dims, n.neighbors = n.neighbors, verbose = F) %>%
    RunTSNE(reduction = reduction_use, dims = 1:dims, num_threads = 8, check_duplicates = F, verbose = F) %>%
    FindNeighbors(reduction = reduction_use, dims = 1:dims, force.recalc = TRUE, verbose = F) %>%
    FindClusters(resolution = resolution, verbose = F)
  
  return(seurat_obj)
}



RunFullSeurat2 <- function(cm, metadata, do.scale = FALSE, do.center = TRUE, doBatch = NULL, var2batch = NULL, samples, resolution = .8, dims = 20, pca_dims = 100, n.neighbors = 30, project, version = 'v4'){
  seurat_obj <- preprocessSeuratObject(cm, project = project, min.cells = 0, min.genes = 0,
                                       scale.factor = 1E5, do.scale = do.scale, do.center = do.center, version = version)
  
  seurat_obj <- AddMetaData(seurat_obj, metadata)
  
  seurat_obj <- RunPCA(
    object = seurat_obj,
    features = VariableFeatures(seurat_obj),
    npcs = pca_dims,
    verbose = F,
    ndims.print  = 1:5,
    nfeatures.print = 5)
  
  if (dims <= 1) {
    dims <- optimizePCA(seurat_obj, dims)
    #message("Using ", dims, " PCs as optimal number!")
  }
  
  if (doBatch == 'harmony'){
    message("Performing batch correction using ", doBatch)
    seurat_obj <- RunHarmony(seurat_obj, group.by.vars = var2batch, theta = 2, max.iter.harmony = 50, plot_convergence = TRUE)
  } else if(doBatch == 'mnn'){
    message("Performing batch correction using ", doBatch)
    seurat_obj.list <- SplitObject(seurat_obj, split.by = var2batch)
    for (i in 1:length(seurat_obj.list)) {
      seurat_obj.list[[i]] <- NormalizeData(seurat_obj.list[[i]], verbose = F)
      seurat_obj.list[[i]] <- FindVariableFeatures(seurat_obj.list[[i]], verbose = F)
    }
    seurat_obj <- RunFastMNN(object.list = seurat_obj.list)
  }
  
  reduction_use <- ifelse(is.null(doBatch), 'pca', doBatch)
  
  seurat_obj <- seurat_obj %>%
    RunUMAP(reduction = reduction_use, dims = 1:dims, n.neighbors = n.neighbors, verbose = F) %>%
    RunTSNE(reduction = reduction_use, dims = 1:dims, num_threads = 8, check_duplicates = F, verbose = F) %>%
    FindNeighbors(reduction = reduction_use, dims = 1:dims, force.recalc = TRUE, verbose = F) %>%
    FindClusters(resolution = resolution, verbose = F)
  
  return(seurat_obj)
}


#RunFullSeurat_v5 <- function(cm, metadata, do.scale = FALSE, do.center = TRUE, doBatch = F, var2batch = NULL, batchMethod = 'harmony', resolution = 1, dims = 20, pca_dims = 100, n.neighbors = 30, project, verbose = FALSE){
#  seurat_obj <- preprocessSeuratObject_v5(cm, project = project, min.cells = 0, min.genes = 0, scale.factor = 1E5, do.scale = do.scale, do.center = do.center, verbose = verbose)
#
#  if (is.data.frame(metadata) | is.matrix(metadata)) {
#    seurat_obj <- AddMetaData(seurat_obj, metadata)
#  } else if (is.vector(metadata)) {
#    seurat_obj <- AddMetaData(seurat_obj, metadata, col.name = 'sample')
#  }
#
#  seurat_obj <- RunPCA(object = seurat_obj, features = VariableFeatures(seurat_obj), npcs = pca_dims, ndims.print  = 1:5, nfeatures.print = 5, verbose = verbose)
#
#  if (dims <= 1) {
#    dims <- optimizePCA(seurat_obj, dims)
#    #message("Using ", dims, " PCs as optimal number!")
#  }
#
#  if (doBatch){
#    seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj@meta.data %>% pull(var2batch))
#    seurat_obj <- suppressWarnings(NormalizeData(seurat_obj, verbose = verbose)) %>%
#      FindVariableFeatures(verbose = verbose) %>%
#      ScaleData(verbose = verbose) %>%
#      RunPCA(npcs = pca_dims, verbose = verbose) %>%
#      FindNeighbors(dims = 1:dims, reduction = "pca", verbose = verbose) %>%
#      FindClusters(resolution = resolution, cluster.name = "unintegrated_clusters", verbose = verbose) %>%
#      RunUMAP(dims = 1:dims, reduction = "pca", reduction.name = "umap.unintegrated", verbose = verbose, return.model = TRUE) %>%
#      RunTSNE(dims = 1:dims, reduction = "pca", reduction.name = "tsne.unintegrated", num_threads = 8, check_duplicates = F, verbose = verbose)
#
#    if (batchMethod != 'all') {
#      message(glue("Performing batch correction using {batchMethod}"))
#      met <- switch(batchMethod, cca = CCAIntegration, rpca = RPCAIntegration, harmony = HarmonyIntegration, mnn = FastMNNIntegration, NA)
#      reduc <- switch(batchMethod, cca = 'integrated.cca', rpca = 'integrated.rpca', harmony = 'harmony', mnn = 'integrated.mnn', NA)
#
#      seurat_obj <- IntegrateLayers(object = seurat_obj, method = met, orig.reduction = "pca", new.reduction = reduc, verbose = verbose)
#      seurat_obj <- seurat_obj %>%
#        FindNeighbors(reduction = reduc, dims = 1:dims, verbose = verbose) %>%
#        FindClusters(resolution = resolution, cluster.name = glue("{batchMethod}_clusters"), verbose = verbose) %>%
#        RunUMAP(reduction = reduc, dims = 1:dims, reduction.name = glue("umap.{batchMethod}"), verbose = verbose, return.model = TRUE) %>%
#        RunTSNE(reduction = reduc, dims = 1:dims, reduction.name = glue("tsne.{batchMethod}"), num_threads = 8, check_duplicates = F, verbose = verbose)
#    } else if (batchMethod == 'all') {
#      integration_methods <- c('cca', 'rpca', 'harmony', 'mnn')
#      for (intMet in integration_methods) {
#        message(glue("Performing batch correction using {intMet}"))
#        met <- switch(intMet, cca = CCAIntegration, rpca = RPCAIntegration, harmony = HarmonyIntegration, mnn = FastMNNIntegration, NA)
#        reduc <- switch(intMet, cca = 'integrated.cca', rpca = 'integrated.rpca', harmony = 'harmony', mnn = 'integrated.mnn', NA)
#
#        seurat_obj <- IntegrateLayers(object = seurat_obj, method = met, orig.reduction = "pca", new.reduction = reduc, verbose = verbose)
#        seurat_obj <- seurat_obj %>%
#          FindNeighbors(reduction = reduc, dims = 1:dims, verbose = verbose) %>%
#          FindClusters(resolution = resolution, cluster.name = glue("{intMet}_clusters"), verbose = verbose) %>%
#          RunUMAP(reduction = reduc, dims = 1:dims, reduction.name = glue("umap.{intMet}"), verbose = verbose, return.model = TRUE) %>%
#          RunTSNE(reduction = reduc, dims = 1:dims, reduction.name = glue("tsne.{intMet}"), num_threads = 8, check_duplicates = F, verbose = verbose)
#      }
#    }
#    seurat_obj <- JoinLayers(seurat_obj)
#  } else {
#    seurat_obj <- seurat_obj %>%
#      FindNeighbors(reduction = 'pca', dims = 1:dims, force.recalc = TRUE, verbose = verbose) %>%
#      FindClusters(resolution = resolution, verbose = verbose) %>%
#      RunUMAP(reduction = 'pca', dims = 1:dims, n.neighbors = n.neighbors, verbose = verbose, return.model = TRUE) %>%
#      RunTSNE(reduction = 'pca', dims = 1:dims, num_threads = 8, check_duplicates = F, verbose = verbose)
#  }
#  return(seurat_obj)
#}

#RunFullSeurat_v5 <- function(cm, metadata, do.scale = FALSE, do.center = TRUE, doBatch = F, var2batch = NULL, batchMethod = 'harmony', resolution = 1, dims = 20, pca_dims = 100, n.neighbors = 30, tsne_perplexity = 30, project, norm.type = 'RNA', mt.pattern = '^MT-', verbose = FALSE){
RunFullSeurat_v5 <- function(cm, metadata, do.scale = FALSE, do.center = TRUE, doBatch = F, var2batch = NULL, batchMethod = 'harmony', resolution = 1, dims = 20, pca_dims = 100, n.neighbors = 30, tsne_perplexity = 30, project, norm.type = 'RNA', verbose = FALSE){
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


## Function for running seurat pipeline
## For subsetting --> rerunning
## Same as Orr's pipeline shared
RunFullSeurat_Immune<- function(cm, samples){
  gcdata<-CreateSeuratObject(counts = cm, min.cells = 0, min.features = 0, project = "")
  gcdata$sampleid<- samples
  nSamples<-table(gcdata$sampleid)
  gcdata<- subset(gcdata, sampleid %in% names(nSamples[nSamples>1]))
  gcdata <- NormalizeData(gcdata, normalization.method = "LogNormalize", scale.factor = 10000)
  gcdata_list<- SplitObject(gcdata, split.by = "sampleid")
  var.genes <- SelectIntegrationFeatures(gcdata_list,
                                         nfeatures = 2000, verbose = TRUE, fvf.nfeatures = 2000,
                                         selection.method = "vst")
  VariableFeatures(gcdata) <- var.genes
  gcdata <- ScaleData(gcdata, split.by = "sampleid", features = VariableFeatures(gcdata),
                      do.center = T, do.scale = F)
  gcdata <- RunPCA(gcdata, features = VariableFeatures(gcdata), npcs = 40)
  gcdata <- FindNeighbors(gcdata, reduction = "pca", dims = 1:20, k.param = 20)
  gcdata <- FindClusters(gcdata, resolution = 1, algorithm = 4, random.seed = 100)  # conda activate /Users/jlabelle/Library/r-miniconda/envs/r-reticulate
  gcdata <- RunUMAP(gcdata, dims = 1:20, reduction = "pca", n.neighbors = 15, min.dist = 0.5, spread = 1, metric = "euclidean", seed.use = 1)
  return(gcdata)
}





## Plot GO results (from runGO)
## Input: list of go_results, 1 for each cell type. Output from runGO.
## Input: n_terms: number of GO terms to plot
plotGO<- function(go_result, n_terms=20){
  library(stringr)
  
  geneSets<-names(go_result)
  allDotPlots<-list()
  for (i in geneSets){
    print(i)
    if (nrow(go_result[[i]]$cluster_summary)>0){
      print("Sig terms present")
      allDotPlots[[i]]<-dotplot(go_result[[i]]$ego,
                                showCategory=n_terms,
                                font.size = 15,
                                title = i,
                                label_format=10) +
        scale_y_discrete(labels=function(x) str_wrap(x, width=60))+
        theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                        color="black", size = 28),
              axis.title = element_text(face = "bold", color="black"),
              axis.text.x = element_text(angle = 45, hjust = 1, color="black"),
              axis.text.y = element_text( color="black",face = "bold"))
    }
  }
  
  allDotPlots<- allDotPlots[!(unlist(lapply(allDotPlots, function(x){is.null(x)})))]
  return(allDotPlots)
}

## function for pseudobulking by variable
pseudobulk_byVariable<- function(cm, meta, variableToGroupBy){
  ## remove groups with only 1 cell
  tmp<- table(meta[,variableToGroupBy])
  group_names<- names(tmp[tmp!=1])
  
  pseudobulk = NULL
  if(any(colnames(cm) != rownames(meta))){print("cm / meta cell names do not match"); break}
  
  for (group in group_names){
    cm_tmp = cm[,meta[,variableToGroupBy] == group]
    pseudobulk = cbind(pseudobulk, rowMeans(cm_tmp))
  }
  pseudobulk<- as.data.frame(pseudobulk)
  colnames(pseudobulk) = group_names
  rownames(pseudobulk) = rownames(cm)
  return(pseudobulk)
}


## Function for heatmap plot, for set of GOI, Using ggplot
## Input count matrix should be pseudobulked counts, or grouped/averaged in some other way
## Column names should be 2 grouping variables, separated by "_". The first half of the column names will be the x axis (subtype usually)
## nameForCounts is just used for labeling legend appropriately
myHeatmap<- function(pseudobulked_cm, max.value=1e6, min.value=-1e6, nameForCounts="TPM", orderSubtypes="None", GOI, facetWrap=TRUE, orderFactors="None"){
  goi_cm<- as.data.frame(t(pseudobulked_cm[GOI,]))
  goi_cm$UniqueID<- rownames(goi_cm)
  goi_cm<- melt(goi_cm); colnames(goi_cm)<- c("UniqueID", "gene", "counts")
  goi_cm$gene<- factor(goi_cm$gene, levels=rev(unique(goi_cm$gene)))
  
  ## Split UniqueID into X axis (first half) and faceting variable (second half)
  tmp<- strsplit(goi_cm$UniqueID, split="_")
  goi_cm$XAxis<- unlist(lapply(tmp, function(x){x[1]}))
  goi_cm$FacetWrap<- unlist(lapply(tmp, function(x){x[2]}))
  
  ## Enforce min/max values
  goi_cm$counts<- ifelse(goi_cm$counts>max.value, max.value, ifelse(goi_cm$counts<min.value, min.value, goi_cm$counts))
  
  ## Order subtypes, if desired
  if(orderSubtypes[1] != "None"){
    goi_cm$XAxis<- factor(goi_cm$XAxis, levels=orderSubtypes)
  }
  
  ## Order factors, if desired
  if(orderFactors[1] != "None"){
    goi_cm$FacetWrap<- factor(goi_cm$FacetWrap, levels=orderFactors)
  }
  
  ## plot
  p<- ggplot(goi_cm, aes(x=XAxis, y= gene, fill=counts))+
    geom_tile()+
    theme_classic()+
    scale_fill_gradient2(high="red", low="yellow", mid="white", name=nameForCounts)+
    scale_x_discrete(position = "top")+
    theme(axis.text.x = element_text(angle=45, hjust=0, color="black", face="bold", size=16),
          axis.text.y = element_text(color="black", face="bold"))+
    xlab("")+ylab("")
  if(facetWrap){
    p<- p + facet_grid(cols=vars(FacetWrap), scale="free_x", space="free_x", switch="x")+
      theme(strip.text.x = element_text(size=16, color="black", face="bold"))
  }
  
  return(p)
}


## Function to find optimal value of PCs in a PCA to use to calculate UMAP
## Input seurat object with PCA already calculated
## csum is the minimum sum that PCs must sum
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


## Function to create enrichment matrix
## Input genes from seurat object [rownames(seurat_obj)]
BuildEnrichmentMatrix <- function(
    genes,
    type = 'GO',
    db = NULL,
    org = 'hsa'
){
  if (is.null(db)){
    db <- FindMSigDB(type, org)
  }
  terms <- names(db)
  enrichment.matrix <- pbapply::pbsapply(db, function(term){
    genes %in% term
  })
  rownames(enrichment.matrix) <- genes
  #enrichment.matrix = enrichment.matrix[, colSums(enrichment.matrix) > 0]
  enrichment.matrix[, colSums(enrichment.matrix) == 0] <- NA
  return(enrichment.matrix)
}


## Function to load all MSigDB information
## Input type that means the DB from MSigDB
FindMSigDB = function(
    type, org
){
  if (org == 'hsa') {
    if (type == 'GO'){
      db = MSigDB$C5[grep('^GO',names(MSigDB$C5))]
    } else if (type == 'HALLMARK'){
      db = MSigDB$H
    } else if (type == 'MOTIF'){
      db = MSigDB$C3[-grep('UNKNOWN|MIR',names(MSigDB$C3))]
    } else if (type == 'PATHWAYS'){
      db = MSigDB$C2[grep('BIOCARTA|REACTOME|KEGG',names(MSigDB$C2))]
    } else if (type == 'BIOCARTA'){
      db = MSigDB$C2[grep('BIOCARTA',names(MSigDB$C2))]
    } else if (type == 'KEGG'){
      db = MSigDB$C2[grep('KEGG',names(MSigDB$C2))]
    } else if (type == 'REACTOME'){
      db = MSigDB$C2[grep('REACTOME',names(MSigDB$C2))]
    } else if (type == 'syngo'){
      db = syngo
    } else {
      db = MSigDB$C5[grep(type, names(MSigDB$C5), ignore.case = TRUE, value = TRUE)]
    }
    return(db)
  } else if (org == 'mm') {
    if (type == 'GO'){
      db = MSigDB$M5[grep('^GO',names(MSigDB$M5))]
    } else if (type == 'HALLMARK'){
      db = MSigDB$MH
    } else if (type == 'MOTIF'){
      db = MSigDB$M3[-grep('UNKNOWN|MIR',names(MSigDB$M3))]
    } else if (type == 'PATHWAYS'){
      db = MSigDB$M2[grep('BIOCARTA|WP_|KEGG',names(MSigDB$M2))]
    } else if (type == 'BIOCARTA'){
      db = MSigDB$M2[grep('BIOCARTA',names(MSigDB$M2))]
    } else if (type == 'WIKI'){
      db = MSigDB$M2[grep('WP_',names(MSigDB$M2))]
    } else if (type == 'REACTOME'){
      db = MSigDB$M2[grep('REACTOME',names(MSigDB$M2))]
    } else {
      db = MSigDB$M5[grep(type, names(MSigDB$M5), ignore.case = TRUE, value = TRUE)]
    }
    return(db)
  }
}


## Function to perform hypergeometric test in a given dataset and DB
## Input geneset as a list of genes
Enrichment = function(
    geneset,
    genes = NULL,
    enrichment.matrix = NULL,
    srt = NULL,
    type = 'GO',
    do.sort = FALSE,
    pval.thresh = 1,
    pval.n = Inf,
    do.print = FALSE){
  
  if (is.null(genes) & is.null(enrichment.matrix)){
    genes = rownames(srt)
    enrichment.matrix = BuildEnrichmentMatrix(genes, type = type)
  }
  if (!is.null(genes) & is.null(enrichment.matrix)){
    enrichment.matrix = BuildEnrichmentMatrix(genes, type = type)
  }
  if (is.null(genes) & !is.null(enrichment.matrix)){
    genes = rownames(enrichment.matrix)
  }
  
  geneset.log = (genes %in% geneset)
  present = geneset.log %*% enrichment.matrix
  pvals = phyper(present, colSums(enrichment.matrix), colSums(1 - enrichment.matrix), sum(geneset.log), lower.tail = FALSE, log.p = FALSE)
  names(pvals) = colnames(enrichment.matrix)
  
  if (do.sort){
    pvals = sort(pvals)
  }
  pvals = pvals[pvals <= pval.thresh]
  if (length(pvals) >= pval.n){
    pvals = sort(pvals)
    pvals = pvals[1:pval.n]
  }
  if (do.print){
    print(pvals)
  }
  return(pvals)
}


## Remove part of a string after a pattern
## @param pattern replacement x
str_remove_after <- function(pattern, replacement, x) {
  gsub(sprintf("\\%s.*",pattern), sprintf("%s",replacement), x)
}


## Remove part of a string before a pattern
## @param pattern replacement x
str_remove_before <- function(pattern, replacement, x) {
  gsub(sprintf(".*\\%s",pattern), sprintf("%s",replacement), x)
}


## Project dataset in a query (dotplot)
projectDataNoScaling <- function(query_cm, query_degs, query_hvgs, query_query_metagene_order, 
                        ref_cm, ref_degs, ref_hvgs, 
                        outFile, cell_type_order, fig_width, fig_height) {
  # Filtering out genes with none expression
  ref_cm <- ref_cm[rowSums(ref_cm) > 1, names(ref_degs)]
  query_cm <- query_cm[rowSums(query_cm) > 0,]
  
  ## Signature genes
  input_genes <- union(ref_hvgs, query_hvgs)
  input_genes <- intersect(input_genes, intersect(rownames(ref_cm), rownames(query_cm)))
  
  # Normalize data 
  ## Tumor
  query_cm <- t(t(query_cm)/colSums(query_cm))*1E4
  query_cm_norm <- log2(query_cm/10+1)
  query_cm_mean <- log2(rowMeans(query_cm)+1)
  query_cm_center <- t(scale(t(query_cm_norm)))
  
  ## Filtering signature genes from reference
  ref_degs <- lapply(ref_degs, function(x) x[x %in% names(query_cm_mean)])
  ref_degs <- ref_degs[unlist(lapply(ref_degs, length)) > 0]
  ref_cm <- ref_cm[, names(ref_degs)]
  
  ## Ref
  ref_cm <- t(t(ref_cm)/colSums(ref_cm))*1E4
  ref_cm_norm <- log2(ref_cm+1)
  ref_cm_mean <- log2(rowMeans(ref_cm)+1)
  ref_cm_center <- t(scale(t(ref_cm_norm)))
  
  # Score metagene programs in ref scRNA-seq data
  message('Scoring query signatures in reference...')
  normal_score <- t(scoreNmfGenes(ref_cm_center, ref_cm_mean, query_degs, verbose = F))
  
  # Score malignant cells with ref DEGs
  message('Scoring reference signatures in query...')
  tumor_score <- scoreNmfGenes(query_cm_center, query_cm_mean, ref_degs, verbose = F)
  tumor_score <- tumor_score[, colnames(normal_score)]
  
  message('Saving plot...')
  
  # Pairwise correlation between ref and epn aggregated cm 
  pairwise_cor <- calculatePairwiseCor(ref_cm_norm, query_cm_norm, input_genes)
  pairwise_cor <- pairwise_cor[, colnames(normal_score)]
  
  ## max-min normalization and melt scoring normal cells by malignant metaprogram
  #tmp <- apply(normal_score, 2, function(x) (x-min(x))/(max(x)-min(x)))
  normal_score_long <- melt(normal_score)
  #normal_score_long <- melt(tmp)
  colnames(normal_score_long) <- c("Cell_type", "Metaprogram", "score")
  
  ## max-min normalization and melt scoring tumor cells by normal DEGs
  #tmp <- apply(tumor_score, 2, function(x) (x-min(x))/(max(x)-min(x)))
  #tumor_score_long <- melt(tmp)
  tumor_score_long <- melt(tumor_score)
  colnames(tumor_score_long) <- c("Cell_type", "Metaprogram", "score")
  
  ## Melt pairwise correlation 
  pairwise_cor_long <- melt(pairwise_cor)
  colnames(pairwise_cor_long) <- c("Cell_type", "Metaprogram", "r")
  
  ## Add score to pairwise-cor
  pairwise_cor_long$normal_score <- normal_score_long$score
  pairwise_cor_long$tumor_score <- tumor_score_long$score
  
  ## Replace "_" with "-"
  pairwise_cor_long$Cell_type <- gsub("_", "-", pairwise_cor_long$Cell_type)
  pairwise_cor_long$Metaprogram <- gsub("_", "-", pairwise_cor_long$Metaprogram)
  
  ## Sort metaprogram order
  pairwise_cor_long$Metaprogram <- factor(pairwise_cor_long$Metaprogram, levels = gsub("_", "-", query_metagene_order))
  pairwise_cor_long$Cell_type <- factor(pairwise_cor_long$Cell_type, levels = cell_type_order)
  
  ## modify all scores < 0 to o (otherwise problem with plotting)
  pairwise_cor_long$normal_score <- ifelse(pairwise_cor_long$normal_score < 0, 0, pairwise_cor_long$normal_score)
  pairwise_cor_long$tumor_score <- ifelse(pairwise_cor_long$tumor_score < 0, 0, pairwise_cor_long$tumor_score)
  
  ## Plot
  pt <- ggplot(pairwise_cor_long, aes(x=Cell_type, y=Metaprogram)) + 
    geom_point(aes(color=tumor_score, size=normal_score, alpha=r)) + scale_size_area(max_size=12) +
    scale_color_gradientn(colors=brewer.pal(n=9, name="Reds")) +
    labs(x = "Malignant metaprogram",
         y = "Normal cell type",
         color = "Expression score\n(tumor cells)", 
         size = "Expression score\n(normal cells)",
         alpha = "Correlation coefficient") +
    theme_classic() + 
    theme(axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 16, hjust = 1, angle = 45, color = 'black'),
          axis.text.y = element_text(size = 16, color = 'black'),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          plot.background = element_rect(fill = alpha("white", 0)),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.5, linetype = 'dashed', colour = "grey"))
  
  ggsave(filename = outFile, width = fig_width, height = fig_height)
  
  return(pt)
}
