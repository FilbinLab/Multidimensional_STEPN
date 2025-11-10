rm(list = ls())

# Load packages -----------------------------------
library(dplyr)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(tidyverse)
library(ggpubr)
library(readxl)
library(qs)
library(SeuratWrappers)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
library(grid)
library(gridExtra)

# Organize environment  -----------------------------------
base_dir <- "/Users/sdaniell/Partners HealthCare Dropbox/Sara Danielli/Filbin lab - sample location/Filbin Lab Active/Sara/Project/Ependymoma/11 - Plots scRNAseq"
data_dir <- file.path(base_dir, "data")
analysis_dir <- file.path(base_dir, 'analysis')
plot_dir <- file.path(analysis_dir, "plots/STEPN/malignant/pseudobulk")
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

resource_dir <- file.path(base_dir, 'scripts/resources')
source(file.path(resource_dir, "Plot_style_v2.R"))
source(file.path(resource_dir, "NMF_helper_function.R"))
source(file.path(resource_dir, "single_cell_preprocessing_helper_functions_CBJr.R"))
source(file.path(resource_dir, "NMF_helper_function.R"))

heatmap_plot <- pc_plot <- list()
hvg_genes <- c(500, 1000, 2000, 3000)

for (gene_nr in hvg_genes) {
  
    # load objects -----------------------------------
    seurat_obj <- qread(file = file.path(base_dir, "data/patient/seurat_obj_ST_malig_annotated.qs"))
    
    # load metadata 
    metadata <- read_xlsx(file.path(base_dir, 'data/metadata/metadata_ST_EPN_Patient.xlsx'))
    metadata <- subset(metadata, metadata$`sc/snRNAseq` %in% c('snRNA-seq', 'scRNA-seq'))
    
  

    # Part 1 - Bulk correlation ----------------------------------------------------
    new_metadata = seurat_obj@meta.data
    labels = as.character(unique(new_metadata$sample))
    list_means = list()
    for(label in labels) {
      cells = rownames(new_metadata[new_metadata$sample == label, ])
      avg_cells = rowMeans(as.matrix(LayerData(seurat_obj, 'data')[, cells]))
      list_means[[label]] <- avg_cells
    }
    avg_bulk = do.call(cbind, list_means)
    
    # try with variable genes
    var_genes <- VariableFeatures(seurat_obj, nfeatures = gene_nr)
    cormM = cor((avg_bulk[var_genes, ]), method="spearman")
    
    
    # annotation
    anno <- data.frame(sample = gsub('^[^-]*-', '', rownames(cormM)),
                       row.names = rownames(cormM))
    matching_indices <- match(anno$sample , metadata$FileName)
    anno$Subtype <- metadata$Subtype[matching_indices]
    
    anno$sample <- NULL
    
    # change rownames to have sampledeID
    matching_indices <- match(rownames(cormM)  , metadata$FileName)
    rownames(cormM) <- metadata$Sample_deID[matching_indices]
    matching_indices <- match(colnames(cormM)  , metadata$FileName)
    
    colnames(cormM) <- metadata$Sample_deID[matching_indices]
    
    
    
    # modify code to add source, location etc
    colAnno <- HeatmapAnnotation(df = anno, which = 'column', 
                                  col = list(Subtype = col_subtype),
                                 simple_anno_size = unit(3, "mm"),
                                 show_annotation_name = T,
                                 annotation_legend_param = list(title_gp = gpar(fontsize = 8),
                                                                title = "State",
                                                                labels_gp = gpar(fontsize = 8)))
    
    
    library(dendextend)
    dend = as.dendrogram(hclust(dist(cormM)), type = "average") # can be used as cluster_rows and columns arg instead of T

    col_fun = colorRamp2(c(0, 0.7, 1), c("#0571B0", "#FCFDFEFF", "#CA0020"))

    
    heatmap_plot[[as.character(gene_nr)]] = Heatmap(cormM, show_column_names = T, show_row_dend = T, 
                 show_column_dend = F, show_row_names = T, 
                 name = "Corr", row_names_gp = gpar(fontsize = 7),
                 column_names_gp = gpar(fontsize = 7),
                 col = col_fun,
                 cluster_rows = dend, row_dend_reorder = F, cluster_columns = dend,
                 column_dend_reorder = F,
                 top_annotation = colAnno,
                 width = ncol(cormM)*unit(2, "mm"), 
                 height = nrow(cormM)*unit(2, "mm"),
                 heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                             title = "Correlation",
                                             labels_gp = gpar(fontsize = 7),
                                             legend_height = unit(2, "cm")))
    
    pdf(file.path(plot_dir, paste0('1_Heatmap_', gene_nr, '.pdf')),
        width = 9, height = 7)
    print(heatmap_plot[[as.character(gene_nr)]])
    dev.off()
    

    
    
    
    
    
    ### Part 2 - PC plot-------------------------------------
    
    ### 1) create pseudobulk matrix
    
    
    # transform into matrix
    cm <- as.matrix(seurat_obj[["RNA"]]$counts)
    
    # extract metadata
    metadata_seurat <- seurat_obj@meta.data
    
    # transform as single cell object
    sce <- SingleCellExperiment(list(counts = cm),
                                colData = metadata_seurat
    )
    
    # first subset the coldata to only have the columns we care about for pseudo-bulking 
    groups <- colData(sce)[, c("sample", "Subtype")]
    
    # create a new SCE object that contains 
    # the pseudo-bulked counts across the provided groups 
    sce <- scuttle::aggregateAcrossCells(sce, 
                                         id = groups)
    
    # column names aren't automatically added to the pseudo-bulked sce, 
    # so let's add them in 
    colnames(sce) <- sce$sample
    
    
    
    
    ## 2) create deseq2
    deseq_object <- DESeq2::DESeqDataSet(sce,
                                         design = ~ Subtype)
    
    # estimate size factors first
    deseq_object <- DESeq2::estimateSizeFactors(deseq_object)
    
    # normalize and log transform to use for visualization
    normalized_object <- DESeq2::vst(deseq_object, 
                                     blind = TRUE)
    
    
    ## 3) PC plot
    pcaData <- plotPCA(normalized_object, intgroup=c("sample", "Subtype"), 
                       ntop = gene_nr,
                       returnData=TRUE)
    
    
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    pc_plot[[as.character(gene_nr)]] <- ggplot(pcaData, aes(PC1, PC2, color=Subtype)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      scale_color_manual(values = col_subtype)+
      theme(plot.title = element_text(hjust=0.5, face="bold"),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),  # Add square border
            axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1, colour="black"),
            axis.text.y = element_text(size=12, colour="black"),
            axis.title=element_text(size=12),
            strip.text = element_text(size = 13, face = "bold")) +
      guides(color = guide_legend(ncol = 1)) + ggtitle(paste0('HVGs = ', gene_nr))
    
    pc_plot[[as.character(gene_nr)]] 
    ggsave(file.path(plot_dir, paste0("2_PCA_", gene_nr, "_hvgs.pdf")), width=5, height=3.5)
    
    
}



# Convert list to grobs and arrange
heatmap_grobs <- lapply(heatmap_plot, function(h) grid.grabExpr(draw(h)))
  
# plot
plot_grid(plotlist = c(pc_plot, heatmap_grobs), nrow =2,  axis = "tb"  )
ggsave(file.path(plot_dir, paste0("Pseudobulk_plots_all.pdf")), width=25, height=9)
 

