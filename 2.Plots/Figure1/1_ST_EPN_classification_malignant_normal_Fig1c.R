# Load packages -----------------------------------
rm(list = ls())

library(data.table)
library(tidyverse)
library(crayon)
library(ape)
library(readxl)
library(qs)
library(cowplot)
library(Seurat)
library(ggpubr)
library(future)
library(SingleCellExperiment)
library(SingleR)
library(celldex)

# Organize environment and create folders -----------------------------------
base_dir <- "/Users/sdaniell/Partners HealthCare Dropbox/Sara Danielli/Filbin lab - sample location/Filbin Lab Active/Sara/Project/Ependymoma/11 - Plots scRNAseq"
resource_dir <- file.path(base_dir, "scripts/resources")
source(file.path(resource_dir, "Plot_style_v2.R"))
source(file.path(resource_dir, "NMF_helper_function.R")) 
source(file.path(resource_dir, "single_cell_preprocessing_helper_functions_CBJr.R"))

data_dir <- file.path(base_dir, "data")
analysis_dir <- file.path(base_dir, 'analysis')
plot_dir <- file.path(analysis_dir, "plots/STEPN/all_cells")
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}


# Define color palettes  -----------------------------------
colors_tumor_calling <- c('grey90', 'lightyellow3', '#3EBCB6FF',  '#BDA14DFF', '#09090CFF', '#D5114EFF')
names(colors_tumor_calling) <- c('Unclassified', 'Malignant', 'Myeloid', 'T-cells', 'Oligodendrocytes', 'Endothelial')


# load in data -----------------------------------
# seurat object
seurat_obj <- qread(file = file.path(base_dir, "data/patient/seurat_obj_ST_normal_malig_fresh_frozen_annotated.qs"))

# Metadata samples
metadata_pat <- read_xlsx(file.path(base_dir, 'data/metadata/metadata_ST_EPN_Patient.xlsx'))


# (1) Check for inferCNV output -----------------------------------

## UMAP By classification
p1 <- DimPlot(object = seurat_obj, reduction = 'umap.harmony', group.by = "cnv", label = FALSE, pt.size = 0.5, shuffle = TRUE, label.size = 0, cols = c( 'red', 'black')) +
  ggtitle(paste0('inferCNV harmony')) 
p1
ggsave(file.path(plot_dir, '1_inferCNV_harmony.pdf'), width = 6, height = 5) 

DimPlot(object = seurat_obj, reduction = 'umap.unintegrated', group.by = "cnv", label = FALSE, pt.size = 0.5, shuffle = TRUE, label.size = 0, cols = c( 'red', 'black')) +
  ggtitle(paste0('inferCNV unintegrated')) 
ggsave(file.path(plot_dir, '2_inferCNV_unintegrated.pdf'), width = 6, height = 5) 

## UMAP By patient subtype
p2 <- DimPlot(object = seurat_obj, reduction = 'umap.harmony', group.by = "Subtype", pt.size = 0.8, shuffle = TRUE,  label = TRUE,label.size = 0,cols = col_subtype) + 
  labs(title = 'ST-EPN', subtitle = 'n = 7,840 cells') + theme_tSNE
p2
ggsave(file.path(plot_dir, '3_subtype_harmony.pdf'), width = 6.5, height = 5) 

DimPlot(object = seurat_obj, reduction = 'umap.unintegrated', group.by = "Subtype", pt.size = 0.8, shuffle = TRUE,  label = TRUE,label.size = 0,cols = col_subtype) + 
  labs(title = 'ST-EPN', subtitle = 'n = 7,840 cells') + theme_tSNE
ggsave(file.path(plot_dir, '4_subtype_unintegrated.pdf'), width = 6.5, height = 5) 

## UMAP By sample
DimPlot(object = seurat_obj, reduction = 'umap.harmony', group.by = "Sample_deID", label = FALSE, pt.size = 0.5, shuffle = TRUE, label.size = 0, cols = col_Sample_deID) +
  ggtitle(paste0('ST-EPN sample harmony')) 
ggsave(file.path(plot_dir, '5_sampleID_harmony.pdf'), width = 8.5, height = 5) 

DimPlot(object = seurat_obj, reduction = 'umap.unintegrated', group.by = "Sample_deID", label = T, pt.size = 0.5, shuffle = TRUE, label.size = 0, cols = col_Sample_deID) +
  ggtitle(paste0('ST-EPN sample unintegrated')) 
ggsave(file.path(plot_dir, '6_sampleID_unintegrated.pdf'), width = 8.5, height = 5) 


## UMAP By classification used for NMF (standard pipeline was used to infer malignant vs non-malignatn cells, therefore results are slightly different)
DimPlot(object = seurat_obj, reduction = 'umap.harmony', group.by = "malignant", label = FALSE, pt.size = 0.5, shuffle = TRUE, label.size = 0) +
  ggtitle(paste0('NMF classification sample harmony')) 
ggsave(file.path(plot_dir, '7_classificationNMF_harmony.pdf'), width = 6.5, height = 5) 

DimPlot(object = seurat_obj, reduction = 'umap.unintegrated', group.by = "malignant", label = T, pt.size = 0.5, shuffle = TRUE, label.size = 0) +
  ggtitle(paste0('NMF classification unintegrated')) 
ggsave(file.path(plot_dir, '8_classificationNMF_unintegrated.pdf'), width = 6.5, height = 5) 



# (2) Check for expression of normal markers -----------------------------------

#Define normal markers 
normal_markers <- c("CD14", "CSF1R", "CD163", ## Myeloids
                    "CD3E", "CD4", ## T-cells
                    "MBP", "PLP1", "TF", ## Oligodendrocytes
                    "VWF", "PECAM1", "CLDN5", ## Endothelial
                    "PDGFRB","ACTA2") ## Mural

## UMAP By clusters
DimPlot(object = seurat_obj, reduction = 'umap.harmony', group.by = "harmony_clusters", label = TRUE, pt.size = 0.5, shuffle = TRUE, 
              label.size = 4) + NoLegend() +
  ggtitle('Seurat clusters') 
ggsave(file.path(plot_dir, '9_UMAP_seurat_clusters_harmony.pdf'), width = 5, height = 5) 

DimPlot(object = seurat_obj, reduction = 'umap.unintegrated', group.by = "unintegrated_clusters", label = TRUE, pt.size = 0.5, shuffle = TRUE, 
        label.size = 4) + NoLegend() +
  ggtitle('Seurat clusters') 
ggsave(file.path(plot_dir, '9_UMAP_seurat_clusters_unintegrated.pdf'), width = 5, height = 5) 

## UMAP By markers
plot <- FeaturePlot(seurat_obj,
            reduction = 'umap.harmony', 
            features = normal_markers,
            cols = c("lightgrey", "red"), ncol = 5, pt.size = 3, combine = F, raster = TRUE)
ggarrange(plotlist = plot, ncol = 6, nrow = 2)
ggsave(file.path(plot_dir, '10_UMAP_markers_harmony.pdf'), width = 20, height = 6) 

plot <- FeaturePlot(seurat_obj,
                  reduction = 'umap.unintegrated', 
                  features = normal_markers,
                  cols = c("lightgrey", "red"), ncol = 5, pt.size = 3, combine = F, raster = TRUE)
ggarrange(plotlist = plot, ncol = 6, nrow = 2)
ggsave(file.path(plot_dir, '11_UMAP_markers_unintegrated.pdf'), width = 20, height = 6) 


## Dotplot by markers
Idents(seurat_obj) <- 'harmony_clusters'
DotPlot(seurat_obj,
        features = normal_markers,
        cols = c("lightgrey", "red")) 
ggsave(file.path(plot_dir, '12_DotPlot_markers_harmony.pdf'), width = 10, height = 4) 



# (3) Use Label transfer to define cell annotations -----------------------------------
seurat_obj$harmony_clusters <- factor(seurat_obj$harmony_clusters, levels = unique(sort(as.numeric( seurat_obj$harmony_clusters ))))

# transform into single cell object 
seurat_obj.sce <- as.SingleCellExperiment(seurat_obj)

# load reference dataset - Human Primary Cell Atlas (Mabbott et al. 2013)
ref <- HumanPrimaryCellAtlasData()

# use SingleR to transfer labels
pred.seurat_obj <- SingleR(test=seurat_obj.sce, 
                           ref = ref, 
                           clusters = seurat_obj$harmony_clusters, 
                           labels = ref$label.main, 
                           assay.type.ref = 1,
                           prune = TRUE)

# inspect results and export (use pruned.labels, which automatically excludes low confidence annotations)
table <- table(pred.seurat_obj$pruned.labels)
write.csv(table, file.path(plot_dir, '13_SingleR_prediction.pdf'))

plotScoreHeatmap(pred.seurat_obj)
ggsave(file.path(plot_dir, '14_SingleR_heatmap.pdf'), width = 10, height = 10) 

plotDeltaDistribution(pred.seurat_obj, ncol = 3)
ggsave(file.path(plot_dir, '15_SingleR_deltas.pdf'), width = 10, height = 10) 


# add predicted labels to Seurat object
Idents(seurat_obj) <- seurat_obj$harmony_clusters 
new.cluster.ids <- pred.seurat_obj$pruned.labels
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj[["SingleR.labels"]] <- Idents(object = seurat_obj)

DimPlot(object = seurat_obj, reduction = 'umap.harmony', group.by = "SingleR.labels", label = TRUE, 
        pt.size = 0.5, raster = FALSE, shuffle = TRUE, label.size = 0) + 
  labs(title = 'SingleR prediction', subtitle = 'n = 7,840 cells') 
ggsave(file.path(plot_dir, '16_SingleR_UMAP_harmony.pdf'), width = 7, height = 5) 

DimPlot(object = seurat_obj, reduction = 'umap.unintegrated', group.by = "SingleR.labels", label = TRUE, 
        pt.size = 0.5, raster = FALSE, shuffle = TRUE, label.size = 0) + 
  labs(title = 'SingleR prediction', subtitle = 'n = 7,840 cells') 
ggsave(file.path(plot_dir, '17_SingleR_UMAP.pdf'), width = 7, height = 5) 



# Manually define clusters as malignant or normal -----------------------------------

# classify  all malignant cells used for NMF as malignant
Idents(seurat_obj) <- 'malignant'
seurat_obj_malignant <- subset(seurat_obj, idents = 'Malignant')
seurat_obj_non_malignant <- subset(seurat_obj, idents = 'Non-malignant')

# extract metadata
seurat_obj_malignant_metadata <- seurat_obj_malignant@meta.data
seurat_obj_non_malignant_metadata <- seurat_obj_non_malignant@meta.data

# classify all malignant cells as malignant
seurat_obj_malignant_metadata$cell_type <- 'Malignant'

# for non-malignant cells, classify based on output of SingleR (unclassified = malignant cells, that were not considered for NMF)
Idents(seurat_obj_non_malignant) <- "harmony_clusters"
new.cluster.ids <- c("Unclassified", # 0
                     "Unclassified", # 1
                     "Unclassified", # 2
                     "Unclassified", # 3
                     "Unclassified", # 4
                     "Unclassified", # 5
                     "Unclassified", # 6
                     "Unclassified", # 7
                     "Unclassified", # 8
                     "Unclassified", # 9
                     "Myeloid", # 10
                     "Unclassified", #11
                     "Unclassified", #12
                     "T-cells", # 13
                     "Endothelial", # 14
                     "Endothelial", # 15
                     "Oligodendrocytes" # 16
)
names(new.cluster.ids) <- levels(seurat_obj_non_malignant)
seurat_obj_non_malignant <- RenameIdents(seurat_obj_non_malignant, new.cluster.ids)
seurat_obj_non_malignant[["cell_type"]] <- Idents(object = seurat_obj_non_malignant)

# extract metadata with cell annotation
seurat_obj_non_malignant_metadata <- seurat_obj_non_malignant@meta.data

# add metadata 
metadata <- rbind(seurat_obj_malignant_metadata, seurat_obj_non_malignant_metadata)
seurat_obj <- AddMetaData(seurat_obj, metadata)

# Plot UMAP with manual classification -----------------------------------
DimPlot(object = seurat_obj, reduction = 'umap.harmony', group.by = "malignant", label = TRUE, cols = c(  'red', 'black'),
        pt.size = 0.5, raster = FALSE, shuffle = TRUE, label.size = 0) +  
  labs(title = 'Manual classification', subtitle = 'n = 7,840 cells') 
ggsave(file.path(plot_dir, '16_UMAP_Manual_classification_harmony.pdf'), width = 6, height = 5) 

DimPlot(object = seurat_obj, reduction = 'umap.unintegrated', group.by = "malignant", label = TRUE, cols =  c(  'red', 'black'),
        pt.size = 0.5, raster = FALSE, shuffle = TRUE, label.size = 0) + 
  labs(title = 'Manual classification', subtitle = 'n = 7,840 cells') 
ggsave(file.path(plot_dir, '17_UMAP_Manual_classification.pdf'), width = 6, height = 5) 

p2 <- DimPlot(object = seurat_obj, reduction = 'umap.harmony', group.by = "cell_type", label = TRUE, cols = colors_tumor_calling,
        pt.size = 0.5, raster = FALSE, shuffle = TRUE, label.size = 4) + 
  labs(title = 'Final cell type classification', subtitle = 'n = 7,840 cells') 
p2
ggsave(file.path(plot_dir, '18_UMAP_Celltype_classification_harmony.pdf'), width = 6.5, height = 5) 


DimPlot(object = seurat_obj, reduction = 'umap.unintegrated', group.by = "cell_type", label = TRUE, cols = colors_tumor_calling,
              pt.size = 0.5, raster = FALSE, shuffle = TRUE, label.size = 4) + 
  labs(title = 'Final cell type classification', subtitle = 'n = 7,840 cells')
ggsave(file.path(plot_dir, '18_UMAP_Celltype_classification.pdf'), width = 6.5, height = 5) 


plot_grid(p1, p2, ncol = 2)
ggsave(file.path(plot_dir, '0_UMAP_ExtendedFigure1.pdf'), width = 13, height = 5) 


DimPlot(object = seurat_obj, reduction = 'umap.unintegrated', group.by = "cell_type", label = TRUE, cols = colors_tumor_calling,
        pt.size = 0.5, raster = FALSE, shuffle = TRUE, label.size = 0) +
  labs(title = 'Final cell type classification', subtitle = 'n = 7,429 cells') 
ggsave(file.path(plot_dir, '19_UMAP_Celltype_classification.pdf'), width = 6, height = 5) 

DimPlot(object = seurat_obj, reduction = 'umap.harmony', group.by = "malignant", label = TRUE, 
        pt.size = 0.5, raster = FALSE, shuffle = TRUE, label.size = 0) + ggtitle('Cell type classification') 
ggsave(file.path(plot_dir, '19_UMAP_Celltype_classification.pdf'), width = 6, height = 5) 

## Dotplot normal markers
Idents(seurat_obj) <- 'cell_type'
levels(seurat_obj) <- c('Malignant', 'Myeloid', 'T-cells', 'Oligodendrocytes', 'Endothelial', 'Unclassified')

DotPlot(seurat_obj,
        features = c("CD14", "SPP1", 'CCL3', ## Myeloids
                     'PYHIN1', 'CD247', 'CD3D', ## T-cells
                     "MOG", "MAG", "MBP", ## Oligodendrocytes
                     "VWF",  "CLDN5" ## Endothelial
                     ),
        dot.scale = 6,
        col.min = 0, 
        #cols = c("grey90", "red3"),
        scale = T, assay = 'RNA'
) +  RotatedAxis()   +
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) 
ggsave(file.path(plot_dir, "20_Dotplot_markers_scaled.pdf"), width=7, height=3)



# Clean up datasets before saving  -----------------------------------

# add metadata with patient info
colnames(metadata_pat) <- gsub('Sample', 'sample', colnames(metadata_pat))
seurat_obj@meta.data <- left_join(seurat_obj@meta.data, metadata_pat)
seurat_obj@meta.data[22:28] <- NULL

# Save datasets -----------------------------------

# bug in Seurat v5 (need to run this, otherwise error subsetting)
rownames(seurat_obj@meta.data) <- Cells(seurat_obj)

# save malignant+normal tumor dataset
qsave(seurat_obj, file.path(base_dir, "data/patient/seurat_obj_ST_normal_malig_annotated.qs"))

# save normal cell dataset
Idents(seurat_obj) = 'malignant'
seurat_obj_normal <- subset(seurat_obj, idents = "Non-malignant")
qsave(seurat_obj_normal, file.path(base_dir, "data/patient/seurat_obj_ST_normal_annotated.qs"))   

