rm(list = ls())

# Load packages -----------------------------------
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
library(writexl)

# Set up environment  -------------------------------------
base_dir <- '/Users/sdaniell/Partners HealthCare Dropbox/Sara Danielli/Sara Danielli/Project/Ependymoma/12 - Xenium'
analysis_dir <- file.path(base_dir, 'analysis/1_preparation')
output_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
plot_dir <- file.path(analysis_dir, 'plot/STEPN')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

data_dir <- "/Users/sdaniell/Partners HealthCare Dropbox/Sara Danielli/Sara Danielli/Project/Ependymoma/11 - Plots scRNAseq/data"

resource_dir <- file.path(base_dir, 'scripts_revisions/resources')
source(file.path(resource_dir, 'single_cell_preprocessing_helper_functions_CBJr.R'))

# Load seurat object with both malignant and normal cells ---------------------
seurat_object <- qread(file = file.path(data_dir, "patient/seurat_obj_ST_normal_malig_annotated.qs"))

# subset to normal cells only
table(seurat_object$malignant)
Idents(seurat_object) <- 'malignant'
seurat_object_normal <- subset(seurat_object, idents = 'Non-malignant')
sort(unique(seurat_object_normal$Sample_deID))

# load seurat object with malignant cells ------------------------------
seurat_obj_mal <- qread(file = file.path(data_dir, "patient/seurat_obj_malignant_annotated2.qs"))


# Create combined annotated object --------------------

# Add malignant cell state to object with all cells 
table(seurat_obj_mal$Metaprogram_noCC)
table(seurat_obj_mal$Metaprogram)
seurat_obj_mal$cell_type <- seurat_obj_mal$Metaprogram_noCC

# concatenate
seurat_object <- merge(seurat_object_normal, seurat_obj_mal)
table(seurat_object$cell_type)
seurat_object <- JoinLayers(seurat_object)

# remove unclassified
Idents(seurat_object) <- 'cell_type'
seurat_object <- subset(seurat_object, idents = 'Unclassified', invert = T)

# reprocess normal+malignant object with RUnFullseurat
metadata <- seurat_object@meta.data

cm <- seurat_object[["RNA"]]$counts
cm_norm <- as.matrix(log2(cm/10+1))
cm_mean <- log2(Matrix::rowMeans(cm)+1)
cm_center <- cm_norm - rowMeans(cm_norm)
seurat_object <- RunFullSeurat_v5(cm = cm, metadata = metadata,  doBatch = F,  project = 'EPN')

# Plot --------------------
colors_tumor_calling <- c('grey90', 'lightyellow3', '#3EBCB6FF',  '#BDA14DFF', '#09090CFF', '#D5114EFF')
names(colors_tumor_calling) <- c('Unclassified', 'Malignant', 'Myeloid', 'T-cells', 'Oligodendrocytes', 'Endothelial')
colors_metaprograms<- c("gray30","#F99E93FF","#9E5E9BFF","#74ADD1FF",'#0F4F8B', "#ACD39EFF","#96410EFF", 'mistyrose1')

names(colors_metaprograms) <- c("Cycling", "Neuroepithelial-like", "Radial-glia-like", 
                                "Embryonic-neuronal-like", "Neuronal-like" ,"Ependymal-like", "MES-like", "Embryonic-like")

colors_pop <- c(colors_tumor_calling, colors_metaprograms)

# UMAP plot
p1 <- DimPlot(seurat_object, reduction = 'umap', group.by = 'cell_type', cols = colors_pop, label.size = 4, label = T)
p2 <- DimPlot(seurat_object, reduction = 'umap', group.by = 'Sample_deID', label.size = 4, label = T)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.3))
ggsave(file.path(plot_dir, '1_UMAP.pdf'), width = 14, height = 4)



# Load marker genes ----------------------------------- 

# Read Xenium base panel  
Xenium_panel <- read_csv(file.path(resource_dir, 'base_panel/Xenium_hBrain_v1_metadata (3).csv')) 
Xenium_panel <- as.data.frame(Xenium_panel)
# convert into list
Xenium_panel <- split(Xenium_panel[, -c(2:6)], f = Xenium_panel$Annotation_Sara)


# Read NMF Marker genes  -----------------------------------   
markers_NMF_200 <- readRDS(file.path(data_dir, 'signatures/nmf_marker_genes_final_annotated_200.rds'))
# remove cycling program
markers_NMF_200 <- markers_NMF_200[!names(markers_NMF_200) =='Cycling']

# select top 
markers_NMF_top <- lapply(markers_NMF_200, head, 100)

#remove RP-, AS-, AL- genes
unique_NMF <- unlist(markers_NMF_top)
unique_NMF <- unique_NMF[!grepl("^(AS|AL|RP|AC0|orf)", unique_NMF)]


names(unique_NMF) <- sapply(unique_NMF, function(element) {
  # Find the name(s) of vectors containing the element
  matches <- names(markers_NMF_top)[sapply(markers_NMF_top, function(vec) element %in% vec)]
  
  # Return the matching name(s) or "none" if no match
  if (length(matches) > 0) {
    paste(matches, collapse = ", ")  # Combine multiple matches
  } else {
    "none"
  }
})

unique_NMF_list <- split(unique_NMF, f = names(unique_NMF))

# only keep genes that are unique to one MP
unique_NMF_list <- unique_NMF_list[names(markers_NMF_top)]
unique_NMF_list <- lapply(unique_NMF_list, unname)
names(unique_NMF_list) <- gsub("Embryonic/neuronal-like", "Embryonic-neuronal-like",  names(unique_NMF_list))


# Export list as df
# Determine the maximum length among all elements in the list
max_length <- max(sapply(unique_NMF_list, length))

# Pad each element in the list with NA to make them of equal length
padded_list <- lapply(unique_NMF_list, function(x) {c(x, rep(NA, max_length - length(x)))})

# Convert padded list to dataframe
my_dataframe <- as.data.frame(do.call(cbind, padded_list))

write_xlsx(my_dataframe, file.path(plot_dir, "0_NMF_marker_genes.xlsx"))




# Identify markers by FindAllMarkers  -----------------------------------   
Idents(seurat_object) <- 'cell_type'
markers <- RunPrestoAll(seurat_object,  only.pos = TRUE, min.cells.feature = 0.25, logfc.threshold = 0.5) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T)

#remove RP-, AS-, AL- genes, 
markers <- markers[!grepl("^(AS|AL|RP|AC0|orf)", markers$gene), ]

markers_list_SS2_all <- split(markers$gene, markers$cluster)
#names(markers_list_SS2_all) <- paste("ss2-", names(markers_list_SS2_all))
# select top 
markers_list_SS2_top <- lapply(markers_list_SS2_all, head, 100)


# Export list as df
# Determine the maximum length among all elements in the list
max_length <- max(sapply(markers_list_SS2_top, length))

# Pad each element in the list with NA to make them of equal length
padded_list <- lapply(markers_list_SS2_top, function(x) {c(x, rep(NA, max_length - length(x)))})

# Convert padded list to dataframe
my_dataframe <- as.data.frame(do.call(cbind, padded_list))

write_xlsx(my_dataframe, file.path(plot_dir, "0_FindAllmarker_genes.xlsx"))




# Check for congruency ----------------------------------------
markers_list_SS2_top_mal <- markers_list_SS2_top[5:11]

common_elements <- lapply(names(markers_list_SS2_top_mal), function(name) {
  if (name %in% names(unique_NMF_list)) {
    intersect(markers_list_SS2_top_mal[[name]], unique_NMF_list[[name]])
  } else {
    NULL  # No match for this name
  }
})
common_elements
names(common_elements) <- names(markers_list_SS2_top_mal)


# Export list as df
# Determine the maximum length among all elements in the list
max_length <- max(sapply(common_elements, length))

# Pad each element in the list with NA to make them of equal length
padded_list <- lapply(common_elements, function(x) {c(x, rep(NA, max_length - length(x)))})

# Convert padded list to dataframe
my_dataframe <- as.data.frame(do.call(cbind, padded_list))

write_xlsx(my_dataframe, file.path(plot_dir, "0_CommonElements_genes.xlsx"))




# Check expresion ----------------------------------------

# base panel
Idents(seurat_object) <- 'cell_type'
DotPlot(seurat_object,
        features = Xenium_panel,
        dot.scale = 6,
        col.min = 0, 
        #cols = c("grey90", "red3"),
        scale = T, assay = 'RNA'
) +  RotatedAxis()   +
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) 
ggsave(file.path(plot_dir, "2_Xenium_base_panel_scaled.pdf"), width=35, height=4)

DotPlot(seurat_object,
        features = Xenium_panel,
        dot.scale = 6,
        col.min = 0, 
        #cols = c("grey90", "red3"),
        scale = F, assay = 'RNA'
) +  RotatedAxis()   +
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) 
ggsave(file.path(plot_dir, "2_Xenium_base_panel_unscaled.pdf"), width=35, height=4)



# NMF genes panel
Idents(seurat_object) <- 'cell_type'
unique_NMF_list <- unique_NMF_list[rev(unique(seurat_object$cell_type)[5:11])]

# select top 
unique_NMF_list <- lapply(unique_NMF_list, head, 17)

DotPlot(seurat_object,
        features = unique_NMF_list,
        dot.scale = 6,
        col.min = 0, 
        #cols = c("grey90", "red3"),
        scale = T, assay = 'RNA'
) +  RotatedAxis()   +
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) 
ggsave(file.path(plot_dir, "3_NMF_markers_scaled.pdf"), width=26, height=4)

DotPlot(seurat_object,
        features = unique_NMF_list,
        dot.scale = 6,
        col.min = 0, 
        #cols = c("grey90", "red3"),
        scale = F, assay = 'RNA'
) +  RotatedAxis()   +
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) 
ggsave(file.path(plot_dir, "3_NMF_markers_unscaled.pdf"), width=26, height=4)




# FindAllMarker genes

# select top 
markers_list_SS2_top <- lapply(markers_list_SS2_top, head, 17)


DotPlot(seurat_object,
        features = markers_list_SS2_top[5:11],
        dot.scale = 6,
        col.min = 0, 
        #cols = c("grey90", "red3"),
        scale = T, assay = 'RNA'
) +  RotatedAxis()   +
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) 
ggsave(file.path(plot_dir, "4_FindAllmarkers_scaled.pdf"), width=40, height=4)

DotPlot(seurat_object,
        features = markers_list_SS2_top[5:11],
        dot.scale = 6,
        col.min = 0, 
        #cols = c("grey90", "red3"),
        scale = F, assay = 'RNA'
) +  RotatedAxis()   +
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) 
ggsave(file.path(plot_dir, "4_FindAllmarkers_unscaled.pdf"), width=26, height=4)
length(unlist(unique_NMF_list))



# common genes 
DotPlot(seurat_object,
        features = common_elements,
        dot.scale = 6,
        col.min = 0, 
        #cols = c("grey90", "red3"),
        scale = T, assay = 'RNA'
) +  RotatedAxis()   +
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) 
ggsave(file.path(plot_dir, "5_common_genes_scaled.pdf"), width=20, height=4)

DotPlot(seurat_object,
        features = common_elements,
        dot.scale = 6,
        col.min = 0, 
        #cols = c("grey90", "red3"),
        scale = F, assay = 'RNA'
) +  RotatedAxis()   +
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) 
ggsave(file.path(plot_dir, "5_common_genes_unscaled.pdf"), width=20, height=4)





# refine final marker gene list -----------------------------------
# take all common genes and add NMF genes for embryonic-neuronal like program
final_panel <- c(common_elements, unique_NMF_list["Embryonic-neuronal-like"])
final_panel[3] <- NULL

final_panel <- final_panel[rev(unique(seurat_object$cell_type)[5:11])]

# final marker genes
DotPlot(seurat_object,
        features = final_panel,
        dot.scale = 6,
        col.min = 0, 
        #cols = c("grey90", "red3"),
        scale = T, assay = 'RNA'
) +  RotatedAxis()   +
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) 
ggsave(file.path(plot_dir, "6_final_panel_scaled.pdf"), width=20, height=4)

DotPlot(seurat_object,
        features = final_panel,
        dot.scale = 6,
        col.min = 0, 
        #cols = c("grey90", "red3"),
        scale = F, assay = 'RNA'
) +  RotatedAxis()   +
  paletteer::scale_colour_paletteer_c("viridis::mako", direction = -1) 
ggsave(file.path(plot_dir, "6_final_panel_unscaled.pdf"), width=20, height=4)





# Export list as df
# Determine the maximum length among all elements in the list
max_length <- max(sapply(final_panel, length))

# Pad each element in the list with NA to make them of equal length
padded_list <- lapply(final_panel, function(x) {c(x, rep(NA, max_length - length(x)))})

# Convert padded list to dataframe
my_dataframe <- as.data.frame(do.call(cbind, padded_list))

write_xlsx(my_dataframe, file.path(plot_dir, "0_final_panel_genes.xlsx"))
