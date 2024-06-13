library(Seurat)
library(tidyverse)
library(cluster)
#library(factoextra)
library(dendextend)
library(weights)
library(ggpubr)
library(matrixStats)
library(readxl)

base_dir <- "/n/scratch/users/s/sad167/EPN/scRNAseq"
resources_dir <- file.path(base_dir, 'scripts/resources')

source(file.path(resources_dir, "single_cell_preprocessing_helper_functions_CBJr.R"))
source(file.path(resources_dir, 'inferCNV_helper_functions.R'))
source(file.path(resources_dir, 'Plotting_helper_functions.R'))

metadata <- read_excel(file.path(base_dir, 'metadata.xlsx')) 

tp <- 'fresh'
metadata_subset <- metadata %>%
  filter(Type == tp)


qc_data <- file.path(base_dir, "analysis/qc/data")
infercnv_dir <- file.path(base_dir, 'analysis/infercnv/fresh')
infercnv_plots <- file.path(base_dir, 'analysis/infercnv/fresh/plots')
if (!dir.exists(infercnv_dir)){dir.create(infercnv_dir, recursive = T)}
if (!dir.exists(infercnv_plots)){dir.create(infercnv_plots, recursive = T)}


cm <- readRDS(file.path(infercnv_dir, "cm_exp_ctrl.rds"))
samples <- readRDS(file.path(infercnv_dir, "samples_exp_ctrl.rds"))
obs_ref <- readRDS(file.path(infercnv_dir, "obs_ref.rds"))
cnv_stat <- readRDS(file.path(infercnv_dir, "call_cnv_w_ctrl.rds"))

meta <- data.frame(sample = samples, type = obs_ref, cnv = cnv_stat)
suppressWarnings(seurat_obj <- RunFullSeurat_v5(cm = cm, metadata = meta,  doBatch = T, var2batch = 'sample', batchMethod = 'harmony', project = 'EPN-fresh'))
saveRDS(seurat_obj, file = file.path(infercnv_dir, "seurat_obj.rds"))

## By cluster
p1 <- DimPlot(object = seurat_obj, group.by = 'seurat_clusters', 
              reduction = 'umap.harmony', label = TRUE, label.box = T, pt.size = 0.1, label.size = 3) +
  ggtitle('Clusters') +
  seurat_theme() +
  NoLegend()

## By samples
p2 <- DimPlot(object = seurat_obj, reduction = 'umap.harmony', group.by = "sample", pt.size = 0.1) +
  ggtitle('Samples') +
  seurat_theme() +
  NoLegend()

## By obs vs ref
p3 <- DimPlot(object = seurat_obj, reduction = 'umap.harmony', group.by = "type", label = TRUE, pt.size = 0.1, label.size = 0, cols = colors_to_use) +
  ggtitle('obs/ref') +
  seurat_theme()

## By CNV
p4 <- DimPlot(object = seurat_obj, reduction = 'umap.harmony', group.by = "cnv", cells = colnames(seurat_obj)[seurat_obj$type == "obs"], cols = colors_to_use, label = TRUE, pt.size = 0.1, label.size = 0) +
  ggtitle('cnv') +
  seurat_theme()

(p1 + p2) / (p3 + p4)
ggsave("umap.pdf", path = infercnv_plots, width = 10, height = 8)

p <- FeaturePlot(seurat_obj,
                 reduction = 'umap.harmony', 
                 features = c("CD14", "CSF1R", ## Myeloids
                              "CD3E", "CD4", ## T-cells
                              "MBP", "PLP1", "TF", "APOD", ## Oligodendrocytes
                              "VWF", ## Endothelial
                              "RGS5"), ## Pericyte
                 cols = c("lightgrey", "red"), ncol = 5, combine = F)
p <- lapply(p, function(x) x + seurat_theme())
ggarrange(plotlist = p, ncol = 5, nrow = 2)
ggsave("normal_marker_genes.pdf", path = infercnv_plots, width = 25, height = 8)


expr_CD14 <- FetchData(seurat_obj, vars = 'CD14')
myeloid_cells_CD14 <- rownames(expr_CD14)[which(expr_CD14$CD14 > quantile(expr_CD14$CD14, 0.95))]
expr_CSF1R <- FetchData(seurat_obj, vars = 'CSF1R')
myeloid_cells_CSF1R <- rownames(expr_CSF1R)[which(expr_CSF1R$CSF1R > quantile(expr_CSF1R$CSF1R, 0.95))]
myeloid_cells <- c(myeloid_cells_CD14, myeloid_cells_CSF1R) %>% unique()
DimPlot(object = seurat_obj, pt.size = 0.1, cells.highlight = myeloid_cells, sizes.highlight = 0.1) +
  seurat_theme()

expr_MBP <- FetchData(seurat_obj, vars = 'MBP')
oligo_cells <- rownames(expr_MBP)[which(expr_MBP$MBP > quantile(expr_MBP$MBP, 0.95))]
DimPlot(object = seurat_obj, pt.size = 0.1, cells.highlight = oligo_cells, sizes.highlight = 0.1) +
  seurat_theme()

expr_VWF <- FetchData(seurat_obj, vars = 'VWF')
endothelial_cells <- rownames(expr_VWF)[which(expr_VWF$VWF > quantile(expr_VWF$VWF, 0.96))]
DimPlot(object = seurat_obj, pt.size = 0.1, cells.highlight = endothelial_cells, sizes.highlight = 0.1) +
  seurat_theme()

expr_RGS5 <- FetchData(seurat_obj, vars = 'RGS5')
pericyte_cells <- rownames(expr_RGS5)[which(expr_RGS5$RGS5 > quantile(expr_RGS5$RGS5, 0.96))]
DimPlot(object = seurat_obj, pt.size = 0.1, cells.highlight = pericyte_cells, sizes.highlight = 0.1) +
  seurat_theme()

normal_cells <- Cells(seurat_obj)[seurat_obj$seurat_clusters  %in% c(2, 4, 3)]
saveRDS(normal_cells, file.path(infercnv_dir, "normal_cells_names.rds"))

seurat_obj$ct <- ifelse(Cells(seurat_obj) %in% normal_cells, 'Normal', 'Malignant')


gene_list <- list(Pericyte = 'RGS5',
                  Oligodendrocyte = c('MBP', 'PLP1', 'TF', 'APOD'),
                  Myeloid = c('CSF1R', 'CD14'),
                  Endothelial = 'VWF')
SCpubr::do_DotPlot(sample = seurat_obj,
                   features = gene_list,
                   group.by = 'ct',
                   flip = FALSE) +
  scale_fill_gradient2(low="blue", mid="white", high="red")
ggsave("dotplot_normal_marker_genes.pdf", path = infercnv_plots, width = 12, height = 4)


## Malignant vs non-malignant cell calling
#############Strategy###############
## All cells in normal clusters and without CNV = Normal
## All cells in samples without cnv are Malignabnt (not used)
## All cells NOT in normal clusters and with CNV = Malignant
## The rest = Low-quality
####################################
# call_tumor <- function(x){
#   ## cells in normal cluster are normal cells
#   if (Cells(seurat_obj)[x] %in% normal_cells & !seurat_obj$cnv[x]){
#     return("Normal")
#   } else if (!(Cells(seurat_obj)[x] %in% normal_cells) & seurat_obj$cnv[x]){
#     return("Malignant")
#   }else{
#     return("Low-quality")
#   }
# }


call_tumor <- function(x){
  ## cells in normal cluster are normal cells
  if (Cells(seurat_obj)[x] %in% normal_cells){
    return("Non-malignant")
  } else{
    return("Malignant")
  }
}
tumor_stat_w_ctrl <- structure(pbapply::pbsapply(seq(length(Cells(seurat_obj))), call_tumor), names = Cells(seurat_obj))

## Remove normal controls
tumor_stat <- tumor_stat_w_ctrl[(seurat_obj$sample != "mg" &
                                   seurat_obj$sample != "od")]
saveRDS(tumor_stat, file = file.path(infercnv_dir, "call_tumor.rds"))

## tSNE colored by malignancy
seurat_obj$tumor <- tumor_stat_w_ctrl
DimPlot(object = seurat_obj, group.by = "tumor", pt.size = 1, cols = c(gg_color_hue(2))) +
  ggtitle('tumor') +
  seurat_theme()
ggsave("tumor.pdf", path = infercnv_plots, width = 6, height = 4)


## Saving results
to_subset <- seurat_obj$type == "obs"
cm_postqc <- GetAssayData(seurat_obj, slot="counts")
cm_postqc <- cm_postqc[, to_subset]
saveRDS(cm_postqc, file.path(infercnv_dir, "cm_postqc.rds"))

to_subset <- seurat_obj$type == "obs" & seurat_obj$tumor == "Malignant"
cm_malignant <- GetAssayData(seurat_obj, slot="counts")
cm_malignant <- cm_malignant[, to_subset]
saveRDS(cm_malignant, file.path(infercnv_dir, "cm_malignant.rds"))


## Preprocessing malignant cells only
seurat_obj <- subset(seurat_obj, cells = colnames(cm_malignant))
metadata <- seurat_obj@meta.data
suppressWarnings(seurat_obj <- RunFullSeurat_v5(cm = cm, metadata = meta,  doBatch = T, var2batch = 'sample', batchMethod = 'harmony', project = 'EPN-fresh'))
saveRDS(seurat_obj, file.path(infercnv_dir, "seurat_obj_malignant.rds"))
