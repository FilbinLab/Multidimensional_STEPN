# Load libraries
library(Seurat)
library(ggplot2)

## Function for UMAP plot with n = 3 columns and many markers without axes
plotMarker <- function(seu_obj, markers){
  p <- FeaturePlot(seu_obj, features = markers, reduction = "umap", combine = FALSE, sort.cell = TRUE)
  for(i in 1:length(p)) {
    p[[i]] <- p[[i]] + NoAxes() + scale_colour_distiller(palette="RdBu")
  }
  cowplot::plot_grid(plotlist = p, ncol=3)
}
