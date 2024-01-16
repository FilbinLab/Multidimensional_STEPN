library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(RColorBrewer)

# Organize environment  -----------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

# color palette
colors_to_use <- structure(c("#F99E93FF","#9E5E9BFF","#74ADD1FF","#ACD39EFF","#96410EFF", '#BDA14DFF', '#3EBCB6FF', '#0169C4FF', '#153460FF' ,'#A56EB6FF'),
                           names = c("Neuroepithelial-like", "Radial glia-like", "NPC-like" ,"Ependymal-like", "Mesenchymal", "T-cell", "Myeloid", "Microglia", "Endothelial", "VLMCs"))


################################################################################
## Neighbors enrichment analysis for malignant programs only
################################################################################
nhood_enrichment_zscore <- read.csv(file.path(base_dir, 'analysis/neighborhood/all_programs/2_neighborhood.csv'), row.names = 1)
nhood_enrichment_zscore[is.na(nhood_enrichment_zscore)] <- 0

ann <- data.frame(rownames(nhood_enrichment_zscore))
colnames(ann) <- c('Metaprogram')
colours <- list('Metaprogram' = colors_to_use)
colAnn <- HeatmapAnnotation(df = ann, which = 'col', col = colours, show_annotation_name = F, show_legend = F)
rowAnn <- HeatmapAnnotation(df = ann, which = 'row', col = colours, show_annotation_name = F, show_legend = F)

breaks <- seq(-50, 80, by= 0.1)
png(file.path(base_dir, 'analysis/neighborhood/all_programs/3_Neighbors_enrichment_average.png'), units = 'in', width = 12, height = 8, res = 600)
    Heatmap(as.matrix(nhood_enrichment_zscore), col = colorRamp2(breaks, colorRampPalette(c("steelblue3", 'white', 'red3'))(length(breaks))), cluster_rows = T, cluster_columns = T, show_column_names = F,
            top_annotation = colAnn, left_annotation = rowAnn, row_names_gp = gpar(fontsize = 10), row_names_side = "left", name = '-',
            heatmap_legend_param = list(legend_direction = "vertical", legend_height = unit(13, "cm")), show_row_dend = F, show_column_dend = F)
    dev.off()
    
    png(file.path(base_dir, 'analysis/neighborhood/all_programs/3_Neighbors_enrichment_average2.png'), units = 'in', width = 12, height = 8, res = 600)
        Heatmap(as.matrix(nhood_enrichment_zscore), col = colorRamp2(breaks, colorRampPalette(c("steelblue3", 'white', 'red3'))(length(breaks))), cluster_rows = T, cluster_columns = T, show_column_names = F,
                top_annotation = colAnn, left_annotation = rowAnn, row_names_gp = gpar(fontsize = 10), row_names_side = "left", name = '-',
                heatmap_legend_param = list(legend_direction = "vertical", legend_height = unit(13, "cm")), show_row_dend = F, show_column_dend = F)
        dev.off()