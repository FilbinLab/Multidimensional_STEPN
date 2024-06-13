library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(RColorBrewer)

# Organize environment  -----------------------------------
base_dir <- file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/12 - Xenium')

resource_dir  <- file.path(base_dir, 'scripts/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))




# Read average neighborhood enrichment z-scores calculated in python  -----------------------------------
nhood_enrichment_zscore <- read.csv(file.path(base_dir, 'analysis/all_Xenium_runs/neighborhood/2_neighborhood.csv'), row.names = 1)
nhood_enrichment_zscore[is.na(nhood_enrichment_zscore)] <- 0

# remove microglia and radial glia-like, as those programs were present only in 1 sample
#rows_to_exclude <- c('Microglia', 'Radial glia-like')
#cols_to_exclude <- c('Microglia', 'Radial.glia.like')

# nhood_enrichment_zscore <- nhood_enrichment_zscore[!(rownames(nhood_enrichment_zscore) %in% rows_to_exclude), 
#                                                    !(colnames(nhood_enrichment_zscore) %in% cols_to_exclude)]

ann <- data.frame(rownames(nhood_enrichment_zscore))
colnames(ann) <- c('Metaprogram')
colours <- list('Metaprogram' = colors_metaprograms_Xenium_grouped)
colAnn <- HeatmapAnnotation(df = ann, which = 'col', col = colours, show_annotation_name = F, show_legend = F)
rowAnn <- HeatmapAnnotation(df = ann, which = 'row', col = colours, show_annotation_name = F, show_legend = F)

# export matrix
write_csv(as.data.frame(nhood_enrichment_zscore), file.path(file.path(base_dir, 'analysis/all_Xenium_runs/neighborhood/6_Neighbors_enrichment_average.csv')))

## Heatmap of correlations
col_fun = colorRamp2(c(-100, 0, 300), c("#0571B0", "#FCFDFEFF", "#CA0020"))

pdf(file.path(base_dir, 'analysis/all_Xenium_runs/neighborhood/6_Neighbors_enrichment_average.pdf'), 
    width = 6, height = 4)
Heatmap(as.matrix(nhood_enrichment_zscore), 
        col = col_fun,
        cluster_rows = T, 
        cluster_columns = T, 
        show_column_names = F,
        top_annotation = colAnn, left_annotation = rowAnn, 
        row_names_gp = gpar(fontsize = 10), row_names_side = "left", name = '-',
        heatmap_legend_param = list(legend_direction = "vertical", legend_height = unit(13, "cm")), show_row_dend = F, show_column_dend = F)
dev.off()
