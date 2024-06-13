# Load packages -----------------------------------
library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(data.table)
library(readxl)
library(dplyr)
library(ggrastr)
library(ComplexHeatmap)
library(cowplot)
library(glue)
library(RColorBrewer)
library(xlsx)
library(ggrepel)


# Organize environment  ----------------------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

resource_dir  <- file.path(base_dir, 'scripts/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))
source(file.path(resource_dir, 'NMF_helper_function.R'))
source(file.path(resource_dir, 'Plotting_helper_functions.R'))

analysis_dir  <- file.path(base_dir, 'analysis/all_Xenium_runs')

plot_dir <- file.path(analysis_dir, 'plots/niche')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

data_dir <- file.path(analysis_dir, 'data/niche')
if (!dir.exists(data_dir)){dir.create(data_dir, recursive = T)}

seurat_obj_dir <- file.path(analysis_dir, 'data/seurat_objects')
if (!dir.exists(seurat_obj_dir)){dir.create(seurat_obj_dir, recursive = T)}

metadata <- read_xlsx(file.path(base_dir, 'SampleIdentifier.xlsx'))

# reorder by sampleID
metadata <- metadata %>% arrange(SampleID)



for (i in seq_along(metadata$Sample)) { 
  
  # read data
  data <- qread(file.path(seurat_obj_dir, paste0('seurat_obj_', metadata$SampleID[i], "_", metadata$SampleName[i], '.qs')))

  # plot niches  ------------------------------------------------
  celltype.plot <- ImageDimPlot(data, group.by = 'group', fov = "fov",  cols = colors_metaprograms_Xenium, border.size = NA, size = 0.5,
                               dark.background = F) + ggtitle("Cell type")
  niche.plot <- ImageDimPlot(data, group.by = "niches", fov = "fov",  cols = col_niches, border.size = NA, size = 0.5,
                           dark.background = F) + ggtitle("Niches")
  celltype.plot | niche.plot
  ggsave(file.path(plot_dir, paste0('1_ImageDimPlot_niche_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=16, height=5)
  
  
  # Save data table with frequency
  table <- data@meta.data %>% group_by(group, niches) %>%
    summarise(n = sum(n()))
  
  # calculate total number of cells per niche   ------------------------------------------------
  table_sum <- table %>% group_by(niches) %>%
    summarise(n = sum(n))
  matching_indices <- match( table$niches, table_sum$niches)
  table$total_cells_per_niche <- table_sum$n[matching_indices]
  
  # calculate frequency of cell type in niche
  table <- table %>% group_by(group, niches) %>%
    summarise(freq = (n/total_cells_per_niche)*100)
  
  # transform into matrix
  cell_frequency <- table %>%
    pivot_wider(names_from = niches, values_from = freq, values_fill = 0)
  cell_frequency2 <- cell_frequency[, !names(cell_frequency) %in% "group"]
  
  cell_frequency_df <- as.matrix(cell_frequency2)
  rownames(cell_frequency_df) <- cell_frequency$group
  write.csv(as.data.frame(cell_frequency_df), row.names = T, file.path(data_dir, paste0('niche_frequency_cell_type_', metadata$SampleID[i],"_", metadata$SampleName[i], '.csv')))
  
  # Heatmap niche   ------------------------------------------------
  
  ann_col <- data.frame(rownames(t(cell_frequency_df)))
  colnames(ann_col) <- c('Niche')
  colours <- list('Niche' = colors_to_use_niches)
  colAnn <- HeatmapAnnotation(df = ann_col, which = 'col', col = colours, show_annotation_name = F, show_legend = F)

  ann_row <- data.frame(colnames(t(cell_frequency_df)))
  colnames(ann_row) <- c('Metaprogram')
  colours <- list('Metaprogram' = colors_to_use_metaprograms)
  rowAnn <- HeatmapAnnotation(df = ann_row, which = 'row', col = colours, show_annotation_name = F, show_legend = F)


  col_fun=circlize::colorRamp2(c(0, 25, 50), c("steelblue3", "white", "red3"))

  pdf(file.path(plot_dir, paste0('2_Heatmap_niche_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')),
      width = 5, height = 4)
  draw(Heatmap(t(t(cell_frequency_df)),
               show_row_dend = T,
               show_column_dend = F,
               cluster_rows = T,
               cluster_columns = F,
               rect_gp = grid::gpar(type = "none"),
               cell_fun = cell_fun,
               col = col_fun,
               column_names_rot = 45,
               show_column_names = T,
               top_annotation = colAnn,
               left_annotation = rowAnn,
               row_names_gp = gpar(fontsize = 10), row_names_side = "left", name = '-',
               heatmap_legend_param = list(legend_direction = "vertical", legend_height = unit(5, "cm")),
               width = ncol(t(t(cell_frequency_df)))*unit(10, "mm"),
               height = nrow(t(t(cell_frequency_df)))*unit(10, "mm")))
  dev.off()


  # Plot niche frequency
  plot_bar(data, data$niches, data$group, colors_metaprograms_Xenium)
  ggsave(file.path(plot_dir, paste0('3_BarPlot_niche_', metadata$SampleID[i], "_", metadata$SampleName[i], '.pdf')), width=6, height=5)


}





# Identify recurrent patterns of niche distribution  ------------------------------------------------

# read frequency of niche files
files <- list.files(data_dir, pattern = 'niche_frequency_cell_type', full.names = T)
samples <- gsub('niche_frequency_cell_type_|.csv', '', basename(files))

tab <- lapply(files, function(x) {
  read.csv(x, row.names = 1) %>%
    rename_all(~glue('Niche_{1:4}')) %>%
    reshape2::melt()
})


results <- list()
for (i in seq_along(files)) {
  results[[i]] <- read.csv(files[i], row.names = 1) %>%
    rename_all(~glue('Niche_{1:4}')) %>%
    rownames_to_column('Program')
  
  if ('VLMCs' %in% results[[i]]$Program & 'Endothelial' %in% results[[i]]$Program) {
    tmp <- data.frame(Program = 'Vessel',
                      Niche_1 = sum(results[[i]][which(results[[i]]$Program == 'VLMCs'),'Niche_1']) + sum(results[[i]][which(results[[i]]$Program == 'Endothelial'),'Niche_1']),
                      Niche_2 = sum(results[[i]][which(results[[i]]$Program == 'VLMCs'),'Niche_2']) + sum(results[[i]][which(results[[i]]$Program == 'Endothelial'),'Niche_2']),
                      Niche_3 = sum(results[[i]][which(results[[i]]$Program == 'VLMCs'),'Niche_3']) + sum(results[[i]][which(results[[i]]$Program == 'Endothelial'),'Niche_3']),
                      Niche_4 = sum(results[[i]][which(results[[i]]$Program == 'VLMCs'),'Niche_4']) + sum(results[[i]][which(results[[i]]$Program == 'Endothelial'),'Niche_4']))
    results[[i]] <- results[[i]][-which(results[[i]]$Program == 'VLMCs' | results[[i]]$Program == 'Endothelial'),]
    results[[i]] <- rbind(results[[i]], tmp)
  } else if ('VLMCs' %in% results[[i]]$Program) {
    results[[i]][which(results[[i]]$Program == 'VLMCs'), 'Program'] <- 'Vessel'
  } else if ('Endothelial' %in% results[[i]]$Program) {
    results[[i]][which(results[[i]]$Program == 'Endothelial'), 'Program'] <- 'Vessel'
  }
  
  if (('Myeloid' %in% results[[i]]$Program) & ('Microglia' %in% results[[i]]$Program) & ('T-cell' %in% results[[i]]$Program)) {
    tmp <- data.frame(Program = 'Immune',
                      Niche_1 = sum(results[[i]][which(results[[i]]$Program == 'Myeloid'),'Niche_1']) + sum(results[[i]][which(results[[i]]$Program == 'Microglia'),'Niche_1']) + sum(results[[i]][which(results[[i]]$Program == 'T-cell'),'Niche_1']),
                      Niche_2 = sum(results[[i]][which(results[[i]]$Program == 'Myeloid'),'Niche_2']) + sum(results[[i]][which(results[[i]]$Program == 'Microglia'),'Niche_2']) + sum(results[[i]][which(results[[i]]$Program == 'T-cell'),'Niche_2']),
                      Niche_3 = sum(results[[i]][which(results[[i]]$Program == 'Myeloid'),'Niche_3']) + sum(results[[i]][which(results[[i]]$Program == 'Microglia'),'Niche_3']) + sum(results[[i]][which(results[[i]]$Program == 'T-cell'),'Niche_3']),
                      Niche_4 = sum(results[[i]][which(results[[i]]$Program == 'Myeloid'),'Niche_4']) + sum(results[[i]][which(results[[i]]$Program == 'Microglia'),'Niche_4']) + sum(results[[i]][which(results[[i]]$Program == 'T-cell'),'Niche_4']))
    results[[i]] <- results[[i]][-which(results[[i]]$Program == 'Myeloid' | results[[i]]$Program == 'Microglia' | results[[i]]$Program == 'T-cell'),]
    results[[i]] <- rbind(results[[i]], tmp)
  } else if (('Myeloid' %in% results[[i]]$Program) & ('Microglia' %in% results[[i]]$Program)) {
    tmp <- data.frame(Program = 'Immune',
                      Niche_1 = sum(results[[i]][which(results[[i]]$Program == 'Myeloid'),'Niche_1']) + sum(results[[i]][which(results[[i]]$Program == 'Microglia'),'Niche_1']),
                      Niche_2 = sum(results[[i]][which(results[[i]]$Program == 'Myeloid'),'Niche_2']) + sum(results[[i]][which(results[[i]]$Program == 'Microglia'),'Niche_2']),
                      Niche_3 = sum(results[[i]][which(results[[i]]$Program == 'Myeloid'),'Niche_3']) + sum(results[[i]][which(results[[i]]$Program == 'Microglia'),'Niche_3']),
                      Niche_4 = sum(results[[i]][which(results[[i]]$Program == 'Myeloid'),'Niche_4']) + sum(results[[i]][which(results[[i]]$Program == 'Microglia'),'Niche_4']))
    results[[i]] <- results[[i]][-which(results[[i]]$Program == 'Myeloid' | results[[i]]$Program == 'Microglia'),]
    results[[i]] <- rbind(results[[i]], tmp)
  } else if (('Myeloid' %in% results[[i]]$Program) & ('T-cell' %in% results[[i]]$Program)) {
    tmp <- data.frame(Program = 'Immune',
                      Niche_1 = sum(results[[i]][which(results[[i]]$Program == 'Myeloid'),'Niche_1']) + sum(results[[i]][which(results[[i]]$Program == 'T-cell'),'Niche_1']),
                      Niche_2 = sum(results[[i]][which(results[[i]]$Program == 'Myeloid'),'Niche_2']) + sum(results[[i]][which(results[[i]]$Program == 'T-cell'),'Niche_2']),
                      Niche_3 = sum(results[[i]][which(results[[i]]$Program == 'Myeloid'),'Niche_3']) + sum(results[[i]][which(results[[i]]$Program == 'T-cell'),'Niche_3']),
                      Niche_4 = sum(results[[i]][which(results[[i]]$Program == 'Myeloid'),'Niche_4']) + sum(results[[i]][which(results[[i]]$Program == 'T-cell'),'Niche_4']))
    results[[i]] <- results[[i]][-which(results[[i]]$Program == 'Myeloid' | results[[i]]$Program == 'T-cell'),]
    results[[i]] <- rbind(results[[i]], tmp)
  } else if (('Microglia' %in% results[[i]]$Program) & ('T-cell' %in% results[[i]]$Program)) {
    tmp <- data.frame(Program = 'Immune',
                      Niche_1 = sum(results[[i]][which(results[[i]]$Program == 'Microglia'),'Niche_1']) + sum(results[[i]][which(results[[i]]$Program == 'T-cell'),'Niche_1']),
                      Niche_2 = sum(results[[i]][which(results[[i]]$Program == 'Microglia'),'Niche_2']) + sum(results[[i]][which(results[[i]]$Program == 'T-cell'),'Niche_2']),
                      Niche_3 = sum(results[[i]][which(results[[i]]$Program == 'Microglia'),'Niche_3']) + sum(results[[i]][which(results[[i]]$Program == 'T-cell'),'Niche_3']),
                      Niche_4 = sum(results[[i]][which(results[[i]]$Program == 'Microglia'),'Niche_4']) + sum(results[[i]][which(results[[i]]$Program == 'T-cell'),'Niche_4']))
    results[[i]] <- results[[i]][-which(results[[i]]$Program == 'Microglia' | results[[i]]$Program == 'T-cell'),]
    results[[i]] <- rbind(results[[i]], tmp)
  } else if ('Myeloid' %in% results[[i]]$Program) {
    results[[i]][which(results[[i]]$Program == 'Myeloid'), 'Program'] <- 'Immune'
  } else if ('Microglia' %in% results[[i]]$Program) {
    results[[i]][which(results[[i]]$Program == 'Microglia'), 'Program'] <- 'Immune'
  } else if ('T-cell' %in% results[[i]]$Program) {
    results[[i]][which(results[[i]]$Program == 'T-cell'), 'Program'] <- 'Immune'
  }
  
  results[[i]] <- results[[i]] %>%
    reshape2::melt() %>%
    mutate(variable = glue('{variable}-{samples[i]}'))
}
results <- do.call(rbind, results)


res <- reshape2::acast(results, variable~Program, value.var = "value")
res[is.na(res)] <- 0


anno <- data.frame(sample = gsub("\\_.*", "", gsub('Niche_[1-4]-', '', rownames(res))),
                   region = unlist(lapply(strsplit(rownames(res), '-'), `[[`, 4)),
                   row.names = rownames(res))

rowAnno <- HeatmapAnnotation(df = anno, which = 'row', 
                             col = list(sample = structure(clusterExperiment::bigPalette[1:length(unique(anno$sample))], names = unique(anno$sample)),
                                        region = structure(clusterExperiment::bigPalette[11:(10+length(unique(anno$region)))], names = unique(anno$region))))

ht <- Heatmap(res, right_annotation = rowAnno, column_names_rot = 45, show_row_names = F)
ht <- draw(ht)

rowOrder <- row_order(ht)


hm_colors <- rev((brewer.pal(n=9, name="RdBu")))
hm_colors <- colorRampPalette(colors = hm_colors)


## Hierarchical clustering of nmf factors with correlation among nmf scores
nmf_factor_hc <- clusterNmfFactors(t(res))
nmf_factor_cor <- nmf_factor_hc$cor_coef[nmf_factor_hc$hc_obj$order, nmf_factor_hc$hc_obj$order]

## Heatmap of correlations
hm_colors <- rev((brewer.pal(n=9, name="RdBu")))
hm_colors <- colorRampPalette(colors = hm_colors)

Heatmap(nmf_factor_cor, col = hm_colors(100), cluster_rows = F, cluster_columns = F, show_column_dend = F,
        show_row_dend = F, show_row_names = T, show_column_names = F, name = 'cor', row_names_gp = gpar(fontsize = 8),
        layer_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", lwd = 0.5, fill = NA))})



## Aggregating programs
nClust <- 6
tmp <- cutree(nmf_factor_hc$hc_obj, nClust)

## Ordering metaprograms order
tmp <- tmp[rownames(nmf_factor_cor)]
tmp <- plyr::mapvalues(tmp, unique(tmp), 1:nClust, warn_missing = FALSE)

## NMF metaprograms
nmf_meta_programs <- lapply(sort(unique(tmp)), function(x) names(tmp)[tmp==x])
names(nmf_meta_programs) <- paste("NMF", sort(unique(tmp)), sep = "_")

anno <- NULL
for (i in 1:length(nmf_meta_programs)) {
  anno <- rbind(anno, data.frame(nmf = names(nmf_meta_programs)[i], sample = nmf_meta_programs[[i]]))
}
rownames(anno) <- anno$sample
anno$sample <- NULL
anno$nmf <- factor(anno$nmf, levels = names(nmf_meta_programs))
annoCol <- list(nmf = structure(colors_to_use[1:length(names(nmf_meta_programs))], names = names(nmf_meta_programs)))

rowAnno <- rowAnnotation(rows = anno_text(rownames(nmf_factor_cor), gp = gpar(fontsize = 8, col = colors_to_use[as.numeric(as.factor(gsub("\\_.*", "", gsub('Niche_[1-4]-', '', rownames(nmf_factor_cor)))))])))
column_ha <- HeatmapAnnotation(df = anno, which = 'column', col = list(nmf = annoCol$nmf), show_legend = F, show_annotation_name = F)
row_ha <- HeatmapAnnotation(df = anno, which = 'row', col = list(nmf = annoCol$nmf), show_annotation_name = F)

ht <- Heatmap(nmf_factor_cor, col = hm_colors(100), cluster_rows = F, cluster_columns = F, show_column_dend = F,
              show_row_dend = F, show_row_names = F, show_column_names = F, name = 'cor',
              right_annotation = rowAnno, top_annotation = column_ha, left_annotation = row_ha,
              row_split = anno$nmf, column_split = anno$nmf, column_title = NULL, row_title = NULL,
              layer_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", lwd = 0.5, fill = NA))})


pdf(file.path(plot_dir, paste0('3_Heatmap_recurrent_niches.pdf')),
    width = 10, height = 7)
print(ht)
dev.off()


# calculate and plot average frequency of each cell type in each recurrent niche ---------------------------
# calculate
niches_programs <- lapply(split(anno, f = anno$nmf), rownames)
mean_programs <- lapply(niches_programs, function(x) {
  results %>% filter(variable %in% x) %>%
    group_by(Program) %>%
    summarise(mean = mean(value))
})
lapply(names(mean_programs), function(x) write.xlsx(mean_programs[[x]], file.path(data_dir, 'recurrent_niche_frequency.xlsx'), sheetName=x, append=TRUE))


# Plot 
piechart_plot <- list()

order <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
           "Neuronal-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
           "Immune", "Vessel", "Oligodendrocytes")

for (i in seq_along(mean_programs)){
  mean_programs[[i]] <- mean_programs[[i]][order(match(mean_programs[[i]]$Program, order)), ]
  
  piechart_plot[[i]] <- ggplot(mean_programs[[i]], aes(x = "", y = mean,  fill = Program)) +
    geom_col() + 
    scale_fill_manual(values=colors_recurrent_metaprograms_Xenium) + 
    coord_polar(theta = "y") + 
    labs(title = names(mean_programs)[i]) +
    theme_void() + NoLegend() + 
    theme(plot.title = element_text(hjust = 0.5))
}

plot_grid(plotlist = piechart_plot, ncol = 1) 
ggsave(file.path(plot_dir, paste0('4_PieChart_recurrent_distribution.pdf')), width=3, height=20)









# Identify recurrent patterns of niche distribution in order vs disordered samples ------------------------------------------------

# split into ordered vs disordered (top 11 vs bottom 12)
disordered <- c('STEPN-10_0010553-Region_4',
                'STEPN-14_0010814-Region_2',
                'STEPN-06_0010619-Region_3',
                'STEPN-14_0010814-Region_1',
                'STEPN-06_0010553-Region_2',
                'STEPN-10_0010619-Region_5',
                'STEPN-16_0010775-Region_1',
                'STEPN-06_0010553-Region_1',
                'STEPN-12_0010575-Region_1',
                'STEPN-12_0010540-Region_4',
                'STEPN-12_0010575-Region_3',
                'STEPN-12_0010575-Region_2')
disordered <- paste(disordered, collapse = "|")
results_disordered <- results[grepl("STEPN-10_0010553-Region_4|STEPN-14_0010814-Region_2|STEPN-06_0010619-Region_3|STEPN-14_0010814-Region_1|STEPN-06_0010553-Region_2|STEPN-10_0010619-Region_5|STEPN-16_0010775-Region_1|STEPN-06_0010553-Region_1|STEPN-12_0010575-Region_1|STEPN-12_0010540-Region_4|STEPN-12_0010575-Region_3|STEPN-12_0010575-Region_2", results$variable), ]

ordered <- c('STEPN-10_0010619-Region_4',
             'STEPN-06_0010619-Region_1',
             'STEPN-06_0010619-Region_2',
             'STEPN-15_0010553-Region_3',
             'STEPN-18_0010498-Region_1',
             'STEPN-17_0010652-Region_4',
             'STEPN-01_0010540-Region_2',
             'STEPN-01_0010540-Region_1',
             'STEPN-01_0010540-Region_3',
             'STEPN-19_0010501-Region_1',
             'STEPN-19_0010501-Region_2')
ordered <- paste(ordered, collapse = "|")
results_ordered <- results[grepl("STEPN-10_0010619-Region_4|STEPN-06_0010619-Region_1|STEPN-06_0010619-Region_2|STEPN-15_0010553-Region_3|STEPN-18_0010498-Region_1|STEPN-17_0010652-Region_4|STEPN-01_0010540-Region_2|STEPN-01_0010540-Region_1|STEPN-01_0010540-Region_3|STEPN-19_0010501-Region_1|STEPN-19_0010501-Region_2", results$variable), ]


# now do correlation for disordered tumors ------------------------------------------------
res <- reshape2::acast(results_disordered, variable~Program, value.var = "value")
res[is.na(res)] <- 0



anno <- data.frame(sample = gsub("\\_.*", "", gsub('Niche_[1-4]-', '', rownames(res))),
                   region = unlist(lapply(strsplit(rownames(res), '-'), `[[`, 4)),
                   row.names = rownames(res))

rowAnno <- HeatmapAnnotation(df = anno, which = 'row', col = list(sample = structure(clusterExperiment::bigPalette[1:length(unique(anno$sample))], names = unique(anno$sample)),
                                                                  region = structure(clusterExperiment::bigPalette[11:(10+length(unique(anno$region)))], names = unique(anno$region))))

ht <- Heatmap(res, right_annotation = rowAnno, column_names_rot = 45, show_row_names = F)
ht <- draw(ht)

rowOrder <- row_order(ht)


hm_colors <- rev((brewer.pal(n=9, name="RdBu")))
hm_colors <- colorRampPalette(colors = hm_colors)


## Hierarchical clustering of nmf factors with correlation among nmf scores
nmf_factor_hc <- clusterNmfFactors(t(res))
nmf_factor_cor <- nmf_factor_hc$cor_coef[nmf_factor_hc$hc_obj$order, nmf_factor_hc$hc_obj$order]

## Heatmap of correlations
hm_colors <- rev((brewer.pal(n=9, name="RdBu")))
hm_colors <- colorRampPalette(colors = hm_colors)

Heatmap(nmf_factor_cor, col = hm_colors(100), cluster_rows = F, cluster_columns = F, show_column_dend = F,
        show_row_dend = F, show_row_names = T, show_column_names = F, name = 'cor', row_names_gp = gpar(fontsize = 8),
        layer_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", lwd = 0.5, fill = NA))})



## Aggregating programs
nClust <- 5
tmp <- cutree(nmf_factor_hc$hc_obj, nClust)

## Ordering metaprograms order
tmp <- tmp[rownames(nmf_factor_cor)]
tmp <- plyr::mapvalues(tmp, unique(tmp), 1:nClust, warn_missing = FALSE)

## NMF metaprograms
nmf_meta_programs <- lapply(sort(unique(tmp)), function(x) names(tmp)[tmp==x])
names(nmf_meta_programs) <- paste("NMF", sort(unique(tmp)), sep = "_")

anno <- NULL
for (i in 1:length(nmf_meta_programs)) {
  anno <- rbind(anno, data.frame(nmf = names(nmf_meta_programs)[i], sample = nmf_meta_programs[[i]]))
}
rownames(anno) <- anno$sample
anno$sample <- NULL
anno$nmf <- factor(anno$nmf, levels = names(nmf_meta_programs))
annoCol <- list(nmf = structure(colors_to_use[1:length(names(nmf_meta_programs))], names = names(nmf_meta_programs)))

rowAnno <- rowAnnotation(rows = anno_text(rownames(nmf_factor_cor), gp = gpar(fontsize = 8, col = colors_to_use[as.numeric(as.factor(gsub("\\_.*", "", gsub('Niche_[1-4]-', '', rownames(nmf_factor_cor)))))])))
column_ha <- HeatmapAnnotation(df = anno, which = 'column', col = list(nmf = annoCol$nmf), show_legend = F, show_annotation_name = F)
row_ha <- HeatmapAnnotation(df = anno, which = 'row', col = list(nmf = annoCol$nmf), show_annotation_name = F)

ht <- Heatmap(nmf_factor_cor, col = hm_colors(100), cluster_rows = F, cluster_columns = F, show_column_dend = F,
              show_row_dend = F, show_row_names = F, show_column_names = F, name = 'cor',
              right_annotation = rowAnno, top_annotation = column_ha, left_annotation = row_ha,
              row_split = anno$nmf, column_split = anno$nmf, column_title = NULL, row_title = NULL,
              layer_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", lwd = 0.5, fill = NA))})


pdf(file.path(plot_dir, paste0('3_Heatmap_recurrent_niches_disordered.pdf')),
    width = 10, height = 7)
print(ht)
dev.off()

# calculate and plot average frequency of each cell type in each recurrent niche 
# calculate for all tumors together 
niches_programs <- lapply(split(anno, f = anno$nmf), rownames)
mean_programs <- lapply(niches_programs, function(x) {
  results %>% filter(variable %in% x) %>%
    group_by(Program) %>%
    summarise(mean = mean(value))
})
lapply(names(mean_programs), function(x) write.xlsx(mean_programs[[x]], file.path(data_dir, 'recurrent_niche_frequency.xlsx'), sheetName=x, append=TRUE))


# Plot 
piechart_plot <- list()

order <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
           "Neuronal-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
           "Immune", "Vessel", "Oligodendrocytes")

for (i in seq_along(mean_programs)){
  mean_programs[[i]] <- mean_programs[[i]][order(match(mean_programs[[i]]$Program, order)), ]
  
  piechart_plot[[i]] <- ggplot(mean_programs[[i]], aes(x = "", y = mean,  fill = Program)) +
    geom_col() + 
    scale_fill_manual(values=colors_recurrent_metaprograms_Xenium) + 
    coord_polar(theta = "y") + 
    labs(title = names(mean_programs)[i]) +
    theme_void() + NoLegend() + 
    theme(plot.title = element_text(hjust = 0.5))
}

plot_grid(plotlist = piechart_plot, ncol = 1) 
ggsave(file.path(plot_dir, paste0('4_PieChart_recurrent_distribution_disordered.pdf')), width=3, height=20)




# now do correlation for ordered tumors ------------------------------------------------
res <- reshape2::acast(results_ordered, variable~Program, value.var = "value")
res[is.na(res)] <- 0



anno <- data.frame(sample = gsub("\\_.*", "", gsub('Niche_[1-4]-', '', rownames(res))),
                   region = unlist(lapply(strsplit(rownames(res), '-'), `[[`, 4)),
                   row.names = rownames(res))

rowAnno <- HeatmapAnnotation(df = anno, which = 'row', col = list(sample = structure(clusterExperiment::bigPalette[1:length(unique(anno$sample))], names = unique(anno$sample)),
                                                                  region = structure(clusterExperiment::bigPalette[11:(10+length(unique(anno$region)))], names = unique(anno$region))))

ht <- Heatmap(res, right_annotation = rowAnno, column_names_rot = 45, show_row_names = F)
ht <- draw(ht)

rowOrder <- row_order(ht)


hm_colors <- rev((brewer.pal(n=9, name="RdBu")))
hm_colors <- colorRampPalette(colors = hm_colors)


## Hierarchical clustering of nmf factors with correlation among nmf scores
nmf_factor_hc <- clusterNmfFactors(t(res))
nmf_factor_cor <- nmf_factor_hc$cor_coef[nmf_factor_hc$hc_obj$order, nmf_factor_hc$hc_obj$order]

## Heatmap of correlations
hm_colors <- rev((brewer.pal(n=9, name="RdBu")))
hm_colors <- colorRampPalette(colors = hm_colors)

Heatmap(nmf_factor_cor, col = hm_colors(100), cluster_rows = F, cluster_columns = F, show_column_dend = F,
        show_row_dend = F, show_row_names = T, show_column_names = F, name = 'cor', row_names_gp = gpar(fontsize = 8),
        layer_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", lwd = 0.5, fill = NA))})



## Aggregating programs
nClust <- 6
tmp <- cutree(nmf_factor_hc$hc_obj, nClust)

## Ordering metaprograms order
tmp <- tmp[rownames(nmf_factor_cor)]
tmp <- plyr::mapvalues(tmp, unique(tmp), 1:nClust, warn_missing = FALSE)

## NMF metaprograms
nmf_meta_programs <- lapply(sort(unique(tmp)), function(x) names(tmp)[tmp==x])
names(nmf_meta_programs) <- paste("NMF", sort(unique(tmp)), sep = "_")

anno <- NULL
for (i in 1:length(nmf_meta_programs)) {
  anno <- rbind(anno, data.frame(nmf = names(nmf_meta_programs)[i], sample = nmf_meta_programs[[i]]))
}
rownames(anno) <- anno$sample
anno$sample <- NULL
anno$nmf <- factor(anno$nmf, levels = names(nmf_meta_programs))
annoCol <- list(nmf = structure(colors_to_use[1:length(names(nmf_meta_programs))], names = names(nmf_meta_programs)))

rowAnno <- rowAnnotation(rows = anno_text(rownames(nmf_factor_cor), gp = gpar(fontsize = 8, col = colors_to_use[as.numeric(as.factor(gsub("\\_.*", "", gsub('Niche_[1-4]-', '', rownames(nmf_factor_cor)))))])))
column_ha <- HeatmapAnnotation(df = anno, which = 'column', col = list(nmf = annoCol$nmf), show_legend = F, show_annotation_name = F)
row_ha <- HeatmapAnnotation(df = anno, which = 'row', col = list(nmf = annoCol$nmf), show_annotation_name = F)

ht <- Heatmap(nmf_factor_cor, col = hm_colors(100), cluster_rows = F, cluster_columns = F, show_column_dend = F,
              show_row_dend = F, show_row_names = F, show_column_names = F, name = 'cor',
              right_annotation = rowAnno, top_annotation = column_ha, left_annotation = row_ha,
              row_split = anno$nmf, column_split = anno$nmf, column_title = NULL, row_title = NULL,
              layer_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", lwd = 0.5, fill = NA))})


pdf(file.path(plot_dir, paste0('3_Heatmap_recurrent_niches_ordered.pdf')),
    width = 10, height = 7)
print(ht)
dev.off()




# calculate and plot average frequency of each cell type in each recurrent niche 
# calculate for all tumors together 
niches_programs <- lapply(split(anno, f = anno$nmf), rownames)
mean_programs <- lapply(niches_programs, function(x) {
  results %>% filter(variable %in% x) %>%
    group_by(Program) %>%
    summarise(mean = mean(value))
})
lapply(names(mean_programs), function(x) write.xlsx(mean_programs[[x]], file.path(data_dir, 'recurrent_niche_frequency.xlsx'), sheetName=x, append=TRUE))


# Plot 
piechart_plot <- list()

order <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
           "Neuronal-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
           "Immune", "Vessel", "Oligodendrocytes")

for (i in seq_along(mean_programs)){
  mean_programs[[i]] <- mean_programs[[i]][order(match(mean_programs[[i]]$Program, order)), ]
  
  piechart_plot[[i]] <- ggplot(mean_programs[[i]], aes(x = "", y = mean,  fill = Program)) +
    geom_col() + 
    scale_fill_manual(values=colors_recurrent_metaprograms_Xenium) + 
    coord_polar(theta = "y") + 
    labs(title = names(mean_programs)[i]) +
    theme_void() + NoLegend() + 
    theme(plot.title = element_text(hjust = 0.5))
}

plot_grid(plotlist = piechart_plot, ncol = 1) 
ggsave(file.path(plot_dir, paste0('4_PieChart_recurrent_distribution_ordered.pdf')), width=3, height=20)







# To plot zoomed images -----------------------------------------------------
i = 22
data <- qread(file.path(seurat_obj_dir, paste0('seurat_obj_', metadata$SampleID[i], "_", metadata$SampleName[i], '.qs')))

ImageDimPlot(data, group.by = 'group',  fov = "fov", cols = colors_metaprograms_Xenium,  border.size = NA, size = 0.5,
             dark.background = T, axes = TRUE, #border.color = "white", border.size = 0.1, coord.fixed = FALSE, nmols = 10000
             ) + 
  #scale_x_continuous(breaks = seq(0,5000,500)) + 
  #scale_y_continuous(breaks = seq(0,5000,500)) + 
  geom_rect(xmin = 2000, xmax = 3000, ymin = 1000, ymax = 2000, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 1500, xmax = 2000, ymin = 4500, ymax = 5500, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 4100, xmax = 4750, ymin = 2100, ymax = 2900, color = 'black', fill = NA, lwd = 0.2) + 
  seurat_theme()
ggsave(file.path(plot_dir, paste0('5_ImageDimPlot_malignant_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=6, height=5)

data[["zoom_1"]] <- Crop(data[["fov"]], x = c(2000, 3000), y = c(1000, 2000), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1", size = 1, border.color = "white", border.size = 0.1, cols = colors_metaprograms_Xenium) + seurat_theme()
ggsave(file.path(plot_dir, paste0('5_Zoom_ImageDimPlot_celltype_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)

data[["zoom_1"]] <- Crop(data[["fov"]], x = c(2000, 3000), y = c(1000, 2000), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1", group.by = 'niches', size = 1, border.color = "white", border.size = 0.1, cols = col_niches) + seurat_theme()
ggsave(file.path(plot_dir, paste0('5_Zoom_ImageDimPlot_niche_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)



i = 16
data <- qread(file.path(seurat_obj_dir, paste0('seurat_obj_', metadata$SampleID[i], "_", metadata$SampleName[i], '.qs')))

ImageDimPlot(data, group.by = 'group',  fov = "fov", cols = colors_metaprograms_Xenium,  border.size = NA, size = 0.5,
             dark.background = T, axes = TRUE, #border.color = "white", border.size = 0.1, coord.fixed = FALSE, nmols = 10000
) + 
  #scale_x_continuous(breaks = seq(0,5000,500)) + 
  #scale_y_continuous(breaks = seq(0,5000,500)) + 
  geom_rect(xmin = 1500, xmax = 2500, ymin = 2000, ymax = 3000, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 1500, xmax = 2000, ymin = 4500, ymax = 5500, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 4100, xmax = 4750, ymin = 2100, ymax = 2900, color = 'black', fill = NA, lwd = 0.2) + 
  seurat_theme()
ggsave(file.path(plot_dir, paste0('6_ImageDimPlot_malignant_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=6, height=5)


data[["zoom_1"]] <- Crop(data[["fov"]], x = c(1500, 2500), y = c(2000, 3000), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1", size = 1, border.color = "white", border.size = 0.1, cols = colors_metaprograms_Xenium) + seurat_theme()
ggsave(file.path(plot_dir, paste0('6_Zoom_ImageDimPlot_celltype_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)

data[["zoom_1"]] <- Crop(data[["fov"]], x = c(1500, 2500), y = c(2000, 3000), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1", group.by = 'niches', size = 1, border.color = "white", border.size = 0.1, cols = col_niches) + seurat_theme()
ggsave(file.path(plot_dir, paste0('6_Zoom_ImageDimPlot_niche_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)




i = 19
data <- qread(file.path(seurat_obj_dir, paste0('seurat_obj_', metadata$SampleID[i], "_", metadata$SampleName[i], '.qs')))

ImageDimPlot(data, group.by = 'group',  fov = "fov", cols = colors_metaprograms_Xenium,  border.size = NA, size = 0.5,
             dark.background = T, axes = TRUE, #border.color = "white", border.size = 0.1, coord.fixed = FALSE, nmols = 10000
) + 
  #scale_x_continuous(breaks = seq(0,5000,500)) + 
  #scale_y_continuous(breaks = seq(0,5000,500)) + 
  geom_rect(xmin = 1300, xmax = 2500, ymin = 1500, ymax = 2500, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 1500, xmax = 2000, ymin = 4500, ymax = 5500, color = 'black', fill = NA, lwd = 0.2) + 
  #geom_rect(xmin = 4100, xmax = 4750, ymin = 2100, ymax = 2900, color = 'black', fill = NA, lwd = 0.2) + 
  seurat_theme()
ggsave(file.path(plot_dir, paste0('6_ImageDimPlot_malignant_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=6, height=5)


data[["zoom_1"]] <- Crop(data[["fov"]], x = c(1300, 2500), y = c(1500, 2500), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1", size = 1, border.color = "white", border.size = 0.1, cols = colors_metaprograms_Xenium) + seurat_theme()
ggsave(file.path(plot_dir, paste0('6_Zoom_ImageDimPlot_celltype_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)

data[["zoom_1"]] <- Crop(data[["fov"]], x = c(1300, 2500), y = c(1500, 2500), coords = "plot")
DefaultBoundary(data[["zoom_1"]]) <- "segmentation"
ImageDimPlot(data, fov = "zoom_1", group.by = 'niches', size = 1, border.color = "white", border.size = 0.1, cols = col_niches) + seurat_theme()
ggsave(file.path(plot_dir, paste0('6_Zoom_ImageDimPlot_niche_', metadata$SampleID[i],"_", metadata$SampleName[i], '.pdf')), width=8, height=5)


