# Load packages -----------------------------------
rm(list = ls())

library(data.table)
library(tidyverse)
library(R.utils)
library(ggpubr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(writexl)
library(tidyverse)
#library(paletteer)
library(readxl)
library(cowplot)
#library(scCustomize)
#library(ComplexHeatmap)
library(qs)
library(ggrastr)

# Organize environment  -----------------------------------
base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/7 - R2 gene expression/ZFTA_vs_YAP1"
analysis_dir <- file.path(base_dir, 'analysis')
resource_dir <- file.path(base_dir, 'scripts/resources')

ref_dir <- file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Developmental datasets/data/Nowakowski_Eze')

plot_dir <- file.path(analysis_dir)

source(file.path(resource_dir, 'Plot_style.R'))

# define color palette
subtype = c("ZFTA-RELA", "ST-YAP1")
col_subtype <- c("#B44672","#B4AF46")
names(col_subtype) <- subtype

# load R2 data  -----------------------------------
# Load gene expression matrix (Z-score)
data <- read_excel(file.path(base_dir, "data/removal_discovery_cohort/ps_avgpres_gse64415geo209_u133p2_box1700152812-datagrabber-.xlsx"))
  # Remove duplicates based on row names
  data <- data[!duplicated(data$gene_id), ]
  # make patient id as rownames
  data <- data %>%
    tibble::column_to_rownames(var = "gene_id")


# Load metadata with patient information
metadata <- read_excel(file.path(base_dir, "data/removal_discovery_cohort/metadata.xlsx"))

# subset R2 dataframe to samples present in metadata (= remove samples that were used in scRNAseq cohort)
data <- data[, intersect(metadata$probeset, colnames(data))]

# Create Seurat object
data_seurat <- CreateSeuratObject(counts = data)
data_seurat[["RNA"]]$data <- data_seurat[["RNA"]]$counts
data_seurat <- AddMetaData(data_seurat, metadata = metadata)

# subset Seurat to patients with YAP1 or RELA fusions only
Idents(data_seurat) = 'molecular_subgroup'
data_seurat <- subset(data_seurat, idents = c('st_epn_yap1', 'st_epn_rela'))

# Rename sample names and populations for final figures -----------------------------------

# Rename sample names
metadata <- data_seurat@meta.data
names_replace <- unique(data_seurat$molecular_subgroup)
new_names <- c("ZFTA-RELA", "ST-YAP1")
# Use dplyr's mutate to substitute specific cells with new names
metadata <-  metadata %>%
  mutate(molecular_subgroup = case_when(
    molecular_subgroup == names_replace[1] ~ new_names[1],
    molecular_subgroup == names_replace[2] ~ new_names[2]
  ))

# reorder based on rownames of Seurat object
metadata <- metadata[order(match(metadata$probeset, colnames(data_seurat))), ]

data_seurat <- AddMetaData(data_seurat, metadata)




# load signatures  -----------------------------------
# load NMF gene list
markers_NMF <- readRDS(file.path(base_dir, 'data/signatures/nmf_marker_genes_final_annotated.rds'))
names(markers_NMF) <- c("Neuroepithelial-like", "MES-like","Ependymal-like","Radial-glia-like","Cycling", "Embryonic-neuronal-like", "Neuronal-like",            "Embryonic-like")

# score dataset -----------------------------------
scores <- c(markers_NMF['Neuroepithelial-like'], markers_NMF['Ependymal-like'], markers_NMF['Neuronal-like'], markers_NMF['Embryonic-like'], markers_NMF['Embryonic-neuronal-like'])
names(scores) <- c('Neuroepithelial-like', 'Ependymal-like', 'Neuronal-like', 'Embryonic-like', 'Embryonic-neuronal-like')
data_seurat <- AddModuleScore(object = data_seurat, assay = 'RNA', features = scores, name = names(scores), seed=5)

# rename metadata names of scores
# identify number of first column with metadata scores
col_start <- length(colnames(data_seurat@meta.data)) - length(names(scores)) +1
# identify number of last column with metadata scores
col_end <- length(colnames(data_seurat@meta.data))
# rename columns with score name
colnames(data_seurat@meta.data)[col_start:col_end] <- names(scores)


# Violin plots of scores -----------------------------------
scores <- c("Neuroepithelial-like",  "Neuronal-like", "Ependymal-like", 'Embryonic-like')
titles <- c("Neuroepithelial-like",  "Neuronal-like", "Ependymal-like", 'Embryonic-like')

p <- VlnPlot(data_seurat, features = scores, group.by = 'molecular_subgroup', combine = FALSE, pt.size=1) 

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() +
    labs (y='Module score, AU', x='', title = titles[[i]]) + 
    scale_fill_manual(values = col_subtype) + 
    geom_boxplot(outlier.shape=NA, width=0.1, fill="white") +
    # stat_compare_means(method = "anova",  
    #                     size = 3, 
    #                     label.y.npc = 0.99,
    #                     label.x.npc = 0.5)+ # Add global p-value
    stat_compare_means(aes(label = after_stat(p.format)),
                        method = 't.test', 
                        size = 5, 
                        label.y.npc = 0.97, 
                        label.x.npc = 0.5) + 
    #scale_y_continuous(limits = c(-1, 2)) + # Set y-axis limits
    NoLegend() +
    theme_vln
}
plot_grid(plotlist = p, ncol=4)
ggsave(file.path(plot_dir, "1_VlnPlot_scores_no_discovery_cohort.pdf"), width=9, height=4)




# Hierarchy plot  -----------------------------------
x1 = "Neuroepithelial-like"
x2 = "Embryonic-like"
y1 = "Neuronal-like"
y2 = "Ependymal-like"

group.by="molecular_subgroup"
sample=data_seurat
variables_to_retrieve <- c(x1, x2, y1,y2)

# And store them as a tibble.
scores <- sample@meta.data[, variables_to_retrieve]
# Shuffle the cells so that we accomplish a random plotting, not sample by sample.

# Compute Y axis values.
d <- apply(scores, 1, function(x){max(x[c(x1, x2)]) - max(x[c(y1, y2)])})

# Compute X axis values.
x <- vapply(seq_along(d), function(x) {
  if (d[x] > 0) {
    d <- log2(abs(scores[x, x1] - scores[x, x2]) + 1)
    ifelse(scores[x, x1] < scores[x, x2], d, -d)
  } else {
    d <- log2(abs(scores[x, y1] - scores[x, y2]) + 1)
    ifelse(scores[x, y1] < scores[x, y2], d, -d)
  }
}, FUN.VALUE = numeric(1))

names(x) <- rownames(scores)

# Define titles for the axis.
x_lab1 <- paste0(y1, "  <---->  ", y2)
x_lab2 <- paste0(x1, "  <---->  ", x2)
y_lab1 <- paste0(y1, "  <---->  ", x1)
y_lab2 <- paste0(x2, "  <---->  ", y2)


# Plot.
df_E <- data.frame(row.names = rownames(scores))
df_E[["set_x"]] <- x
df_E[["set_y"]] <- d
df_E[["group.by"]] <- sample@meta.data[, group.by]
df_E[["subtype"]] <- data_seurat$molecular_subgroup

colors.use=c("#B44672","#B47846","#46B478","#46B4AF","#4682B4","#B4AF46")
names(colors.use) = c("ZFTA-RELA","ZFTA-Cluster 1","ZFTA-Cluster 2","ZFTA-Cluster 3","ZFTA-Cluster 4","ST-YAP1")


samples = unique(df_E$group.by)
samples = c("ZFTA-RELA","ZFTA-Cluster 1","ZFTA-Cluster 2","ZFTA-Cluster 3","ZFTA-Cluster 4", "ST-YAP1")

for (i in 1:length(samples)) {
  tmp <- df_E$group.by==samples[i]
  tmp <- gsub("TRUE",samples[i],tmp)
  tmp <- gsub("FALSE","Other",tmp)
  df_E <- cbind(df_E,tmp)
}
subtypes <- c("ZFTA-RELA","ZFTA-Cluster 1","ZFTA-Cluster 2","ZFTA-Cluster 3","ZFTA-Cluster 4", "ST-YAP1")
colnames(df_E)[5:10] <- subtypes

ggplot(df_E,aes(x=set_x,
                y=set_y,
                col=group.by))+
  geom_point(size=3)+
  scale_color_manual(values=colors.use)+
  scale_y_continuous(limits = c(-1.5, 1.5)) + # Set y-axis limits
  scale_y_continuous(limits = c(-1.5, 1.5)) + # Set y-axis limits
  theme_cellstate_plot + 
  labs(x = 'Relative meta-module score', y = 'Relative meta-module score', title = 'ST-EPN cellular hierarchy', subtitle = 'n = 6,933 cells') +
  geom_vline(xintercept = 0, linetype="dotted") + geom_hline(yintercept = 0, linetype="dotted") 
ggsave(file.path(plot_dir, "2_CellState_plot_no_discovery_cohort_nolegend.pdf"), width=5.0, height=3.5)

