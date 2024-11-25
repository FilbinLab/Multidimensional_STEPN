# Load packages -----------------------------------
library(Seurat)
library(tidyverse)
library(glue)
library(gplots)
library(qs)
library(reshape)


# Organize environment  -----------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

coherence_dir <- file.path(base_dir, 'analysis/6_coherence/data')
cellid_dir <- file.path(base_dir, 'analysis/3_program_annotation/data')

plot_dir  <- file.path(base_dir, 'analysis/6_coherence/plots')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}


resource_dir  <- file.path(base_dir, 'scripts_revisions/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))
source(file.path(resource_dir, "xenium_preprocessing_helper_functions_SD_v2.R"))
source(file.path(resource_dir, 'color_palette.R'))




## Read metadata ------------------------------------------------------
metadata <- read_xlsx(file.path(base_dir, 'scripts_revisions/SampleIdentifier.xlsx'))
metadata <- metadata %>% filter(!SampleName %in% "STEPN14_Region_2")

# reorder by sampleID
metadata <- metadata %>% arrange(SampleName)

FileName <- metadata$SampleName

order_metaprograms <- c("Cycling",  "Embryonic-like", "Neuroepithelial-like", "Radial-glia-like", 
                        "Embryonic-neuronal-like", 
                        "Neuronal-like" ,"Ependymal-like", "MES-like", 
                        "T-cells", "Myeloid",  "Endothelial",  "Oligodendrocytes", 'Astrocyte', 'Neurons', 'Unknown')


# Box plots by MP and sample ----------------------

plot_coherence_score(FileName, cellid_dir, coherence_dir, plot_dir, colors_metaprograms_Xenium, col_sampling) 



# Grid plots ----------------------


# Read annotated datasets
SampleName_arg = "STEPN19_Region_2"

data <- qread(file.path(seurat_obj_dir, paste0('seurat_obj_', SampleName_arg, '.qs')))
print(paste0('Successful reading of ', SampleName_arg))

adata <- LayerData(data, 'counts', 'SCT')
metadata <- data.frame(cell_id = data@images$fov$centroids@cells,
                       x_centroid = data@images$fov$centroids@coords[,1],
                       y_centroid = data@images$fov$centroids@coords[,2])


## Pre-processing data
preprocessing <- prepareData(adata, metadata)

## Normalizing data
normalized_matrix <- norm_data(t(as.matrix(preprocessing$cell_feature_matrix)))
normalized_matrix <- sqrt(normalized_matrix) + sqrt(normalized_matrix + 1)
normalized_matrix <- norm_data(normalized_matrix)

metaprogram <- plyr::mapvalues(rownames(normalized_matrix), from = names(data$cell_type), to = as.character(unname(data$cell_type)), warn_missing = F)

## Gridding
gridding <- grid_spatial(norm_data = t(normalized_matrix), spatial = preprocessing$spatial, 
                         variable = metaprogram, nbins = round(sqrt(dim(data@assays$SCT)[2])))
print(paste0('nbins = ', round(sqrt(dim(data@assays$SCT)[2]))))

coherence <- coherence_score(grid_df = gridding, variable = 'Metaprogram')

ggplot(gridding, aes(x = y, y = x, color = Metaprogram)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = colors_metaprograms_Xenium) +
  theme_void() +
  guides(colour = guide_legend(override.aes = list(size = 4)))
ggsave(file.path(plot_dir, paste0('1_Gridding', SampleName_arg, '.pdf')), width = 9, height = 8)


mat <- gridding %>% mutate(x = as.character(x)) %>% mutate(y = as.character(y))
mat <- reshape2::acast(mat, y~x, value.var = 'Metaprogram') %>% as.data.frame()
#rownames(mat) <- 1:nrow(mat)
#colnames(mat) <- 1:ncol(mat)


# read spatial results
spatial_coherence <- qread(file.path(coherence_dir, paste0('results_', SampleName_arg, '.qs')))
spatial_coherence$results_df

# Summing all matrices
summed_matrix <- Reduce(`+`, spatial_coherence$results_df)

# View the result
print(summed_matrix)

#change x/y values with actual positions
rownames(summed_matrix) <- rownames(mat) 
colnames(summed_matrix) <- colnames(mat) 

summed_matrix <- reshape2::melt(summed_matrix)

df <- as.data.frame(as.table(summed_matrix))  # Convert matrix to long format
colnames(df) <- c("x", "y", "value")  # Rename columns

# Convert x and y to numeric if necessary
df$x <- as.numeric(as.character(df$x))
df$y <- as.numeric(as.character(df$y))
df$value <- as.numeric(as.character(df$value))


# create vector for color density
col <- c(
  "white", "#E0F3DBFF", "#CCEBC5FF", "#A8DDB5FF", 
  "#7BCCC4FF", "#4EB3D3FF", "#2B8CBEFF", "#0868ACFF", "#084081FF"
)
names(col) <- 0:8
col


# Plot density of spatial cohernece score
ggplot(df, aes(x = x, y = y, color = factor(value))) +
  geom_point(size = 0.2) +
  scale_color_manual(values = col) +
  theme_void() +
  guides(colour = guide_legend(override.aes = list(size = 4)))
ggsave(file.path(plot_dir, paste0('1b_Gridding_coherence', SampleName_arg, '.pdf')), width = 9, height = 8)


