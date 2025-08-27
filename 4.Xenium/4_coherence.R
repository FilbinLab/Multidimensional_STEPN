# this file takes a single input, the sample_name (such as 'STEPN47_Region_2'). One can simply create a variable 
# by running something like: "sample_name <- 'STEPN47_Region_2'", and run this script line by line, skipping the 
# first part where the sample_name is read from supplied arguments, or can execute this script from the terminal

# reading argument supplied if ran from terminal or slurm job
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

# optionally, one can manually run this file line by line, when the value of sample_name is initiated here only
# sample_name <- 'STEPN47_Region_2'

# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/FLXenium/functions/')

# source the ependymoma project specific functions
source('~/ependymoma/xenium/scripts_revisions/resources/epn_functions.R')

# get the vars for epn analyses
all_vars <- SetUpEpendymomaGlobalVars(sample_name)
coherence_dir <- all_vars$coherence_dir
data_dir <- glue('{coherence_dir}/data')
plot_dir <- glue('{coherence_dir}/plots')

# Read annotated datasets
data <- qread(all_vars$annotation_file)

# get the cell feature matrix from data
cell_feature_matrix <- LayerData(data, 'counts', 'SCT')

# NOTE: while extracting the x and y coordinates of centroids of cells, we have to be super careful, and reverse 
# them to maintain x coordinates being horizontal (in object's image information, x is vertical, hence we need to reverse it)
metadata <- data.frame(cell_id = data@images$fov$centroids@cells,
                       y_centroid = data@images$fov$centroids@coords[,1],
                       x_centroid = data@images$fov$centroids@coords[,2])

# Pre-processing data (filtering some cells and genes and col renaming in metadata)
preprocessing <- prepareData(cell_feature_matrix, cells_info = metadata)

# Since most of the slices are almost square, using same number of bins for rows and cols
nbins <- sqrt(dim(data@assays$SCT)[2])
nbinsrowcol <- c(nbins, nbins)

# now, using the nbinsrowcol, make the grid with appropriate number of rows and cols 
# make a grid over the cells, and return a df containing each grid-square's x and y coordinates, and the assigned 
# most common celltype (from those beneath a square of the grid) for each square.
gridding <- grid_spatial(norm_data = data@assays$SCT$data, spatial = preprocessing$spatial,
                         variable = data$cell_type, nbinsrow = round(nbinsrowcol[1]), nbinscol = round(nbinsrowcol[2]))

# compute the coherence, meaning finding out the number of neighbors for each square in the grid with same associated celltype as itself
coherence <- coherence_score(grid_df = gridding, variable = 'Metaprogram', empty_neighbors_increase_coherence = FALSE)

# Save output 
qsave(coherence, glue('{data_dir}/results_{sample_name}.qs'))

# Plot grids with density results
plot <- gridding_coherence_density(gridding, coherence, sample_name, color_coherence_density, colors_metaprograms_Xenium)
ggsave(glue('{plot_dir}/{sample_name}.pdf'), width = 18, height = 8)


