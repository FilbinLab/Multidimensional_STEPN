# code for main figure 3 panel b

# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/Multidimensional_STEPN/4.Xenium/resources/FLXenium/')

# source the ependymoma project specific functions
source('~/Multidimensional_STEPN/4.Xenium/resources/epn_functions.R')

# get the vars for epn analyses
sample_name = 'STEPN19_Region_1'
all_vars <- SetUpEpendymomaGlobalVars(sample_name)
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

# Since most of the slices are almost square, always using same number of bins for rows and cols
nbins <- sqrt(dim(data@assays$SCT)[2])
nbinsrowcol <- c(nbins, nbins)

# now, using the nbinsrowcol, make the grid with appropriate number of rows and cols 
# make a grid over the cells, and return a df containing each grid-square's x and y coordinates, and the assigned 
# most common celltype (from those beneath a square of the grid) for each square.
gridding <- grid_spatial(norm_data = data@assays$SCT$data, spatial = preprocessing$spatial,
                         variable = data$cell_type, nbinsrow = round(nbinsrowcol[1]), nbinscol = round(nbinsrowcol[2]))

# compute the coherence, meaning finding out the number of neighbors for each square in the grid with same associated celltype as itself
coherence <- coherence_score(grid_df = gridding, variable = 'Metaprogram', empty_neighbors_increase_coherence = FALSE)

# make the coherence gridding plot
plot <- gridding_coherence_density(gridding, coherence, sample_name, color_coherence_density, colors_metaprograms_Xenium)







