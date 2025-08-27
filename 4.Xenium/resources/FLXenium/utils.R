#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~BASIC UTILITY FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This file contains all the functions which I have most frequently accessed directly (during first few months of FLXenium).

# RunScriptBySubmittingJobs(script = '~/ependymoma/xenium/scripts_revisions/other_analyses/running_scvi_on_all_samples.py', vector_of_inputs = 'temp', cores = 4, ram = 40, time = '08:00', gpu = TRUE, sbatch_location = '~/temp/scvi_run.sbatch', script_title = 'scvi_run')
# RunScriptBySubmittingJobs(script = '~/hope/other_analyses/preprocess_with_subset_N20192360_1A_Region_3.R', vector_of_inputs = 'temp', cores = 4, ram = 100, time = '02:00', gpu = FALSE, sbatch_location = '~/temp/prep_N20192360_1A_Region_3.sbatch', script_title = 'prep_N20192360_1A_Region_3')
# RunScriptBySubmittingJobs(script = '~/dmg_costanza/4_program_annotation.R', vector_of_inputs = heavy_samples, cores = 4, ram = 40, time = '01:00', gpu = FALSE, sbatch_location = '~/temp/annot_dmg_heavy.sbatch', script_title = 'annot_dmg_heavy')
# RunScriptBySubmittingJobs(script = '~/G34RV/patient_samples/3_program_annotation.R', vector_of_inputs = heavy_samples, cores = 4, ram = 40, time = '01:00', gpu = FALSE, sbatch_location = '~/temp/annot_dmg_heavy.sbatch', script_title = 'annot_dmg_heavy')

# run a script for a bunch of input values (originally intended to be sample_names), and submit corresponding jobs.
# time is in format HH:MM. ram is in GBs. This function can be used for just small, simple analyses, and actual, 
# finalized steps, in which case, the default value of sbatch_location should probably be modified.
# NOTE: since originally it was designed for looping through a vector of sample_names, the script should ideally take just one input (sample_name)
RunScriptBySubmittingJobs <- function(script, vector_of_inputs, cores = 4, ram = 40, time = '01:00', gpu = FALSE,
                                      sbatch_location = '~/temp/temp.sbatch', script_title = NULL, log_dir = NULL){
  library(glue)
  
  # according to script extension, assign value to variable script_language. Currently, only python and R are accepted
  if (tools::file_ext(script) == 'py'){ 
    script_language = 'python'
  }else if (tools::file_ext(script) == 'R'){ 
    script_language = 'R'
  }else{
    stop('Incorrect extension of file!')
  }
  
  # by default, we log in the location of the script
  if (is.null(log_dir)) log_dir <- dirname(sbatch_location)
  # by default, we have the script title as the name of the script file without extension
  if (is.null(script_title)) script_title <- tools::file_path_sans_ext(basename(script))
  
  # first, create the sbatch script
  CreateSbatchScript(script, cores, ram, time, gpu, sbatch_location, script_language)
  
  final_command <- ''
  # then, submit the set of commands to run the sbatch script on each input argument
  for (input in vector_of_inputs){
    # if currently in the last iteration, last_command = T else F. This is used to add "; \" at the end of commands so that when we paste in terminal, they execute properly
    # if (input == vector_of_inputs[length(vector_of_inputs)]) last_command <- TRUE else last_command <- FALSE
    command <- SubmitSbatchForInput(input, sbatch_location, script_title = script_title, log_dir = log_dir)
    final_command <- paste0(final_command, command)
  }
  # print the concatenated command
  print(final_command)
  
  print('NOTE: Please close the script now, as it would not be a good idea to make changes into it after submitting the jobs!')
}

# Create a sbatchscript at specified location with the given arguments. Some of the assumptions made: 
# i) the script takes exactly one argument. ii) Gonna load the same, standard set of modules generally 
# done for Rscripts. ram is in GBs, time is in format HH:MM.
# Designed to be called from RunScriptBySubmittingJobs().
CreateSbatchScript <- function(script, cores, ram, time, gpu, sbatch_location, script_language){
  library(glue)
  # create the vector storing the information to be written in the sbatch_script file
  if (script_language == 'python'){
    # decide according to value of supplied argument gpu, if partition gonna be short or gpu
    if (gpu == FALSE){
      slurm_script <- c(
        "#!/bin/bash",
        "#SBATCH -N 1",
        glue("#SBATCH -c {cores}"),
        glue("#SBATCH -t 0-{time}"),
        glue("#SBATCH --mem={ram}GB"),
        "#SBATCH -p short",
        "#SBATCH --job-name %j",
        "#SBATCH -o %j.out",
        "#SBATCH -e %j.err",
        "module purge",
        "module load gcc/14.2.0 conda/miniforge3",
        "conda activate /home/shk490/.conda/envs/sc",
        paste("python", script, "${1}")
      )
    }else{
      slurm_script <- c(
        "#!/bin/bash",
        "#SBATCH -N 1",
        glue("#SBATCH -c {cores}"),
        glue("#SBATCH -t 0-{time}"),
        glue("#SBATCH --mem={ram}GB"),
        "#SBATCH -p gpu",
        "#SBATCH --gres=gpu:1",
        "#SBATCH --job-name %j",
        "#SBATCH -o %j.out",
        "#SBATCH -e %j.err",
        "module purge",
        "module load gcc/14.2.0 conda/miniforge3",
        "conda activate /home/shk490/.conda/envs/sc",
        paste("python", script, "${1}")
      )
    }
  }else if (script_language == 'R'){
    # ensure that if entered this section, then gpu was FALSE, because we generally don't use gpu in r code.
    if(gpu != FALSE) stop('Function CreateSbatchScript currently not configured to run R scripts on GPU!')
    slurm_script <- c(
      "#!/bin/bash",
      "#SBATCH -N 1",
      glue("#SBATCH -c {cores}"),
      glue("#SBATCH -t 0-{time}"),
      glue("#SBATCH --mem={ram}GB"),
      "#SBATCH -p short",
      "#SBATCH --job-name %j",
      "#SBATCH -o %j.out",
      "#SBATCH -e %j.err",
      "module purge",
      "module load gcc/14.2.0 R/4.4.2 python/3.13.1 proj gdal udunits geos imageMagick/7.1.1 hdf5/1.14.5",
      paste("Rscript", script, "${1}")
    )
  }else{
    stop ('Incorrect script language supplied!')
  }
  # now write the script at the appropriate location
  writeLines(slurm_script, sbatch_location)
}

# run a shell command running the sbatch file at sbatch_location, with the one and only argument: 'input'
# script_title is just used for naming the log files etc.
SubmitSbatchForInput <- function(input, sbatch_location, script_title, log_dir){
  library(glue)
  
  # make the command string (last command is just a boolean telling that this is the last command in the loop, so that we don't need to add ': \' at the end of the command)
  command <- glue('sbatch \\
                  -o {log_dir}/{input}_{script_title}.out \\
                  -e {log_dir}/{input}_{script_title}.err \\
                  -J {input}_{script_title} \\
                  {sbatch_location} {input}; ')
  # run the command
  # system(command)
  return (command) # just printing right now, so that can copy and paste in terminal, because not able to run sbatch from within Rstudio for some reason.
}

# Generate the embeddings of all the functions in a specified directory. 
# It saves the embeddings and metadata_df in the out_dir.
GenerateFunctionEmbeddings <- function(files_location = '/home/shk490/FLXenium/functions', 
        out_dir = '/home/shk490/FLXenium/search_utility', flxenium_location = '/home/shk490/FLXenium'){
  # store in a var, the path to script to run
  script_to_run <- glue('{flxenium_location}/search_utility/generate_embeddings.py')
  
  # source the python file to make its functions available in current R env
  print('Sourcing the generate_embeddings.py file...')
  source_python(script_to_run)
  
  # finally, supply the appropriate arguments to the generate_and_save_embeddings function, and generate 
  # the embeddings and metadata_df files
  # I use a fixed model_path, because there seems to be an issue when sentence-transformers tries to download the model when this is run from reticulate
  generate_and_save_embeddings(files_location, out_dir, model_path = '/home/shk490/.cache/huggingface/hub/models--sentence-transformers--all-MiniLM-L6-v2/snapshots/c9745ed1d9f207416be6d2e6f8de32d1f16199bf')
}

# Find a particular function using a description of it.
# embeddings dir is where the embeddings and metadata_df are stored. flxenium_location determines 
# location of the python script query.py (where the actual query function is there)
ff <- function(query, top_k = 5, embeddings_dir = '/home/shk490/FLXenium/search_utility', flxenium_location = '/home/shk490/FLXenium'){
  library(reticulate)
  # query = 'find what is the size of various folders'
  
  # store in a var, the path to script to run
  script_to_run <- glue('{flxenium_location}/search_utility/query.py')
  
  # source the python file to make its functions available in current R env
  # print('Sourcing the query.py file...')
  source_python(script_to_run)
  
  find_function(query, top_k, embeddings_dir, model_path = '/home/shk490/.cache/huggingface/hub/models--sentence-transformers--all-MiniLM-L6-v2/snapshots/c9745ed1d9f207416be6d2e6f8de32d1f16199bf')
}

# detach all packages loaded. This is useful when we want to test any piece of code and want to start clean.
DetachAll <- function(){
  invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
}

# ReadVectorFromConfig function reads vectors from the config file. It gracefully handles when empty the vectors are empty in the 
# file, which is required in some cases.
ReadVectorFromConfig <- function(config, vector_name){
  # warning('Consider using the updated and more general function ReadVectorFromYAML!')
  if(is.null(config$filters[[vector_name]])){
    return (character(0))
  }else{
    return (unlist(strsplit(config$filters[[vector_name]], split = ",")))
  }
}

# function to add the expression values (default assay) of a vector of genes to the meta.data of a seurat obj
AddExpressionToMetadata <- function(data_fn, genes_of_interest){
  # simply get the expression_df using FetchData function
  expression_df <- FetchData(data_fn, vars = genes_of_interest)
  # cbind the expression values of the genes to meta.data
  data_fn@meta.data <- cbind(data_fn@meta.data, expression_df)
  # return the seurat object
  return (data_fn)
}

# function to get a dataframe with marker genes' information for all the identities in a seurat 
# object (colbinded df, and 1 vs all comparison)
TopNMarkersDF <- function(data_fn, ident_col = 'seurat_clusters', top_n_markers = 20){
  # ensure that the Idents is right of data_fn
  Idents(data_fn) <- ident_col
  
  # create a dataframe in which we will store the information for markers for all the identities of data
  df <- data.frame(row.names = seq(from = 1, to = top_n_markers, by = 1))
  # loop through each unique identity value, and get the top n markers for it, and append to our df
  for (identity in unique(data_fn@meta.data[,ident_col])){
    # get the markers dataframe
    markers <- FindMarkers(data_fn, ident.1 = identity)
    # make another column which will store the genes (which are currently in index), and reset the index
    markers[,'genes'] <- rownames(markers)
    rownames(markers) <- NULL
    # rename the columns in markers df to have the current identity suffixed
    colnames(markers) <- paste0(colnames(markers), '_', identity)
    # append to our final df
    df <- cbind(df, markers[1:top_n_markers, ])
  }
  # return the final df
  return (df)
}

# Subset the cells in a seurat object to a particular percentage. Usually useful for performing 
# some tests on big seurat objects.
SubsetProportionCells <- function(data, percentage){
  # subset cells for faster processing
  subset_cells <- sample(Cells(data), percentage*ncol(data), replace = FALSE)
  data <- subset(data, cells = subset_cells) 
  return (data)
}

# remove specific genes from the seurat object itself. Doing this after SCTransform has been done is a bit 
# tricky because of seurat's intricacies. Hence, a dedicated function for it.
RemoveSpecificGenesPostSCT <- function(data, features_to_remove){
  
  # get the features to keep by subtracting the features_to_remove from features of data
  features_to_keep <- setdiff(rownames(data), features_to_remove)
  
  # complicated way to subset seurat.obj after it has SCT assay
  tmp_SCT_features_attributes <- slot(object = data[['SCT']], name = "SCTModel.list")[[1]]@feature.attributes 
  tmp_SCT_features_attributes <- tmp_SCT_features_attributes[features_to_keep, ]
  slot(object = data[['SCT']], name = "SCTModel.list")[[1]]@feature.attributes <- tmp_SCT_features_attributes
  data <- subset(data, features = features_to_keep)
  return (data)
}

# function to take a df, limit to the rows identifying a particular sample, and compute the proportions of 
# celltypes, with option to limit to specific celltypes only. It returns a df.
GetCelltypeProportions <- function(df_fn, sample_name_vec, identifier_col, celltype_col, limit_to_celltypes = NA){
  
  # to make the below code more readable, will just rename the appropriate cols in df with standard names
  df_fn$cell_type <- df_fn[,celltype_col]
  df_fn$sample_identifier <- df_fn[,identifier_col]
  
  # first subset the df to just the rows with sample_name as identifier
  df_fn <- df_fn %>%
    select(c(sample_identifier, cell_type)) %>%
    filter(sample_identifier %in% sample_name_vec)
  
  # if we also limit to certain celltypes then we need to filter the celltypes too
  if (!is.na(limit_to_celltypes[1])){
    df_fn <- df_fn %>% filter(cell_type %in% limit_to_celltypes)
  }
  
  # now get the percentage values (composition) of each celltype
  df_fn <- df_fn %>%
    group_by(cell_type) %>%
    summarise(counts = n()) %>%
    mutate(percentage = counts/sum(counts))
  
  return (df_fn)    
}

# gracefully read ref_projection (depending on if it is the actual object or the path to it) 
ReadRefProjection <- function(ref_projection){
  # check if ref_projection supplied is the actual seurat object or just the string and 
  # accordingly load seurat object for label transfer 
  if (inherits(ref_projection, 'character')){ # the class is ("glue" "character") when we do glue, hence, using inherits function is more robust
    seurat_object <- qread(ref_projection)
  }else if(class(ref_projection) == 'Seurat'){
    seurat_object <- ref_projection
  }else{
    warning('supplied ref_projection is neither a string (location of the object) nor the actual object!')
  }
  return (seurat_object)
}

# Mean center the entries of a df, assuming all its cols are numerical. If center_by = 'col', then from each col, 
# subtract its mean. If center_by = 'row', then from each row, subtract its mean.
MeanCenterDF <- function(df, center_by){
  
  if (center_by == 'row'){
    df_means <- matrix(rowMeans(df), nrow = nrow(df), ncol = ncol(df), byrow = FALSE)
    mean_centered_df <- df - df_means
  }else if(center_by == 'col'){
    df_means <- matrix(colMeans(df), nrow = nrow(df), ncol = ncol(df), byrow = TRUE)
    mean_centered_df <- df - df_means
  }else{
    stop("Incorrect value for center_by argument specified. It should be one of 'row' or 'col'!")
  }
  return (mean_centered_df)
}

# Neaten the predictions df, which is obtained after calling the TransferData function from seurat. This function removes 
# the prediction.score.max and predicted.id cols and renames the colnames according to Filbin lab standard.
NeatenPredictionsDF <- function(predictions_df){
  # drop unnecessary cols
  predictions_df$prediction.score.max <- NULL 
  predictions_df$predicted.id = NULL          
  # rename the columns of predictions_df to not have prediction.score at their start
  colnames(predictions_df) <- gsub('prediction.score.', '', colnames(predictions_df))
  colnames(predictions_df) <- gsub('\\.', '-', colnames(predictions_df))
  return (predictions_df)
}

# A more general version of the AddClusterHomogeneityInfo function. It allows left joining on a seurat object. 
# Have to do that in a special way because seurat has a bug in doing that in a standard way.
# NOTE: Ensure that in the df, the cell_ids are stored in a col named 'cell_id'
LeftJoinOnData <- function(data, df){
  
  # Need following complicated mechanism of adding to metadata because otherwise the plots dont have any data due to a bug in seurat
  # add cell_id col in data
  data$cell_id <- rownames(data@meta.data)
  # make new meta.data var storing the information from data
  meta.data <- data@meta.data[,'cell_id']
  # do a left join to get the appropriate information added to each cell
  meta.data <- data@meta.data %>%
    left_join(df, by = 'cell_id')
  # finally, add the meta.data variable's information to data
  data <- AddMetaData(data, meta.data)
  
  return (data)
}

# min-max scale a vector and returned the scaled vector
MinMaxScaleVector <- function(vec, target_min = 0, target_max = 1){
  # vec = c(1,2,34,5)
  vector_min = min(vec)
  vector_max = max(vec)
  scaled_vec <- (vec-vector_min)*(target_max-target_min)/(vector_max-vector_min) + target_min
  return (scaled_vec)
}

# min-max scale a df (ensure that all values in the df are numeric) and return the df with
# values in the range (min, max)
MinMaxScaleDF <- function(df, target_min = 0, target_max = 1){
  # target_min = -100
  # target_max = 100
  df_min <- min(df)
  df_max <- max(df)
  scaled_df <-  target_min + (df-df_min)*(target_max-target_min)/(df_max-df_min)
  return (scaled_df)
}

# get the bounding box coordinates (xmin, xmax, ymin, ymax) vector from seurat object (xenium). 
# NOTE: Just like the standard I have assumed in all the projects, the x axis is horizontal and y axis is vertical.
GetBoundingBox <- function(data){
  # get the bounding box
  bounding_box <- data@images$fov$centroids@bbox
  
  # get the actual coordinates from the bounding_box. Also reverse x and y, because they are in Image format in 
  # data (x is vertical and y is horizontal). Hence, by reversing we get it in the form we want: x being horizontal 
  # and y vertical.
  ymin <- bounding_box['x', 'min']
  ymax <- bounding_box['x', 'max']
  xmin <- bounding_box['y', 'min']
  xmax <- bounding_box['y', 'max']
  return (c(xmin, xmax, ymin, ymax))
}

# return appropriate width and height for imagedimplot, given a sample, and area we decide the plots should take
# also have some computation for computing the appropriate point size for imagedimplot.
# Returns: plot_width, plot_height, point_size for imagedimplot
ComputeImageDimPlotWidthHeightSize <- function(data, area = 100){
  message(blue('Ensure that the ImageDimPlot is made without Legend. As that is how the ratio calibration was done for imagedimplot.'))
  
  # first get the bounding box for the object
  bbox_coors <- GetBoundingBox(data)
  x_span <- bbox_coors[2] - bbox_coors[1]
  y_span <- bbox_coors[4] - bbox_coors[3]
  # get the appropriate width and height of the plot if the area was as supplied
  width = SplitNumberInAspectRatio(area, x_span, y_span)[1]
  height = SplitNumberInAspectRatio(area, x_span, y_span)[2]
  
  # now, compute the appropriate dot size
  total_area <- ComputeTotalAreaOfCells(data)
  # Use the empirical equation I have found to work well for different density vs size argument in ImageDimPlot to get the 
  # appropriate point size. NOTE: below equation was calibrated for area=100 (default value of argument of this fn). 
  # For other area values, the equation may look different. Plot is saved on my desmos.
  pt_size = -1.3*1e-7*total_area + 1.99
  
  # also, we clip pt_size to a minimum of 0.4. 
  if (pt_size < 0.4) pt_size = 0.4
  
  return (c(width, height, pt_size))
}

# get the summed area of all the cells in a sample
ComputeTotalAreaOfCells <- function(data){
  list_of_areas <- lapply(data@images$fov$segmentations@polygons, function(x){x@area})
  total_area_of_cells <- sum(unlist(list_of_areas))
  return(total_area_of_cells)
}

# get the number of cells by area of a xenium sample. Originally made to compute appropriate dot size in imagedimplot
ComputeDensityOfSample <- function(data){
  # first, get the bounding box
  bbox_coors <- GetBoundingBox(data)
  x_span <- bbox_coors[2] - bbox_coors[1]
  y_span <- bbox_coors[4] - bbox_coors[3]
  
  # first get the area of all cells
  # area <- ComputeTotalAreaOfCells(data)
  
  # then, get number of cells
  num_cells <- dim(data)[2]
  
  density = num_cells/(x_span * y_span)
  return (density)
}

# get the sizes of the xenium folders of the samples listed in the provided metadata. 
# returns a df with cols: SampleName, size
# rawdatapath_prefix is the string to be added before the values in RawDataPath col in metadata, to get the absolute addresses of the xenium folders
# NOTE: ensure that the sample names are stored in a col named SampleName in metadata, and location in col named RawDataPath
GetSizesOfXeniumFolders <- function(metadata, rawdatapath_prefix){
  
  # rawdatapath_prefix <- '/n/scratch/users/s/shk490/ependymoma/xenium/data/raw_data'
  
  # get the vector of the individual locations of all the SampleNames
  individual_locations <- paste(rawdatapath_prefix, metadata$RawDataPath, sep = '/')
  
  sizes = c()
  pb <- progress_bar$new(format = "Obtaining sizes of xenium samples on disk [:bar] :percent", total = length(individual_locations), clear = FALSE)
  for (location in individual_locations){
    # get the whole du command output
    whole_du_output <- system(glue('du -s {location}'), intern = TRUE) # we are not having -h flag because we want all the sizes in same unit and later, we will divide by appropriate number to get them all in MBs
    # get the size part from it
    size_character <- str_split_1(whole_du_output, pattern = '\t')[1]
    # make it into numeric, and convert into MBs for more readability
    size_numeric <- round(as.numeric(size_character)/1024)
    # append
    sizes <- c(sizes, size_numeric)
    pb$tick()
  }
  # add the sizes as new col to metadata and return just the appropriate subset of metadata
  metadata$sizes <- sizes
  df <- metadata %>% select(SampleName, sizes)
  return (df)
}

# save the heatmap in specified directory with specified width and height. This function is for the plots which are made 
# using ComplexHeatmap, as can't use ggsave on them. This is the ggsave analogue for such plots. It should be used after 
# calling the Heatmap() function.
SaveHeatmap <- function(plot, out_dir, width, height){
  # create the pdf file in which we will add the plot information
  pdf(out_dir, width, height)
  
  # It works even with below draw command commented until we are using rstudio, because it calls the show() command when we make the Heatmap() constructor. So that needs to be kept in mind. I am commenting it for removing the necessity of supplying the plot object.
  # draw the plot on the device which we created using the pdf command
  draw(plot)
  
  # turn off the graphic device
  dev.off()
}

# cluster nmf factors, which is a common task we have in our pipeline. nmf_basis can atleast be a matrix, whose 
# columns we want to compute correlation between
# I haven't looked into this function in detail, just copied and pasted from EPN helper functions.
clusterNmfFactors <- function(nmf_basis, cor_method="pearson"){
  nmf_dist = 1-cor(nmf_basis, method=cor_method)
  hc_nmf = hclust(as.dist(nmf_dist), method="ward.D2")
  result = list()
  result[["cor_coef"]] = 1-nmf_dist
  result[["dist"]] = nmf_dist
  result[["hc_obj"]] = hc_nmf
  return(result)
}

# make the hires imagedimplot of a xenium seurat object (so that the default boundary is just changed in the function 
# scope and outside it remains non-changed)
MakeHiresImageDimPlot <- function(data, group_by_col, colors, axes = FALSE){
  DefaultBoundary(data[['fov']]) <- 'segmentation'
  plot <- ImageDimPlot(data, group.by = group_by_col, dark.background = F, cols = colors, border.size = 0, axes = axes)
  return(plot)
}

# make the hires imagefeatureplot of a xenium seurat object (so that the default boundary is just changed in the function 
# scope and outside it remains non-changed)
MakeHiresImageFeaturePlot <- function(data, plot_raw_counts = T, ...){
  if (plot_raw_counts != T) warning('Raw counts are not being plotted, but counts from DefaultAssay. It has been found that when SCT values are being plotted, values below 0 are not shown!')
  
  if (plot_raw_counts == T){
    DefaultAssay(data) <- 'Xenium'
  }
  
  DefaultBoundary(data[['fov']]) <- 'segmentation'
  plot <- ImageFeaturePlot(data, dark.background = F, border.size = 0, ...)
  return(plot)
}

# make a rasterized imagefeatureplot, as there is not inbuilt seurat's feature of making the plot rasterized
MakeRasterizedImageFeaturePlot <- function(data, features, axes = FALSE, dpi = 200, size = 0.3){
  library(ggrastr) 
  plot <- ImageFeaturePlot(data, features = features, dark.background = F, border.size = 0, axes = axes, size = size)
  plot <- rasterize(plot, layers = 'Point', dpi = dpi)
  return (plot)
}

# make a rasterized imagedimplot, as there is not inbuilt seurat's feature of making the plot rasterized
MakeRasterizedImageDimPlot <- function(data, group.by = 'cell_type', cols = NULL, axes = FALSE, dpi = 200, size = 0.3){
  library(ggrastr) 
  plot <- ImageDimPlot(data, group.by = group.by, cols = cols, dark.background = F, border.size = 0, axes = axes, size = size)
  plot <- rasterize(plot, layers = 'Point', dpi = dpi)
  return (plot)
}

# wrapper around Seurat's imagefeatureplot, to show the counts in imagefeatureplot instead of normalized 
# values. Showing normalized (sct) values also leads to negative values not being shown
ImageFeaturePlotCounts <- function(data, ...){
  DefaultAssay(data) <- 'Xenium'
  plot <- ImageFeaturePlot(data, ...)
  return (plot)
}

# run ggsave('~/plot.pdf')
sp <- function(width = NULL, height = NULL){
  if (is.null(width)){
    width = dev.size()[1]
  }
  if (is.null(height)){
    height = dev.size()[2]
  }
  ggsave('~/plot.pdf', width = width, height = height)
}

# add two new cols to object: nCount_Xenium_log, nFeature_Xenium_log, which help in visualizing the 
# areas with low quality cells in samples. This function also shows the violinplot on the available device, and
# calling ggsave just after this function will save the violinplots.
AddLoggedCountInformation <- function(data, count_colname = 'nCount_Xenium_log', feature_colname = 'nFeature_Xenium_log'){
  
  data@meta.data[, count_colname] <- log1p(data@meta.data$nCount_Xenium)
  data@meta.data[, feature_colname] <- log1p(data@meta.data$nFeature_Xenium)
  # also make the violin plots showing distribution of the values to get idea of quality
  show(VlnPlot(data, features = c(count_colname, feature_colname), ncol = 2, pt.size = 0, group.by = "orig.ident"))
  return (data)
}

# copy the summary htmls from the information in a metadata file (and rename the metadata files 
# appropriately), and save in specified folder
# data_dir_prefix specifies what is to be prepended to the values in data_dir_col to get absolute paths to xenium folders of samples
CopyAndRenameHTMLs <- function(metadata, sample_name_col, data_dir_prefix, data_dir_col, out_dir){
  
  # ensure that out_dir exists
  if (!dir.exists(out_dir)){
    stop(glue('out_dir {out_dir} doesnt exist. Make it manually first. This is to ensure user is paying attn.'))
  }
  
  # loop through each entry in metadata, and copy the corresponding html file
  for (i in seq(1, nrow(metadata))){
    # get current value of sample_name
    sample_name = metadata[i, sample_name_col]
    
    # get current location of html_file
    html_file <- glue('{data_dir_prefix}/{metadata[i, data_dir_col]}/analysis_summary.html')
    
    # now copy the html dir to the specified location
    system(glue('cp {html_file} {out_dir}/{sample_name}.html'))
  }
  return
}

# on an object, which has been obtained by combining multiple samples, perform sctransform per sample, and then merge and
# harmony and other preprocessing steps. This is generally a useful step applied on single cell reference objects.
# Reference: https://github.com/satijalab/seurat/issues/4896#issuecomment-894360415
# source_colname is the column in sc identifying which sample each cell came from
# even if object has SCT precomputed or whatever, this function by default extracts counts matrix from RNA assay's counts slot.
# features is vec of genes which we will keep in the normalized object if they are variable in atleast one of the subsets (acc to source_colname)
# of sc (generally the genes available in xenium data for which this object will be used as reference). If NULL, then default behaviour of sctransform is used
RunSCTransformPerSampleThenPreprocess <- function(sc, source_colname, residual_features = NULL, run_harmony = TRUE){
  message(blue('Splitting object and running SCTransform on each sub-object individually...'))
  sc.list <- SplitObject(sc, split.by = source_colname)
  sc.list <- lapply(X = sc.list, 
                                 FUN = SCTransform, 
                                 method = "glmGamPoi", 
                                 return.only.var.genes = FALSE,
                                 residual.features = residual_features,
                                  )
  # this basically yields the union of all the features which were found to be variable in atleast one sc.list entry
  var.features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 3000) 
  sc.sct <- merge(x = sc.list[[1]], y = sc.list[2:length(sc.list)], merge.data=TRUE) # here, in sc[['SCT']]@scale.data (which is what RunPCA takes information from), merge will yield just the intersection of features in the list. Hence, we need to recompute scaling below.
  VariableFeatures(sc.sct) <- var.features
  # we do rescaling if we supply residual features only (because only then there are very few features found to be highly variable in all samples, and we need more which we get by rescaling)
  if (!is.null(residual_features)){
    sc.sct <- ScaleData(sc.sct) # now, we will have all the genes in scale.data
  }
  
  message(blue('Running PCA, harmony if run_harmony=TRUE, UMAP, Clustering...'))
  if (run_harmony == TRUE){
    sc.sct <- RunPCA(sc.sct, verbose = TRUE)
    sc.sct <- harmony::RunHarmony(sc.sct, assay.use="SCT", group.by.vars = source_colname)
    sc.sct <- RunUMAP(sc.sct, reduction = "harmony", dims = 1:30)
    sc.sct <- FindNeighbors(sc.sct, reduction = "harmony", dims = 1:30) %>% FindClusters()
  } else{
    sc.sct <- RunPCA(sc.sct, verbose = TRUE)
    sc.sct <- RunUMAP(sc.sct, reduction = "pca", dims = 1:30)
    sc.sct <- FindNeighbors(sc.sct, reduction = "pca", dims = 1:30) %>% FindClusters()
  }
  
  message(blue('Displaying UMAP showing integration quality...'))
  show(DimPlot(sc.sct, group.by = source_colname))
  return (sc.sct)
}

# version with scvi integration instead of harmony of func RunSCTransformPerSampleThenHarmony.
# NOTE: this function should be run only when there was no RETICULATE_PYTHON env variable supplied 
# (in the .Renviron file for example), which this Rstudio session was started.
RunSCTransformPerSampleThenSCVI <- function(sc, source_colname, residual_features = NULL){
  message(blue('Splitting object and running SCTransform on each sub-object individually...'))
  sc.list <- SplitObject(sc, split.by = source_colname)
  sc.list <- lapply(X = sc.list,
                    FUN = SCTransform,
                    method = "glmGamPoi",
                    return.only.var.genes = FALSE,
                    residual.features = residual_features,
  )
  # this basically yields the union of all the features which were found to be variable in atleast one sc.list entry
  var.features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 3000)
  sc.sct <- merge(x = sc.list[[1]], y = sc.list[2:length(sc.list)], merge.data=TRUE) # here, in sc[['SCT']]@scale.data (which is what RunPCA takes information from), merge will yield just the intersection of features in the list. Hence, we need to recompute scaling below.
  VariableFeatures(sc.sct) <- var.features
  # we do rescaling if we supply residual features only (because only then there are very few features found to be highly variable in all samples, and we need more which we get by rescaling)
  if (!is.null(residual_features)){
    sc.sct <- ScaleData(sc.sct) # now, we will have all the genes in scale.data
  }
  library(sceasy)
  adata <- convertFormat(sc.sct, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
  anndata::write_h5ad(adata, "/n/data1/dfci/pedonc/filbin/lab/users/shk490/g34rv/analysis/1_preparation/Chen_2020/shashank_objects/scvi_test.h5ad")
  
  # So what I am thinking is that I can just run sctransform sample-wise on the object here, convert to adata,
  # and then save it. Then, in a python job, which has gpu access, I will run scvi on that object, and 
  # visualize the UMAP.
  
  # # run setup_anndata
  # scvi$model$SCVI$setup_anndata(adata)
  # 
  # # create the model
  # model = scvi$model$SCVI(adata)
  # 
  # # train the model
  # model$train()
  # 
  # # to specify the number of epochs when training:
  # # model$train(max_epochs = as.integer(400))
  # 
  # 
  # 
  # message(blue('Running PCA, harmony, UMAP, Clustering...'))
  # sc.sct <- RunPCA(sc.sct, verbose = TRUE)
  # 
  # qsave(sc.sct, '~/sc_sct.qs')
  # 
  # sc.sct_integrated <- IntegrateLayers(temp, method = scVIIntegration, orig.reduction = "pca",
  #                 new.reduction = 'reduc',
  #                 conda_env = '/home/shk490/.conda/envs/scvi',
  #                 num_workers = 2,
  #                 verbose = T)
  # /home/shk490/.conda/envs/scvi
  # 
  # # reticulate use appropriate conda env
  # reticulate::use_python('/home/cao385/envs/scvi/bin/python')
  # use_condaenv()
  # use_condaenv("/home/cao385/envs/scvi", conda = '', required = TRUE)
  # 
  # scvi <- import("scvi", convert = FALSE)
  # 
  # 
  # sc.sct <- harmony::RunHarmony(sc.sct, assay.use="SCT", group.by.vars = source_colname)
  # sc.sct <- RunUMAP(sc.sct, reduction = "harmony", dims = 1:30)
  # sc.sct <- FindNeighbors(sc.sct, reduction = "harmony", dims = 1:30) %>% FindClusters()
  # message(blue('Displaying UMAP showing integration quality...'))
  # p1 <- DimPlot(sc.sct, group.by = 'ident', label = T) + NoLegend()
  # p2 <- DimPlot(sc.sct, group.by = source_colname)
  # show(p1+p2)
  # 
  # return (sc.sct)
}

# run basic spatial data preprocessing steps. Original use case is to be run just after a spatial data object is 
# subsetted for some reason (like spatial subset or malignant subset).
# NOTE: this function is supposed to be used with a xenium data's object (because it uses all the genes in doing pca), 
# and might not work if a single cell object (with thousands of genes) is supplied
RunBasicSCTAndPreprocessSpatial <- function(data){
  # run sctransform
  data <- SCTransform(data, assay = 'Xenium', vars.to.regress = c('nFeature_Xenium', 'nCount_Xenium'))
  # perform pca with 50 npcs
  data <- RunPCA(data, npcs = 50, features = rownames(data), verbose = FALSE)
  # run umap with numpcs coming from optimizePCA func
  data <- RunUMAP(data, dims = 1:optimizePCA(data, 0.8), verbose = FALSE)
  # run find neighbors with numpcs coming from optimizePCA func
  data <- FindNeighbors(data, reduction = "pca", dims = 1:optimizePCA(data, 0.8), verbose = FALSE)
  # find clusters using fixed resolution of 0.8
  data <- FindClusters(data, resolution = 0.8, verbose = FALSE) # the resolution here is not important for the use case of subsetting to mal cells during annotation and annotating them
  # show the new embedding
  show(DimPlot(data))
  return (data)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Same as PlotCelltypeCompositionsOld, just takes in cellid_dir instead of seurat_objects_vec
# make the celltype composition barplot with all the samples specified. Currently made at the level of sample names only. 
# If later need to include being able to also make sample-id averaged plot, will have to make modifications.
# normal celltypes are used to determine which celltypes will be plotted below the x axis.
PlotCelltypeCompositions <- function(sample_names, cellid_dir, colors, normal_celltypes = c('Myeloid', 'Astrocyte', 'Neurons', 'Oligodendrocytes', 'T-cells', 'Endothelial'),
                                     average_by_sample_id = FALSE, sample_ids = NULL){
  
  # ensure that if average_by_sample_id is TRUE, sample_ids is supplied 
  if (average_by_sample_id == TRUE & is.null(sample_ids)) stop('If average_by_sample_id is TRUE, sample_ids should be supplied!')
  
  # get the metaprogram proportions of celltypes in each sample (not averaging by sampleid)
  mp_props <- GetMetaprogramProportions(sample_names, cellid_dir, average_by_sample_id, sample_ids)
  
  # make the proportion values of the normal celltypes to be negative
  mp_props <- mp_props %>% mutate(proportions_normal_mal = ifelse(group %in% normal_celltypes, -proportions, proportions))
  
  # now make the plot
  plot <- ggplot(mp_props, aes(x = Identifier, y = proportions_normal_mal, fill = group)) + 
    geom_col(color = 'black', size = 0.3) + 
    scale_fill_manual(values = colors) +
    theme_classic() + 
    ylab('Celltype Proportions') + 
    xlab('Sample Name') + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  return (plot)
}

# PlotCelltypeCompositions function takes in a list of annotated seurat objects, an identifier for each of them, 
# indicating their names for the classification we are interested in, and the title for that classification, and the 
# celltype col. It then makes a barplot showing the celltype compositions in all the seurat objects. 
PlotCelltypeCompositionsOld <- function(seurat_objects_vec, identifiers_vec, identifier_colname, colors, celltype_colname = 'cell_type', return_fig = FALSE, output_file = NULL, width = 5, height = 5){
  
  # some checks about arguments
  if (return_fig == FALSE & is.null(output_file)) stop('Either return_fig should be true or output_file should be supplied!')
  
  # create an empty df, where we will store the meta.data for our objects
  df <- data.frame()
  # add the appropriate cols to all the seurat objects and gradually append to the df
  for (i in seq_along(seurat_objects_vec)){
    seurat_object <- seurat_objects_vec[[i]]
    seurat_object@meta.data[,identifier_colname] <- identifiers_vec[i]
    df <- rbind(df, seurat_object@meta.data)
  }
  
  # we dont need to do manual grouping (to get counts of each celltype in each sample). stat = 'count' automatically does the grouping for us.
  plot <- ggplot(df, aes(fill = cell_type, x = .data[[identifier_colname]])) + 
    geom_bar(position = 'fill', stat = 'count') + 
    scale_fill_manual(values = colors) +
    theme_minimal() + 
    ggtitle(glue('Celltype compositions for different {identifier_colname}s'))
  
  # save or return the plot
  if (return_fig == TRUE){
    return (plot)
  }else{
    ggsave(output_file, width = width, height = height)
  }
}

# make splitted UMAP plot of the single cell mapping scores for all the different celltypes. It assumes that the 
# annotation scores are in data@meta.data, and the columns ending in '_sc' store the annotation scores for a 
# particular celltype. Also makes the dimplot of celltype predictions in the specified col.
PlotScMappingAnnotationScores <- function(data, colors, annotation_col_to_plot = 'predicted_label_snRNAseq'){
  # get the columns which have scores to be plotted 
  sc_annotation_score_colnames <- colnames(data@meta.data)[grepl('_sc$', colnames(data@meta.data))]
  # make the featureplot
  p1 <- FeaturePlot(data, features = sc_annotation_score_colnames)
  p2 <- DimPlot(data, group.by = annotation_col_to_plot, cols = colors)
  # return the plot
  return (p1+p2)
}

# make splitted UMAP plot of the spatial panel annotation scores for all the different celltypes. It assumes that the 
# annotation scores are in data@meta.data, and the columns ending in '_UCell' store the annotation scores for a 
# particular celltype. Also makes the dimplot of celltype predictions in the specified col.
PlotSpatialPanelAnnotationScores <- function(data, colors, annotation_col_to_plot = 'predicted_label_UCell'){
  # get the columns which have scores to be plotted 
  spatial_annotation_score_colnames <- colnames(data@meta.data)[grepl('_UCell$', colnames(data@meta.data))]
  # the col 'predicted_label_UCell' also ends in '_UCell', but we don't want it. Hence, remove it from the vector
  spatial_annotation_score_colnames <- spatial_annotation_score_colnames[-which(spatial_annotation_score_colnames == 'predicted_label_UCell')]
  # make the featureplot
  p1 <- FeaturePlot(data, features = spatial_annotation_score_colnames)
  p2 <- DimPlot(data, group.by = annotation_col_to_plot, cols = colors)
  # return the plot
  return (p1+p2)
}

# make barplot visualizing proportions of celltypes in each of the categories identified by a particular col of data
# column is a column in data@meta.data with categorical values
# NOTE: celltype annotations should be stored in 'cell_type' col
CelltypePropsByCol <- function(data, column, colors, order_metaprograms = NULL){
  
  # if order_metaprograms is NULL, then just sort in ascending order
  if (is.null(order_metaprograms)) order_metaprograms <- sort(unique(as.character(df$cell_type)))
  
  df <-  data@meta.data
  df$cell_type <- factor(df$cell_type, levels = order_metaprograms)
  
  # for easy plotting, just make another column in data, with same information as the supplied column variable
  data$column <- data@meta.data[, column]
  
  # Plot niche frequency
  plot <- ggplot(df, aes(column, fill = factor(cell_type, levels = order_metaprograms))) +
    scale_fill_manual(values = colors) + 
    geom_bar(position = "fill", color = "black") +
    guides(fill = guide_legend(title = "Cell type")) +
    labs(y = 'Proportion, %', x = '') + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1, colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 12)) 
  
  return (plot)
}

# add a proportions col in a df with three cols: one col according to which we group rows (like cluster), another 
# column whose proportions (like celltype) in each group need to be computed, and another column storing counts, 
# using which we compute the proportions. 
# grouping_col, class_col, and counts_col are names of the cols in df
GroupAndAddProportionsCol <- function(df, grouping_col, class_col, counts_col){
  # df <- pie_chart_df
  # grouping_col = 'cluster'
  # class_col = 'group'
  # counts_col = 'celltype_counts_in_cluster'
  
  # store in a new df, the total counts of each group (unique value in grouping_col)
  df_total_counts_group <- df %>% group_by(!!sym(grouping_col)) %>% summarise(total_counts_group = sum(!!sym(counts_col)))
  
  # perform a left join so that we have information for total_counts_group too in the df
  df <- df %>% left_join(df_total_counts_group, by = grouping_col)
  
  # now compute the proportions 
  df <- df %>% mutate(proportions = !!sym(counts_col)/total_counts_group)
  
  return (df)
}

# get a df storing the mean expression of each gene in the spatial data (and same gene's mean 
# expression in sc data), for purpose of qc on the genes
GetGeneExpressionMeansForXeniumAndSc <- function(data, sc_obj){
  # make df storing information from spatial data about the gene_names, mean_counts, rank (accoding to descending order)
  xen_means <- data.frame(
    mean_counts = rowMeans(data[["Xenium"]]$counts),
    Genes = rownames(data[["Xenium"]]$counts)) %>%
    arrange(desc(mean_counts)) %>%
    mutate(Rank = 1:n())
  # make df storing information from sc data about the gene_names, mean_counts, rank (accoding to descending order)
  sc_means <- data.frame(
    mean_counts = rowMeans(sc_obj[["RNA"]]$counts),
    Genes = rownames(sc_obj[["RNA"]]$counts)) %>%
    arrange(desc(mean_counts)) %>%
    mutate(Rank = 1:n())
  
  # join mean counts per cell by gene name (and keep only xenium genes)
  merged_means <- merge(xen_means, sc_means, by.x = "Genes", by.y = "Genes", all.x = TRUE)
  return(merged_means)
}

# take a plot and return its non-legend plot and legend plots separately, so that they can be saved then
SplitPlotIntoWithAndWithoutLegend <- function(plot){
  plot_without_legend <- plot + NoLegend() + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
  legend <- as_ggplot(ggpubr::get_legend(plot)) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
  return (list(plot_without_legend, legend))
}

# create a customized legend
# labels_named_vector holds labels and their colors to be made (just like a color palette)
# legend_marks_shape indicates if the legend marks (the shape that is colored), are circles, squares, etc. Valid values = 'square' or 'circle'
# legend_marks_bg doesnt have an effect if legend_marks_shape is 'circle'
CreateCustomLegend <- function(legend_title, labels_named_vector, legend_marks_shape = 'square', legend_marks_bg = 'black'){
  
  # create a dummy df
  df <- data.frame(labels = factor(names(labels_named_vector), levels = names(labels_named_vector)))
  
  if (legend_marks_shape == 'circle'){
    plot <- ggplot(df, aes(color = labels, x = 1, y = 1)) + 
      geom_point(size = 5) + 
      scale_color_manual(values = labels_named_vector, name = legend_title) + 
      theme_classic() + 
      theme(legend.text = element_text(size = 11))
  }else if (legend_marks_shape == 'square'){
    plot <- ggplot(df, aes(fill = labels, x = 1)) + 
      geom_bar(color = legend_marks_bg) + 
      scale_fill_manual(values = labels_named_vector, name = legend_title) + 
      theme(legend.text = element_text(size = 11))
  }else{
    stop("Please specify correct value of legend_marks_shape. Can be 'circle' or 'square'")
  }
  # now get the legend of the plot
  plot_split <- SplitPlotIntoWithAndWithoutLegend(plot)
  return (plot_split[[2]])
}

# make the plot showing expression of all genes in single cell data and in spatial data, so that  
# possibly we can remove the genes which are highly different expressionwise
# Ensure that marker_genes (panel df) has genes and their source in cols: Genes, Source
MakeScatterplotGEXSpatialSC <- function(data, sample_name, ref_projection, marker_genes){
  sc_obj <- qread(ref_projection)
  genes_expr <- GetGeneExpressionMeansForXeniumAndSc(data, sc_obj)
  
  # also join the information from the panel about the customness of the genes
  genes_info <- marker_genes %>% select(Genes, Source)
  genes_expr <- genes_expr %>% left_join(genes_info, by = 'Genes')
  # finally make the plot
  library(ggpmisc)
  # get a subset of genes_expr with extreme values of mean_counts.x and mean_counts.y so that can plot their text
  genes_expr_subset <- genes_expr %>% filter(!is.na(mean_counts.y)) %>% # remove the genes present in sc but not in panel
    filter(mean_counts.x < quantile(mean_counts.x, 0.05) | mean_counts.x > quantile(mean_counts.x, 0.95) | 
             mean_counts.y < quantile(mean_counts.y, 0.05) | mean_counts.y > quantile(mean_counts.y, 0.95))
  plot <- ggplot(genes_expr, aes(x = mean_counts.x, y = mean_counts.y, color = Source)) + 
    geom_point(size = 0.5) + 
    stat_poly_eq() + 
    geom_smooth(method = 'lm', color = 'black') + 
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
    scale_x_log10(limits = c(0.0001, 1000000)) + scale_y_log10(limits = c(0.0001, 1000000)) +
    ylab('SC Reference Mean Expression') + xlab('Xenium Mean Expression') + ggtitle(glue('GEX Correlation in {sample_name}')) + 
    theme_classic() + 
    geom_text(data = genes_expr_subset, aes(label = Genes), size = 2)
  return (plot)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~COHERENCE FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# normalize the data in cell_feature_matrix. Was originally made for coherence (Sara's version), but probably can 
# have other uses too later, so keeping it.
norm_data <- function(cell_feature_matrix) {
  counts_per_cell <- rowSums(cell_feature_matrix)
  counts_greater_than_zero <- counts_per_cell[counts_per_cell > 0]
  after <- median(counts_greater_than_zero)
  counts_per_cell <- counts_per_cell / after
  mat <- cell_feature_matrix/counts_per_cell
  return(mat)
}

# get the cell_id csv files of supplied sample_names, combine them, and return proportions of celltypes.
# Uniqueness of this function: it combines celltypes from all the samples, hence will have one proportion 
# value for each celltype and not different values according to some identifier
CombineAndReturnMetaprogramProportions <- function(sample_names, cellid_dir){
  # sample_names <- c("STEPN58_Region_1", "STEPN58_Region_2")
  df <- data.frame()
  # loop through all the samples and read their cellID csv file and rbind to df
  pb <- progress_bar$new(total = length(sample_names))
  for (sample_name in sample_names){
    # read the cellID csv for current sample
    cellID <- read.csv(glue('{cellid_dir}/cell_ID_{sample_name}.csv'))
    # bind to the dataframe
    df <- rbind(df, cellID)
    pb$tick()
  }
  # get the proportions of the celltypes
  df <- df %>% group_by(group) %>% summarise(counts = n()) %>% 
    mutate(proportions = (100*counts)/sum(counts))
  return (df)
}

# load the information about the coherences and proportions (and possibly other pieces of information for each 
# celltype in each sample) of celltypes in the different samples. This df is long form, with cols like: 
# Identifier, celltype, coherence, proportion, pixel_proportion, etc.
GetCoherencesAndProportionsLongForm <- function(sample_names, cellid_dir, coherence_dir, average_by_sample_id = FALSE, sample_ids = NULL){
  # ensure that if average_by_sample_id is TRUE, sample_ids is supplied 
  if (average_by_sample_id == TRUE & is.null(sample_ids)) stop('If average_by_sample_id is TRUE, sample_ids should be supplied!')
  
  # start with metaprogram_proportions because it is already in long form when returned from the function GetMetaprogramProportions
  metaprogram_proportions <- GetMetaprogramProportions(sample_names, cellid_dir, average_by_sample_id, sample_ids)
  # now get the coherence information
  coherences <- GetCoherencesLongForm(sample_names, coherence_dir, average_by_sample_id, sample_ids)
  # just do a join between metaprogram_proportions and coherences using the Identifier col
  final_df <- metaprogram_proportions %>% left_join(coherences, by = c('group', 'Identifier'))
  return (final_df)
}

# get the coherence information from samples in long form. 
# While originally implementing, loading the coherence information from the grid version (EPN project), and averaging 
# the mps coherences to get the sample coherence
GetCoherencesLongForm <- function(sample_names, coherence_dir, average_by_sample_id, sample_ids){
  # ensure that if average_by_sample_id is TRUE, sample_ids is supplied 
  if (average_by_sample_id == TRUE & is.null(sample_ids)) stop('If average_by_sample_id is TRUE, sample_ids should be supplied!')
  
  # we will loop through all the samples and get the information for all of their celltypes in different rows and will also
  # have a col for sample_id. At the very end, according to value of average_by_sample_id, we will group and summarise.
  coherences_df <- data.frame()
  for (i in seq_along(sample_names)){
    sample_name <- sample_names[i]
    sample_id <- sample_ids[i]
    # sample_name = 'STEPN58_Region_1'
    # get the df containing information for each celltype in different rows for this sample.
    sample_coherence_df <- GetCoherenceOfSampleLongForm(sample_name, coherence_dir)
    # add the information about the current sample's sample_id in another col
    sample_coherence_df$SampleID <- sample_id
    coherences_df <- rbind(coherences_df, sample_coherence_df)
  }
  
  # now, according to if we are averaging by sample id, average the coherences_df or not
  if (average_by_sample_id == TRUE){
    # remove the information of sample names as we are not gonna deal with that resolution 
    coherences_df$SampleName <- NULL
    # group by sample_id and celltype and get the average of the values in the other cols
    coherences_df <- coherences_df %>% group_by(SampleID, group) %>% summarise(across(where(is.numeric), mean), .groups = 'drop')
    # now rename the sample_id col to identifier
    coherences_df$Identifier <- coherences_df$SampleID
    coherences_df$SampleID <- NULL
  }else{
    # just rename the sample_name col as Identifier
    coherences_df$Identifier <- coherences_df$SampleName
    coherences_df$SampleName <- NULL
  }
  
  # return the coherences_df
  return (coherences_df)
}

# get the df containing information for each celltype in different rows for a sample.
GetCoherenceOfSampleLongForm <- function(sample_name, coherence_dir){
  # NOTE: this function is designed according to the idea that the relevant coherence information of a sample will 
  # only contain numerical values and not character type information. Hence, we average in the function 
  # GetCoherencesLongForm which takes data from this function according to this assumption.
  
  # originally, when making this function, using output of the grid coherence version. But possibly will change this 
  # in the future, which would mean most of this functions contents will change.
  
  # sample_name = 'STEPN12_Region_3'
  results <- qread(glue('{coherence_dir}/results_{sample_name}.qs'))
  # just getting the metaprogram-wise coherence information and calling it df as in this df only, will add the other information for this sample's celltypes
  df <- results$coherence_score_program
  colnames(df) <- c('group', 'coherence')
  # add sample_name information
  df$SampleName <- sample_name
  # add total sample_coherence information
  df$sample_coherence <- results$coherence_score
  # add pixel counts information
  df$pixel_counts <- results$counts$n
  # add pixel_proportions_information
  df$pixel_proportion <- df$pixel_counts/sum(df$pixel_counts)
  
  return (df)
}

# get the concatenated nhood zscore df (long form) of the specified samples. If average_by_sample_id == T, 
# the zscores for the available celltypes are averaged for the samples belonging to same sample_id
GetConcatenatedNhoodZscores <- function(sample_names, nhood_dir, average_by_sample_id = FALSE, sample_ids = NULL){
  library(progress)
  
  # ensure that if average_by_sample_id is TRUE, sample_ids is supplied 
  if (average_by_sample_id == TRUE & is.null(sample_ids)) stop('If average_by_sample_id is TRUE, sample_ids should be supplied!')
  
  # we loop through all sample_names, and concatenate their nhood zscore values for each celltype pair (which 
  # are already saved in long form)
  nhood_zscores <- data.frame()
  for (i in seq_along(sample_names)){
    # i = 1
    sample_name = sample_names[i]
    nhood_zscore <- read.csv(glue('{nhood_dir}/data/zscores_{sample_name}.csv'))[, -1] # dropping first col because it has reduntant information (about row number)
    nhood_zscore$celltype_pair <- paste0(nhood_zscore$celltype1, '_', nhood_zscore$celltype2)
    nhood_zscore$SampleName <- sample_name
    # if average_by_sample_id is True, then also attach the information about sample_ids
    if (average_by_sample_id == TRUE){
      sample_id = sample_ids[i]
      nhood_zscore$SampleID <- sample_id
    }
    nhood_zscores <- rbind(nhood_zscores, nhood_zscore)
  }
  # now, if average_by_sample_id was true, we will average the values by SampleID col, else not
  if (average_by_sample_id == TRUE){
    nhood_zscores <- nhood_zscores %>% group_by(SampleID, celltype_pair) %>% summarise(zscore = mean(zscore), .groups = 'drop')
    # now, recreate the celltype1 and celltype2 cols in nhood_zscores
    nhood_zscores$celltype1 <- str_split_i(nhood_zscores$celltype_pair, '_', 1) # 1 specifies which of the two splitted strings to obtain
    nhood_zscores$celltype2 <- str_split_i(nhood_zscores$celltype_pair, '_', 2) 
  }
  return (nhood_zscores)
}

# read the nhood count values (from corresponding files saved in nhood_dir) and concatenate and return the df. Does 
# same job as GetConcatenatedNhoodZscores, but reads the count values instead of the zscore values
GetConcatenatedNhoodCounts <- function(sample_names, nhood_dir, average_by_sample_id = FALSE, sample_ids = NULL){
  library(progress)
  
  # ensure that if average_by_sample_id is TRUE, sample_ids is supplied 
  if (average_by_sample_id == TRUE & is.null(sample_ids)) stop('If average_by_sample_id is TRUE, sample_ids should be supplied!')
  
  # we loop through all sample_names, and concatenate their nhood zscore values for each celltype pair (which 
  # are already saved in long form)
  nhood_counts <- data.frame()
  for (i in seq_along(sample_names)){
    # i = 1
    sample_name = sample_names[i]
    nhood_count <- read.csv(glue('{nhood_dir}/data/nhood_counts_{sample_name}.csv'))[, -1] # dropping first col because it has reduntant information (about row number)
    nhood_count$celltype_pair <- paste0(nhood_count$celltype1, '_', nhood_count$celltype2)
    nhood_count$SampleName <- sample_name
    # if average_by_sample_id is True, then also attach the information about sample_ids
    if (average_by_sample_id == TRUE){
      sample_id = sample_ids[i]
      nhood_count$SampleID <- sample_id
    }
    nhood_counts <- rbind(nhood_counts, nhood_count)
  }
  # now, if average_by_sample_id was true, we will average the values by SampleID col, else not
  if (average_by_sample_id == TRUE){
    nhood_counts <- nhood_counts %>% group_by(SampleID, celltype_pair) %>% summarise(nhood_count = mean(nhood_count), .groups = 'drop')
    # now, recreate the celltype1 and celltype2 cols in nhood_counts
    nhood_counts$celltype1 <- str_split_i(nhood_counts$celltype_pair, '_', 1) # 1 specifies which of the two splitted strings to obtain
    nhood_counts$celltype2 <- str_split_i(nhood_counts$celltype_pair, '_', 2) 
  }
  return (nhood_counts)
}

# get the combined, long df of metaprogram proportions of the supplied list of Identifier (sample_names or sample_ids), 
# by looking at the directory where the cell_id csvs are stored. 
# average_by_sample_id indicates if we want to return the average coherences averaged by replicates. This was required 
# in EPN project. sample_ids is a vector of length same as sample_names which holds the ids for each sample.
GetMetaprogramProportions <- function(sample_names, cellid_dir, average_by_sample_id = FALSE, sample_ids = NULL){
  library(progress)
  
  # ensure that if average_by_sample_id is TRUE, sample_ids is supplied 
  if (average_by_sample_id == TRUE & is.null(sample_ids)) stop('If average_by_sample_id is TRUE, sample_ids should be supplied!')
  
  df <- data.frame()
  # if we average by sample id, the sample_names to pass to CombineAndReturnMetaprogramProportions will be different, hence if-else
  if (average_by_sample_id == TRUE){
    # find the unique sample_ids in the supplied sample_ids vector and loop through them
    unique_sample_ids <- unique(sample_ids)
    pb <- progress_bar$new(format = "Loading metaprogram proportions for all sample_ids [:bar] :percent", total = length(unique_sample_ids), clear = FALSE)
    for (sample_id in unique_sample_ids){
      # find the sample_names for current sample_id and get the combined MP proportions
      sample_names_for_current_id <- sample_names[which(sample_ids == sample_id)]
      df_current_id <- CombineAndReturnMetaprogramProportions(sample_names_for_current_id, cellid_dir)
      # df_current_id <- CombineAndReturnMetaprogramProportionsTemp(sample_names_for_current_id, cellid_dir)
      # add identifier for current sample_id
      df_current_id$Identifier <- sample_id
      # append to final df
      df <- rbind(df, df_current_id)
      pb$tick()
    }
  }else{
    # if not averaging replicates, then simply loop through the supplied sample_names and obtain their individual Proportions
    pb <- progress_bar$new(format = "Loading metaprogram proportions for all sample_names [:bar] :percent", total = length(unique(sample_names)), clear = FALSE)
    for (sample_name in sample_names){
      df_current_sample <- CombineAndReturnMetaprogramProportions(sample_name, cellid_dir)
      df_current_sample$Identifier <- sample_name
      df <- rbind(df, df_current_sample)
      pb$tick()
    }
  }
  return (df)
}

# read and concatenate the cell_id csvs of the supplied samples, and average according to sample_id if 
# average_by_sample_id == TRUE
GetConcatenatedCellIdDF <- function(sample_names, cellid_dir, average_by_sample_id = FALSE, sample_ids = NULL){
  # ensure that if average_by_sample_id is TRUE, sample_ids is supplied 
  if (average_by_sample_id == TRUE & is.null(sample_ids)) stop('If average_by_sample_id is TRUE, sample_ids should be supplied!')
  
  # loop through each sample_name, read the cell_id df, add the Identifier col appropriately (according to 
  # average_by_sample_id argument), and append to a final df
  concatenated_cell_ids <- data.frame()
  pb <- progress_bar$new(format = 'Reading and concatenating cell_ids... [:bar] :percent', total = length(sample_names))
  for (i in seq_along(sample_names)){
    sample_name <- sample_names[i]
    # read the csv
    cellid_df <- read.csv(glue('{cellid_dir}/cell_ID_{sample_name}.csv'))
    # add the Identifier col appropriately
    if (average_by_sample_id == TRUE) cellid_df$Identifier <- sample_ids[i]
    else cellid_df$Identifier <- sample_names[i]
    # append to the final df
    concatenated_cell_ids <- rbind(concatenated_cell_ids, cellid_df)
    pb$tick()
  }
  # return the concatenated cellids df, with the Identifier col added
  return (concatenated_cell_ids)
}

# new version of the function, which loads the new, cell-level coherences
GetAverageSampleCoherences <- function(sample_names, coherence_dir, average_by_sample_id = FALSE, sample_ids = NULL){
  
  # ensure that if average_by_sample_id is TRUE, sample_ids is supplied 
  if (average_by_sample_id == TRUE & is.null(sample_ids)) stop('If average_by_sample_id is TRUE, sample_ids should be supplied!')
  
  # read the coherence list of each sample and append into df with identifier col called Identifier
  average_coherence <- data.frame()
  # not working on the metaprogram level coherences right now.
  for (sample_name in sample_names){
    # sample_name = 'STEPN12_Region_3'
    results <- read.csv(glue('{coherence_dir}/{sample_name}.csv'))
    # append the whole sample's coherence into the average_coherence dataframe
    average_coherence <- rbind(average_coherence, data.frame(Identifier = sample_name, coherence = mean(results$scores)))
  }
  
  # if we want to average the coherence values by sample id (average over the replicates), then add a new column 
  # for sample_ids and groupby. 
  if (average_by_sample_id == TRUE){
    average_coherence$Identifier <- sample_ids
    average_coherence <- average_coherence %>% group_by(Identifier) %>% summarise(across(where(is.numeric), mean)) # where(is.numeric) like functions are selection_helpers and are designed to be used only with functions like across(), select(), rename()
  }
  return (average_coherence)
}

# plot linear regression plots showing coherence of all the celltypes, with coherence of a specific celltype on 
# the x axis.
CelltypeCoherenceCelltypeCoherenceLR <- function(coherence_celltype, metadata, coherence_dir, average_by_sample_id = TRUE, colors){
  # coherence_celltype = c('MES-like')
  library(reshape)
  # since we don't want discrepancies between metadata and sample_names, we are extracting sample_names from metadata only. 
  # Hence, one needs to ensure that metadata has only those entries which need to be worked on
  sample_names <- metadata$SampleName
  
  # get the df storing coherence values for each SampleName or SampleID. colnames: Identifier, coherence
  average_spatial_coherence_df <- GetAverageSampleCoherencesGridVersion(sample_names, coherence_dir, average_by_sample_id, sample_ids = metadata$SampleID, fill_na = NA)
  # rename the hyphens with underscore in colnames of average_spatial_coherence_df, because that will make the below steps much easier
  colnames(average_spatial_coherence_df) <- gsub('-', '_', colnames(average_spatial_coherence_df))
  
  # get the vector of celltypes from columns in average_spatial_coherence_df of the form 'coherence_'
  celltypes <- colnames(average_spatial_coherence_df)[grepl('coherence_', colnames(average_spatial_coherence_df))]
  celltypes <- gsub('coherence_', '', celltypes)
  
  # Also change hyphens into underscores on the supplied argument coherence_celltype
  coherence_celltype <- gsub('-', '_', coherence_celltype)
  
  # loop through each celltype, do linear regression, make the plot and store in list
  plots_list <- list()
  for (celltype in celltypes){
    # celltype <- 'MES_like'
    # perform linear regression and get results
    result <- PerformLinearRegression(average_spatial_coherence_df, col1 = glue('coherence_{coherence_celltype}'), 
                                      col2 = glue('coherence_{celltype}'), adjust_pval_by = length(celltypes))
    plots_list[[celltype]] <- PlotLinearRegression(average_spatial_coherence_df, result, color = colors[celltype])
  }
  plot <- patchwork::wrap_plots(plots_list)
  return (plot)
}

# test relation between the coherence of specific celltypes with proportion of specific celltypes. Originally 
# written for EPN project.
CelltypeCoherenceCelltypeProportionLR <- function(coherence_celltypes, proportion_celltypes, metadata, cellid_dir, coherence_dir, average_by_sample_id){
  # coherence_celltypes = c('MES-like')
  # proportion_celltypes = c('Endothelial')
  library(reshape)
  # since we don't want discrepancies between metadata and sample_names, we are extracting sample_names from metadata only. 
  # Hence, one needs to ensure that metadata has only those entries which need to be worked on
  sample_names <- metadata$SampleName
  
  # read proportion of each cell type in each tumor 
  message(green('-------------------------------------------------------'))
  message(blue('[1/3] Reading metaprogram proportions'))
  
  # get the long df which has proportions of the different metaprograms in each sample_name/sample_id. (cols = Identifier, group, counts, proportions)
  metaprogram_proportion <- GetMetaprogramProportions(sample_names = sample_names, cellid_dir, average_by_sample_id, sample_ids = metadata$SampleID)
  # limiting to just the relevant columns
  metaprogram_proportion <- metaprogram_proportion[ , c("Identifier", "group", "proportions")]
  # ensure that every combination of Identifier and group is present in the data, and if they were not present till now, their proportions will be 0
  metaprogram_proportion_complete <- metaprogram_proportion %>% 
    complete(Identifier, group, fill = list(proportions = 0))
  # cast the just created df into wide form so that we have Identifier x group shape of df.
  metaprogram_proportion_wide <- cast(metaprogram_proportion_complete, Identifier~group, value = 'proportions')
  
  message(green('-------------------------------------------------------'))
  message(blue('[2/3] Get spatial coherence scores'))
  
  # get the df storing coherence values for each SampleName or SampleID. colnames: Identifier, coherence
  average_spatial_coherence_df <- GetAverageSampleCoherencesGridVersion(sample_names, coherence_dir, average_by_sample_id, sample_ids = metadata$SampleID)
  # replace the zeros by NA, because the zeros most likely originate from samples without that celltype
  average_spatial_coherence_df[average_spatial_coherence_df == 0] <- NA
  
  message(green('-------------------------------------------------------'))
  message(blue('[3/3] Organize data for performing Linear Regression'))
  
  # combine the information about celltypes proportions into average_spatial_coherence_df
  linear_regression_df <- average_spatial_coherence_df %>% 
    left_join(metaprogram_proportion_wide, by = 'Identifier')
  # make the cols in linear_regression_df to have dots instead of spaces because as.formula function which we 
  # use below doesn't handle columns with hyphens well
  colnames(linear_regression_df) <- gsub('-', '.', colnames(linear_regression_df))
  # Also do the same on the supplied arguments coherence_celltypes and proportion_celltypes
  coherence_celltypes <- gsub('-', '.', coherence_celltypes)
  proportion_celltypes <- gsub('-', '.', proportion_celltypes)
  # loop through each combination, do linear regression, make the plot and store in list
  plots_list <- list()
  # get the number of tests we are performing. We will adjust the p-val by this number.
  adjust_pval_by = length(proportion_celltypes) * length(coherence_celltypes)
  for (proportion_celltype in proportion_celltypes){
    for (coherence_celltype in coherence_celltypes){
      # coherence_celltype <- 'Neurons'
      # proportion_celltype <- 'MES.like'
      # perform linear regression and get results
      coherence_col <- paste0('coherence_', coherence_celltype) # need this new string because the coherence information is in col named in the form: coherence_celltype
      result <- PerformLinearRegression(linear_regression_df, proportion_celltype, coherence_col, adjust_pval_by)
      plots_list[[paste0(coherence_celltype, '~', proportion_celltype)]] <- PlotLinearRegression(linear_regression_df, result, x_label = glue('Proportion of {proportion_celltype}'))
    }
  }
  plot <- patchwork::wrap_plots(plots_list)
  return (plot)
}

# perform linear regression on a df using the formula col1 ~ col2. We do this commonly, hence making a function for it.
# Return a df storing the relevant results we need. col1 is on the x axis, col2 on y by default (in the plots that are made from result of this).
# adjust_pval_by holds the number of tests being performed, so that it can multiply that number to the p-value (bonferroni correction)
PerformLinearRegression <- function(linear_regression_df, col1, col2, adjust_pval_by = 1){
  
  # corner case: when col1 == col2, then create a duplicate col2 in df so that lm doesn't throw 
  # error due to predictor and response cols being same. There is difference in creating a new col 
  # with '_duplicate' suffixed and also in creating the result df, because in that, we use the col1 
  # for both x and y instead of col1 for x and col2 for y
  if (col1 == col2){
    # create duplicate col with '_duplicate' as suffix
    linear_regression_df[, glue('{col2}_duplicate')] <- linear_regression_df[, col2]
    # change the value of variable col2
    col2 <- glue('{col2}_duplicate')
    # Dynamically create the formula. Make quoted because else as.formula() function splits string open at the hyphen 
    formula_name <- paste(col1, "~", col2)
    formula <- as.formula(formula_name)
    print(glue('Using formula \"{formula_name}\" for LR. Ensure its as intended.'))
    # Fit the linear model
    model <- lm(formula, data = linear_regression_df)
    # Extract coefficients summary
    summary_model <- summary(model)
    coefficients <- summary_model$coefficients
    # Extract R-squared and p-value
    r_squared <- summary_model$r.squared
    p_value <- coefficients[2, "Pr(>|t|)"]  # p-value of the predictor term
    # also obtain the pearson correlation because the reviewer asked for it (in EPN)
    correlation <- cor(x = linear_regression_df[, col2], y = linear_regression_df[, col1])
    # Append results to dataframe, but if col1 and col2 were same and hence we added a new col above, 
    # then don't use the added col's name for x/y in result. Because the linear_regression_df outside 
    # this function doesn't have that column and hence we face problems in plotting
    result <- data.frame(x = col1, y = col1, Estimate = coefficients[2, "Estimate"], 
                         StdError = coefficients[2, "Std. Error"], tValue = coefficients[2, "t value"], pValue = p_value, adj_pValue = p_value*adjust_pval_by, Rsquared = r_squared, correlation = correlation[1])
  }else{
    # Dynamically create the formula. Make quoted because else as.formula() function splits string open at the hyphen 
    formula_name <- paste(col1, "~", col2)
    formula <- as.formula(formula_name)
    print(glue('Using formula \"{formula_name}\" for LR. Ensure its as intended.'))
    # Fit the linear model
    model <- lm(formula, data = linear_regression_df)
    # Extract coefficients summary
    summary_model <- summary(model)
    coefficients <- summary_model$coefficients
    # Extract R-squared and p-value
    r_squared <- summary_model$r.squared
    p_value <- coefficients[2, "Pr(>|t|)"]  # p-value of the predictor term
    # also obtain the pearson correlation because the reviewer asked for it (in EPN)
    correlation <- cor(x = linear_regression_df[, col2], y = linear_regression_df[, col1])
    # Append results to dataframe
    result <- data.frame(x = col1, y = col2, Estimate = coefficients[2, "Estimate"], 
                         StdError = coefficients[2, "Std. Error"], tValue = coefficients[2, "t value"], pValue = p_value, adj_pValue = p_value*adjust_pval_by, Rsquared = r_squared, correlation = correlation[1])
  }
  return (result)
}

# use the linear regression result with cols: x, y, Estimate, StdError, tValue, pValue, MES.like to make plots and return plots as specified
# if color_col is specified, then color should be a color_palette denoting what values in the df are gonna be used for the dots coloring
# if color_col is not specified, then color should be a color value which is the color by which all the dots will be colored for that specific LR plot. If no color is specified in this case, then black is used
PlotLinearRegression <- function(linear_regression_df, result, x_label = NULL, y_label = NULL, color = NULL, color_col = NULL){
  # add new cols for relevant aesthetics in linear_regression_df
  linear_regression_df$x_axis <- linear_regression_df[[result$x]]
  linear_regression_df$y_axis <- linear_regression_df[[result$y]]
  # linear_regression_df$y_axis <- linear_regression_df[[result$x]]
  # get the x_label and y_label to use in plot according to what user has specified
  if (is.null(x_label)) x_label <- result$x
  if (is.null(y_label)) y_label <- result$y
  
  # if color is NULL, just making an ad-hoc color palette for coloring the plot
  if (is.null(color_col)){
    # if color also not specified, then just use black for dots
    if(is.null(color)) color <- 'black'
    linear_regression_df$color <- 'dot_color'
    # create a named vector holding the color
    functions_color_palette <- c(color)
    names(functions_color_palette) <- 'dot_color'
  } else{
    # if supplied a color pallette
    linear_regression_df$color <- linear_regression_df[[color_col]]
    functions_color_palette <- color
  }
  
  # make the plot
  plot <- ggplot(linear_regression_df, aes(x = x_axis, y = y_axis, color = color)) +
    labs(x = x_label, y = y_label) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black') +
    geom_point(size = 4) +
    scale_color_manual(values = functions_color_palette) + 
    theme_classic() + 
    theme(legend.position = 'none') +
    labs(subtitle = paste0('r = ', round(result$correlation, 2), ', R² = ', round(result$Rsquared, 2), ', pval: ', signif(result$pValue, 3), ', adj. pval: ', signif(result$adj_pValue, 3)))
  
  return (plot)
}

# Made because wanted to answer if presence of MES affects coherence of other things in the sample, in EPN project.
GetAverageSampleCoherencesAfterRemovingCelltype <- function(celltype_to_remove, sample_names, coherence_dir, average_by_sample_id = FALSE, sample_ids = NULL){
  # celltype_to_remove <- 'MES-like'
  # ensure that if average_by_sample_id is TRUE, sample_ids is supplied 
  if (average_by_sample_id == TRUE & is.null(sample_ids)) stop('If average_by_sample_id is TRUE, sample_ids should be supplied!')
  
  # read the coherence list of each sample and append into df with identifier col called Identifier
  average_coherence <- data.frame()
  for (sample_name in sample_names){
    # sample_name = 'STEPN01_Region_3'
    results <- qread(glue('{coherence_dir}/results_{sample_name}.qs'))
    # now get the average coherence after removing the pixels from provided celltype
    cell_types <- names(results$results_df)
    total_coherence_celltype_list <- lapply(results$results_df, sum) # list storing celltype:total_coherence
    npixels_celltype_df <- results$counts
    total_pixels_with_cells <- sum(npixels_celltype_df$n)
    num_pixels_of_celltype_of_interest <- npixels_celltype_df[['n']][npixels_celltype_df[['Metaprogram']] == celltype_to_remove]
    total_pixels_with_cells_after_removing_celltype_to_remove <- total_pixels_with_cells-num_pixels_of_celltype_of_interest
    total_coherence <- sum(as.numeric(total_coherence_celltype_list))
    total_coherence_after_removing_celltype_to_remove <- total_coherence - total_coherence_celltype_list[[celltype_to_remove]]
    new_coherence <- total_coherence_after_removing_celltype_to_remove/total_pixels_with_cells_after_removing_celltype_to_remove
    average_coherence <- rbind(average_coherence, data.frame(Identifier = sample_name, coherence = new_coherence))
  }
  
  # if we want to average the coherence values by sample id (average over the replicates), then groupby accordingly. 
  if (average_by_sample_id == TRUE){
    average_coherence$Identifier <- sample_ids
    average_coherence <- average_coherence %>% group_by(Identifier) %>% summarise(across(where(is.numeric), mean)) # where(is.numeric) like functions are selection_helpers and are designed to be used only with functions like across(), select(), rename()
  }
  return (average_coherence)
}

# make linear regression plots similar to OverallCoherenceCelltypeProportionLinearRegression, but have 
# coherence of a mp instead of its proportion on the x axis.
OverallCoherenceCelltypeCoherenceLinearRegression <- function(metadata, coherence_dir, average_by_sample_id, colors){
  library(reshape)
  # since we don't want discrepancies between metadata and sample_names, we are extracting sample_names from metadata only. 
  # Hence, one needs to ensure that metadata has only those entries which need to be worked on
  sample_names <- metadata$SampleName
  
  # get the df storing coherence values for each SampleName or SampleID. colnames: Identifier, coherence
  average_spatial_coherence_df <- GetAverageSampleCoherencesGridVersion(sample_names, coherence_dir, average_by_sample_id, sample_ids = metadata$SampleID, fill_na = NA)
  # rename the hyphens with underscore in colnames of average_spatial_coherence_df, because that will make the below steps much easier
  colnames(average_spatial_coherence_df) <- gsub('-', '_', colnames(average_spatial_coherence_df))
  # calculate scaled score (get values to 0-1 range and add as a new col)
  average_spatial_coherence_df <- average_spatial_coherence_df %>% mutate(scaled_spatial_coherence = MinMaxScaleVector(coherence))
  
  # get the vector of celltypes from columns in average_spatial_coherence_df of the form 'coherence_'
  celltypes <- colnames(average_spatial_coherence_df)[grepl('coherence_', colnames(average_spatial_coherence_df))]
  celltypes <- gsub('coherence_', '', celltypes)
  
  # loop through each celltype and perform linear regression between the coherence of that celltype in the 
  # samples and the overall scaled_spatial_coherence of the sample
  plots_list <- list()
  for (celltype in celltypes){
    # celltype = 'Myeloid'
    results <- PerformLinearRegression(average_spatial_coherence_df, glue('coherence_{celltype}'), 'scaled_spatial_coherence', adjust_pval_by = length(celltypes))
    # now make the linear regression plot using the results df
    plots_list[[celltype]] <- PlotLinearRegression(average_spatial_coherence_df, results, color = colors[celltype])
  }
  plot <- patchwork::wrap_plots(plots_list)
  return (plot)
}

# make boxplots comparing *celltype-coherence* across different parameters. Difference between this and 
# CoherenceComparisonBoxplots is that CoherenceComparisonBoxplots has overall coherence on y-axis, and this
# has for one MP. covariate is a column name from metadata, against which the comparison will be made.
MPLevelCoherenceComparisonBoxplots <- function(metadata, covariate, coherence_dir, average_by_sample_id, colors){
  # Ensure that metadata has only those entries which need to be worked on
  sample_names <- metadata$SampleName
  # Ensure that covariate has exactly 2 unique values in metadata
  if (length(unique(metadata[[covariate]])) != 2) stop('covariate col in metadata should have exactly 2 unique values!')
  
  message(green('-------------------------------------------------------'))
  message(blue('[1/3] Get spatial coherence scores'))
  
  # get the df storing coherence values for each SampleName or SampleID. colnames: Identifier, coherence
  average_spatial_coherence_df <- GetAverageSampleCoherencesGridVersion(sample_names, coherence_dir, average_by_sample_id, sample_ids = metadata$SampleID, fill_na = NA)
  # limit to the cols which have relevant information for this function
  average_spatial_coherence_df <- average_spatial_coherence_df %>% select(Identifier, starts_with('coherence_'))

  message(green('-------------------------------------------------------'))
  message(blue('[3/3] Add additional metadata for each sample'))
  
  # first add a 'Identifier' col in metadata, according to value of average_by_sample_id
  if (average_by_sample_id == TRUE) metadata$Identifier <- metadata$SampleID else metadata$Identifier <- metadata$SampleName
  # limit metadata to just the cols of interest
  metadata <- metadata %>% select(Identifier, .data[[covariate]])
  # perform left_join on metadata and the coherence df
  df <- metadata %>% left_join(average_spatial_coherence_df, by = 'Identifier')
  # just get the unique rows in df. There can be duplicates because when we have average average_by_sample_id == TRUE,
  # we have multiple rows in df for each SampleID, which have essentially same value
  df <- unique(df)
  
  # now, make the df into long form, so that we don't have the coherence of different mps in different columns
  # but in one so that we can plot boxplots separated by celltype
  df <- df %>% pivot_longer(cols = starts_with('coherence_'), names_to = 'celltype', values_to = 'coherence')
  # change the names of entries in df$celltype from something like coherence_Myeloid to Myeloid
  df$celltype <- gsub('coherence_', '', df$celltype)
  
  # compute the df storing pval information for each comparison (coherence of a celltype vs covariate)
  pval_df <- ComputePvalsForEachComparison(df, category_col = 'celltype', values_col = 'coherence', covariate_col = covariate)
  # just to make the user aware, print the information about which test was used.
  print('printing the pval_df for the user to be aware of the statistical tests used.')
  print(pval_df)
  # perform a left join on df with pval_df
  df <- df %>% left_join(pval_df, by = 'celltype')
  
  # combine the pval information in df into a new col (in col named subplot_title)
  df$subplot_title <- paste0(df$celltype, ' (p-val: ', signif(df$pval, 3), ', adj.p-val: ', signif(df$adj_pval, 3), ')')
  
  # make the boxplots for the different celltypes (y = coherence, x = covariate values)
  plot <- ggplot(df, aes(x = .data[[covariate]], y = coherence, fill = celltype)) + 
    geom_boxplot() + 
    geom_jitter(height = 0.1, width = 0.3) + 
    facet_wrap(~subplot_title) + 
    scale_fill_manual(values = colors) + 
    theme_minimal()
  
  return (plot)
}

# compute the pval corresponding to each value in a category col in a df. The pval is computed 
# for the comparison of the values (in values_col) for the different unique entries in covariate col.
# Eg- category_col = celltype, values_col = celltype_coherence, covariate_col = Source (primary/recurrence).
# In which case, the pvals will indicate how much a particular celltype's coherence varies in primary vs recurrence.
# Returns a df with cols: [category, pval, adj_pval, test_used].
# Assumes only 2 unique values in covariate_col
ComputePvalsForEachComparison <- function(df, category_col, values_col, covariate_col, test_to_use = 'wilcoxon'){
  
  # initialize the df which we will return
  pval_df <- data.frame()
  
  # get the unique values in category_col
  categories <- unique(df[[category_col]])
  # get unique covariate values. We are assuming same covariate values (2 in number) to 
  # be compared for all the categories 
  covariate_values <- unique(df[[covariate_col]])
  
  # loop through each category in categories and add their entry in pval_df
  for (category in categories){
    # category = 'Embryonic-neuronal-like'
    # get the subset of df for this category
    df_subset <- df[df[, category_col] == category, ]
    # remove the rows with na for the current df_subset. This is required for our checking mechanism for <3 datapoints to work.
    df_subset <- df_subset %>% drop_na() 
    
    # initialize the test variable, which takes value either 'none' (when not enough values) or 
    # value of test_to_use variable.
    test = test_to_use
    
    # loop through each covariate value, and see if for any value of covariate_value for current 
    # category, we have less than 3 datapoints. In that case, just have the pval as NA
    for (covariate_value in covariate_values){
      if (sum(df_subset[, covariate_col] == covariate_value) < 2){ 
        test = 'none'
        pval <- NA
        # create the new_row df which we will append to pval_df
        new_row <- data.frame(category_col = category, pval = pval, adj_pval = pval*length(categories), test_used = test)
        # append the new row to pval_df
        pval_df <- rbind(pval_df, new_row)
        break
      }
    }
    # if test is none, meaning less than 3 datapoints for some covariate, move to next category
    if (test == 'none') next
    
    # according to the supplied value of test, get the pval
    if (test == 't_test'){
      pval <- GetTTestPval(df_subset, covariate_col, values_col)
    }else if (test == 'wilcoxon'){
      pval <- GetWilcoxonPval(df_subset, covariate_col, values_col)
    }
    # create the new_row df which we will append to pval_df
    new_row <- data.frame(category_col = category, pval = pval, adj_pval = pval*length(categories), test_used = test)
    # append the new row to pval_df
    pval_df <- rbind(pval_df, new_row)
  }
  # change the names of the category_col of pval_df to what category_col variable holds
  colnames(pval_df)[1] <- category_col
  return (pval_df)
}

# get the t test p value for a comparison in a df, where the information for the categories 
# is stored in covariate_col and the actual values stored in values_col
GetTTestPval <- function(df, covariate_col, values_col){
  
  # get the 2 unique values in covariate_col of df
  covariate_1 <- unique(df[[covariate_col]])[1]
  covariate_2 <- unique(df[[covariate_col]])[2]
  
  # get the two vectors for the 2 covariate_values (holding values from values_col of df) 
  vector_1 <- df[df[[covariate_col]] == covariate_1, ][[values_col]]
  vector_2 <- df[df[[covariate_col]] == covariate_2, ][[values_col]]
  
  # perform the test and get the values. This part can be changed for different applications. 
  # Originally writing with focus on EPN project
  t_test_result <- t.test(vector_1, vector_2, alternative = 'two.sided', paired = F, var.equal = T)
  
  return (t_test_result$p.value)
}

# get the wilcox test (rank sum test) p value for a comparison in a df, where the information 
# for the categories is stored in covariate_col and the actual values stored in values_col
GetWilcoxonPval <- function(df, covariate_col, values_col){
  # get the 2 unique values in covariate_col of df
  covariate_1 <- unique(df[[covariate_col]])[1]
  covariate_2 <- unique(df[[covariate_col]])[2]
  
  # get the two vectors for the 2 covariate_values (holding values from values_col of df) 
  vector_1 <- df[df[[covariate_col]] == covariate_1, ][[values_col]]
  vector_2 <- df[df[[covariate_col]] == covariate_2, ][[values_col]]
  
  # perform the test and get the values. This part can be changed for different applications. 
  # Originally writing with focus on EPN project
  wilcox_test_result <- wilcox.test(vector_1, vector_2, alternative = 'two.sided', paired = F)
  
  return (wilcox_test_result$p.value)
}

# make the plots comparing coherence values across different parameters
# covariates is a column name from metadata, against which the comparison will be made.
CoherenceComparisonBoxplots <- function(metadata, covariate, coherence_dir, average_by_sample_id, colors){
  
  # Ensure that metadata has only those entries which need to be worked on
  sample_names <- metadata$SampleName
  
  message(green('-------------------------------------------------------'))
  message(blue('[1/3] Get spatial coherence scores'))
  
  # get the df storing coherence values for each SampleName or SampleID. colnames: Identifier, coherence
  average_spatial_coherence_df <- GetAverageSampleCoherencesGridVersion(sample_names, coherence_dir, average_by_sample_id, sample_ids = metadata$SampleID)
  # calculate scaled score (get values to 0-1 range and add as a new col)
  average_spatial_coherence_df <- average_spatial_coherence_df %>% mutate(scaled_spatial_coherence = MinMaxScaleVector(coherence))
  
  message(green('-------------------------------------------------------'))
  message(blue('[3/3] Add additional metadata for each sample'))
  
  # first add a 'Identifier' col in metadata, according to value of average_by_sample_id
  if (average_by_sample_id == TRUE) metadata$Identifier <- metadata$SampleID else metadata$Identifier <- metadata$SampleName
  # limit metadata to just the cols of interest
  metadata <- metadata %>% select(Identifier, .data[[covariate]])
  # perform left_join on metadata and the coherence df
  df <- metadata %>% left_join(average_spatial_coherence_df, by = 'Identifier')
  # just get the unique rows in df. There can be duplicates because when we have average average_by_sample_id == TRUE,
  # we have multiple rows in df for each SampleID, which have essentially same value
  df <- unique(df)
  
  # perform the shapiro test on scaled_spatial_coherence for every value of the covariate
  for (covariate_value in unique(df[[covariate]])){
    # covariate_value <- 'Primary'
    df_filtered <- df %>% filter(.data[[covariate]] == covariate_value)
    shapiro_result <- shapiro.test(df_filtered$scaled_spatial_coherence)
    print(glue('shapiro normality test p-value for covariate value \"{covariate_value}\": {shapiro_result$p.value}'))
  }
  # perform the bartlett test to check for equalness of variance of the scaled_spatial_coherence in each of covariate values
  bartlett_test_list = list()
  for (covariate_value in unique(df[[covariate]])){
    # covariate_value <- 'Primary'
    df_filtered <- df %>% filter(.data[[covariate]] == covariate_value)
    bartlett_test_list[[covariate_value]] <- df_filtered$scaled_spatial_coherence
  }
  bartlett_test_result <- bartlett.test(bartlett_test_list)
  print(glue('bartlett variance-similarity test p-value for given covariate \"{covariate}\": {bartlett_test_result$p.value}'))
  
  # now make the boxplot of interest. Depending on the p-values of shapiro and bartlett tests, method of stat_compare_means would change.
  # have explicit computation of p-value
  pval <- t.test((df %>% filter(Source == 'Primary'))$scaled_spatial_coherence,
                 (df %>% filter(Source == 'Recurrence'))$scaled_spatial_coherence,
                 alternative = "two.sided", paired = FALSE, var.equal = TRUE)$p.value
  
  plot <- ggplot(df, aes(x = .data[[covariate]], y = scaled_spatial_coherence, fill = .data[[covariate]])) + 
    geom_boxplot() + 
    geom_jitter(height = 0.1, width = 0.3) + 
    theme_minimal() + 
    ylab('Spatial Coherence Score') + 
    theme(plot.subtitle = element_text(hjust = 0.5)) + 
    labs(subtitle = glue('p-value: {signif(pval, 3)}')) + 
    scale_fill_manual(values = colors)
  
  return (plot)
}

# get df with number of cells in each sample in sample_names. Originally made because I wanted to find out if the coherence 
# has some dependence on number of cells in the sample
GetNumCells <- function(sample_names, cellid_dir){
  df <- data.frame()
  # loop through all the samples and read their cellID csv file and rbind to df
  for (sample_name in sample_names){
    # read the cellID csv for current sample
    cellID <- read.csv(glue('{cellid_dir}/cell_ID_{sample_name}.csv'))
    # bind to the dataframe
    df <- rbind(df, data.frame(Identifier = sample_name, num_cells = nrow(cellID)))
  }
  return (df)
}

# get the celltypes of all the cells near each cell (determined by threshold), and the vector storing coherence scores for each cell too.
# returns a vector of coherence scores and a list with names as the cell_indices, and entries as vectors of celltypes of cells which were near 
# to the cell with that cell_index. Originally made for a specific type of coherence computation, but possibly can have other applications in the future too.
# NOTE: Ensure that the celltypes are stored in data in col named 'cell_type'
ComputeCoherenceCellLevel <- function(data, distance_threshold, percent_mode){
  
  cell_names <- rownames(data@meta.data)
  near_cells_celltypes_list <- list() # each entry looks like cell_index:neighbors_celltypes_vector
  coherence_scores <- c() # vector in which we will store the coherence scores of all the cells in order
  # we compute distance_threshold_sq so that we don't have to compute the square roots of all distances inside the loop. Hence, faster this way as only one square computation
  distance_threshold_sq <- distance_threshold*distance_threshold
  pb <- progress_bar$new(format = "Finding neighbors for each cell [:bar] :percent", total = length(cell_names), clear = FALSE)
  for (center_cell_index in 1:length(cell_names)){
    # compute the differences of each cells x and y coordinates from x and y coordinates of center cell  
    deltas <- sweep(data@images$fov@boundaries$centroids@coords, 2, data@images$fov@boundaries$centroids@coords[center_cell_index, ], FUN = '-')
    deltas_sq <- deltas*deltas
    dist_sqs <- rowSums(deltas_sq) # this stores the squares of the distances of each cell from current cell
    
    # now, get the celltypes of the cells which are 
    celltypes_of_cells_near <- data$cell_type[dist_sqs < distance_threshold_sq]
    near_cells_celltypes_list[[center_cell_index]] <- celltypes_of_cells_near
    
    # also get the coherence score of current cell
    coherence_score <- ComputeCoherenceScoreFromNearCelltypesVec(data$cell_type[center_cell_index], celltypes_of_cells_near, percent_mode)
    coherence_scores <- c(coherence_scores, coherence_score)
    pb$tick()
  }
  return (list('coherence_scores' = coherence_scores, 'near_cells_celltypes_list' = near_cells_celltypes_list))
}

# compute the percentage of nearby cells which are of same category as celltype (intended to be celltype of center cell)
# celltype is a string with the celltype of the cell in focus. celltypes_of_cells_near is a vector storing 
# the celltypes of the nearby cells. 
# percent_mode indicates if coherence is a proportion celltype similar out of how many cells near or just the raw number
# returns just the score value (int if percent mode is FALSE, else decimal between 0 and 1)
ComputeCoherenceScoreFromNearCelltypesVec <- function(celltype, celltypes_of_cells_near, percent_mode){
  # score is the number of cells out of current cell's neighbors which are of same celltype as itself
  if (length(celltypes_of_cells_near) == 1) return (0) # in the edge case of no neighbor (identified by length of vector = 1, meaning itself), score is zero
  score <- sum(celltypes_of_cells_near == celltype) - 1 # we subtract 1 because the cell itself is always included in the vector of neighbors, but we don't want to include it.
  num_neighbors <- length(celltypes_of_cells_near) - 1 # since the cell itself is also technically not a neighbor
  if (percent_mode == TRUE){
    return (score/num_neighbors)
  }else{
    return (score)
  }
}



