# Here are the functions which deal with doing something with the actual spatial coordinates of the file (like subsetting according 
# to spatial coordinates etc.). These were originally written for the G34RV project, but are useful for other projects too.

# function to return the coordinates of the centroid of a cell, given its name and the data object
GetCellCentroid <- function(data, cell_name){
  return (data@images$fov@boundaries$centroids[cell_name]@coords)
}

# make the image feature plot with rectangle
ImageFeaturePlotWithRect <- function(data, xmin, xmax, ymin, ymax, features, axes, dark.background = F){
  if (dark.background == T){
    plot <- ImageFeaturePlot(data, axes = axes, dark.background = T, features = features) + 
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, color = 'white')
  }else{
    plot <- ImageFeaturePlot(data, axes = axes, dark.background = F, features = features) + 
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, color = 'black')
  }
  return (plot)
}

# function to take a data and show where in it is the region we selected (by the xmin, xmax, ymin and ymax coordinates).
# NOTE: Our xmin and xmax are horizontal coordinates only. 
ImageDimPlotWithRect <- function(data, xmin, xmax, ymin, ymax, colors = NULL, axes = TRUE){
  if (is.null(colors)){
    plot <- ImageDimPlot(data, axes = axes, dark.background = FALSE) + 
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, color = 'black')
  }else{
    plot <- ImageDimPlot(data, axes = axes, dark.background = FALSE, cols = colors) + 
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, color = 'black')
  }
  return (plot)
}

# same as ImageDimPlotWithSquare, but center is a list of coordinates here. Each entry of the list is a vector 
# of length 2, storing the x and y coordinates. It assumes that the colors are in the global environment.
ImageDimPlotWithSquares <- function(data, centers, size, col){
  # if a vector was supplied, convert it into list, as below code assumes that centers is a list
  if (class(centers) == 'numeric'){centers <- list(centers)}
  
  # convert the list of centers to dataframe and add relevant columns
  df <- as.data.frame(do.call(rbind, centers))
  df[, 'xmin'] <- df[,'V1'] - size
  df[, 'xmax'] <- df[,'V1'] + size
  df[, 'ymin'] <- df[,'V2'] - size
  df[, 'ymax'] <- df[,'V2'] + size
  
  # make the plot, with the information about the base plot and the added rectangles
  plot <- ImageDimPlot(data, axes = TRUE, dark.background = FALSE) +
    geom_rect(data = df, mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, color = col, inherit.aes = FALSE)
  return (plot)
}

# function to compute the distance between two coordinates (during making, the class of coordinates is assumed to be 
# "matrix", eg: [1,] 6032.142 4383.116)
DistanceBetweenCoords <- function(coord1, coord2){
  x1 <- coord1[1] 
  x2 <- coord2[1]
  y1 <- coord1[2] 
  y2 <- coord2[2]
  
  distsq <- (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)
  return (sqrt(distsq))
}

# another version of DistanceBetweenCoords which takes 4 numbers as input instead of two coords (class matrix). 
DistanceBetweenCoordsGeneral <- function(x1, x2, y1, y2){
  distsq <- (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)
  return (sqrt(distsq))
}

# version of DistanceBetweenCoordsGeneral which doesnt do square root for speed.
DistanceSquareBetweenCoords <- function(x1, x2, y1, y2){
  distsq <- (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)
  return (distsq)
}

# function to check if two cells are in vicinity of each other, according to a supplied threshold
InVicinity <- function(data, sender_cell, receiver_cell, distance_thresh = 30){
  # get the coordinates of the two cells
  sender_centroid <- GetCellCentroid(data, sender_cell)
  receiver_centroid <- GetCellCentroid(data, receiver_cell)
  
  dist <- DistanceBetweenCoords(sender_centroid, receiver_centroid)
  if (dist < distance_thresh){
    return (TRUE)
  }else{
    return (FALSE)
  }
}

# get subsets of data (small squares), around a center point. Initially intended to visualize the locations of 
# interest (size around 20-30), but can be used for more purposes by modifying the size argument.
# NOTE: Always, the x coordinates are assumed to be horizontal and y coordinates to be vertical (even though in the 
# imagedimplots x axis label is 'y', for us its x coordinate). This is done for uniformity in subsetting and plotting.
SubsetToSquare <- function(data, center, size){ # size is distance to boundary from centre
  # get the xmin, xmax, ymin, ymax values
  xmin = center[1] - size
  xmax = center[1] + size
  ymin = center[2] - size
  ymax = center[2] + size
  
  print(glue('Coordinates of rectangle for area in focus: {xmin}, {xmax}, {ymin}, {ymax}'))
  
  # crop the data
  cropped_data <- Crop(data[['fov']], x = c(xmin, xmax), y = c(ymin, ymax), coords = 'plot') 
  data_sub <- subset(data, cells = Cells(cropped_data))
  
  return (data_sub)
}

# same as SubsetToSquare, just takes xmin, xmax, ymin, ymax as inputs
# NOTE: Always, the x coordinates are assumed to be horizontal and y coordinates to be vertical (even though in the 
# imagedimplots x axis label is 'y', for us its x coordinate). This is done for uniformity in subsetting and plotting.
SubsetToRect <- function(data, xmin, xmax, ymin, ymax){
  # crop the data
  cropped_data <- Crop(data[['fov']], x = c(xmin, xmax), y = c(ymin, ymax), coords = 'plot') 
  data_sub <- subset(data, cells = Cells(cropped_data))
  return (data_sub)
}

# function to sample a given number of coordinates in the specified area.
# returns a list with each entry a point (vector containing two values, the x and y coordinates)
SamplePointsWithinCoors <- function(xmin, xmax, ymin, ymax, size){
  
  # some assertions
  stopifnot(xmin < xmax)
  stopifnot(ymin < ymax)
  
  # get the x and y coordinates of all points
  x_coors <- sample(seq(from = xmin, to = xmax, by = 1), size = size)
  y_coors <- sample(seq(from = ymin, to = ymax, by = 1), size = size)
  
  # create a list of the points
  points <- list()
  for (i in seq(from = 1, to = size, by = 1)){
    points[[i]] <- c(x_coors[i], y_coors[i])
  }
  return (points)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~ FUNCTIONS SPECIFIC TO THE SENDER-RECEIVER PART OF GABAERGIC NICHES ~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# add a new column named 'sender_or_receiver', storing if cell is sender, receiver or none. While 
# creating it, have only in mind, usage in G34RV project.
AddSenderReceiverInfo <- function(data, sender_genes, receiver_genes){
  # create a new col storing boolean according to yes/no expression of sender genes
  data <- AddColCellsExpressingAnyOfGenes(data, sender_genes, new_col_name = 'sender')
  
  # create a new col storing boolean according to yes/no expression of receiver genes
  data <- AddColCellsExpressingAnyOfGenes(data, receiver_genes, new_col_name = 'receiver')
  
  # using the just created cols, add information into a new col in data. 
  data$sender_or_receiver = 'none'
  data$sender_or_receiver[data$sender & data$receiver] = 'sender_and_receiver'
  data$sender_or_receiver[data$sender & !data$receiver] = 'sender'
  data$sender_or_receiver[!data$sender & data$receiver] = 'receiver'
  
  # just rearrange the order of the levels, because the default, alphabetical order isn't the best
  data$sender_or_receiver <- factor(data$sender_or_receiver, levels = c('sender', 'receiver', 'sender_and_receiver', 'none'))
  
  # remove the 'sender' and 'receiver' columns
  data$sender <- NULL
  data$receiver <- NULL
  
  # return the data object
  return (data)
}

# identify the cells which have nonzero expression of any of the genes from the given vector, add the new 
# boolean column indicating that information, and return the updated data object
AddColCellsExpressingAnyOfGenes <- function(data, genes, new_col_name){
  # make a new col with appropriate name, and initialize with all False
  data@meta.data[, new_col_name] <- FALSE
  # for the cells which express any of the given genes, make the value in that col as True
  data@meta.data[colSums(data@assays$Xenium$counts[genes, , drop = FALSE]) > 0, new_col_name] <- TRUE
  return (data)
}

# function to get the df storing celltype of sender and receiver celltypes appropriately. Assumes celltypes 
# to be stored in col named cell_type, and receiver/sender information to be stored in sender_or_receiver
GetInteractionPairsFromLOI <- function(data, locations_of_interest, distance_threshold){
  # now, for all the locations in locations_of_interest, check the celltypes of the cells composing the area around 
  # the location, with square radius equal to distance_threshold  
  interaction_pairs_df <- data.frame()
  for (location_index in seq_along(locations_of_interest[,'x'])){
    # get the corresponding data subset
    data_subset <- SubsetData(data, centre = locations_of_interest[location_index,], size = distance_threshold)
    # get the sender and receiver celltypes in this data_subset
    sender_celltype <- unname(data_subset$cell_type[data_subset$sender_or_receiver == 'sender'][1]) # we take the first index, to counter the small number of cases where there will be multiple senders in the small area
    receiver_celltype <- data_subset$cell_type[data_subset$sender_or_receiver == 'receiver']
    # if receiver celltype is none (meaning receiver_celltype is empty factor), then skip this iteration
    if (length(receiver_celltype) == 0){next}
    # now append to the dataframe of interaction pairs, the interaction pairs for the current location
    interaction_pairs_df <- rbind(interaction_pairs_df, data.frame(sender_celltype = sender_celltype, receiver_celltype = receiver_celltype))
    print(glue('iteration {location_index} done...'))
  }
  return(interaction_pairs_df)
}

# function to visualize a specific location of interest, with the cells colored by their receiver or 
# sender-ness. It assumes the sender/receiver information to be stored in col named sender_or_receiver
VisualizeLocationOfInterest <- function(data, centre, distance_threshold){
  # get a small area focusing on one specific spot where interaction is happening (keeping equal to distance 
  # threshold because then we will know exactly which cells were considered as direct neighbors)
  data_subset <- SubsetData(data, centre = centre, size = distance_threshold)
  # make the default boundary as segmentation so that in imagedimplots, we get cell boundaries instead of cell centroids
  DefaultBoundary(data_subset[['fov']]) <- "segmentation"
  # finally make the imagedimplot
  return (ImageDimPlot(data_subset, group.by = 'sender_or_receiver'))
}

# function to get the centroids of those cells which express sender gene and are also within a 
# distance threshold to another cell expressing a receiver gene
FindLocationsOfInterest <- function(data, sender_cells, receiver_cells, distance_threshold = 15){
  locations_of_interest <- data.frame()
  # loop through all combinations of sender and receiver cells
  for (sender_cell in sender_cells){
    print(glue('processing sender cell {sender_cell}...'))
    in_vicinity <- FALSE
    for (receiver_cell in receiver_cells){
      if (in_vicinity){next} # if for the given sender, we have already found a nearby receiver, then don't need to proceed further as the coordinate has already been added in df
      # check if they are in vicinity
      in_vicinity <- InVicinity(data, sender_cell, receiver_cell, distance_thresh = distance_threshold)
      # if they were in vicinity, then append the centroid of the sender cell in locations_of_interest
      if (in_vicinity){
        sender_centroid <- GetCellCentroid(data, sender_cell)
        locations_of_interest <- rbind(locations_of_interest, data.frame(x = sender_centroid[1], y = sender_centroid[2]))
      }
    }
  }
  return (locations_of_interest)
}

# get factors of a number (argument called area) in the provided aspect ratio. 
# This function was originally made for the coherence step, to do the gridding according to 
# dimensions of the tissue, while keeping the number of grid-pixels approximately equal to num-cells.
# I originally made this func for coherence bin deciding for EPN project.
SplitNumberInAspectRatio <- function(area, width, height){
  # get the coefficient to which the provided width and height need to be multiplied
  coeff <- area/(width*height)
  
  # now we just need to multiply the provided width and height to sqrt of the obtained coefficient to get 
  # the right values which multiply to be approximately = area, and are in required ratio
  nrow = height*sqrt(coeff)
  ncol = width*sqrt(coeff)
  
  return (c(ncol, nrow))
}








