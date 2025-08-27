# in this file, I load the required libraries, and have the functions to set up the global variables

# Load necessary libraries. User should have them installed 
library(tidyverse)
library(glue)
library(qs)
library(Seurat)
library(readxl)
library(patchwork) # keeping as required because it is used quite widely in the projects
library(crayon) # required for colorful status update messages
library(yaml)
library(UCell) 
library(progress) # for making progress bars
library(paletteer)
library(RColorBrewer)

# OTHER, OPTIONAL LIBRARIES
# library(SeuratWrappers)
# library(mclust)
# library(progress)
# library(ComplexHeatmap)
# library(speckle) # for the propeller function, used for computing p-vals for composition comparisons
# library(ggpubr) # used for getting legends from plot etc
# library(parallel)
# library(cowplot)
# library(circlize)
# library(reticulate) # for being able to use FindFunction() function, which runs python code for running embedding based search
# library(ggrastr) # for rasterizing ImageDimPlots and ImageFeaturePlots which don't have a rasterize functionality in seurat
options(future.globals.maxSize = 20 * 1024^3) # 20gb

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PROJECT SET UP FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# source the script with project specific functions for ependymoma project
LoadEPNFunctions <- function(){
  source('/home/shk490/ependymoma/xenium/scripts_revisions/epn_functions.R')
}







