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
annotation_file <- all_vars$annotation_file
niche_dir <- all_vars$niche_dir

# read in the seurat object
data <- qread(annotation_file)

# Perform niche analysis
data <- NicheAnalysis(data, niches_min = 4, niches_max = 9, n_neighbors = 20)

# save the information about the niches in csv 
data$cell_id <- rownames(data@meta.data)
data$group <- data$cell_type
df <- data@meta.data[, c('cell_id', 'group', 'niche_4', 'niche_5', 'niche_6', 'niche_7', 'niche_8', 'niche_9')]
write_csv(df, glue('{niche_dir}/{sample_name}.csv'))

# save the actual object too, just in case
qsave(data, glue('{niche_dir}/{sample_name}.qs'))


