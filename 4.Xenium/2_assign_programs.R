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

# get all the global variable values from all_vars 
preprocessing_file <- all_vars$preprocessing_file
ref_projection <- all_vars$ref_projection
marker_genes_list <- all_vars$marker_genes_list
malignant_celltypes <- all_vars$malignant_celltypes
manual_annotation_yaml <- all_vars$manual_annotation_yaml
colors <- all_vars$colors
plot_dir <- glue('{all_vars$annotation_dir}/plots')
data_dir <- glue('{all_vars$annotation_dir}/data')

# read the data object
data <- qread(preprocessing_file)

# run while testing, for fast iteration speed
# data <- SubsetProportionCells(data, 0.1)

# Assign labels and save resultant data and plots
XeniumCellAssignment(data, ref_projection, sample_name, marker_genes_list, malignant_celltypes,
                     manual_annotation_yaml, colors, data_dir, plot_dir)







