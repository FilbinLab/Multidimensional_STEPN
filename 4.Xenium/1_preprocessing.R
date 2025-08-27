# this file takes a single input, the sample_name (such as 'STEPN47_Region_2'). One can simply create a variable 
# by running something like: "sample_name <- 'STEPN47_Region_2'", and run this script line by line, skipping the 
# first part where the sample_name is read from supplied arguments, or can execute this script from the terminal 
# using something like: $ Rscript 1_preprocessing.R STEPN47_Region_2

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
sample_raw_data_dir <- all_vars$sample_raw_data_dir

# the location where the preprocessed file will be stored. For convenience, we 
# are creating the directory if it doesn't exist already
preprocessing_dir <- all_vars$preprocessing_dir
if (!dir.exists(preprocessing_dir)) {dir.create(preprocessing_dir, recursive = T)}

# Process sample
RunFullXeniumEPN(smp = sample_name, 
                 rawDir = sample_raw_data_dir, 
                 outFile = all_vars$preprocessing_file)



