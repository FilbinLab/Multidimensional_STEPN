args <- commandArgs(trailingOnly = TRUE)
print(args)

SampleName_arg <- args[1]

library(xlsx)
options(future.globals.maxSize = 6000 * 1024^2)

# Organize environment  ----------------------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium')

analysis_dir <- file.path(base_dir, 'analysis/2_preprocessing')

resource_dir  <- file.path(base_dir, 'scripts_revisions/resources')
source(file.path(resource_dir, "xenium_preprocessing_helper_functions_SD_v2.R"))

seurat_obj_dir_ref <- file.path(base_dir, 'analysis/1_preparation/data')

# read metadata ---------------------------------------------------------
metadata <- read_xlsx(file.path(base_dir, 'scripts_revisions/SampleIdentifier.xlsx'))
metadata <- metadata %>% filter(SampleName == SampleName_arg)
sampleDirs <- file.path(base_dir, 'data/raw_data', metadata$RawDataPath)

# List folder path ------------------------------------------------
analysis_dir_o <- file.path(analysis_dir,  SampleName_arg)
if (!dir.exists(analysis_dir_o)) {dir.create(analysis_dir_o, recursive = T)}


# Process samples ------------------------------------------------
lapply(sampleDirs, function(x) {
  RunFullXenium(smp = SampleName_arg,
                rawDir = sampleDirs,
                outDir = analysis_dir_o)
})
