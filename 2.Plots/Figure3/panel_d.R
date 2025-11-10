# code for main figure 3 panel d

# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/Multidimensional_STEPN/4.Xenium/resources/FLXenium/')

# source the ependymoma project specific functions
source('~/Multidimensional_STEPN/4.Xenium/resources/epn_functions.R')

# get the vars for epn analyses
all_vars <- SetUpEpendymomaGlobalVarsGeneral()
cellid_dir <- all_vars$cellid_dir
order_metaprograms <- all_vars$order_metaprograms
malignant_celltypes <- all_vars$malignant_celltypes
# subset metadata to just zfta-rela samples
metadata <- all_vars$metadata
metadata <- metadata %>% filter(Subtype == 'ZFTA-RELA')
sample_names <- metadata$SampleName
sample_ids <- metadata$SampleID
coherence_data_dir <- glue('{all_vars$coherence_dir}/data')

### make the plot
OverallCoherenceCelltypeProportionLinearRegression(
  metadata, cellid_dir, coherence_data_dir, average_by_sample_id = TRUE,
  colors = colors_metaprograms_Xenium_dot_separator, transform_proportions = T)