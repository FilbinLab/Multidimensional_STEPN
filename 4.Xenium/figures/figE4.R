# code for extended figure 4

# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/FLXenium/functions/')

# source the ependymoma project specific functions
source('~/ependymoma/xenium/scripts_revisions/resources/epn_functions.R')

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
zfta_reference <- qread(all_vars$zfta_reference)

### panel b

PlotCompositionComparisonSpatialSC(zfta_reference, metadata, cellid_dir)

### panel d

PlotXeniumMetaprogramsCoherence(
  metadata, cellid_dir, coherence_data_dir, order_metaprograms, colors_metaprograms_Xenium, 
  malignant_celltypes, average_by_sample_id = FALSE, order_by_coherence = TRUE)

### panel e

OverallCoherenceCelltypeProportionLinearRegression(
  metadata, cellid_dir, coherence_data_dir, average_by_sample_id = TRUE,
  colors = colors_metaprograms_Xenium_dot_separator)

### panel f

OverallCoherenceCelltypeProportionLinearRegressionAfterRemoval(
  celltype_to_remove = 'MES-like', metadata, cellid_dir, coherence_data_dir, average_by_sample_id = TRUE)













