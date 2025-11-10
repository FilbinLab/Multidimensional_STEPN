# code for extended figure 6 panel b

# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/Multidimensional_STEPN/4.Xenium/resources/FLXenium/')

# source the ependymoma project specific functions
source('~/Multidimensional_STEPN/4.Xenium/resources/epn_functions.R')

# get the vars for epn analyses
all_vars <- SetUpEpendymomaGlobalVarsGeneral()
niche_dir <- glue('{all_vars$niche_dir}/data')
# subset metadata to just zfta-rela samples
metadata <- all_vars$metadata
metadata <- metadata %>% filter(Subtype == 'ZFTA-RELA')

# make the niche heatmap, corresponding pie-charts showing composition of celltypes in it, and save it in 
# plot_dir (change value of plot_dir argument according to requirement)
HeatmapNicheCorrelation(metadata, niche_dir, niche_resolution = 6, n_clusters = 6, plot_dir = '~/Multidimensional_STEPN/4.Xenium/figures/results', plot_suffix = '')




