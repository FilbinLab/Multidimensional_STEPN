# code for main figure 2 panel c

# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/Multidimensional_STEPN/4.Xenium/resources/FLXenium/')

# source the ependymoma project specific functions
source('~/Multidimensional_STEPN/4.Xenium/resources/epn_functions.R')

# get the vars for epn analyses
all_vars <- SetUpEpendymomaGlobalVarsGeneral()
sc <- qread(all_vars$sc_mal_path) # read the sc object with all mal cells from all samples and subtypes

# make the plot
DimPlot(sc, group.by = 'Metaprogram_noCC', cols = colors_metaprograms_Xenium)
