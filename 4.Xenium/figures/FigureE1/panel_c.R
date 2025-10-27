# code for extended figure 1 panel c

# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/Multidimensional_STEPN/4.Xenium/resources/FLXenium/')

# source the ependymoma project specific functions
source('~/Multidimensional_STEPN/4.Xenium/resources/epn_functions.R')

# get the vars for epn analyses
all_vars <- SetUpEpendymomaGlobalVarsGeneral()
# read the single cell object with information about mal cells from all samples from all subtypes 
sc <- qread(all_vars$sc_full_path)

# declare the markers and make the featureplots

myeloid_markers <- c('CD14', 'CSF1R')
FeaturePlot(sc, features = myeloid_markers)

tcell_markers <- c('CD3E', 'CD4')
FeaturePlot(sc, features = tcell_markers)

oligo_markers <- c('MBP', 'PLP1')
FeaturePlot(sc, features = oligo_markers)

endo_markers <- c('VWF', 'CLDN5')
FeaturePlot(sc, features = endo_markers)
