# code for extended figure 3 panel a

# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/Multidimensional_STEPN/4.Xenium/resources/FLXenium/')

# source the ependymoma project specific functions
source('~/Multidimensional_STEPN/4.Xenium/resources/epn_functions.R')

# get the vars for epn analyses
all_vars <- SetUpEpendymomaGlobalVarsGeneral()
# read the single cell object with information about mal cells from all samples from all subtypes 
sc <- qread(all_vars$sc_mal_path)




# make the plot
DimPlot(sc, group.by = 'Sample_deID', cols = col_sampleid_all, pt.size = 1) + theme_void() + NoLegend() + ggtitle(NULL)
ggsave('~/plot.pdf', width = 4.8, height = 4.3)


