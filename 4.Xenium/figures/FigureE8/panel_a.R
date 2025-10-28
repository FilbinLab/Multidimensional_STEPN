# code for extended figure 8 panel a

# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/Multidimensional_STEPN/4.Xenium/resources/FLXenium/')

# source the ependymoma project specific functions
source('~/Multidimensional_STEPN/4.Xenium/resources/epn_functions.R')

# get the vars for epn analyses
all_vars <- SetUpEpendymomaGlobalVarsGeneral()
# read the single cell object with information about mal cells from all samples from all subtypes 
sc <- qread(all_vars$sc_mal_path)

# subset sc to just ZFTA-RELA and remove Embryonic-like cells as they are just 9 out of 2854 cells
Idents(sc) <- 'Subtype'
sc <- subset(sc, idents = 'ZFTA-RELA')
Idents(sc) <- 'Metaprogram_noCC'
sc <- subset(sc, idents = 'Embryonic-like', invert = TRUE)
# make the Metaprogram_noCC col into factor for correct order of MPs in dotplot
sc$Metaprogram_noCC <- factor(sc$Metaprogram_noCC, levels <- rev(c('Neuroepithelial-like', 'Radial-glia-like', 'Embryonic-neuronal-like', 'Neuronal-like', 'Ependymal-like', 'MES/Hypoxia')))
Idents(sc) <- 'Metaprogram_noCC' # apparently after making it factor, need to reassign as Ident to be able to have the mps in correct order

# since dotplot function extracts the information from data layer but we want counts information, we copy into data, the counts information.
DotPlot(sc, features = rev(c('CCDC40', 'MAP3K19', 'DNER', 'SOX10', 'TUBB3', 'DCX', 'CENPF', 'ARL6IP1', 'VIM', 'S100B', 'NES')),   
        col.min = 0, col.max = 1.8) + 
  scale_color_gradientn(colors = color_dotplot) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ylab(NULL) + xlab(NULL) + NoLegend()
ggsave('~/plot.pdf', width = 4.5, height = 2.3)
