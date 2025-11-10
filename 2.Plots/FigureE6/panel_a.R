# code for extended figure 6 panel a

# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/Multidimensional_STEPN/4.Xenium/resources/FLXenium/')

# source the ependymoma project specific functions
source('~/Multidimensional_STEPN/4.Xenium/resources/epn_functions.R')

# get the vars for epn analyses
sample_name = 'STEPN12_Region_3'
all_vars <- SetUpEpendymomaGlobalVars(sample_name)

# read the data
data <- qread(all_vars$niche_file)
DefaultAssay(data) <- 'SCT'

### working on the first area of interest
# get the two sets of coordinates 
x = 2090
y = 2910
size = 300
# subset the main data object to the coordinates
data_subset <- SubsetToSquare(data, c(x, y), size = size)
data_subset$niche_6 <- factor(data_subset$niche_6)
data_subset$niche_9 <- factor(data_subset$niche_9)
data_subset$niche_12 <- factor(data_subset$niche_12)
# make the hires image plot
p1 <- MakeHiresImageDimPlot(data_subset, group_by_col = 'cell_type', colors = colors_metaprograms_Xenium, axes = FALSE)
p2 <- MakeHiresImageDimPlot(data_subset, group_by_col = 'niche_6', colors = colors_niches, axes = FALSE)
p1+p2
ggsave('~/plot.pdf', width = 10, height = 6)

