# code for extended figure 5

# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/FLXenium/functions/')

# source the ependymoma project specific functions
source('~/ependymoma/xenium/scripts_revisions/resources/epn_functions.R')

all_vars <- SetUpEpendymomaGlobalVarsGeneral()
metadata <- all_vars$metadata
metadata <- metadata %>% filter(Subtype == 'ZFTA-RELA')
sample_names <- metadata$SampleName

# loop through all samples and append their imagedimplot into plots_list
library(ggrastr) # for rasterizing plots which keeps the size of image low
plots_list <- list()
for (i in seq_along(sample_names)){
  sample_name = sample_names[i]
  all_vars <- SetUpEpendymomaGlobalVars(sample_name)
  data <- qread(all_vars$annotation_file)
  
  # make the plot and rasterize and add to list
  plot <- ImageDimPlot(data, fov = 'fov', group.by = 'cell_type', cols = all_vars$colors, dark.background = F) + NoLegend() + ggtitle(sample_name)
  rasterized_plot <- rasterize(plot, layers = 'Point', dpi = 100)
  plots_list[[i]] <- rasterized_plot
  print(glue('{i}/{length(sample_names)} done...'))
}

# now that we have the list of plots, simply use patchwork to wrap the plots and make a combined plot object
final_plot <- patchwork::wrap_plots(plots_list, ncol = 5)






