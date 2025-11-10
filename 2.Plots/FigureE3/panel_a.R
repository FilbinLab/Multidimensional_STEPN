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

# make the subtype col into factor type so that I can have the order of subtypes appropriately in the plot
sc$Subtype <- factor(sc$Subtype, levels = c('ZFTA-RELA', 'ZFTA-Cluster 1', 'ZFTA-Cluster 2', 'ZFTA-Cluster 3', 'ZFTA-Cluster 4', "ST-YAP1"))

# visualize as boxplots for each mp
df <- sc@meta.data %>% select(Sample_deID, Metaprogram, Subtype) %>% 
  group_by(Sample_deID, Metaprogram) %>% 
  mutate(counts = n()) %>% unique() %>% ungroup() %>% 
  group_by(Sample_deID) %>% mutate(sample_counts = sum(counts)) %>% 
  mutate(proportion = (counts)/sample_counts)
# loop through each mp and make boxplot corresponding to it
plots_list <- list()
for (i in seq_along(unique(df$Metaprogram))){
  mp = unique(df$Metaprogram)[i]
  # mp = 'Embryonic-neuronal-like'
  df_subset <- df %>% filter(Metaprogram == mp)
  plot <- ggplot(df_subset, aes(x = Subtype, fill = Metaprogram, y = proportion)) + 
    geom_boxplot(outliers = F) + 
    geom_jitter(size = 1) +
    theme_classic() + 
    scale_fill_manual(values = colors_metaprograms_Xenium) + 
    ylab(NULL) + xlab(NULL) + NoLegend() + ggtitle(mp) + 
    theme(plot.title = element_text(hjust = 0.5, size = 11), 
          axis.text = element_text(size = 10), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  # scale_y_continuous(limits = c(0, 1))
  plots_list[[i]] <- plot
  # ggsave(glue('plots/boxplot_{mp}_proportion_comparison_nolg.pdf'), width = 2, height = 2.6)
}
patchwork::wrap_plots(plots_list, ncol = 4)
ggsave(glue('plots/boxplot_celltypes_proportion_comparison2_nolg.pdf'), width = 7.6, height = 4.9)

# now do the statistic on the above plot
library(speckle)
results <- propeller(clusters = sc$Metaprogram, sample = sc$Sample_deID, group = sc$Subtype)
# apparently when there are more than one groups, and anova is done, there is no col indicating the clusters. Hence add that col before saving
results$cluster <- rownames(results)
results <- results[, c(ncol(results), 1:(ncol(results)-1))]

