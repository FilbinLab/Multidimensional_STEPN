# code for extended figure 1 panel g

# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/Multidimensional_STEPN/4.Xenium/resources/FLXenium/')

# source the ependymoma project specific functions
source('~/Multidimensional_STEPN/4.Xenium/resources/epn_functions.R')

# get the vars for epn analyses
all_vars <- SetUpEpendymomaGlobalVarsGeneral()

# read the sc object and reference
sc <- qread(all_vars$sc_mal_path)
ref <- qread(all_vars$sc_reference)

# rename the 'MGE inhibitory neuron lineage' program to just 'MGE inhibitory', because shorter
ref$cell_type[ref$cell_type == "MGE inhibitory neuron lineage"] <- 'MGE inhibitory'
# the order in which we want the metaprograms in the plot
order_programs <- c("Neuroepithelial", "Radial glia", "Early radial glia", "Newborn/maturing neuron",
                    "MGE inhibitory", "Neuronal", "IPCs", "Mesenchymal", "OPC", "Astrocyte", "Choroid")
rev_order_programs <- c("Choroid", "Astrocyte", "OPC", "Mesenchymal", "IPCs", "Neuronal", "MGE inhibitory",
                        "Newborn/maturing neuron", "Early radial glia", "Radial glia", "Neuroepithelial")

# get the df which has the information relevant to plotting
df_long <- ref@meta.data %>% select(cell_type, Age_in_Weeks)
rownames(df_long) <- NULL
# get counts for each combination of timepoint and celltype and make the df into wide form as easy to normalize then
df <- df_long %>% group_by(Age_in_Weeks, cell_type) %>% summarise(counts = n(), .groups = 'drop')
df <- df %>% pivot_wider(names_from = Age_in_Weeks, values_from = counts, values_fill = 0)
# make the numeric cols of df, normalized
for (col in colnames(df)){
  if (is.numeric(df[[col]])){
    df[[col]] <- df[[col]]/sum(df[[col]])
  }
}
# make back into long form
df <- df %>% pivot_longer(!cell_type, names_to = 'Age_in_Weeks', values_to = 'proportion')

# make sure the Age in Weeks col is in factor form (so that in plot, ordering is correct)
df$Age_in_Weeks <- factor(as.numeric(df$Age_in_Weeks))
# make sure the celltypes are also in the correct order
df$cell_type <- factor(df$cell_type, levels = rev_order_programs)
# make the plot
plot <- ggplot(df, aes(y = cell_type, x = Age_in_Weeks, colour = cell_type)) + 
  geom_point(aes(size = proportion, alpha = proportion == 0)) + 
  scale_color_manual(values = colors_sc_ref, guide = 'none') + 
  theme_classic() + 
  scale_alpha_manual(values = c(1,0)) + 
  guides(alpha = 'none') + 
  scale_size_continuous(breaks = c(0.05, 0.25, 0.50, 0.75, 1.00)) + 
  labs(size = 'Proportion') + 
  xlab('Age in Weeks') + 
  ylab('Cell Type') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
plot + NoLegend()
ggsave('~/dotplot_celltype_vs_age_classic.pdf', width = 7.8, height = 2.5)
