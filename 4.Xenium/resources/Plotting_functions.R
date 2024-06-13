options(future.globals.maxSize = +Inf)

# Color palettes -------------------------------------

colors_metaprograms_Xenium <- c("gray30","#F99E93FF","#9E5E9BFF",
                                '#0F4F8B',"#ACD39EFF","#96410EFF",'grey80', 
                                '#FFF087FF',  'turquoise3', 'turquoise2', 'violetred3', 'violetred2', '#000000FF')

names(colors_metaprograms_Xenium) <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
                                       "Neuronal-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
                                       "T-cell", "Myeloid", "Microglia", "Endothelial", "VLMCs",  "Oligodendrocytes")


colors_to_use_metaprograms <- structure(colors_metaprograms_Xenium,
                           names = names(colors_metaprograms_Xenium))


colors_metaprograms_Xenium_grouped <- c("gray30","#F99E93FF","#9E5E9BFF",
                                          '#0F4F8B',"#ACD39EFF","#96410EFF",'grey80', 
                                          'turquoise3', 'violetred3',  '#000000FF')

names(colors_metaprograms_Xenium_grouped) <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
                                                 "Neuronal-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
                                                 "Immune", "Vascular", "Oligodendrocytes")



colors_recurrent_metaprograms_Xenium <- c("gray30","#F99E93FF","#9E5E9BFF",
                                '#0F4F8B',"#ACD39EFF","#96410EFF",'grey80', 
                                  'turquoise3', 'violetred3',  '#000000FF')

names(colors_recurrent_metaprograms_Xenium) <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
                                       "Neuronal-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
                                       "Immune", "Vessel", "Oligodendrocytes")






colors_groups_barplot <- c(colors_metaprograms_Xenium,
                           'grey80','grey80','grey80','grey80','grey80','grey80','grey80')
names(colors_groups_barplot) <- c(names(colors_metaprograms_Xenium),
                                  '12', '13', '9', '10', '11', '14', '15')


col_niches <- c('#58148e', '#15a2a2', '#ea515f','#d0a03b','#0466c8')

colors_to_use_niches <- structure(col_niches,
                                        names = c('1', '2', '3', '4', '5'))


col_normal_malignant <- c('#FFC72CFF', '#582C83FF')
names(col_normal_malignant) <- c('Non-malignant', 'Malignant')


# Barplot function
plot_bar <- function(seurat_obj, x_var, y_var, colors){
  ggplot(seurat_obj@meta.data, aes(x_var, fill = y_var)) +
    scale_fill_manual(values = colors) + 
    geom_bar(position = "fill", color="black") +
    labs (y='Proportion', x='') + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=14, angle = 90, vjust = 0.5, hjust=1, colour="black"),
          axis.text.y = element_text(size=14, colour="black"),
          axis.title=element_text(size=14),
          legend.text = element_text(size = 14),
          legend.title = element_blank()) 
}

seurat_theme <- function(){
  theme_bw() +
    theme(panel.background = element_rect(colour = "black", size=0.1),
          plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
          axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}



cell_fun = function(j, i, x, y, width, height, fill) {
  grid::grid.rect(x = x, y = y, width = width *0.99, 
                  height = height *0.99,
                  gp = grid::gpar(col = "grey", 
                                  fill = fill, lty = 1, lwd = 0.5))
}
