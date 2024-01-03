# Color palettes -------------------------------------
colors_groups <- c("gray30","#F99E93FF","#9E5E9BFF","#74ADD1FF","#ACD39EFF","#96410EFF", 'grey80',
                   '#b81702', '#370617',  '#ffba08', '#2fb5c7',  '#ff0072')
names(colors_groups) <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
                          "NPC-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
                          "Immune", "Microglia", "Endothelial", "Neurons", "VLMCs")

col_niches <- c('#58148e', '#15a2a2', '#ea515f','#d0a03b','#0466c8')


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
