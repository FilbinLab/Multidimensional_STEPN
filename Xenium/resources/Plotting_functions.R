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
