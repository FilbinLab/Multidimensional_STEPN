options(future.globals.maxSize = +Inf)

# Color palettes -------------------------------------
colors_groups <- c("gray30","#F99E93FF","#9E5E9BFF","#74ADD1FF","#ACD39EFF","#96410EFF", 'grey80',
                   '#BDA14DFF', '#3EBCB6FF', '#0169C4FF', '#153460FF', '#D5114EFF' ,'#A56EB6FF' ,'#4B1C57FF')

names(colors_groups) <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
                          "NPC-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
                          "T-cell", "Myeloid", "Microglia", "Endothelial", "Neurons", "VLMCs", "Oligodendrocytes")

colors_immune <- c(rep('grey95', 7), '#BDA14DFF', '#3EBCB6FF','#0169C4FF', rep('grey95', 4))

names(colors_immune) <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
                          "NPC-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
                          "T-cell", "Myeloid", "Microglia", "Endothelial", "Neurons", "VLMCs", "Oligodendrocytes")


colors_mesen <- c(rep('grey95', 5),"#96410EFF", rep('grey95', 8))

names(colors_mesen) <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
                          "NPC-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
                          "T-cell", "Myeloid", "Microglia", "Endothelial", "Neurons", "VLMCs", "Oligodendrocytes")


colors_vascular <- c(rep('grey95', 10), '#153460FF', 'grey95' ,'#A56EB6FF' ,'grey95')

names(colors_vascular) <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
                          "NPC-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
                          "T-cell", "Myeloid", "Microglia", "Endothelial", "Neurons", "VLMCs", "Oligodendrocytes")



colors_groups_barplot <- c("gray30","#F99E93FF","#9E5E9BFF","#74ADD1FF","#ACD39EFF","#96410EFF", 'grey80',
                           '#BDA14DFF', '#3EBCB6FF', '#0169C4FF', '#153460FF', '#D5114EFF' ,'#A56EB6FF' ,'#4B1C57FF', 
                   'grey80','grey80','grey80','grey80','grey80','grey80','grey80')
names(colors_groups_barplot) <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
                                  "NPC-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
                                  "T-cell", "Myeloid", "Microglia", "Endothelial", "Neurons", "VLMCs", "Oligodendrocytes",
                          '12', '13', '9', '10', '11', '14', '15')


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