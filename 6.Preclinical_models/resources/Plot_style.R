library(dplyr)
library(Seurat)
library(ggplot2)

## load signatures for scoring 
## Read ZFTA::RELA patient scores
NMFmarkers <- readRDS(file.path(resource_dir, "ZFTA.gene.set.merged.rds"))
## Combine neuronal signatures
NMFmarkers$Neurons <- c(NMFmarkers$ZFTA_NPC.like.1, NMFmarkers$ZFTA_NPC.like.3,NMFmarkers$ZFTA_NPC.like.2)
## Create list with all markers
ref.gene.EPN <- list(NMFmarkers$ZFTA_MES_like, NMFmarkers$ZFTA_Radial.glia.like, NMFmarkers$ZFTA_Ependymal.like, NMFmarkers$ZFTA_Early_Glial_Progenitor, NMFmarkers$ZFTA_Cycling, NMFmarkers$Neurons)
names(ref.gene.EPN) <- c('Mesenchymal', 'Radial_glia', 'Ependymal', 'Early_progenitor', 'Cycling', 'NPC')

## Color scales

col_line <- paletteer::paletteer_d("vapoRwave::avanti")

paletteer::paletteer_d("vapoRwave::avanti")



## Style plots
theme_vln <- theme(panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.background = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                  axis.text.y = element_text(size=12, colour="black"),
                  axis.title=element_text(size=12),
                  plot.title = element_text(size=12, face="bold")) 


theme_cellstate_plot <-  theme(panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               plot.background = element_blank(),
                               panel.background = element_blank(),
                               axis.line = element_line(colour = "black"),
                               axis.title=element_text(size=12),
                               axis.text.x = element_text(size=12, vjust = 0.5, colour="black"),
                               axis.text.y = element_text(size=12, colour="black"),
                               plot.title = element_text(hjust = 0.5, angle = 0, size = 14, face = "bold", vjust = 1),
                               plot.subtitle = element_text(hjust = 0.5),
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 12, face = "bold"),
                               legend.box.background = element_rect(color="black", size=1),
                               legend.key = element_rect(fill = 'white')) 

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
          axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
          axis.text.y = element_text(size=12, colour="black"),
          axis.title=element_text(size=12)) 
}



