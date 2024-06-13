library(dplyr)
library(Seurat)
library(ggplot2)

## Style plots

# Barplot function
theme_ggplot = theme(legend.position = "none",
      plot.title = element_text(hjust=0.5, face="bold"),
      panel.border = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1, colour="black"),
      axis.text.y = element_text(size=12, colour="black"),
      axis.title=element_text(size=12),
      strip.text = element_text(size = 13, face = "bold"))

theme_ggplot_legend = theme(plot.title = element_text(hjust=0.5, face="bold"),
                     panel.border = element_blank(),
                     plot.background = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1, colour="black"),
                     axis.text.y = element_text(size=12, colour="black"),
                     axis.title=element_text(size=12),
                     strip.text = element_text(size = 13, face = "bold"))
