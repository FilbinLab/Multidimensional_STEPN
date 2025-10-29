library(dplyr)
library(Seurat)
library(ggplot2)

## Style plots

# color palettes
col_sample = c('#370617',
          '#e01e37',
          '#f6cacc',
          '#ffba08',
          '#ffa200',
          '#d4d700',
          '#55a630',
          '#8be8d7',
          '#2fb5c7',
          '#0377a8',
          '#002855',
          '#a564d3',
          '#ff5ca5',
          '#ffb9d8',
          '#bc1f66',
          '#dcb9a1')

col_cycling <- c('red', 'black')


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

## Style plots

theme_vln <- theme(panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.background = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   axis.text.x = element_text(size=14, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                   axis.text.y = element_text(size=14, colour="black"),
                   axis.title=element_text(size=14),
                   plot.title = element_text(size=14, face="bold")) 


theme_cellstate_plot <-  theme(panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               plot.background = element_blank(),
                               panel.background = element_blank(),
                               axis.line = element_line(colour = "black"),
                               axis.title=element_text(size=12),
                               axis.text.x = element_text(size=12, vjust = 0.5, colour="black"),
                               axis.text.y = element_text(size=12, colour="black"),
                               plot.title = element_text(hjust = 0.5, angle = 0, size = 12, face = "bold", vjust = 1),
                               plot.subtitle = element_text(hjust = 0.5),
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 12, face = "bold"),
                               legend.box.background = element_rect(color="black", size=1),
                               legend.key = element_rect(fill = 'white')) 

theme_tSNE <-  theme( plot.title = element_text(hjust = 0.5, angle = 0, size = 14, face = "bold", vjust = 1),
                      plot.subtitle = element_text(hjust = 0.5, size = 12),
                      legend.text = element_text(size = 12),
                      legend.title = element_text(size = 12, face = "bold")) 

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

## Customized ggplot2 theme for plotting
seurat_theme <- function(){
  theme_bw() +
    theme(panel.background = element_rect(colour = "black", size=0.1),
          plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
          axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}


## Piechart theme
piechart_theme <- theme(panel.border = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        plot.background = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.title=element_blank(),
                        #legend text 
                        legend.text=element_text(size=12),
                        #style facet text
                        #strip.text.x = element_text(size=12, face="bold"),
                        #strip.text.y = element_text(size=12, face="bold"),
                        #strip.background = element_rect(fill="grey80"),
                        # space facet
                        panel.spacing = unit(2, "lines"))


# Violin plot
gg_violinplot_style <- function(df, fill_var) {
  ggplot(df, aes(x=Sampling, y=freq, fill = fill_var)) +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
          axis.text.y = element_text(size=12, colour="black"),
          axis.title=element_text(size=12),
          #legend text 
          legend.text=element_text(size=12),
          #style facet text
          strip.text.x = element_text(size=12, face="bold", vjust = 0),
          strip.text.y = element_text(size=12, face="bold", vjust = 0),
          strip.background = element_rect(fill=NA, size=1.5),
          # space facet
          panel.spacing = unit(1, "lines")) +
    scale_y_continuous(labels = scales::percent) 
}

