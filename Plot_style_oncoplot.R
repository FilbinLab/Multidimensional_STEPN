# theme oncoplot
theme_oncoplot_no_legend <- theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), 
                                    axis.line = element_blank(),
                                    axis.title.x=element_blank(),
                                    axis.text.x=element_blank(),
                                    #axis.title.y = element_blank(),
                                    axis.title.y = element_text(angle = 0, vjust = 0.5),  # Set angle to 0 for horizontal title
                                    axis.text.y = element_blank(),
                                    axis.ticks=element_blank(),
                                    legend.position = "none",
                                    plot.margin = unit(c(0, 0, 0, 0), "cm")) 

theme_oncoplot_legend <- theme(panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_blank(), 
                               axis.line = element_blank(),
                               axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.title.y = element_blank(),
                               axis.text.y = element_blank(),
                               axis.ticks=element_blank(),
                               legend.key.size = unit(0.3, "cm"),
                               legend.text = element_text(size = 10),
                               legend.title = element_text(size = 10),
                               plot.margin = unit(c(0, 0, 0, 0), "cm")) 
