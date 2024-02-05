# Load packages -----------------------------------
#library(Seurat)
library(ggplot2)
library(tidyverse)
#library(qs)
#library(ggpubr)
library(data.table)
library(readxl)
library(dplyr)
library(plyr)
library(ggpubr)

# Organize environment  -----------------------------------
base_dir <- file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/12 - Xenium')

data_dir <- file.path(base_dir, 'analysis/order_disorder/data')

output_dir <- file.path(base_dir, 'analysis/order_disorder/correlation')
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}


samples <- c('BT2126', '7EP41', '3EP8', '7EP1',
             'BT775', 'BT1717', 'BT2169', 'BT1804',
             '11EP22',  'BT1743')
names(samples) <- c('BT2126', '7EP41', '3EP8', '7EP1',
                    'BT775', 'BT1717', 'BT2169', 'BT1804',
                    '11EP22',  'BT1743')

colors_metaprograms_ZFTA <- c("gray30","#F99E93FF","#9E5E9BFF","#74ADD1FF","#ACD39EFF","#96410EFF", 'grey80', '#BDA14DFF', '#3EBCB6FF', '#0169C4FF', '#153460FF', '#D5114EFF' ,'#A56EB6FF' ,'#4B1C57FF', '#09090CFF')

names(colors_metaprograms_ZFTA) <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
                                     "NPC-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
                                     "T-cell", "Myeloid", "Microglia", "Endothelial", "Neurons", "VLMCs", "Oligodendrocytes")

col_sample = c('#e01e37',
              # '#f6cacc',
               '#ffba08',
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

# read metaprogram proportion
proportion <- read_csv(file.path(data_dir, 'Xenium_Metaprogram_proportion.csv'))
colnames(proportion)[1] <- 'sample'
colnames(proportion)[2] <- 'metaprogram'
proportion <- as.data.frame(proportion)


# read data (from Carlos)
spatial_coherence <- readRDS(file.path(data_dir, 'average_per_sample.rds'))

# convert to df
spatial_coherence_df <- ldply(spatial_coherence, data.frame)
colnames(spatial_coherence_df)[1] <- 'sample'
colnames(spatial_coherence_df)[3] <- 'spatial_coherence'
spatial_coherence_df$spatial_coherence <- as.numeric(spatial_coherence_df$spatial_coherence)
spatial_coherence_df <- as.data.frame(spatial_coherence_df)

# Add new columns containing metaprogram proportion
spatial_coherence_df <- left_join(spatial_coherence_df, proportion, by = c("metaprogram", "sample"))

# Add information about malignant or non-malignant
spatial_coherence_df$malignant <- ifelse(spatial_coherence_df$metaprogram %in% c('NPC-like', 'Ependymal-like', 'Neuroepithelial-like',
                                           'Mesenchymal', 'Radial glia-like'), "Malignant", "Non-malignant")

# remove columns for which we have only 1 data point available
# spatial_coherence_df <- spatial_coherence_df %>% 
#   dplyr::filter(metaprogram != "Unassigned" & metaprogram != "Radial glia-like" & metaprogram != "Oligodendrocytes" & metaprogram != "Microglia")



# calculate normalized score
spatial_coherence_df$normalized_spatial_coherence <- spatial_coherence_df$spatial_coherence*spatial_coherence_df$Frequency

# calculate scaled score
min <- min(spatial_coherence_df$spatial_coherence)
max <- max(spatial_coherence_df$spatial_coherence)
range <- max - min

spatial_coherence_df$scaled_spatial_coherence <- (spatial_coherence_df$spatial_coherence - min)/range

# Calculate average and SEM
spatial_coherence_df2 <- spatial_coherence_df %>%
  dplyr::group_by(malignant)  %>%
  dplyr::summarise(average = mean(scaled_spatial_coherence),
                  sem = sd(scaled_spatial_coherence) / sqrt(n())
  )

spatial_coherence_df2

# plot decreasing order metaprogram
spatial_coherence_df$metaprogram <- reorder(spatial_coherence_df$metaprogram, spatial_coherence_df$scaled_spatial_coherence)

ggplot(spatial_coherence_df, aes(x = metaprogram, y = scaled_spatial_coherence, color = sample, size = Frequency)) +
  geom_point() +
  labs(x = "Metaprogram",
       y = "Spatial coherence") +
  scale_color_manual(values = col_sample) +
  stat_boxplot(geom = "errorbar", color = "black") +
  stat_summary(fun = mean, colour = "black", geom = "point", shape = 8, size = 4, show.legend = FALSE) + 
  theme_minimal() + theme(panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          plot.background = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"),
                          axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                          axis.text.y = element_text(size=12, colour="black"),
                          axis.title=element_text(size=12),
                          plot.title = element_text(size=10, face="bold")) 
ggsave(file.path(output_dir, "1_spatial_coherence_by_metaprogram.pdf"), width=5, height=4)


# plot decreasing order sample
spatial_coherence_df$sample <- reorder(spatial_coherence_df$sample, spatial_coherence_df$scaled_spatial_coherence)

ggplot(spatial_coherence_df, aes(x = sample, y = scaled_spatial_coherence, color = metaprogram, size = Frequency)) +
  geom_point() +
  labs(x = "Sample",
       y = "Spatial coherence score") +
  scale_color_manual(values = colors_metaprograms_ZFTA) +
  stat_boxplot(geom = "errorbar", color = "black") +
  stat_summary(fun = mean, colour = "black", geom = "point", shape = 8, size = 4, show.legend = FALSE) + 
  theme_minimal() + theme(panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          plot.background = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"),
                          axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                          axis.text.y = element_text(size=12, colour="black"),
                          axis.title=element_text(size=12),
                          plot.title = element_text(size=10, face="bold")) 
ggsave(file.path(output_dir, "1_spatial_coherence_by_sample.pdf"), width=5, height=4)


# plot malignant vs non-malignant coherence
spatial_coherence_df2 <- spatial_coherence_df %>%
  dplyr::group_by(malignant, sample)  %>%
  dplyr::summarise(scaled_spatial_coherence = mean(scaled_spatial_coherence),
                   Frequency = mean(Frequency))

ggplot(spatial_coherence_df2, aes(x = malignant, y = scaled_spatial_coherence, color = sample, size = Frequency)) +
  geom_point() +
  stat_boxplot(geom = "errorbar", color = "black") +
  stat_summary(fun = mean, colour = "black", geom = "point", shape = 8, size = 4, show.legend = FALSE) + 
  labs(x = "Metaprogram",
       y = "Spatial coherence") +
  scale_color_manual(values = col_sample) +
  #geom_boxplot(width = 0.07) +
  theme_minimal() + theme(panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          plot.background = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"),
                          axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                          axis.text.y = element_text(size=12, colour="black"),
                          axis.title=element_text(size=12),
                          plot.title = element_text(size=10, face="bold")) 
ggsave(file.path(output_dir, "3_spatial_coherence_by_malignant.pdf"), width=3, height=4.5)



spatial_coherence_df2 <- spatial_coherence_df %>%
  dplyr::group_by(malignant, metaprogram) %>%
  dplyr::summarize(
    scaled_spatial_coherence = mean(scaled_spatial_coherence),
    Frequency = mean(Frequency)) 
write_csv(spatial_coherence_df2, file.path(output_dir, 'malignant_non_malignant.csv'))

ggplot(spatial_coherence_df2, aes(x = malignant, y = scaled_spatial_coherence, color = metaprogram, size = Frequency)) +
  geom_point() + 
  stat_boxplot(geom = "errorbar", color = "black") +
  stat_summary(fun = mean, colour = "black", geom = "point", shape = 8, size = 4, show.legend = FALSE) + 
labs(x = "Metaprogram",
       y = "Spatial coherence") +
  scale_color_manual(values = colors_metaprograms_ZFTA) +
  #geom_boxplot(width = 0.07) +
  theme_minimal() + theme(panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          plot.background = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"),
                          axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                          axis.text.y = element_text(size=12, colour="black"),
                          axis.title=element_text(size=12),
                          plot.title = element_text(size=10, face="bold")) 
ggsave(file.path(output_dir, "3_spatial_coherence_by_malignant_metaprogram.pdf"), width=3.5, height=4.5)







