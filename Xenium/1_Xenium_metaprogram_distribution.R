# Load packages -----------------------------------
#library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
#library(ggpubr)
library(data.table)
library(readxl)
library(dplyr)
library(ggrastr)
library(plyr)
library(tidyverse)

# Organize environment  -----------------------------------
base_dir <- file.path('/n/scratch/users/s/sad167/EPN/Xenium/analysis')

resource_dir  <- file.path('/n/scratch/users/s/sad167/EPN/Ependymoma2023/Xenium/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))

data_files <- c('20231020__200939__BT2126_BT1745/data/individual/cell_ID_0010652-Region_4.csv',
              
              '20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010575-Region_1.csv',
              '20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010575-Region_2.csv',
              '20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010575-Region_3.csv',
              '20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010619-Region_1.csv',
              '20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010619-Region_2.csv',
              '20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010619-Region_3.csv',
              '20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010619-Region_4.csv',
              '20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010619-Region_5.csv',
              
              '20231107__203958__BT1717_BT775/data/individual/cell_ID_0010501-Region_1.csv',
              '20231107__203958__BT1717_BT775/data/individual/cell_ID_0010501-Region_2.csv',
              #'20231107__203958__BT1717_BT775/data/individual/cell_ID_0010814-Region_1.csv',
              #'20231107__203958__BT1717_BT775/data/individual/cell_ID_0010814-Region_2.csv',
              
              '20231109__203408__BT1804_BT2169/data/individual/cell_ID_0010498-Region_1.csv', 
              '20231109__203408__BT1804_BT2169/data/individual/cell_ID_0010775-Region_1.csv', 
              
              '20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/individual/cell_ID_0010540-Region_1.csv',
              '20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/individual/cell_ID_0010540-Region_2.csv',
              '20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/individual/cell_ID_0010540-Region_3.csv',
              '20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/individual/cell_ID_0010540-Region_4.csv',
              '20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/individual/cell_ID_0010553-Region_1.csv',
              '20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/individual/cell_ID_0010553-Region_2.csv',
              '20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/individual/cell_ID_0010553-Region_3.csv',
              '20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/individual/cell_ID_0010553-Region_4.csv')

metadata <- list()
proportions <- list()
for (i in seq_along(data_files)) {
  metadata[i] <- read.table(file.path(base_dir, data_files[i]), header = TRUE, sep = ',', row.names = "cell_id")
  #metadata[i] <- fread(file.path(base_dir, data_files[i]), data.table = FALSE)
  proportions[[i]] <- prop.table(table(metadata[i]))
}

names(metadata) <- c('0010652-Region_4', '0010575-Region_1', '0010575-Region_2', '0010575-Region_3', 
                     '0010619-Region_1', '0010619-Region_2', '0010619-Region_3', '0010619-Region_4', '0010619-Region_5', 
                     '0010501-Region_1', '0010501-Region_2', 
                     #'0010814-Region_1', '0010814-Region_2',
                     '0010498-Region_1', '0010775-Region_1', 
                     '0010540-Region_1', '0010540-Region_2', '0010540-Region_3', '0010540-Region_4', 
                     '0010553-Region_1', '0010553-Region_2', '0010553-Region_3', '0010553-Region_4')
# names(proportions) <- c('0010652-Region_4', '0010575-Region_1', '0010575-Region_2', '0010575-Region_3', 
#                      '0010619-Region_1', '0010619-Region_2', '0010619-Region_3', '0010619-Region_4', '0010619-Region_5', 
#                      '0010501-Region_1', '0010501-Region_2', '0010814-Region_1', '0010814-Region_2',
#                      '0010498-Region_1', '0010775-Region_1', 
#                      '0010540-Region_1', '0010540-Region_2', '0010540-Region_3', '0010540-Region_4', 
#                      '0010553-Region_1', '0010553-Region_2', '0010553-Region_3', '0010553-Region_4')

names(proportions) <- c('BT2126_Region_4', '7EP41_Region_1', '7EP41_Region_2', '7EP41_Region_3', 
                        '3EP8_Region_1.1', '3EP8_Region_2.1', '3EP8_Region_3', '7EP1_Region_4.1', '7EP1_Region_5', 
                        'BT775_Region_1', 'BT775_Region_2', 
                        #'BT1717_Region_1', 'BT1717_Region_2',
                        'BT2169_Region_1', 'BT1804_Region_1', 
                        '11EP22_Region_1', '11EP22_Region_2', '11EP22_Region_3', '7EP41_Region_4', 
                        '3EP8_Region_1.2', '3EP8_Region_2.2', 'BT1743_Region_3', '7EP1_Region_4.2')

# bind rows together and transform into df
metaprogram_proportion <- bind_rows(proportions, .id = "Xenium_Region")
metaprogram_proportion <- as.data.frame(metaprogram_proportion)

# change df structure for plotting
df_melted <- metaprogram_proportion %>%
  gather(key = "Metaprogram", value = "Frequency", -Xenium_Region)
#write_csv(df_melted, file.path(base_dir, 'Xenium_Metaprogram_proportion.csv'))

# reorder so that it is the same as in the paper
df_melted$Metaprogram <- factor(df_melted$Metaprogram, 
                                levels = c('Neuroepithelial-like', "Radial glia-like", 'NPC-like',
                                           "Ependymal-like", "Mesenchymal", 'Myeloid', 'T-cell', 'Endothelial', 'VLMCs', 'Microglia',
                                           'Oligodendrocytes', 'Neurons', 'Unassigned', '13', '9', '10', '11', '14', '12', '15'))


# Plot a stacked bar plot
ggplot(df_melted, aes(x = Xenium_Region, y = Frequency, fill = Metaprogram)) +
  scale_fill_manual(values = colors_groups_barplot) +
  geom_bar(stat = "identity", position = "fill", color="black") +
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
ggsave(file.path(base_dir, 'Metaprogram_proportion.pdf'), width=10, height=5)


# export average proportion
metaprogram_proportion_summary <- metaprogram_proportion %>%
  summarise_all(mean,  na.rm = TRUE)

metaprogram_proportion_summary <- t(metaprogram_proportion_summary)
metaprogram_proportion_summary <- as.data.frame(metaprogram_proportion_summary)

metaprogram_proportion_summary <- metaprogram_proportion_summary[ ,-c(1, 11:17)]

write_csv(metaprogram_proportion_summary, file.path(base_dir, 'Xenium_Metaprogram_proportion_summary.csv'))

## Recreate plot only with 1 region per sample ---------------------------

data_files <- c('20231020__200939__BT2126_BT1745/data/individual/cell_ID_0010652-Region_4.csv',
                '20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010575-Region_3.csv',
                '20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010619-Region_1.csv',
                '20231102__215055__7EP1_7EP41_3EP8/data/individual/cell_ID_0010619-Region_5.csv',
                '20231107__203958__BT1717_BT775/data/individual/cell_ID_0010501-Region_2.csv',
                '20231109__203408__BT1804_BT2169/data/individual/cell_ID_0010498-Region_1.csv', 
                '20231109__203408__BT1804_BT2169/data/individual/cell_ID_0010775-Region_1.csv', 
                '20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/individual/cell_ID_0010540-Region_1.csv',
                '20231208__193657__3EP8_BT1743_7EP1_11EP22_7EP41/data/individual/cell_ID_0010553-Region_3.csv')

metadata <- list()
proportions <- list()
for (i in seq_along(data_files)) {
  metadata[i] <- read.table(file.path(base_dir, data_files[i]), header = TRUE, sep = ',', row.names = "cell_id")
  #metadata[i] <- fread(file.path(base_dir, data_files[i]), data.table = FALSE)
  proportions[[i]] <- prop.table(table(metadata[i]))
}

names(metadata) <- c('BT2126_Region_4', '7EP41_Region_3', '3EP8_Region_1', '7EP1_Region_5',
                     'BT775_Region_2', 'BT2169_Region_1', 'BT1804_Region_1',
                     '11EP22_Region_1',  'BT1743_Region_3')

names(proportions) <- c('BT2126_Region_4', '7EP41_Region_3', '3EP8_Region_1', '7EP1_Region_5',
                     'BT775_Region_2', 'BT2169_Region_1', 'BT1804_Region_1',
                     '11EP22_Region_1',  'BT1743_Region_3')


# bind rows together and transform into df
metaprogram_proportion <- bind_rows(proportions, .id = "Xenium_Region")
metaprogram_proportion <- as.data.frame(metaprogram_proportion)

# change df structure for plotting
df_melted <- metaprogram_proportion %>%
  gather(key = "Metaprogram", value = "Frequency", -Xenium_Region)
#write_csv(df_melted, file.path(base_dir, 'Xenium_Metaprogram_proportion.csv'))

# reorder so that it is the same as in the paper
df_melted$Metaprogram <- factor(df_melted$Metaprogram, 
                                levels = c('Neuroepithelial-like', "Radial glia-like", 'NPC-like',
                                           "Ependymal-like", "Mesenchymal", 'Myeloid', 'T-cell', 'Endothelial', 'VLMCs', 'Microglia',
                                           'Oligodendrocytes', 'Neurons', 'Unassigned', '13', '9', '10', '11', '14', '12', '15'))


# Plot a stacked bar plot
ggplot(df_melted, aes(x = Xenium_Region, y = Frequency, fill = Metaprogram)) +
  scale_fill_manual(values = colors_groups_barplot) +
  geom_bar(stat = "identity", position = "fill", color="black") +
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
ggsave(file.path(base_dir, 'Metaprogram_proportion_selected.pdf'), width=6, height=5)

