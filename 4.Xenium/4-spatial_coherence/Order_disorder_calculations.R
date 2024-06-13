# Load packages -----------------------------------
#library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(reshape)
library(data.table)
library(readxl)
library(dplyr)
library(plyr)
library(ggpubr)

# Organize environment  -----------------------------------
base_dir <- file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/12 - Xenium')

data_dir <- file.path(base_dir, 'analysis/all_Xenium_runs/data')

output_dir <- file.path(base_dir, 'analysis/order_disorder/correlation')
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}

resource_dir <- file.path(base_dir, 'scripts/resources')
source(file.path(resource_dir, 'Plotting_functions.R'))


colors_metaprograms_ZFTA <- c("gray30","#F99E93FF","#9E5E9BFF",
                                '#0F4F8B',"#ACD39EFF","#96410EFF",'grey80', 
                                '#FFF087FF',  'turquoise3', 'turquoise2', 'violetred3', 'violetred2', '#000000FF')

names(colors_metaprograms_ZFTA) <- c("Cycling", "Neuroepithelial-like", "Radial glia-like", 
                                       "Neuronal-like" ,"Ependymal-like", "Mesenchymal", "Unassigned", 
                                       "T-cell", "Myeloid", "Microglia", "Endothelial", "VLMCs",  "Oligodendrocyte")

col_sample = c('#e01e37',
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
               '#bc1f66')


# read metadata
metadata <- read_xlsx(file.path(base_dir, 'scripts/SampleIdentifier.xlsx'))
# reorder by sampleID
metadata <- metadata %>% arrange(SampleID)

metadata$SampleNameComplete <- paste0(metadata$SampleID, "_", metadata$SampleName)

FileName <- metadata$SampleNameComplete

# read  proportion of each cell type in each tumor
proportion <- read_csv(file.path(data_dir, 'metaprogram_frequency.csv'))
proportion <- proportion[, -which(names(proportion) == "...1")]
colnames(proportion)[2] <- 'metaprogram'
proportion$SampleNameComplete <- paste0(proportion$SampleID, "_", proportion$Xenium_Region)
proportion <- as.data.frame(proportion)


# read spatial coherence score data (from Carlos)
spatial_coherence <- qread(file.path(base_dir, 'analysis/order_disorder/data/results.qs'))

# extract global scores
average_spatial_coherence_df <- c()
average_spatial_coherence_program <- list()

    for (i in seq_along(FileName)) {
      average_spatial_coherence_df[i] <-spatial_coherence[[i]]$coherence_score
      average_spatial_coherence_program[[i]] <- spatial_coherence[[i]]$coherence_score_program
    }

names(average_spatial_coherence_df) = names(average_spatial_coherence_program) = names(spatial_coherence)

# reorganize and export average_spatial_coherence_df ------------------------------------------
average_spatial_coherence_df <- as.data.frame(average_spatial_coherence_df)
average_spatial_coherence_df$sample <- rownames(average_spatial_coherence_df)
colnames(average_spatial_coherence_df)[1] <- 'average_spatial_coherence'

 # calculate scaled score
 min <- min(average_spatial_coherence_df$average_spatial_coherence)
 max <- max(average_spatial_coherence_df$average_spatial_coherence)
 range <- max - min
 average_spatial_coherence_df$scaled_spatial_coherence <- (average_spatial_coherence_df$average_spatial_coherence - min)/range
 
 # calculate average score per sample
 average_spatial_coherence_df$SampleID <- metadata$SampleID
 average_spatial_coherence2 <- average_spatial_coherence_df %>% 
   dplyr::group_by(SampleID) %>%
   dplyr::summarise(mean = mean(scaled_spatial_coherence))

   # calculate scaled score
   min <- min(average_spatial_coherence2$mean)
   max <- max(average_spatial_coherence2$mean)
   range <- max - min
   average_spatial_coherence2$scaled_spatial_coherence <- (average_spatial_coherence2$mean - min)/range
 
 
 # export dataframe
 write_csv(as.data.frame(average_spatial_coherence_df), file.path(output_dir, 'Spatial_coherence_score_by_sample.csv'))
 
 # export dataframe
 write_csv(as.data.frame(average_spatial_coherence2), file.path(output_dir, 'Spatial_coherence_score_averge_by_sample.csv'))
 
 
 # add information on sample ID
 matching_indices <- match(average_spatial_coherence2$SampleID, metadata$SampleID)
 average_spatial_coherence2$Source <- metadata$Source[matching_indices]
 
 # reorganize dataframe with metaprogram proportion
 proportion2 <- proportion[ , c("SampleNameComplete", "metaprogram", "freq")]
 
 proportion2 <- proportion2 %>% 
   complete(nesting(SampleNameComplete, metaprogram), fill = list(freq = 0))  
 
 proportion3 <- cast(proportion2, SampleNameComplete~metaprogram, mean)
 
 # export dataframe
 write_csv(as.data.frame(proportion3), file.path(output_dir, 'Program_frequency.csv'))
 
 
# reorganize and export average_spatial_coherence_program (per program) ------------------------------------------
 average_spatial_coherence_program

 # add column with sample name
 for (i in seq_along(average_spatial_coherence_program)) {
   average_spatial_coherence_program[[i]]$SampleName <- names(average_spatial_coherence_program)[i]
 }
 
 # merge information
 average_spatial_coherence_program_df <- average_spatial_coherence_program[[1]]
 
 for (i in 2:length(average_spatial_coherence_program)) {
   average_spatial_coherence_program_df <- rbind(average_spatial_coherence_program_df, average_spatial_coherence_program[[i]])
 }
 
 # group all immune and all vascular cells
 average_spatial_coherence_program_df$metaprogram <- gsub("Microglia|T-cell|Myeloid", "Immune", average_spatial_coherence_program_df$metaprogram)
 average_spatial_coherence_program_df$metaprogram <- gsub("Endothelial|VLMCs", "Vascular", average_spatial_coherence_program_df$metaprogram)
 
 # remove oligos, radial glia (present only in 1 sample) and unassigned cells 
 average_spatial_coherence_program_df <- average_spatial_coherence_program_df[average_spatial_coherence_program_df$metaprogram != "Oligodendrocytes", ] 
 average_spatial_coherence_program_df <- average_spatial_coherence_program_df[average_spatial_coherence_program_df$metaprogram != "Unassigned", ] 
 average_spatial_coherence_program_df <- average_spatial_coherence_program_df[average_spatial_coherence_program_df$metaprogram != "Radial glia-like", ] 
 
 # calculate scaled score
 min <- min(average_spatial_coherence_program_df$average)
 max <- max(average_spatial_coherence_program_df$average)
 range <- max - min
 average_spatial_coherence_program_df$scaled_spatial_coherence <- (average_spatial_coherence_program_df$average - min)/range
 
 # add information on sample ID
 matching_indices <- match(average_spatial_coherence_program_df$SampleName, metadata$SampleNameComplete)
 average_spatial_coherence_program_df$SampleID <- metadata$SampleID[matching_indices]
 average_spatial_coherence_program_df$Source <- metadata$Source[matching_indices]

 
# Plot scores ------------------------------------------

ggplot(average_spatial_coherence_program_df, aes(x = reorder(metaprogram, scaled_spatial_coherence), y = scaled_spatial_coherence)) +
  #geom_violin(width = 1, color = 'black') + 
  geom_boxplot(fatten = NULL, outlier.shape = NA, width = 0.5, color = 'black') +
  geom_point(aes(group = SampleName, fill = SampleID), size = 2, shape = 21, stroke = 0, position = position_dodge(0.2)) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.4, size = 1, linetype = "solid") + 
  scale_fill_manual(values = col_sample) +
  labs(x = "Metaprogram", y = "Spatial coherence score") +
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
ggsave(file.path(output_dir, "1_spatial_coherence_by_metaprogram.pdf"), width=6, height=5)




# primary vs relapse
ggplot(average_spatial_coherence2, aes(x = Source, y = scaled_spatial_coherence)) +
  geom_violin(width = 1, color = 'black') + 
  #geom_boxplot(fatten = NULL, outlier.shape = NA, width = 0.7) + 
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.7, size = 1, linetype = "solid") + 
  geom_point(aes(group = SampleID, fill = SampleID), size = 3, shape = 21, stroke = 0, position = position_dodge(0.2)) +
  geom_line(aes(group = SampleID), position = position_dodge(0.2)) + 
  labs(x = "Source", y = "Spatial coherence score") +
  scale_fill_manual(values = col_sample) +
  theme_minimal() + theme(panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          plot.background = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"),
                          axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                          axis.text.y = element_text(size=12, colour="black"),
                          axis.title=element_text(size=12),
                          plot.title = element_text(size=10, face="bold")) + 
  ggpubr::stat_compare_means(method = "t.test", size = 4)
ggsave(file.path(output_dir, "3_spatial_coherence_by_source.pdf"), width=3.5, height=4)

