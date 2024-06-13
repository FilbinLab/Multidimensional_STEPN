# Load packages -----------------------------------
rm(list = ls())

library(here)
library(tidyr)
library(dplyr)
library(readr)
library(glue)
library(tibble)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(ggpubr)

# Organize environment  -----------------------------------
base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/11 - Plots scRNAseq"
resources_dir <- file.path(base_dir, 'scripts/resources')

source(file.path(resources_dir, 'Plot_style_oncoplot.R'))

# create data directory for storing output
analysis_dir <- file.path(base_dir, 'analysis')
plot_dir <- file.path(analysis_dir, 'plots/oncoplot')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

# Read metadata file with patient information -----------------------------------
oncoprint_input_tumors <- read_excel(file.path(base_dir, 'data/metadata/metadata_ST_EPN_Patient.xlsx')) 

# Color codes  -----------------------------------
col_gender <- c('F' = 'lightpink2', 'M' = '#CBD5E8', 'NA' = 'grey90')
col_sampling <- c('Primary' = '#f7c297', 'Recurrence' = '#ffecb8', 'NA' = 'grey90')
col_spatial <- c('Fresh frozen in-situ'= 'grey60', 'FFPE in-situ' = 'black')
col_scRNAseq <- c('scRNA-seq'= 'lightskyblue2', 'snRNA-seq' = 'mediumpurple1')
col_fusion <- c(paletteer::paletteer_d("miscpalettes::pastel")[c(5:11, 4)])
names(col_fusion) <- unique(oncoprint_input_tumors$Fusion)
col_subtypes = c('ZFTA-RELA' = "#B44672",'ZFTA-Cluster 1' = "#B47846", 'ZFTA-Cluster 2' = "#46B478",
                 'ZFTA-Cluster 3' = "#46B4AF", 'ZFTA-Cluster 4' = "#4682B4",'ST-YAP1' = "#B4AF46")

col_patientID <- c('#370617',
                   '#e01e37',
                   '#f6cacc',
                   '#ff0000',
                   '#da5552',
                   '#fec89a',
                   '#ffd7ba',
                   '#ffba08',
                   '#f48c06',
                   '#ffea00',
                   '#ffa200',
                   '#ff6000',
                   '#eeef20',
                   '#d4d700',
                   '#aacc00',
                   '#55a630',
                   '#2b9348',
                   '#007f5f',
                   '#b9fbc0',
                   '#10451d',
                   '#27a300',
                   '#b4fadc',
                   '#8be8d7',
                   '#63d4cc',
                   '#51ccd1',
                   '#2fb5c7',
                   '#0377a8',
                   '#0466c8',
                   '#023e7d',
                   '#002855',
                   '#a564d3',
                   '#bf99f2',
                   '#60308c',
                   '#511f73',
                   '#ff0072',
                   '#ff2e8c',
                   '#ff5ca5',
                   '#ffa2cb',
                   '#ffb9d8',
                   '#bc1f66',
                   '#e41b60',
                   '#dcb9a1',
                   '#ebcfbc',
                   '#d0a03b')

col_patientID <- sample(col_patientID)
names(col_patientID) <- unique(oncoprint_input_tumors$PatientID)


# (1) Oncoplot patients -----------------------------------

# Read metadata file with patient information -----------------------------------
oncoprint_input_tumors <- read_excel(file.path(base_dir, 'data/metadata/metadata_ST_EPN_Patient.xlsx')) 

#order dataframe based on patient subtypes
desired_order <- c('ZFTA-RELA' ,'ZFTA-Cluster 1' , 'ZFTA-Cluster 2', 'ZFTA-Cluster 3', 'ZFTA-Cluster 4', 'ST-YAP1')
oncoprint_input_tumors$Subtype <- factor(oncoprint_input_tumors$Subtype, levels = desired_order)

oncoprint_input_tumors <- oncoprint_input_tumors %>%
  arrange(Subtype, PatientID, Sample_deID) 
oncoprint_input_tumors


# lock in factor level name to make sure patients are plotted in same order as dataframe
oncoprint_input_tumors$Sample_deID <- factor(oncoprint_input_tumors$Sample_deID, levels = oncoprint_input_tumors$Sample_deID)


# Create bar plots with patient information -----------------------------------

# create individual bar plots
p1 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = Subtype), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_subtypes , na.value = "grey95", guide = guide_legend(ncol = 2))  + 
  theme_oncoplot_no_legend + 
  labs(y = 'Subtype')

p1_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = Subtype), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_subtypes , na.value = "grey95", guide = guide_legend(ncol = 2))  + 
  theme_oncoplot_legend + 
  labs(y = 'Subtype')
legend_p1 <- get_legend(p1_leg)

p1b <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = as.numeric(Score_Subclass)), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"), na.value = "gray90") + 
  theme_oncoplot_no_legend + 
  labs(y = 'Methylation score')

p1b_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = as.numeric(Score_Subclass)), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"), na.value = "gray90") + 
  theme_oncoplot_legend + 
  labs(y = 'Methylation score')
legend_p1b <- get_legend(p1b_leg)

p2 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = Gender), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_gender, na.value = "grey95", guide = guide_legend(ncol = 1))  +
  theme_oncoplot_no_legend + 
  labs(y = 'Gender')

p2_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = Gender), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_gender, na.value = "grey95", guide = guide_legend(ncol = 1))  +
  theme_oncoplot_legend + 
  labs(y = 'Gender')
legend_p2 <- get_legend(p2_leg)


p3 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = Age),  width = 0.9, height = 0.9) +
  geom_text(aes(label = Age), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c(rep('white', 18)), na.value = "white")  +
  theme_oncoplot_no_legend + 
  labs(y = 'Age')

p3_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = Age),  width = 0.9, height = 0.9) +
  geom_text(aes(label = Age), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c(rep('white', 18)), na.value = "white")  +
  theme_oncoplot_legend + 
  labs(y = 'Age')
legend_p3 <- get_legend(p3_leg)


p4 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = Sampling), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_sampling, na.value = "grey95", guide = guide_legend(ncol = 1)) +  
  theme_oncoplot_no_legend + 
  labs(y = 'Sampling')

p4_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = Sampling), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_sampling, na.value = "grey95", guide = guide_legend(ncol = 1)) +  
  theme_oncoplot_legend + 
  labs(y = 'Sampling')
legend_p4 <- get_legend(p4_leg)


p5 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = Fusion), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_fusion, na.value = "grey95", guide = guide_legend(ncol = 2)) + 
  theme_oncoplot_no_legend + 
  labs(y = 'Fusion')

p5_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = Fusion), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_fusion, na.value = "grey95", guide = guide_legend(ncol = 2)) + 
  theme_oncoplot_legend + 
  labs(y = 'Fusion')
legend_p5 <- get_legend(p5_leg)

p6 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill =  `sc/snRNAseq`), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_scRNAseq, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_no_legend + 
  labs(y = 'sc/snRNA-seq')

p6_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill =  `sc/snRNAseq`), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_scRNAseq, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_legend + 
  labs(y = 'sc/snRNA-seq')
legend_p6 <- get_legend(p6_leg)

p7 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = Spatial), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_spatial, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_no_legend + 
  labs(y = 'Spatial') +
  theme(axis.text.x = element_text(color = "black", size = 11, angle = 90, hjust = 1, vjust = 0.5))


p7_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = Spatial), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_spatial, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_legend + 
  labs(y = 'Spatial') 
legend_p7 <- get_legend(p7_leg)

p8 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = PatientID), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_patientID, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_no_legend + 
  labs(y = 'Patient ID') #+
  #theme(axis.text.x = element_text(color = "black", size = 11, angle = 90, hjust = 1, vjust = 0.5))

p8_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample_deID, y = 1)) +
  geom_tile(aes(fill = PatientID), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_patientID, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_legend + 
  labs(y = 'Patient ID') 
legend_p8 <- get_legend(p8_leg)



# Print oncoplot -----------------------------------
plot_grid(p1, p1b, p2, p3, p4, p5, p6, p7, p8,
          align = "v", axis = "l", ncol = 1,
          rel_heights = c(rep(0.2, 3), 0.18, rep(0.2, 4), 0.87))
ggsave(file.path(plot_dir, "Oncoprint.pdf"), width=10, height=2.9, dpi=300)

plot_grid(p1, p1b, p2, p3, p4, p5, p6, p7, p8,
          align = "v", axis = "l", ncol = 1,
          rel_heights = c(rep(0.2, 3), 0.18, rep(0.2, 3), 0.87, 0.2))
ggsave(file.path(plot_dir, "Oncoprint.pdf"), width=10, height=2.9, dpi=300)


# print legends
plot_grid(as_ggplot(legend_p1),
          as_ggplot(legend_p1b),
          as_ggplot(legend_p2),
          as_ggplot(legend_p4),
          as_ggplot(legend_p5),
          as_ggplot(legend_p6),
          as_ggplot(legend_p7),
          axis = 'tb',
          ncol = 7,
          rel_widths = c(0.4, 0.1, 0.15, 0.15, 0.4, 0.2, 0.2))
ggsave(file.path(plot_dir, "Legend_Oncoprint.pdf"), width=18, height=1.5, dpi=300)

  






# (2) Oncoplot Cell lines -----------------------------------

# Read metadata file with patient information -----------------------------------
oncoprint_input_tumors <- read_excel(file.path(base_dir, 'data/metadata/metadata_ST_EPN_CL.xlsx')) 

#order dataframe based on patient subtypes
oncoprint_input_tumors <- oncoprint_input_tumors %>%
  arrange(Sample) 
oncoprint_input_tumors

#transform age into character variable (otherwise issues with plotting)
oncoprint_input_tumors$Age <- as.character(oncoprint_input_tumors$Age)

# lock in factor level name to make sure patients are plotted in same order as dataframe
oncoprint_input_tumors$Sample <- factor(oncoprint_input_tumors$Sample, levels = oncoprint_input_tumors$Sample)


# Create bar plots with cell line information -----------------------------------

# create individual bar plots
p1 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Subtype), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_subtypes , na.value = "grey95", guide = guide_legend(ncol = 2))  + 
  theme_oncoplot_no_legend + 
  theme(axis.title.y = element_blank())

p1_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Subtype), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_subtypes , na.value = "grey95", guide = guide_legend(ncol = 2))  + 
  theme_oncoplot_legend 
legend_p1 <- get_legend(p1_leg)

p1b <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = as.numeric(Score_Subclass)), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"), na.value = "gray90") + 
  theme_oncoplot_no_legend + 
  theme(axis.title.y = element_blank())

p1b_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = as.numeric(Score_Subclass)), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"), na.value = "gray90") + 
  theme_oncoplot_legend 
legend_p1b <- get_legend(p1b_leg)


p2 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Gender), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_gender, na.value = "grey95", guide = guide_legend(ncol = 1))  +
  theme_oncoplot_no_legend + 
  theme(axis.title.y = element_blank())

p2_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Gender), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_gender, na.value = "grey95", guide = guide_legend(ncol = 1))  +
  theme_oncoplot_legend 
legend_p2 <- get_legend(p2_leg)


p3 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Age),  width = 0.9, height = 0.9) +
  geom_text(aes(label = Age), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = col_gender, na.value = "white")  +
  theme_oncoplot_no_legend + 
  theme(axis.title.y = element_blank())

p3_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Age),  width = 0.9, height = 0.9) +
  geom_text(aes(label = Age), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = col_gender, na.value = "white")  +
  theme_oncoplot_legend 
legend_p3 <- get_legend(p3_leg)


p4 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Sampling), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_sampling, na.value = "grey95", guide = guide_legend(ncol = 1)) +  
  theme_oncoplot_no_legend + 
  theme(axis.title.y = element_blank())

p4_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Sampling), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_sampling, na.value = "grey95", guide = guide_legend(ncol = 1)) +  
  theme_oncoplot_legend 
legend_p4 <- get_legend(p4_leg)


p5 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Fusion), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_fusion, na.value = "grey95", guide = guide_legend(ncol = 2)) + 
  theme_oncoplot_no_legend + 
  theme(axis.title.y = element_blank())

p5_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Fusion), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_fusion, na.value = "grey95", guide = guide_legend(ncol = 2)) + 
  theme_oncoplot_legend  
legend_p5 <- get_legend(p5_leg)


p6 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill =  `sc/snRNAseq`), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_scRNAseq, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_no_legend  + 
  theme(axis.title.y = element_blank())

p6_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill =  `sc/snRNAseq`), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_scRNAseq, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_legend 
legend_p6 <- get_legend(p6_leg)

p7 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Spatial), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_spatial, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_no_legend  + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black", size = 11, angle = 90, hjust = 1, vjust = 0.5))


p7_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Spatial), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_spatial, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_legend 
legend_p7 <- get_legend(p7_leg)

p8 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = PatientID), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = c('grey10', 'grey50', 'grey75'), na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_no_legend +
  theme(axis.title.y = element_blank()) #+
#theme(axis.text.x = element_text(color = "black", size = 11, angle = 90, hjust = 1, vjust = 0.5))

p8_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = PatientID), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values =  c('grey10', 'grey50', 'grey75'), na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_legend   + 
  theme(axis.title.y = element_blank())
legend_p8 <- get_legend(p8_leg)




# Print oncoplot -----------------------------------
plot_grid(p1, p1b, p2, p3, p4, p5, p6, p7, p8,
          align = "v", axis = "l", ncol = 1,
          rel_heights = c(rep(0.2, 3), 0.15, rep(0.2, 3), 0.65, 0.2))
ggsave(file.path(plot_dir, "Oncoprint_CL.pdf"), width=0.8, height=3, dpi=300)


# (3) Oncoplot PDXs -----------------------------------

# Read metadata file with patient information -----------------------------------
oncoprint_input_tumors <- read_excel(file.path(base_dir, 'data/metadata/metadata_ST_EPN_PDX.xlsx')) 

#order dataframe based on patient subtypes
oncoprint_input_tumors <- oncoprint_input_tumors %>%
  arrange(Sample) 
oncoprint_input_tumors

#transform age into character variable (otherwise issues with plotting)
oncoprint_input_tumors$Age <- as.character(oncoprint_input_tumors$Age)

# lock in factor level name to make sure patients are plotted in same order as dataframe
oncoprint_input_tumors$Sample <- factor(oncoprint_input_tumors$Sample, levels = oncoprint_input_tumors$Sample)


# Create bar plots with PDX information -----------------------------------

# create individual bar plots
p1 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Subtype), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_subtypes , na.value = "grey95", guide = guide_legend(ncol = 2))  + 
  theme_oncoplot_no_legend + 
  theme(axis.title.y = element_blank())

p1_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Subtype), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_subtypes , na.value = "grey95", guide = guide_legend(ncol = 2))  + 
  theme_oncoplot_legend 
legend_p1 <- get_legend(p1_leg)


p1b <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = as.numeric(Score_Subclass)), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"), na.value = "gray90") + 
  theme_oncoplot_no_legend + 
  theme(axis.title.y = element_blank())

p1b_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = as.numeric(Score_Subclass)), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"), na.value = "gray90") + 
  theme_oncoplot_legend 
legend_p1b <- get_legend(p1b_leg)



p2 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Gender), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_gender, na.value = "grey95", guide = guide_legend(ncol = 1))  +
  theme_oncoplot_no_legend + 
  theme(axis.title.y = element_blank())

p2_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Gender), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_gender, na.value = "grey95", guide = guide_legend(ncol = 1))  +
  theme_oncoplot_legend 
legend_p2 <- get_legend(p2_leg)


p3 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Age),  width = 0.9, height = 0.9) +
  geom_text(aes(label = Age), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = col_gender, na.value = "white")  +
  theme_oncoplot_no_legend + 
  theme(axis.title.y = element_blank())

p3_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Age),  width = 0.9, height = 0.9) +
  geom_text(aes(label = Age), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = col_gender, na.value = "white")  +
  theme_oncoplot_legend 
legend_p3 <- get_legend(p3_leg)


p4 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Sampling), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_sampling, na.value = "grey95", guide = guide_legend(ncol = 1)) +  
  theme_oncoplot_no_legend + 
  theme(axis.title.y = element_blank())

p4_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Sampling), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_sampling, na.value = "grey95", guide = guide_legend(ncol = 1)) +  
  theme_oncoplot_legend 
legend_p4 <- get_legend(p4_leg)


p5 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Fusion), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_fusion, na.value = "grey95", guide = guide_legend(ncol = 2)) + 
  theme_oncoplot_no_legend + 
  theme(axis.title.y = element_blank())

p5_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Fusion), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_fusion, na.value = "grey95", guide = guide_legend(ncol = 2)) + 
  theme_oncoplot_legend  
legend_p5 <- get_legend(p5_leg)


p6 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill =  `sc/snRNAseq`), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_scRNAseq, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_no_legend  + 
  theme(axis.title.y = element_blank())

p6_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill =  `sc/snRNAseq`), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_scRNAseq, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_legend 
legend_p6 <- get_legend(p6_leg)

p7 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Spatial), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_spatial, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_no_legend  + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color = "black", size = 11, angle = 90, hjust = 1, vjust = 0.5))


p7_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = Spatial), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = col_spatial, na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_legend 
legend_p7 <- get_legend(p7_leg)

p8 <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = PatientID), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = c('grey10', 'grey50', 'grey75'), na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_no_legend +
  theme(axis.title.y = element_blank()) #+
#theme(axis.text.x = element_text(color = "black", size = 11, angle = 90, hjust = 1, vjust = 0.5))

p8_leg <- oncoprint_input_tumors %>%
  ggplot(aes(x = Sample, y = 1)) +
  geom_tile(aes(fill = PatientID), colour = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values =  c('grey10', 'grey50', 'grey75'), na.value = "grey95", guide = guide_legend(ncol = 1)) + 
  theme_oncoplot_legend   + 
  theme(axis.title.y = element_blank())
legend_p8 <- get_legend(p8_leg)





# Print oncoplot -----------------------------------
plot_grid(p1,  p1b, p2, p3, p4, p5, p6, p7, p8,
          align = "v", axis = "l", ncol = 1,
          rel_heights = c(rep(0.2, 3), 0.15, rep(0.2, 3), 0.65, 0.2))
  ggsave(file.path(plot_dir, "Oncoprint_PDXs.pdf"), width=0.8, height=3, dpi=300)








