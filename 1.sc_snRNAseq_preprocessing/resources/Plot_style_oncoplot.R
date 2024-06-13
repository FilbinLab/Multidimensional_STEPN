base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Project_HOPE/SARA"

# Read metadata file with patient information -----------------------------------
oncoprint_input_tumors <- read_excel(file.path(base_dir, 'metadata/metadata_10X_smartseq.xlsx')) 

# remove samples with no goodQC cells left postQC
badqc <- c('841169', 'HBAD-D2')
oncoprint_input_tumors <- oncoprint_input_tumors[!oncoprint_input_tumors$FileName %in% badqc,]

# remove samples without primary/relapse information (from 10X cohort)
remove <- c(NA)
oncoprint_input_tumors <- oncoprint_input_tumors[!oncoprint_input_tumors$Source_simplified %in% remove,]


# color palettes
pal_patientID = c('#370617',
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

# Order dataframe -----------------------------------

#order dataframe based on patient subtypes
oncoprint_input_tumors <- oncoprint_input_tumors %>%
  arrange(desc(Technology), Location_new, PatientID, Source) 
oncoprint_input_tumors

# lock in factor level name to make sure patients are plotted in same order as dataframe
oncoprint_input_tumors$SampleName <- factor(oncoprint_input_tumors$SampleName, levels = oncoprint_input_tumors$SampleName)

# transform age into numeric
oncoprint_input_tumors$Age <- as.numeric(oncoprint_input_tumors$Age)

# Color codes  -----------------------------------
col_gender <- c('F' = '#F08080FF', 'M' = '#87CEEBFF', 'NA' = 'grey90')
col_source_simplified <- c('Primary' = '#CB74ADFF', 'Recurrence' = '#FBE4C6FF', 'NA' = 'grey90')
col_scRNAseq <- c('10X'= 'grey60', 'Smart-seq2' = 'grey30')
col_diagnosis <- c('HGG' = '#009474FF', 'DMG' = '#EFDDCFFF', 'GBM' = '#DCBE9BFF', "protoplasmatic Astroblastoma" = '#11C2B5FF',
                   "Fibrillary astrocytoma" = '#40E0D0FF', 'H3 WT' = '#B0986CFF')
col_location <- c("Hemispheric" = '#79AD41FF', "Posterior Fossa" = '#34B6C6FF', "Midline" = '#4063A3FF')

col_patientID <- pal_patientID[1:38]
names(col_patientID) <- unique(oncoprint_input_tumors$PatientID, 'NA' = 'grey90')


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
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 12, face = 'bold'),
                               plot.margin = unit(c(0, 0, 0, 0), "cm")) 
