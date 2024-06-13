base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Project_HOPE/SARA"
metadata <- read_excel(file.path(base_dir, 'metadata/metadata.xlsx'))
FileName <- unique(metadata$FileName)

colors_tumor_normal <- c('#d8b365', '#5ab4ac')

colors_samples <- as.vector(paletteer::paletteer_c("scico::roma", n = length(FileName)))
names(colors_samples) <- FileName

col_signatures <- as.vector(rev(paletteer::paletteer_d("ghibli::LaputaMedium")))[1:6]
names(col_signatures) <- c('OPC-like', 'MES-like', 'MES-like-AC-like', 'AC-like', 'Cycling', 'NPC-like')

color_populations_immune <- paletteer::paletteer_d("MoMAColors::VanGogh")[2:6]

color_cell_type_all_malignant <-  paletteer::paletteer_d("MoMAColors::VanGogh")

color_cell_type <- c(col_signatures, color_populations_immune)

col_cycling_prop <- c('grey20', 'grey90')
names(col_cycling_prop) <- c('Cycling', 'Non-cycling')

col_gender <- c('F' = '#F08080FF', 'M' = '#87CEEBFF', 'NA' = 'grey90')
col_source_simplified <- c('Primary' = '#FBE4C6FF', 'Recurrence' = '#CB74ADFF', 'NA' = 'grey90')
col_scRNAseq <- c('10X'= 'grey60', 'Smart-seq2' = 'grey30')
col_diagnosis <- c('HGG' = '#009474FF', 'DMG' = '#EFDDCFFF', 'GBM' = '#DCBE9BFF', "protoplasmatic Astroblastoma" = '#11C2B5FF',
                   "Fibrillary astrocytoma" = '#40E0D0FF', 'H3 WT' = '#B0986CFF')
col_location <- c("Hemispheric" = '#79AD41FF', "Posterior Fossa" = '#34B6C6FF', "Midline" = '#4063A3FF')


pal_a = c('#370617',
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


col_patientID <- pal_a[1:length(unique(metadata$PatientID))]
names(col_patientID) <- unique(metadata$PatientID)

