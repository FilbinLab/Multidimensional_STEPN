col_subtype = c('ZFTA-RELA' = "#B44672",'ZFTA-Cluster 1' = "#B47846", 'ZFTA-Cluster 2' = "#46B478",
                'ZFTA-Cluster 3' = "#46B4AF", 'ZFTA-Cluster 4' = "#4682B4",'ST-YAP1' = "#B4AF46")

colors_metaprograms_Xenium <- c("gray50","#F99E93FF","#9E5E9BFF","#74ADD1FF",'#0F4F8B', "#ACD39EFF","#96410EFF", "#96410EFF", 'mistyrose1',
                                'grey90', 'grey90',
                                '#FFF087FF',  '#F47942FF', 'violetred3', '#AECEBFFF', '#30A0A4FF', '#55A19EFF')

names(colors_metaprograms_Xenium) <- c("Cycling", "Neuroepithelial-like", "Radial-glia-like", 
                                       "Embryonic-neuronal-like", "Neuronal-like" ,"Ependymal-like", "MES-like", "MES/Hypoxia", "Embryonic-like",
                                       "Unknown", 'other',
                                       "T-cells", "Myeloid",  "Endothelial",  "Oligodendrocytes", 'Astrocyte', 'Neurons')

colors_metaprograms_Xenium_dot_separator <- colors_metaprograms_Xenium
names(colors_metaprograms_Xenium_dot_separator) <- gsub('-', '.', names(colors_metaprograms_Xenium_dot_separator))

colors_metaprograms_Xenium_underscore_separator <- colors_metaprograms_Xenium
names(colors_metaprograms_Xenium_underscore_separator) <- gsub('-', '_', names(colors_metaprograms_Xenium_underscore_separator))

# mes_coherence_category (the two broad categories in which mes-like cells' coherence belonged) colors
colors_mes_coh_cat <- c('#E1C16E', '#7B3F00')
names(colors_mes_coh_cat) <- c('0', '1')

# colors_niches <- as.vector(paletteer::paletteer_d("colRoz::c_decresii"))
colors_niches <- c("#DCA761FF", "#C6C16DFF", "#8B9C94FF", "#628CA5FF", 
                   "#5A6C7AFF", "#514F5CFF", "#f6cacc", "#f7c297", "#ffecb8", "#BDCDFF", '#F47942FF', "#1C8356")
names(colors_niches) <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12')

colors_niches_mps <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "white")
# colors_niches_mps <- c('#B57FA4', '#35B800', '#E84746', '#52A8BC', '#370617', '#EFE185')
names(colors_niches_mps) <- c('A', 'B', 'C', 'D', 'E', 'F', 'NA')

# colors for the spatial clusters obtained using cellcharter
colors_spatial_clusters <- c('0' = "#5A5156", '1' = "#E4E1E3", '2' = "#FE00FA", '3' = "#3283FE", '4' = "#FEAF16",
                             '5' = "#B00068", '6' = "#90AD1C", '7' = "#AA0DFE", '8' = "#F8A19F", '9' = "#325A9B",
                             '10' = "#C4451C", '11' = "#1C8356", '12' = "#85660D", '13' = "#FBE426", '14' = "#1CBE4F",
                             '15' = "#FA0087", '16' = "#F7E1A0", '17' = "#C075A6", '18' = "#782AB6", '19' = "#AAF400",
                             '20' = "#BDCDFF", '21' = "#822E1C", '22' = "#B5EFB5", '23' = "#7ED7D1", '24' = "#1C7F93",
                             '25' = "#66B0FF")

# code to get a certain number of equally spaced colors
# library(scales)
# show_col(hue_pal()(100))
# show_col(pal_brewer("div")(5))
# library(Polychrome)
# show_col(dark.colors(n = 24))
# show_col(kelly.colors(n = 22))
# show_col(alphabet.colors(n = 26))
# show_col(sky.colors(n = 24))
# show_col(palette36.colors(n = 36))
# show_col(green.armytage.colors(n = 26))
# show_col(glasbey.colors(n = 32))

# Palette_Function	Max Colors	Intended Use
# kelly.colors(n = 22)	22	Kelly’s 22-color palette
# glasbey.colors(n = 32)	32	Glasbey’s 32-color palette
# green.armytage.colors(n = 26)	26	Armytage’s green-based palette
# palette36.colors(n = 36)	36	A 36-color palette
# alphabet.colors(n = 26)	26	Alphabet palette (letters a–z)
# light.colors(n = 24)	24	24 lighter shades
# dark.colors(n = 24)	24	24 darker shades
# sky.colors(n = 24)	24	Sky-inspired pastels

col_sampling <- c('Primary' = '#f7c297', 'Recurrence' = '#ffecb8', 'NA' = 'grey90')

col_patientID <-  c('#370617',
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

# all_vars <- SetUpEpendymomaGlobalVarsGeneral()
# metadata <- all_vars$metadata
# metadata = metadata %>% filter(Subtype != 'ZFTA-RELA')
# sample_ids = metadata$SampleID
# sort(unique(sample_ids))

col_sampleid_zfta <- c("STEPN-01" = '#95191d', 
                       "STEPN-06" = '#16396c', 
                       "STEPN-10" = '#0daa4f', 
                       "STEPN-12" = '#762213', 
                       "STEPN-14" = '#da8e29', 
                       "STEPN-15" = '#0f8745', 
                       "STEPN-16" = '#008865', 
                       "STEPN-17" = '#9f1b1e', 
                       "STEPN-19" = '#26813e', 
                       "STEPN-45" = '#025786', 
                       "STEPN-46" = '#801215', 
                       "STEPN-47" = '#a21c24', 
                       "STEPN-48" = '#0d3d69', 
                       "STEPN-49" = '#93181c', 
                       "STEPN-50" = '#0f6c35', 
                       "STEPN-51" = '#a21c29', 
                       "STEPN-52" = '#871518', 
                       "STEPN-53" = '#b82690', 
                       "STEPN-54" = '#009383', 
                       "STEPN-55" = '#009a4f', 
                       "STEPN-56" = '#62893c', 
                       "STEPN-57" = '#b94d25', 
                       "STEPN-58" = '#066d36')

col_sampleid_all <- c("STEPN-01" = '#95191d', 
                       "STEPN-02" = '#2E91E5',
                         "STEPN-03" = "#E15F99",
                         "STEPN-04" = "#1CA71C",
                         "STEPN-05" = "#FB0D0D",
                       "STEPN-06" = '#16396c', 
                       "STEPN-07" = "#DA16FF",
                         "STEPN-08" = "#222A2A",
                         "STEPN-09" = "#B68100",
                       "STEPN-10" = '#0daa4f',
                       "STEPN-11" = "#750D86",
                       "STEPN-12" = '#762213',
                       "STEPN-13" = "#EB663B",
                       "STEPN-14" = '#da8e29', 
                       "STEPN-15" = '#0f8745', 
                       "STEPN-16" = '#008865', 
                       "STEPN-17" = '#9f1b1e', 
                       "STEPN-19" = '#26813e',
                       "STEPN-20" = "#511CFB",
                       "STEPN-45" = '#025786', 
                       "STEPN-46" = '#801215', 
                       "STEPN-47" = '#a21c24', 
                       "STEPN-48" = '#0d3d69', 
                       "STEPN-49" = '#93181c', 
                       "STEPN-50" = '#0f6c35', 
                       "STEPN-51" = '#a21c29', 
                       "STEPN-52" = '#871518', 
                       "STEPN-53" = '#b82690', 
                       "STEPN-54" = '#009383', 
                       "STEPN-55" = '#009a4f', 
                       "STEPN-56" = '#62893c', 
                       "STEPN-57" = '#b94d25', 
                       "STEPN-58" = '#066d36',
                       "STEPN-21" = "#00A08B",
                       "STEPN-22" = '#e01e37',
                         "STEPN-23" = "#FC0080",
                         "STEPN-24" = '#f6cacc',
                         "STEPN-18" = '#370617',
                         "STEPN-25" = "#778AAE",
                         "STEPN-26" = "#862A16",
                         "STEPN-27" = "#A777F1",
                         "STEPN-28" = "#620042",
                         "STEPN-59" = '#a564d3',
                         "STEPN-29" = "#DA60CA",
                         "STEPN-30" = "#6C4516",
                         "STEPN-31" = "#0D2A63",
                         "STEPN-32" = "#AF0038",
                         "STEPN-33" = '#55a630',
                         "STEPN-34" = '#8be8d7',
                         "STEPN-35" = "#875692",
                         "STEPN-36" = "#f38400",
                         "STEPN-37" = "#a1caf1",
                         "STEPN-38" = "#be0032",
                         "STEPN-39" = '#2fb5c7',
                         "STEPN-40" = "#848482",
                         "STEPN-60" = '#ff5ca5',
                         "STEPN-61" = '#ffb9d8',
                         "STEPN-41" = "#0067a5",
                         "STEPN-42" = "#f99379",
                         "STEPN-43" = '#0377a8',
                         "STEPN-44" = '#002855'
                         )

col_sampleid_noncanonical = c("STEPN-18" = '#370617', 
                              "STEPN-22" = '#e01e37', 
                              "STEPN-24" = '#f6cacc', 
                              "STEPN-33" = '#55a630', 
                              "STEPN-34" = '#8be8d7', 
                              "STEPN-39" = '#2fb5c7', 
                              "STEPN-43" = '#0377a8', 
                              "STEPN-44" = '#002855', 
                              "STEPN-59" = '#a564d3', 
                              "STEPN-60" = '#ff5ca5', 
                              "STEPN-61" = '#ffb9d8')

# the colors for the metaniches
# cluster_colors <- c(
#   1 = '#addbc7',
#   2 = '#858bc4',
#   3 = '#ebdaea',
#   4 = '#dce57e',
#   5 = '#f9cb85',
#   6 = '#a1bfe4',
# )


# create vector for color density
color_coherence_density <- c(
  'white', '#FFFFD9FF', "#E0F3DBFF", "#CCEBC5FF", "#A8DDB5FF", 
  "#7BCCC4FF", "#4EB3D3FF", "#2B8CBEFF", "#0868ACFF", "#084081FF"
)
names(color_coherence_density) <- 0:8

############  continuous palette ############
color_coherence_density_continuous <- c(
  'white', '#FFFFD9FF', "#E0F3DBFF", "#CCEBC5FF", "#A8DDB5FF", 
  "#7BCCC4FF", "#4EB3D3FF", "#2B8CBEFF", "#0868ACFF", "#084081FF"
)
# names(color_coherence_density_continuous) <- seq(from = 0, to = 1, length = length(color_coherence_density_continuous))

# get a function which takes a number as input and interpolates the provided colors vector into those many intermediate colors.
# this function can be used as input to scale_color_gradientn() for the colors argument
# pal_fn <- colorRampPalette(color_coherence_density_continuous)

# color_coherence_density_continuous_palette <- scale_color_gradientn(colors = pal_fn(100))

# seq(from = 0, to = 1, length = 100)

############ discrete palette but without white ############

color_coherence_density_without_white <- c(
  '#FFFFD9FF', "#E0F3DBFF", "#CCEBC5FF", "#A8DDB5FF", 
  "#7BCCC4FF", "#4EB3D3FF", "#2B8CBEFF", "#0868ACFF", "#084081FF"
)

### color palette for dotplot

color_dotplot <- c('#d5eadb', '#67c3ac', '#30acac', '#318fa8', '#3180a4', '#3a5d9d', '#403e7d', '#36264e', '#372854')









