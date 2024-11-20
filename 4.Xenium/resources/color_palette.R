library(paletteer)

colors_metaprograms_Xenium <- c("gray50","#F99E93FF","#9E5E9BFF","#74ADD1FF",'#0F4F8B', "#ACD39EFF","#96410EFF", 'mistyrose1',
                                'grey90', 
                                '#FFF087FF',  '#F47942FF', 'violetred3', '#AECEBFFF', '#30A0A4FF', '#55A19EFF')

names(colors_metaprograms_Xenium) <- c("Cycling", "Neuroepithelial-like", "Radial-glia-like", 
                                       "Embryonic-neuronal-like", "Neuronal-like" ,"Ependymal-like", "MES-like", "Embryonic-like",
                                       "Unknown", 
                                       "T-cells", "Myeloid",  "Endothelial",  "Oligodendrocytes", 'Astrocyte', 'Neurons')


colors_niches <- paletteer::paletteer_d("colRoz::c_decresii")
