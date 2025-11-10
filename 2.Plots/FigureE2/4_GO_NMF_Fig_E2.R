# Load packages -----------------------------------
library(Seurat)
library(tidyverse)
library(enrichR)
library(patchwork)
library(glue)
library(Matrix)
library(DropletUtils)
library(SeuratWrappers)
library(dorothea)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(qs)
library(enrichplot)
library(writexl)
library(paletteer)


# Organize environment  -----------------------------------

base_dir <- "/Users/sdaniell/Partners HealthCare Dropbox/Sara Danielli/Filbin lab - sample location/Filbin Lab Active/Sara/Project/Ependymoma/11 - Plots scRNAseq"

genelist_dir <- file.path(base_dir, "analysis/plots/STEPN/malignant")

plot_dir <- file.path(base_dir, 'analysis/plots/GO_NMF')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

resource_dir <- file.path(base_dir, 'scripts/resources')

source(file.path(resource_dir,'Plotting_helper_functions.R'))


# load seurat obj  -----------------------------------
seurat_obj <- qread(file = file.path(base_dir, "data/patient/seurat_obj_malignant_annotated2.qs"))

# load NMF gene list -----------------------------------
marker_list <- readRDS(file.path(base_dir, 'data/signatures/nmf_marker_genes_final_annotated.rds'))
names(marker_list) <- c("Neuroepithelial-like", "MES-Hypoxia","Ependymal-like","Radial-glia-like","Cycling", "Embryonic-neuronal-like", "Neuronal-like", "Embryonic-like")

# select top 30
marker_list <- lapply(marker_list, head, 30)
marker_list



# Run GO with NMF marker genes (top30) -----------------------------------
plotList <- list()
go_table <- list()
gse_list <- list()
programs <- names(marker_list)

for (program in programs) {
  tryCatch({
    markers <- marker_list[[program]]
    
    gse_list[[program]] <- enrichGO(gene = markers, 
                                    ont = "ALL", 
                                    keyType = "SYMBOL", 
                                    minGSSize = 3, 
                                    maxGSSize = 800,
                                    pvalueCutoff = 0.05, 
                                    OrgDb = org.Hs.eg.db, 
                                    pAdjustMethod = "BH")
    gse_list[[program]]@result <- gse_list[[program]]@result %>% arrange(p.adjust)
    plotList[[program]] <- dotplot(gse_list[[program]], showCategory = 10) + ggtitle(program) + theme(plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1))
    go_table[[program]] <- gse_list[[program]]@result
  }, error = function(e) {
    cat("Error occurred for program:", program, "\n")
    cat("Error message:", conditionMessage(e), "\n")
    # You can choose to skip this program or handle the error differently
  })
}

# export  -----------------------------------

# export plots
pt <- wrap_plots(plotList[c(1:length(plotList))])
ggsave(plot = pt, file.path(plot_dir, '1_enrich_GO_top30_NMF_markers.pdf'), width = 20, height = 15)

# export df
write_xlsx(go_table, file.path(plot_dir, '2_enrich_GO_top30_NMF_markers_source_data.xlsx'))


# export plot for paper
gse_list_df <- lapply(gse_list, as.data.frame)

# Add column with dataframe names
df_list_with_names <- Map(function(df, name) {
  df$name <- name
  return(df)
}, gse_list_df, names(gse_list_df))

# concatenate datasets
concatenated_df <- do.call(rbind, df_list_with_names)

# order MPs
order_MP <- c("Cycling", "MES-Hypoxia", "Embryonic-like", "Neuroepithelial-like", "Radial-glia-like", 
              "Embryonic-neuronal-like", "Neuronal-like"  ,"Ependymal-like")


concatenated_df$name <- factor(concatenated_df$name, levels = order_MP)

concatenated_df_subset <- concatenated_df %>% 
  group_by(name) %>% 
  arrange(p.adjust, .by_group = T) %>%
  slice(1:5)

concatenated_df_subset$name <- factor(concatenated_df_subset$name, levels = order_MP)


# Dot plot ---------------------------------------------------
concatenated_df_subset$Description <- factor(concatenated_df_subset$Description, levels = unique(concatenated_df_subset$Description))


ggplot(concatenated_df_subset, aes(x = Description, y = name, size = Count, color = -log(p.adjust))) +
  geom_point() +
  labs(x = "GO terms", y = "Metaprogram", color = "p.adjust") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1), 
        axis.text.y = element_text(color = "black"),  # Keep axis text for both axes
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=12)) +
  scale_color_paletteer_c("ggthemes::Blue-Teal") +
  geom_hline(yintercept = seq_along(unique(concatenated_df_subset$name)) - 0.5, 
             color = "grey50", linetype = "dashed") 
ggsave(file.path(plot_dir, '3_enrich_GO_top30_NMF_markers_paper.pdf'), width = 8, height = 4.5)




