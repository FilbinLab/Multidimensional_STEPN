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


# Organize environment  -----------------------------------

#base_dir <- '/Volumes/Sara_PhD/scRNAseq_data'
base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/11 - Plots scRNAseq"

genelist_dir <- file.path(base_dir, "analysis/plots/STEPN/malignant")

plot_dir <- file.path(base_dir, 'analysis/plots/GO')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

resource_dir <- file.path(base_dir, 'scripts/resources')

source(file.path(resource_dir,'Plotting_helper_functions.R'))



# load Seurat object -----------------------------------
seurat_obj <- qread(file = file.path(base_dir, "data/patient/seurat_obj_ST_combined_malig_annotated.qs"))
Idents(seurat_obj) <- 'Metaprogram'

# Run GSEA and GO -----------------------------------
plotList <- list()
gse_table <- list()
tumor_programs <- unique(seurat_obj$Metaprogram)
programs <- unique(seurat_obj$Metaprogram)

for (program in programs) {
  tryCatch({
  markers <- FindMarkers(seurat_obj, 
                         ident.1 = program, 
                         logfc.threshold = 0.5, 
                         min.pct = 0.15,
                         only.pos = F) %>% filter(p_val_adj < 0.05)
  
  markers <- markers %>% arrange(-avg_log2FC) %>% rownames_to_column('gene')
  markers <- markers %>% pull(avg_log2FC, gene)
  
  gse <- gseGO(geneList = markers, 
               ont = "ALL", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "BH")
  gse@result <- gse@result %>% filter(NES > 0)
  plotList[[program]] <- dotplot(gse, x ='NES', showCategory = 20) + ggtitle(program) + theme(plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1))
  gse_table[[program]] <- gse@result
  }, error = function(e) {
    cat("Error occurred for program:", program, "\n")
    cat("Error message:", conditionMessage(e), "\n")
    # You can choose to skip this program or handle the error differently
  })
}

pt <- wrap_plots(plotList[c(1:length(plotList))])
ggsave(plot = pt, file.path(plot_dir, 'enrich_GSEA.pdf'), width = 25, height = 20)
write_xlsx(gse_table, file.path(plot_dir, '1_GSEA.xlsx'))


plotList <- list()
go_table <- list()
gse_list <- list()
programs <- unique(seurat_obj$Metaprogram)
for (program in programs) {
  tryCatch({
  markers <- FindMarkers(seurat_obj, 
                         ident.1 = program, 
                         logfc.threshold = 0.5, 
                         min.pct = 0.15,
                         only.pos = F) %>% filter(p_val_adj < 0.05)
  markers <- markers %>% arrange(-avg_log2FC) %>% rownames_to_column('gene')
  markers <- markers %>% top_n(500, wt = avg_log2FC) %>% pull(gene)
  
  gse_list[[program]] <- enrichGO(gene = markers, 
                  ont = "ALL", 
                  keyType = "SYMBOL", 
                  minGSSize = 3, 
                  maxGSSize = 800,
                  pvalueCutoff = 0.05, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")
  gse@result <- gse@result %>% arrange(p.adjust)
  plotList[[program]] <- dotplot(gse, showCategory = 20) + ggtitle(program) + theme(plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1))
  go_table[[program]] <- gse@result
  }, error = function(e) {
    cat("Error occurred for program:", program, "\n")
    cat("Error message:", conditionMessage(e), "\n")
    # You can choose to skip this program or handle the error differently
  })
}

pt <- wrap_plots(plotList[c(1:length(plotList))])
ggsave(plot = pt, file.path(plot_dir, 'enrich_GO.pdf'), width = 30, height = 25)
write_xlsx(go_table, file.path(plot_dir, '2_GO.xlsx'))


# Plot GO terms for paper -------------------------------

gse_list_df <- lapply(gse_list, as.data.frame)

# Add column with dataframe names
df_list_with_names <- Map(function(df, name) {
  df$name <- name
  return(df)
}, gse_list_df, programs)

# concatenate datasets
concatenated_df <- do.call(rbind, df_list_with_names)


concatenated_df_subset <- concatenated_df %>% 
  group_by(name) %>% 
  arrange(p.adjust, .by_group = T) %>%
  slice(1:5)

# Dot plot
concatenated_df_subset$Description <- factor(concatenated_df_subset$Description, levels = unique(concatenated_df_subset$Description))

ggplot(concatenated_df_subset, aes(x = name, y = Description, size = Count, color = -log(p.adjust))) +
  geom_point() +
  #scale_size_continuous(range = c(1, 10)) +  # Adjust the range of dot sizes as needed
 # facet_grid(. ~ name) +
  labs(x = "Description", y = "Count", color = "p.adjust") +
  theme_minimal()  + # Customize the theme if needed+  # Customize the theme if needed
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1), 
        axis.text.y = element_text(color = "black"), # Rotate x-axis labels by 45 degrees
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.background = element_blank())  + 
 scale_color_paletteer_c("ggthemes::Blue-Teal", direction = -1)
ggsave(file.path(plot_dir, 'GO_top5_paper.pdf'), width = 6, height = 8)

