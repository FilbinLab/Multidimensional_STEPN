################################################################################
## Projection for hCortex
################################################################################



################################################################################
## Projection for earlyDev
################################################################################



################################################################################
## DotPlot for hCortex + earlyDev references
################################################################################






















################################################################################
############################## earlyDev by AgePCW ##############################
################################################################################
ref <- readRDS('/n/scratch/users/c/cao385/Ependymoma/dj83/data/Normal_ref_data/Nowakowski_earlydev_natureneuro2021/seurat objects/seurat_obj_Nowakowski_earlyDev.rds')
table(test$Cell.Type)

ref$Age_PCW <- ref$Age
ref$Age_PCW <- gsub(12, 4, ref$Age_PCW)
ref$Age_PCW <- gsub(13, 4, ref$Age_PCW)
ref$Age_PCW <- gsub(14, 5, ref$Age_PCW)
ref$Age_PCW <- gsub(15, 6, ref$Age_PCW)
ref$Age_PCW <- gsub(19, 7, ref$Age_PCW)
ref$Age_PCW <- gsub(20, 8, ref$Age_PCW)
ref$Age_PCW <- gsub(22, 8, ref$Age_PCW)



## Reference
ref_degs <- readRDS("/n/scratch/users/c/cao385/Ependymoma/dj83/data/Normal_ref_data/Nowakowski_earlydev_natureneuro2021/Projection RDS/ref_degs/ref_DE_genes_by_PCWAge_top50.rds")
ref_hvgs <- readRDS("/n/scratch/users/c/cao385/Ependymoma/dj83/data/Normal_ref_data/Nowakowski_earlydev_natureneuro2021/Projection RDS/hvgs/hvg.rds")
ref_cm <- readRDS("/n/scratch/users/c/cao385/Ependymoma/dj83/data/Normal_ref_data/Nowakowski_earlydev_natureneuro2021/Projection RDS/pseudobulk/pseudobulk_by_PCWAge.rds")
ref_degs <- ref_degs[order(match(names(ref_degs),colnames(ref_cm)))]


## Tumor
pb_degs <- qread(file.path(analysis_dir, "DEGs/ST_samples_DE_list.qs"))
pb_hvgs <- qread(file.path(analysis_dir, "Projection RDS/ST_hvgs.qs"))
pb_cm <- qread(file.path(analysis_dir, "Projection RDS/ST_samples_pseudobulk.qs"))


## Signature genes
input_genes = union(ref_hvgs, pb_hvgs)
input_genes = intersect(input_genes, intersect(rownames(ref_cm), rownames(pb_cm)))



# Normalize data
## Ref
ref_cm_list <- NormCenter(ref_cm)
ref_cm_norm <- ref_cm_list$norm_data
ref_cm_center <- ref_cm_list$center_data
ref_cm_mean = log2(rowMeans(ref_cm)+1)

## Tumor
pb_cm_norm = log2(pb_cm+1)
pb_cm_mean = log2(rowMeans(pb_cm)+1)
pb_cm_center = t(scale(t(pb_cm_norm)))
pb_cm_norm <- as.matrix(pb_cm_norm)


# Score metagene programs in ref scRNA-seq data 
normal_score = scoreNmfGenes(ref_cm_center, ref_cm_mean, pb_degs)
normal_score = t(normal_score)

tumor_score = scoreNmfGenes(pb_cm_center, pb_cm_mean, ref_degs)
colnames(tumor_score) = gsub("-","_",colnames(tumor_score))
tumor_score = tumor_score[rownames(normal_score), colnames(normal_score)]


# Pairwise correlation between ref and epn aggregated cm 
pairwise_cor = calculatePairwiseCor(ref_cm_norm, pb_cm_norm, input_genes)
colnames(pairwise_cor) = gsub("-","_",colnames(pairwise_cor))
pairwise_cor = pairwise_cor[rownames(normal_score), colnames(normal_score)]
##saveRDS(pairwise_cor, file=paste0(analysis_dir, ref_df, "_pairwise_cor.rds"))


# Dotplot

## Prepare data
## max-min normalization and melt scoring normal cells by malignant metaprogram
tmp = apply(normal_score, 2, function(x) (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)))
normal_score_long = melt(tmp)
colnames(normal_score_long) = c("Cell_type", "Subtype", "score")

## max-min normalization and melt scoring tumor cells by normal DEGs
tmp = apply(tumor_score, 2, function(x) (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)))
tumor_score_long = melt(tmp)
colnames(tumor_score_long) = c("Cell_type", "Subtype", "score")

## Melt pairwise correlation 
pairwise_cor_long = melt(pairwise_cor)
colnames(pairwise_cor_long) = c("Cell_type", "Subtype", "r")

## Add score to pairwise-cor
pairwise_cor_long$normal_score = normal_score_long$score
pairwise_cor_long$tumor_score = tumor_score_long$score

## Replace "_" with "-"
#pairwise_cor_long$Cell_type = gsub("_", "-", pairwise_cor_long$Cell_type)
#pairwise_cor_long$Metaprogram = gsub("_", "-", pairwise_cor_long$Metaprogram)

pairwise_cor_long[is.na(pairwise_cor_long)] = 0




## Plot
pairwise_cor_long2 <- pairwise_cor_long

ref_order_trimmed = c("PCW4", "PCW5", "PCW6", "PCW7", "PCW8")


sample_order <- c("3EP11","4EP49","3EP8","BT1743","4EP51","7EP35","7EP41","11EP22","3EP54","3EP67","7EP1","11EP8","7EP9","16EP8","4EP53","MUV43R1","MUV43R2","MUV43R4","WEPN9","4EP44","9EP47","9EP35","BT268","I128034","I128034R1","I128034R2","1230717","1599417","2192118","RD050618","RD141019","RD251119","4717EP17","FR0671","I050420","4EP46","9EP45","11EP21","3EP29")


p = ggplot(pairwise_cor_long2, aes(x=factor(Subtype,sample_order), y=factor(Cell_type,ref_order_trimmed))) +
  geom_point(aes(color=tumor_score,size=normal_score,alpha=r)) + scale_size_area(max_size=12) +
  scale_color_gradientn(colors=brewer.pal(n=9, name="Reds")) +
  labs(title = paste0("ST-EPNs"," "),
       x = "",
       y = "",
       color = "Expression score\n(tumor cells)", 
       size = "Expression score\n(normal cells)",
       alpha="Correlation Coefficient") +
  theme_classic() + 
  theme(plot.title = element_text(size=24, hjust = 0.5),
        axis.title = element_text(size=24),
        axis.text.x=element_text(size=20, hjust=1,angle=45),
        axis.text.y=element_text(size=20),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'dashed', colour = "grey"))+coord_flip()
p + theme(plot.margin = margin(10,10,10,40))

ggsave("/n/scratch/users/c/cao385/ST_samples_earlyDev_projection_AgePCW.png", width=9.5, height=18, dpi = 1200)






View(as.data.frame(ref_degs$ref_degs))

View(as.data.frame(ref_degs))



data <- qread('/n/scratch/users/c/cao385/Linnarsson_2022/seurat.qs')
DimPlot(data, reduction = "umap.full", group.by = 'CellClass', cols = clusterExperiment::bigPalette, alpha = 0.1)

## DEGs
Idents(data) <- data$CellClass
markers <- RunPrestoAll(data, only.pos = T, logfc.threshold = 0.5) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T)
qsave(markers, "/n/scratch/users/c/cao385/Linnarsson_2022/markers_CellClass.qs")


####################### Linnarsson (aiming neural-crest) #######################
## Reference
ref_degs <- qread("/n/scratch/users/c/cao385/Linnarsson_2022/markers_CellClass.qs")
ref_degs <- ref_degs %>% group_by(cluster) %>% top_n(50, wt = avg_log2FC)
ref_degs <- lapply(split(ref_degs, f = ref_degs$cluster), function(x) x$gene)
ref_hvgs <- qread("/n/scratch/users/c/cao385/Linnarsson_2022/hvgs.qs")
ref_cm <- qread("/n/scratch/users/c/cao385/Linnarsson_2022/agg_cm_mean_CellClass.qs")
ref_degs <- ref_degs[order(match(names(ref_degs),colnames(ref_cm)))]


## Tumor
pb_degs <- qread(file.path(analysis_dir, "DEGs/ST_samples_DE_list.qs"))
pb_hvgs <- qread(file.path(analysis_dir, "Projection RDS/ST_hvgs.qs"))
pb_cm <- qread(file.path(analysis_dir, "Projection RDS/ST_samples_pseudobulk.qs"))


## Signature genes
input_genes = union(ref_hvgs, pb_hvgs)
input_genes = intersect(input_genes, intersect(rownames(ref_cm), rownames(pb_cm)))



# Normalize data
## Ref
ref_cm_list <- NormCenter(ref_cm)
ref_cm_norm <- ref_cm_list$norm_data
ref_cm_center <- ref_cm_list$center_data
ref_cm_mean = log2(rowMeans(ref_cm)+1)

## Tumor
pb_cm_norm = log2(pb_cm+1)
pb_cm_mean = log2(rowMeans(pb_cm)+1)
pb_cm_center = t(scale(t(pb_cm_norm)))
pb_cm_norm <- as.matrix(pb_cm_norm)


# Score metagene programs in ref scRNA-seq data 
normal_score = scoreNmfGenes(ref_cm_center, ref_cm_mean, pb_degs)
normal_score = t(normal_score)

tumor_score = scoreNmfGenes(pb_cm_center, pb_cm_mean, ref_degs)
colnames(tumor_score) = gsub("-","_",colnames(tumor_score))
tumor_score = tumor_score[rownames(normal_score), colnames(normal_score)]


# Pairwise correlation between ref and epn aggregated cm 
pairwise_cor = calculatePairwiseCor(ref_cm_norm, pb_cm_norm, input_genes)
colnames(pairwise_cor) = gsub("-","_",colnames(pairwise_cor))
pairwise_cor = pairwise_cor[rownames(normal_score), colnames(normal_score)]
##saveRDS(pairwise_cor, file=paste0(analysis_dir, ref_df, "_pairwise_cor.rds"))


# Dotplot

## Prepare data
## max-min normalization and melt scoring normal cells by malignant metaprogram
tmp = apply(normal_score, 2, function(x) (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)))
normal_score_long = melt(tmp)
colnames(normal_score_long) = c("Cell_type", "Subtype", "score")

## max-min normalization and melt scoring tumor cells by normal DEGs
tmp = apply(tumor_score, 2, function(x) (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)))
tumor_score_long = melt(tmp)
colnames(tumor_score_long) = c("Cell_type", "Subtype", "score")

## Melt pairwise correlation 
pairwise_cor_long = melt(pairwise_cor)
colnames(pairwise_cor_long) = c("Cell_type", "Subtype", "r")

## Add score to pairwise-cor
pairwise_cor_long$normal_score = normal_score_long$score
pairwise_cor_long$tumor_score = tumor_score_long$score

## Replace "_" with "-"
#pairwise_cor_long$Cell_type = gsub("_", "-", pairwise_cor_long$Cell_type)
#pairwise_cor_long$Metaprogram = gsub("_", "-", pairwise_cor_long$Metaprogram)

pairwise_cor_long[is.na(pairwise_cor_long)] = 0




## Plot
pairwise_cor_long2 <- pairwise_cor_long

ref_order_trimmed <- unique(pairwise_cor_long2$Cell_type)

sample_order <- rev(c("3EP11","4EP49","3EP8","BT1743","4EP51","7EP35","7EP41","11EP22","3EP54","3EP67","7EP1","11EP8","7EP9","16EP8","4EP53","MUV43R1","MUV43R2","MUV43R4","WEPN9","4EP44","9EP47","9EP35","BT268","I128034","I128034R1","I128034R2","1230717","1599417","2192118","RD050618","RD141019","RD251119","4717EP17","FR0671","I050420","4EP46","9EP45","11EP21","3EP29"))


pt2 = ggplot(pairwise_cor_long2, aes(x=factor(Subtype,sample_order), y=factor(Cell_type,ref_order_trimmed))) +
  geom_point(aes(color=tumor_score,size=normal_score,alpha=r)) + scale_size_area(max_size=12) +
  scale_color_gradientn(colors=brewer.pal(n=9, name="Reds")) +
  labs(title = paste0("ST-EPNs"," "),
       x = "",
       y = "",
       color = "Expression score\n(tumor cells)", 
       size = "Expression score\n(normal cells)",
       alpha="Correlation Coefficient") +
  theme_classic() + 
  theme(plot.title = element_text(size=24, hjust = 0.5),
        axis.title = element_text(size=24),
        axis.text.x=element_text(size=20, hjust=1,angle=45),
        axis.text.y=element_text(size=20),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'dashed', colour = "grey"))+coord_flip()
pt2 <- pt2 + theme(plot.margin = margin(10,10,10,40))
ggsave(plot = pt2, "/n/scratch/users/c/cao385/Ependymoma/dj83/New/ST_samples_Linnarsson_Crest.png", width = 9.5, height = 18, dpi = 1200)




test <- qread('/n/scratch/users/c/cao385/seurat.qs')
DimPlot(test, group.by = 'CellType', cols = clusterExperiment::bigPalette, reduction = 'umap.harmony', split.by = 'CellType')


test <- qread("/n/scratch/users/c/cao385/Ependymoma/dj83/data/Normal_ref_data/Nowakowski_earlydev_natureneuro2021/seurat objects/seurat_obj_Nowakowski_earlyDev_updated.qs")
DimPlot(test, group.by = 'Cell.Type', cols = clusterExperiment::bigPalette, reduction = 'umap', split.by = 'Cell.Type')