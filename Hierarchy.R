library(qs)
library(Seurat)
library(tidyverse)
library(stringr)
library(patchwork)
library(glue)

base_dir <- "/n/scratch/users/c/cao385/Ependymoma/dj83"
script_dir <- "helper_scripts"
source(file.path(base_dir,script_dir, "single_cell_preprocessing_helper_functions.R"))
source(file.path(base_dir, script_dir, "NMF_helper_function.R"))

seurat_dir <- file.path(base_dir, "data/seurat_objects")
ZFTA_dir <- file.path(base_dir, "data/ZFTA_frozen/")
fig_dir <- file.path(base_dir, "analysis/ST_all/Figures/Hierarchy Plots/")
analysis_dir <- file.path(base_dir, "analysis/ST_all/")
ref_dir <- file.path(base_dir, "data/Normal_ref_data/Kriegstein_Human_Cortex/")
ref <- qread(file.path(ref_dir,"seurat objects/seurat_obj_hCortex_Kriegstein.qs"))

#ref2.genes <- readRDS('/n/scratch/users/c/cao385/Ependymoma/dj83/data/Normal_ref_data/Nowakowski_earlydev_natureneuro2021/Projection RDS/ref_degs/ref_DE_genes_by_CellTypes_allMarkers.rds')
#ref2.genes <- ref2.genes[!str_detect(ref2.genes$gene, "RP11"), ]
#top.markers2 <- ref2.genes %>%
#  group_by(cluster) %>%
#  top_n(95, avg_log2FC) %>%
#  arrange(avg_log2FC, .by_group = TRUE)
#top.markers2 <- lapply(split(top.markers2, f = top.markers2$cluster), function(x) x$gene)
top.markers2 <- qread(glue('{base_dir}/data/Normal_ref_data/Nowakowski_earlydev_natureneuro2021/EarlyDev_markers.qs'))

ref.genes <- read.csv(file.path(ref_dir, "Nowakowski2017_markers.csv"))
ref.genes <- ref.genes[ref.genes$avg_diff > 0,]
ref.genes <- ref.genes[!str_detect(ref.genes$gene, "RP11"), ]

top.markers <- ref.genes %>%
  group_by(cluster) %>%
  top_n(95, avg_diff) %>%
  arrange(avg_diff, .by_group = TRUE)

sample.clusters <- unique(top.markers$cluster)
sample.gene.list <- list()

for (i in 1:length(sample.clusters)) {
  tmp <- top.markers[top.markers$cluster==sample.clusters[i],] %>% 
    arrange(-avg_diff)
  tmp_genes <- c(tmp$gene)
  sample.gene.list <- append(sample.gene.list,list(tmp_genes))
}

names(sample.gene.list) <- sample.clusters

sample.gene.list$Neurons <- unique(c(top.markers2$Neuronal[1:10], sample.gene.list$`EN-PFC1`[1:10], sample.gene.list$`EN-PFC2`[1:10], sample.gene.list$`EN-PFC3`[1:10], sample.gene.list$`EN-V1-1`[1:10], sample.gene.list$`EN-V1-2`[1:10], sample.gene.list$`EN-V1-3`[1:10], sample.gene.list$`IN-CTX-CGE1`[1:10], sample.gene.list$`IN-CTX-CGE2`[1:10], sample.gene.list$`IN-CTX-MGE1`[1:10], sample.gene.list$`IN-CTX-MGE2`[1:10], sample.gene.list$`nEN-early1`[1:10], sample.gene.list$`nEN-early2`[1:10], sample.gene.list$`nEN-late`[1:10], sample.gene.list$nIN1[1:10], sample.gene.list$nIN2[1:10], sample.gene.list$nIN3[1:10], sample.gene.list$nIN4[1:10], sample.gene.list$nIN5[1:10]))


x <- sample.gene.list$Neurons

#new_list <- list(sample.gene.list$`RG-early`,sample.gene.list$Choroid, x)
#names(new_list) <- c("RG-early","Ependymal","Neurons")

#sample.gene.list$`RG-early`, top.markers2$Neuroepithelial

new_list <- list(unique(c(top.markers2$Neuroepithelial[1:50], sample.gene.list$`RG-early`[1:50])), sample.gene.list$Choroid, x)
names(new_list) <- c("Neuroepithelial","Ependymal","Neurons")

#qsave(new_list,paste0(ref_dir,"Nowakowski2017_markers_list.qs"))
ST <- qread(file.path(seurat_dir,"seurat_obj_ST_all_malig_filtered.qs"))

ST@meta.data <- ST@meta.data[,1:16]
ST <- AddModuleScore(
  object = ST,
  features = new_list,
  ctrl = 5,
  name = 'program_features'
)
colnames(ST@meta.data)[17:19] <- names(new_list)

#cm_norm <- as.matrix(log2(LayerData(ST, 'counts')/10+1))
#cm_mean <- log2(Matrix::rowMeans(LayerData(ST, 'counts'))+1)
#cm_center <- cm_norm - rowMeans(cm_norm)
#nmf_score <- t(scoreNmfGenes(cm_center, cm_mean, new_list))
#ST <- AddMetaData(ST, nmf_score)



x1 <- "Ependymal"
x2 <- "Neurons"
y1 <- "Neuroepithelial"
#y1 = "RG-early"
group.by <- "subtype"
sample <- ST


variables_to_retrieve <- c(x1, x2, y1, group.by)

# And store them as a tibble.
scores <- sample@meta.data[, variables_to_retrieve]
scores[["cell"]] <- rownames(scores)
# Shuffle the cells so that we accomplish a random plotting, not sample by sample.
scores <- tidyr::tibble(scores)

# Compute the scores for the X axis.
x <- pbapply::pbvapply(seq_len(nrow(scores)), function(x) {
  score_1 <- scores[x, x1] + stats::runif(1, min=0, max=0.15)
  score_2 <- scores[x, x2] + stats::runif(1, min=0, max=0.15)
  d <- max(score_1, score_2, na.rm = TRUE)
  output <- ifelse(score_1 > score_2, d, -d)
  return(output)
}, FUN.VALUE = numeric(1))


# Compute the scores for the Y axis.
y <- pbapply::pbvapply(seq_len(nrow(scores)), function(x) {
  score_1 <- scores[x, x1] + stats::runif(1, min=0, max=0.15)
  score_2 <- scores[x, x2] + stats::runif(1, min=0, max=0.15)
  d <- max(score_1, score_2, na.rm = TRUE)
  output <- as.data.frame(scores)[x, y1] - d
  return(output)
}, FUN.VALUE = numeric(1))

names(x) <- scores[["cell"]]
names(y) <- scores[["cell"]]

# Define titles.
lab_x <- paste0(x2, "  <---->  ", x1)
lab_y <- 'Neuroepithelial/RG-early'

df <- data.frame("set_x" = x, "set_y" = y, "group.by" = scores[[group.by]])

samples = unique(df$group.by)

for (i in 1:length(samples)) {
  tmp <- df$group.by==samples[i]
  tmp <- gsub("TRUE",samples[i],tmp)
  tmp <- gsub("FALSE","Other",tmp)
  df <- cbind(df,tmp)
}
subtypes <- c("ZFTA_RELA","ZFTA_Cluster1","ZFTA_Cluster2","ZFTA_Cluster3","ZFTA_Cluster4","ST_YAP1")
colnames(df)[4:9] <- subtypes

# Define limits of plots
lim <- max(abs(x))
lim_x <- c(-lim, lim)
lim <- max(abs(y))
lim_y <- c(-lim, lim)

#Define colors 
subtypes <- c("ZFTA_RELA","ZFTA_Cluster1","ZFTA_Cluster2","ZFTA_Cluster3","ZFTA_Cluster4","ST_YAP1")
colors.use=c("#B44672","#B47846","#46B478","#46B4AF","#4682B4","#B4AF46")
df$group.by <- factor(df$group.by, levels = unique(df$group.by))


ggplot(df,aes(x = set_x, y = set_y, col = group.by))+
  geom_point(size = 2)+
  scale_color_manual(values=colors.use)+
  xlim(-1,1)+
  ylim(lim_y)+
  theme_classic()+
  labs(x=lab_x, y=lab_y)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, family="Helvetica"),
        axis.title = element_text(size = 15, family="Helvetica",face="bold"), axis.text.x = element_text(size=10,hjust = 0.5, family="Helvetica",face="bold"), axis.text.y=element_text(size=10,family="Helvetica",face="bold"),legend.text = element_text(family="Helvetica", size=12,face="bold")) +
  theme(legend.text=element_text(family="Helvetica"))+
  theme(legend.title=element_blank())+
  theme(legend.direction="horizontal",legend.position="bottom")
ggsave("/n/scratch/users/c/cao385/Ependymoma/dj83/New/ST_hierarchy_subtype.pdf", width = 10.5, height = 10, dpi = 500)




colors.use=c("#B44672","#B47846","#46B478","#46B4AF","#4682B4","#B4AF46")

p1 <- ggplot(df %>% arrange(ZFTA_RELA))+
  geom_point(aes(x=set_x, y=set_y,color=ZFTA_RELA),size=1)+
  scale_colour_manual(values=c("grey90","#B44672"), name="", labels=c("Other","ZFTA-RELA"))+
  xlim(-1,1)+
  ylim(lim_y)+
  theme_classic()+
  labs(x=lab_x, y=lab_y)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, family="Helvetica"),
        axis.title = element_text(size = 15, family="Helvetica",face="bold"), axis.text.x = element_text(size=10,hjust = 0.5, family="Helvetica",face="bold"), axis.text.y=element_text(size=10,family="Helvetica",face="bold"),legend.text = element_text(family="Helvetica", size=12,face="bold")) +
  theme(legend.text=element_text(family="Helvetica"))+
  theme(legend.title=element_blank())+
  theme(legend.direction="horizontal",legend.position="bottom")
ggsave(plot = p1, "/n/scratch/users/c/cao385/Ependymoma/dj83/New/ZFTA-RELA_highlight_hierarchy.pdf",width=10.5,height=10,dpi=500)

p2 <- ggplot(df %>% arrange(ZFTA_Cluster1))+
  geom_point(aes(x=set_x, y=set_y,color=ZFTA_Cluster1),size=1)+
  scale_colour_manual(values=c("grey90","#B47846"), name="", labels=c("Other","ZFTA-Cluster 1"))+
  xlim(-1,1)+
  ylim(lim_y)+
  theme_classic()+
  labs(x=lab_x, y=lab_y)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, family="Helvetica"),
        axis.title = element_text(size = 15, family="Helvetica",face="bold"), axis.text.x = element_text(size=10,hjust = 0.5, family="Helvetica",face="bold"), axis.text.y=element_text(size=10,family="Helvetica",face="bold"),legend.text = element_text(family="Helvetica", size=12,face="bold")) +
  theme(legend.text=element_text(family="Helvetica"))+
  theme(legend.title=element_blank())+
  theme(legend.direction="horizontal",legend.position="bottom")
ggsave(plot = p2, "/n/scratch/users/c/cao385/Ependymoma/dj83/New/ZFTA-Cluster 1_highlight_hierarchy.pdf",width=10.5,height=10,dpi=500)

p3 <- ggplot(df %>% arrange(ZFTA_Cluster2))+
  geom_point(aes(x=set_x, y=set_y,color=ZFTA_Cluster2),size=1)+
  scale_colour_manual(values=c("grey90","#46B478"), name="", labels=c("Other","ZFTA-Cluster 2"))+
  xlim(-1,1.0)+
  ylim(lim_y)+
  theme_classic()+
  labs(x=lab_x, y=lab_y)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, family="Helvetica"),
        axis.title = element_text(size = 15, family="Helvetica",face="bold"), axis.text.x = element_text(size=10,hjust = 0.5, family="Helvetica",face="bold"), axis.text.y=element_text(size=10,family="Helvetica",face="bold"),legend.text = element_text(family="Helvetica", size=12,face="bold")) +
  theme(legend.text=element_text(family="Helvetica"))+
  theme(legend.title=element_blank())+
  theme(legend.direction="horizontal",legend.position="bottom")
ggsave(plot = p3, "/n/scratch/users/c/cao385/Ependymoma/dj83/New/ZFTA-Cluster 2_highlight_hierarchy.pdf",width=10.5,height=10,dpi=500)

p4 <- ggplot(df %>% arrange(ZFTA_Cluster3))+
  geom_point(aes(x=set_x, y=set_y,color=ZFTA_Cluster3),size=1)+
  scale_colour_manual(values=c("grey90","#46B4AF"), name="", labels=c("Other","ZFTA-Cluster 3"))+
  xlim(-1,1)+
  ylim(lim_y)+
  theme_classic()+
  labs(x=lab_x, y=lab_y)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, family="Helvetica"),
        axis.title = element_text(size = 15, family="Helvetica",face="bold"), axis.text.x = element_text(size=10,hjust = 0.5, family="Helvetica",face="bold"), axis.text.y=element_text(size=10,family="Helvetica",face="bold"),legend.text = element_text(family="Helvetica", size=12,face="bold")) +
  theme(legend.text=element_text(family="Helvetica"))+
  theme(legend.title=element_blank())+
  theme(legend.direction="horizontal",legend.position="bottom")
ggsave(plot = p4, "/n/scratch/users/c/cao385/Ependymoma/dj83/New/ZFTA-Cluster 3_highlight_hierarchy.pdf", width = 10.5, height = 10, dpi = 500)

p5 <- ggplot(df %>% arrange(ZFTA_Cluster4))+
  geom_point(aes(x=set_x, y=set_y,color=ZFTA_Cluster4),size=1)+
  scale_colour_manual(values=c("grey90","#4682B4"), name="", labels=c("Other","ZFTA-Cluster 4"))+
  xlim(-1,1)+
  ylim(lim_y)+
  theme_classic()+
  labs(x=lab_x, y=lab_y)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, family="Helvetica"),
        axis.title = element_text(size = 15, family="Helvetica",face="bold"), axis.text.x = element_text(size=10,hjust = 0.5, family="Helvetica",face="bold"), axis.text.y=element_text(size=10,family="Helvetica",face="bold"),legend.text = element_text(family="Helvetica", size=12,face="bold")) +
  theme(legend.text=element_text(family="Helvetica"))+
  theme(legend.title=element_blank())+
  theme(legend.direction="horizontal",legend.position="bottom")
ggsave(plot = p5, "/n/scratch/users/c/cao385/Ependymoma/dj83/New/ZFTA-Cluster 4_highlight_hierarchy.pdf", width = 10.5, height = 10, dpi = 500)

p6 <- ggplot(df %>% arrange(ST_YAP1))+
  geom_point(aes(x=set_x, y=set_y,color=ST_YAP1),size=1)+
  scale_colour_manual(values=c("grey90","#B4AF46"), name="", labels=c("Other","ST-YAP1"))+
  xlim(-1,1)+
  ylim(lim_y)+
  theme_classic()+
  labs(x=lab_x, y=lab_y)+
  theme(plot.title = element_text(hjust = 0.5, size = 20, family="Helvetica"),
        axis.title = element_text(size = 15, family="Helvetica",face="bold"), axis.text.x = element_text(size=10,hjust = 0.5, family="Helvetica",face="bold"), axis.text.y=element_text(size=10,family="Helvetica",face="bold"),legend.text = element_text(family="Helvetica", size=12,face="bold")) +
  theme(legend.text=element_text(family="Helvetica"))+
  theme(legend.title=element_blank())+
  theme(legend.direction="horizontal",legend.position="bottom")
ggsave(plot = p6, "/n/scratch/users/c/cao385/Ependymoma/dj83/New/ST-YAP1_highlight_hierarchy.pdf", width = 10.5, height = 10, dpi = 500)


wrap_plots(p1, p2, p3, p4, p5, p6)
ggsave("/n/scratch/users/c/cao385/Ependymoma/dj83/New/ST_hierarchy.pdf", width = 15, height = 10, dpi = 500)
