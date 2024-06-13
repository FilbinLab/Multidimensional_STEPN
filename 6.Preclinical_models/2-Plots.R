library(ggplot2)
library(cowplot)
library(tidyverse)

colors_to_use <- clusterExperiment::bigPalette

base_dir <- "/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/6 - scRNAseq models"

metadata <- read_excel(file.path(base_dir, 'metadata.xlsx')) 

FileName <- unique(metadata$FileName)
Project <- unique(metadata$Project)
Type <- unique(metadata$Type)


files <- sprintf('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/6 - scRNAseq models/analysis/qc/data/individual/%s-qc_raw_data.rds', unique(metadata$FileName))
fnames <- gsub('-qc_raw_data.rds', '', basename(files))

qc_files <- lapply(files, readRDS) %>% 
  setNames(fnames) %>% 
  bind_rows(.id = "tumor")

qc <- qc_files %>% 
  group_by(sample, pass_qc) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  left_join(qc_files[,c('tumor', 'sample'),], by = 'sample', multiple = "first")


qc$sample <- factor(qc$sample, levels = qc %>% filter(pass_qc == 'TRUE') %>% arrange(-freq) %>% pull(sample) %>% as.character())


h1 <- ggplot(qc) +
  geom_bar(mapping = aes(x = sample, y = freq, fill = pass_qc), stat = "identity", width = 1) +
  guides(fill = guide_legend(ncol = 2, title = 'Good quality')) +
  labs(x = 'Sample', y = 'Percentage') +
  theme(axis.text.x.bottom = element_blank(), 
        axis.ticks = element_blank(),
        panel.spacing.x = unit(0, "mm"),
        axis.title.x = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white", linewidth = 0.5, linetype = "solid"),
        plot.margin = unit(c(0.1, 0, 0, 0), "cm"), 
        legend.key.size = unit(0.5,"line")) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = c('#0085A6', '#F9C960')) +
  facet_grid(.~tumor, scales = "free", space = "free")
h2 <- ggplot(qc) +
  geom_bar(mapping = aes(x = sample, y = 1, fill = tumor), stat = "identity", width = 1) +
  scale_fill_manual(values = colors_to_use[-2]) +
  theme_void() +
  theme(panel.spacing.x = unit(0, "mm"), 
        strip.background = element_blank(),
        strip.text.x = element_blank(), 
        legend.key.size = unit(0.5,"line")) +
  facet_grid(.~tumor, scales = "free", space = "free")
legend <- plot_grid(get_legend(h2), get_legend(h1), ncol = 1)
h1 <- h1 + theme(legend.position = "none")
h2 <- h2 + theme(legend.position = "none")
plot <- plot_grid(h2, h1, align = "v", ncol = 1, axis = "tb", rel_heights = c(0.5, 15))
plot_grid(plot, legend, nrow = 1, rel_widths = c(10, 1.5))

ggsave('qc.pdf', path = '/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/6 - scRNAseq models/analysis/qc', width = 8, height = 4)
ggsave('qc_legend.pdf', path = '/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Project/Ependymoma/6 - scRNAseq models/analysis/qc', width = 8, height = 8)
