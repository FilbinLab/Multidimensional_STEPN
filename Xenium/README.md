## Xenium 

- Daeun's metadata path with programs annotation: `/labs/mfilbin/homes/djeong/dj83/analysis/Xenium/samples_metadata`
- Parameters for Jupyter session in O2:
  * *Modules to be preloaded*: `gcc/9.2.0 python/3.9.14 cmake/3.14.1 R/4.3.1`
  * *Jupyter Environment*: `source /home/cao385/envs/jupytervenv/bin/activate`
  * *Total Job Memory in GB for the job*: Usually `48` is enough
- To run niches analysis using Seurat package (click [here](https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#mouse-brain-10x-genomics-xenium-in-situ:) for more informations):
```
data <- qread('/n/scratch3/users/c/cao385/Xenium/results/20231115__203753__BT1873_BT1733/0010697-Region_1/0010697-Region_1.qs')
Idents(data) <- data$SCT_snn_res.0.1
data <- RenameIdents(data, '0' = 'AC-like', '1' = 'Ependymal-like-1', '2' = 'Glio_angiogenesis', '3' = 'Tcells', '4' = 'Neurons', '5' = 'Normal')
data$Annotation <- ifelse(Idents(data) %in% c('Tcells', 'Normal'), 'Immune', 'Malignant')

data <- BuildNicheAssay(object = data, fov = "fov", group.by = "Annotation", niches.k = 5, neighbors.k = 30)
celltype.plot <- ImageDimPlot(data, group.by = "Annotation", size = 0.8, cols = colors_to_use[1:2], dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(data, group.by = "niches", size = 0.8, dark.background = F) + ggtitle("Niches") + scale_fill_manual(values = colors_to_use)
celltype.plot | niche.plot

t(table(data$Annotation, data$niches))
```
- Neighbor analysis in `Neighbors-Daeun.nbconvert.ipynb` script.
