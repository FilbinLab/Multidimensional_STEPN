library(ggplot2)

## Make a bar graph to summarize proportion of clusters/programs/etc in each sample
## @para: x, x-axis data
## @para: y, y-axis data
## @para: x_order: order of x-axis variables
## @para: y_order: order of y-axis variables
## @para: col_names: column names for df for plotting
## @para: x_var: name of x-axis variable to plot
## @para: y_var: name of y-axis variable to plot
## @para: fill_var: name of variable for filling the bar graph
## @para: plot_title: title of the plot
## @para: x_title: title of x-axis
## @para: y_title: title of y-axis
## @para: bar_position: stack vs dodge
## The rest parameters are sizes of different labels
plotProportion <- function(x, y, x_order, y_order, col_names,
                           x_var, y_var, fill_var, colors,
                           plot_title="", x_title="", y_title="",
                           bar_position = "stack",
                           title_size=32, axis_title_size=28,
                           x_text_size=20, x_text_angle=45, x_text_hjust=1,
                           y_text_size=24, legend_title_size=28, legend_text_size=24){
    crosstab = table(x, y)
    crosstab = crosstab/rowSums(crosstab)
    crosstab = crosstab[x_order, y_order]
    crosstab = data.frame(crosstab)

    colnames(crosstab) = col_names
    ggplot(crosstab, aes_string(x=x_var, y=y_var, fill=fill_var)) +
        geom_bar(stat="identity", color="black", position = bar_position) +
        scale_fill_manual(values = colors) +
        ggtitle(plot_title) + xlab(x_title) + ylab(y_title) +
        theme_classic() + theme(plot.title = element_text(hjust=0.5, face="bold", size=title_size),
                                axis.title = element_text(face="bold", size=axis_title_size),
                                axis.text.x = element_text(angle=x_text_angle, hjust=x_text_hjust, size=x_text_size),
                                axis.text.y = element_text(size=y_text_size),
                                legend.title = element_text(face="bold", size=legend_title_size),
                                legend.text = element_text(size=legend_text_size))
}

## Make a grid of violin plots for gene expressions for each cell subpopulation
## @para: cm, count matrix (log transformed value)
## @para: annotation, cell annotation of each cell
## @para: genes, genes to plot
## @para: drug, name of the drug (only used for image name)
## @para: out, output directory
## The rest label sizes and image sizes
plotGridViolin <- function(cm, annotation, genes, drug, out = seurat_fig_folder,
                           y_axis_title = 24, x_axis_text = 24, strip_text = 12,
                           out_width = 16, out_height = 12){
  tmp = data.frame(cm[,genes])
  tmp$metaprogram = annotation
  tmp = melt(tmp, id.vars = "metaprogram")
  colnames(tmp) = c("metaprogram", "gene", "expression")

  ggplot(tmp, aes(x=metaprogram, y=expression, fill = metaprogram)) +
    geom_violin(scale = "width") +
    xlab("") + ylab("Log2 expression\n") +
    scale_fill_manual(values = color_scheme) +
    facet_grid(gene~.) +
    theme(axis.title.y = element_text(size = y_axis_title, face = "bold"),
          axis.text.x = element_text(size = x_axis_text, face = "bold", angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          strip.text = element_text(size=strip_text, face="bold"),
          axis.ticks.y = element_blank(),
          plot.margin = margin(10,10,10,70),
          legend.position = "none")
  ggsave(paste0("PF_", drug, "_targets_violin.png"),
         path = out, width = out_width, height = out_height)
}

## Plot dotplot with p value and jaccard index
plotProgramOverlap <- function(a, b, a_name, b_name, x_lab, y_lab,
                               color_lab = "-log10(p-value)", size_lab = "Jaccard Index",
                               title_size=32, axis_title_size=28,
                               x_text_size=20, x_text_angle=45, x_text_hjust=1,
                               y_text_size=24, legend_title_size=28, legend_text_size=24){
    ## Compute pairwise fisher exact pvalue and jaccard index
    bg_genes = unique(c(unlist(a), unlist(b)))
    pvalue_res = NULL
    jaccard_res = NULL

    for (prog1 in a){
        tmp = NULL
        tmp2 = NULL
        for (prog2 in b){
            tmp = c(tmp, fisher_test(prog1, prog2, bg_genes))
            tmp2 = c(tmp2, jaccard_index(prog1, prog2))
        }
        pvalue_res = rbind(pvalue_res, tmp)
        jaccard_res = rbind(jaccard_res, tmp2)
    }

    rownames(pvalue_res) = names(a); colnames(pvalue_res) = names(b)
    rownames(jaccard_res) = names(a); colnames(jaccard_res) = names(b)

    ## Combine p value and jaccard index in long format
    pvalue_res = melt(pvalue_res)
    jaccard_res = melt(jaccard_res)
    df = cbind.data.frame(pvalue_res$Var1, pvalue_res$Var2, pvalue_res$value, jaccard_res$value)
    colnames(df) = c(a_name, b_name, "P_value", "Jaccard_Index")
    df$P_value = -log10(df$P_value)

    ## Plot dotplot
    ggplot(df, aes_string(x=a_name, y=b_name)) +
        geom_point(aes(color=P_value, size=Jaccard_Index)) +
        scale_size_area(max_size=20) +
        scale_color_gradientn(colors=(brewer.pal(n=9, name="YlOrRd"))) +
        xlab(x_lab) + ylab(y_lab) +
        labs(color=color_lab, size=size_lab) +
        theme_classic() +
        theme(plot.title = element_text(hjust=0.5, face="bold", size=title_size),
              axis.title = element_text(face="bold", size=axis_title_size),
              axis.text.x = element_text(angle=x_text_angle, hjust=x_text_hjust, size=x_text_size),
              axis.text.y = element_text(size=y_text_size),
              legend.title = element_text(face="bold", size=legend_title_size),
              legend.text = element_text(size=legend_text_size),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(size = 0.5, linetype = 'dashed', colour = "grey"),
              plot.margin = margin(10,10,10,10))
}

## Heatmap for metaprogram scores
## @para: nmf_score, metaprogram score (row = cells, col = scores of each metaprogram)
## @para: annotation, metaprogram annotation of each cell
## @para: metagene_order, order of metaprograms to plot
## @para: cutoff, upper and lower cutoff to trim values, default = 2.5
## @para: col_pal, color pallete to use, default = "RdBu"
## @para: col_range, range of color to use from the color pallete, default = 1:9
## @para: out_dir, name of the output directory
## @para: out_name, name of the output figure file, default = "nmf_cell_score_sorted.png"
## @para: font_size, font size of labels of heatmap, default = 16
plotMetaScore <- function(nmf_score, annotation, metagene_order,
                          cutoff=2.5, col_pal="RdBu", col_range=1:9,
                          out_dir = nmf_fig_folder,
                          out_name = "nmf_cell_score_sorted.png",
                          font_size = 16){

    ## Sort cells based on expressed program and make a side bar
    cell_order = annotation
    names(cell_order) = rownames(nmf_score)
    #cell_order = factor(cell_order, levels = metagene_order)
    #cell_order = cell_order[order(cell_order)]

    ## Sort cells and metaprograms (cells: expressed program; metaprograms: specified order)
    nmf_score_sorted = nmf_score[names(cell_order), metagene_order]
    nmf_score_sorted = as.matrix(nmf_score_sorted)

    ## Trim values
    nmf_score_sorted_heatmap = ifelse(nmf_score_sorted > cutoff, cutoff,
                                      ifelse(nmf_score_sorted < -cutoff, -cutoff,
                                             nmf_score_sorted))

    ## heatmap for NMF scores
    hm_colors = rev(brewer.pal(n=9, name=col_pal))[col_range]
    hm_colors = colorRampPalette(colors = hm_colors)
    png(filename = file.path(out_dir, out_name), units = 'in', width = 15, height = 8, res = 300)
    aheatmap(t(nmf_score_sorted_heatmap),
             color = hm_colors(100),
             Rowv = NA, Colv = NA, labCol = NA,
             annCol = cell_order, fontsize = font_size)
    dev.off()
}

## Heatmap for metaprogram scores
## @para: cm_center, centered count matrix (row = genes, col = cells)
## @para: nmf_score, metaprogram score (row = cells, col = scores of each metaprogram)
## @para: annotation, metaprogram annotation of each cell
## @para: metagene_order, order of metaprograms to plot
## @para: cutoff, upper and lower cutoff to trim values, default = 2.5
## @para: col_pal, color pallete to use, default = "RdBu"
## @para: col_range, range of color to use from the color pallete, default = 1:9
## @para: out_dir, name of the output directory
## @para: out_name, name of the output figure file, default = "nmf_cell_score_sorted.png"
## @para: font_size, font size of labels of heatmap, default = 16
plotMetaGeneExpr <- function(cm_center, nmf_score, annotation, metagene_order,
                             nmf_genes = nmf_marker_genes_final,
                             cutoff=2.5, col_pal="RdBu", col_range=1:9,
                             out_dir = nmf_fig_folder,
                             out_name = "nmf_marker_gene_expr.png",
                             font_size = 16){

    cell_order = annotation
    names(cell_order) = rownames(nmf_score)
    cell_order = factor(cell_order, levels = metagene_order)
    cell_order = cell_order[order(cell_order)]

    ## heatmap for expressions of genes in merged NMF factors
    hm_colors = rev(brewer.pal(n=9, name=col_pal))[col_range]
    hm_colors = colorRampPalette(colors = hm_colors)
    cm_center = as.matrix(cm_center)
    cm_heatmap = ifelse(cm_center > cutoff, cutoff,
                        ifelse(cm_center < -cutoff, -cutoff, cm_center))

    ## Only keep genes that are also in this dataset
    nmf_marker_genes = unlist(nmf_genes[metagene_order])
    nmf_marker_genes = nmf_marker_genes[nmf_marker_genes %in% rownames(cm_center)]

    ## Prepare row annotations
    nmf_marker_gene_vec = gsub('[0-9]+', '', names(nmf_marker_genes))

    ## Plot
    png(filename = file.path(out_dir, out_name), units = 'in', width = 15, height = 8, res = 300)
    aheatmap(cm_heatmap[nmf_marker_genes, names(cell_order)],
             color = hm_colors(100), labCol = NA,
             Rowv = NA, Colv = NA,
             annRow = factor(nmf_marker_gene_vec), annCol = cell_order)
    dev.off()
}


## Run interactive shiny plot
## Input is the dataframe, pre-plotted ggplot, and the names of x/y axes
RunShiny<- function(df, myPlot, x, y){
  shinyApp(
    ui = basicPage(
      downloadButton('download',"Download the data"),
      fluidRow(
        column(width = 4,
               plotOutput("plot", height=300,
                          click = "plot_click",  # Equiv, to click=clickOpts(id="plot_click")
                          hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
                          brush = brushOpts(id = "plot_brush")
               ),
               h4("Clicked points"),
               tableOutput("plot_clickedpoints"),
               h4("Brushed points"),
               tableOutput("plot_brushedpoints")
        ),
        column(width = 4,
               verbatimTextOutput("plot_clickinfo"),
               verbatimTextOutput("plot_hoverinfo")
        ),
        column(width = 4,
               wellPanel(actionButton("newplot", "New plot")),
               verbatimTextOutput("plot_brushinfo")
        )
      )
    ),
    server = function(input, output, session) {
      data <- reactive({
        df
      })
      output$plot <- renderPlot({
        myPlot
      })
      output$plot_clickedpoints <- renderTable({
        # For base graphics, we need to specify columns, though for ggplot2,
        # it's usually not necessary.
        res <- nearPoints(data(), input$plot_click, x, y)
        if (nrow(res) == 0)
          return()
        res
      })
      output$plot_brushedpoints <- renderTable({
        res <- brushedPoints(data(), input$plot_brush, x, y)
        if (nrow(res) == 0)
          return()
        res
      }

      )
      output$download <- downloadHandler(
        filename = function(){"thename.csv"},
        content = function(fname){
          write.csv(thedata(), fname)
        } )
    }
  )
}


## Customized ggplot2 theme for plotting
seurat_theme <- function(){
  theme_bw() +
    theme(panel.background = element_rect(colour = "black", size=0.1),
          plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
          axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}


## Colors for plotting
# clusterExperiment::bigPalette
colors_to_use <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3", "#bd18ea", "cyan", "#B2DF8A", "#FB9A99",
                   "deeppink4", "#00B3FFFF", "#CAB2D6", "#FFFF99", "#05188a", "#CCFF00FF", "cornflowerblue", "#f4cc03", "black", "blueviolet", "#4d0776",
                   "maroon3", "blue", "#E5D8BD", "cadetblue4", "#e5a25a", "lightblue1", "#F781BF", "#FC8D62", "#8DA0CB", "#E78AC3", "green3",
                   "#E7298A", "burlywood3", "#A6D854", "firebrick", "#FFFFCC", "mediumpurple", "#1B9E77", "#FFD92F", "deepskyblue4", "yellow3", "#00FFB2FF",
                   "#FDBF6F", "#FDCDAC", "gold3", "#F4CAE4", "#E6F5C9", "#FF00E6FF", "#7570B3", "goldenrod", "#85848f", "lightpink3", "olivedrab", "cadetblue3")


## Generates density plot for specific gene or genes
## Input is the processed Seurat object, the reduction type and features names
d.plot <- function(sobj, reduction = 'tsne', slot = 'data', features) {
  rd <- ifelse(reduction == 'tsne', 'tSNE', 'UMAP')

  cell_embeddings <- Embeddings(sobj, reduction)

  exp_data <- FetchData(sobj, vars = features, slot = slot)

  res <- apply(exp_data, 2, Nebulosa:::calculate_density, cell_embeddings, "wkde", 1)
  res <- apply(res, 2, function(x) scales::rescale(x, to = c(0,1)))

  if (ncol(res) == 1) {
    p <- ggplot(data.frame(cell_embeddings, feature = res[,1])) +
      aes_string(paste0(rd, '_1'), paste0(rd, '_2'), color = "feature") +
      geom_point(shape = 16, size = 1) +
      labs(x = paste0(rd, '-1'), y = paste0(rd, '-2'), title = features, color = guide_legend('Density')) +
      seurat_theme() +
      scale_color_viridis_c(option = 'inferno')
  } else {
    plotList <- list()
    for (i in 1:ncol(res)) {
      plotList[[i]] <- ggplot(data.frame(cell_embeddings, feature = res[,i])) +
        aes_string(paste0(rd, '_1'), paste0(rd, '_2'), color = "feature") +
        geom_point(shape = 16, size = 1) +
        labs(x = paste0(rd, '-1'), y = paste0(rd, '-2'), title = colnames(res)[i], color = guide_legend('Density')) +
        seurat_theme() +
        scale_color_viridis_c(option = 'inferno')
    }
    p <- ggarrange(plotlist = plotList, common.legend = T, legend = 'right')
  }
  return(p)
}


## Monocle plot with scaled expression
## Input is the cds object and the gene or genes name
plotMonocle <- function(cds, gene) {
  if (sum(gene %in% rownames(cds)) == 0) {
    stop('None gene found in dataset')
  }

  if (sum(gene %in% rownames(cds)) != length(gene)) {
    gene <- gene[gene %in% rownames(cds)]
  }

  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(monocle))
  suppressPackageStartupMessages(library(scales))
  suppressPackageStartupMessages(library(ggpubr))
  suppressPackageStartupMessages(library(viridis))

  return_rotation_mat <- function(theta) {
    theta <- theta/180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),
           nrow = 2)
  }
  monocle_theme_opts <- function() {
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
      theme(panel.border = element_blank()) +
      theme(axis.line.x = element_line(size=0.25, color="black")) +
      theme(axis.line.y = element_line(size=0.25, color="black")) +
      theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
      theme(panel.background = element_rect(fill='white')) +
      theme(legend.key=element_blank())
  }


  tmp <- cds@assayData$exprs[gene, ]
  if (length(gene) == 1) {
    tmp <- rescale(tmp, to = c(-2,2))
    cds[[gene]] <- tmp
  } else {
    tmp <- apply(tmp, 1, function(x) rescale(x, to = c(-2,2)))
    for (i in 1:ncol(tmp)) {
      cds[[colnames(tmp)[i]]] <- tmp[,i]
    }
  }

  pt <- plot_cell_trajectory(cds, color_by = gene[1])

  reduced_dim_coords <- reducedDimK(cds)
  ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>%
    select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
    mutate(sample_name = rownames(.), sample_state = rownames(.))
  dp_mst <- minSpanningTree(cds)

  edge_df <- dp_mst %>% igraph::as_data_frame() %>% select_(source = "from", target = "to") %>% left_join(ica_space_df %>% select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>% left_join(ica_space_df %>% select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")

  rot_mat <- return_rotation_mat(0)
  cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
  cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
  edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
  edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)


  data_df <- pt$data

  g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))

  g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1",
                                   y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1",
                                   yend = "target_prin_graph_dim_2"), size = 0.75,
                        linetype = "solid", na.rm = TRUE, data = edge_df)

  mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
  branch_point_df <- ica_space_df %>% dplyr::slice(match(mst_branch_nodes, sample_name)) %>% mutate(branch_point_idx = seq_len(n()))
  g <- g + geom_point(aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"), size = 5, na.rm = TRUE, branch_point_df) +
    geom_text(aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2", label = "branch_point_idx"), size = 4, color = "white", na.rm = TRUE, branch_point_df)

  g <- g + monocle_theme_opts() + xlab("Component 1") + ylab("Component 2") +
    theme(legend.position = "top", legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill = "white"))

  plotlist <- list()
  for (i in 1:length(gene)) {
    plotlist[[i]] <- g + geom_point(data = data_df[which(data_df[[gene[i]]] < 0), ], aes_string(color = paste0('`', gene[i], '`')), size = I(1), na.rm = TRUE) +
      geom_point(data = data_df[which(data_df[[gene[i]]] > 0), ], aes_string(color = paste0('`', gene[i], '`')), size = I(1.5), na.rm = TRUE) +
      scale_color_viridis(option = 'C', discrete = F, end = 0.9) + ggtitle(gene[i]) +
      theme(plot.title = element_text(hjust = 0.5)) + labs(color = "")
  }

  pt2 <- ggarrange(plotlist = plotlist, common.legend = T)

  return(pt2)
}


## Plot proportions for a liger object
## Input is a liger object
plot_integrated_clusters <- function (ligerObj) {
  ## take a liger object, plot distributions over samples
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)

  count_table <- table(as.character(ligerObj@clusters), ligerObj@cell.data$dataset)
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)

  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)

  sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)

  colnames(melt_mtx)[2] <- "dataset"

  p1 <- ggplot(cluster_size, aes(y = cluster,x = value)) +
    geom_bar(position = "dodge", stat = "identity",fill = "grey60") +
    theme_bw() +
    scale_x_log10() +
    labs(x = "Cells per cluster, log10 scale", y = "")
  p2 <- ggplot(melt_mtx, aes(x = cluster,y = value, fill = dataset)) +
    geom_bar(position = "fill", stat = "identity") +
    theme_bw() +
    coord_flip() +
    scale_fill_manual(values = colors_to_use) +
    labs(x = "Factor number", y = "Fraction of cells in each dataset") +
    theme(legend.position = "top") +
    scale_y_continuous(labels = scales::percent, expand = c(0,0))

  p2 + p1 + plot_layout(widths = c(3,1))
}


vln_custom <- function(object, feature, y_int = NULL, pt_size = -1) {
  p1 <- VlnPlot(object = object, features = feature, pt.size = pt_size, cols = "gray80") +
    NoLegend() + ylab(NULL) + xlab(NULL) +
    theme(title = element_text(size = 8, face = "plain"),
          axis.text.x = element_blank())

  if (!is.null(y_int)) p1 + geom_hline(yintercept = y_int, color = "red")
  else p1

}


vln_fun <- function(object, clust_df, criterion) {
  # plot the labels at a value slightly below the max, to ensure they're shown
  # within plot limits
  clust_df$y <- max(object[[criterion]]) * 0.95
  Seurat::VlnPlot(object, criterion, cols = object@misc$colours, pt.size = 0.2) +
    theme(legend.position = "none") +
    geom_text(data = clust_df, aes(label = N, x = ident, y = y)) +
    xlab("Cluster")
}


scatter_custom <- function(object, feature, var1 = 'nCount_RNA', var2 = 'nFeature_RNA', filt = "pre") {
  object@meta.data %>%
    ggplot(aes_string(x=var1, y=var2, color=feature)) +
    geom_point(size = 0.8) +
    scale_colour_manual(values = colors_to_use) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    ggtitle(feature) +
    theme(plot.title = element_text(hjust = 0.5))
}





# Barplot function
plot_bar <- function(seurat_obj, x_var, y_var, colors){
  ggplot(seurat_obj@meta.data, aes(x_var, fill = y_var)) +
    scale_fill_manual(values = colors) + 
    geom_bar(position = "fill", color="black") +
    labs (y='Proportion', x='') + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
          axis.text.y = element_text(size=12, colour="black"),
          axis.title=element_text(size=12)) 
}

