#' global.R

# Copyright (C) Carlos Biagi Jr
#
# Tis is a free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# This software is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, see <http://www.gnu.org/licenses/>.

options(shiny.usecairo=TRUE)
dpi <- 96

################################################################################
## Load packages
################################################################################
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(glue)
  library(tidyverse)
})


################################################################################
#> Ependymoma
################################################################################
dataset_reductions <- readRDS('data/reductions.rds')
nCount_RNA <- dataset_reductions$nCount_RNA
dataset_reductions <- dataset_reductions |> dplyr::select(-nCount_RNA)

counts <- readRDS('data/counts.rds')
hierarchy <- readRDS('data/hierarchy.rds')

colors_Sample <- c('STEPN-01' = "#5e60ce", 'STEPN-02' = "#511f73", 'STEPN-03' = "#0377a8", 'STEPN-04' = "#aacc00", 'STEPN-05' = "#bf99f2", 'STEPN-06' = "#f6cacc", 'STEPN-07' = "#2fb5c7", 'STEPN-08' = "#fec89a", 'STEPN-09' = "#f15bb5", 'STEPN-10' = "#8ecae6", 'STEPN-11' = "#ffa200", 'STEPN-12' = "#f48c06", 'STEPN-13' = "#ff5ca5", 'STEPN-14' = "#8be8d7", 'STEPN-15' = "#a564d3", 'STEPN-18' = "#00bbf9", 'STEPN-19' = "#ef476f", 'STEPN-20' = "#ffba08", 'STEPN-21' = "#ff6000", 'STEPN-22' = "#b4fadc", 'STEPN-23' = "#002855", 'STEPN-24' = "#06d6a0", 'STEPN-25' = "#6d597a", 'STEPN-26' = "#ff0000", 'STEPN-27' = "#ffb9d8", 'STEPN-28' = "#eeef20", 'STEPN-29' = "#ffa2cb", 'STEPN-30' = "#10451d", 'STEPN-31' = "#bc1f66", 'STEPN-32' = "#073b4c", 'STEPN-33' = "#ffea00", 'STEPN-34' = "#007f5f", 'STEPN-35' = "#370617", 'STEPN-36' = "#b9fbc0", 'STEPN-37' = "#023e7d", 'STEPN-38' = "#ff0072", 'STEPN-39' = "#9b5de5", 'STEPN-40' = "#ebcfbc", 'STEPN-41' = "#da5552", 'STEPN-42' = "#cdb4db", 'STEPN-43' = "#51ccd1", 'STEPN-44' = "#dcb9a1")
colors_Metaprograms <- c("Cycling" = "gray30", "Neuroepithelial-like" = "#F99E93FF", "Radial-glia-like" = "#9E5E9BFF", "Embryonic-neuronal-like" = "#74ADD1FF", "Neuronal-like" = '#0F4F8B' ,"Ependymal-like" = "#ACD39EFF", "MES/Hypoxia" = "#96410EFF", "Embryonic-like" = "mistyrose1")
colors_Subtype <- c('ZFTA-RELA' = "#B44672",'ZFTA-Cluster 1' = "#B47846", 'ZFTA-Cluster 2' = "#46B478", 'ZFTA-Cluster 3' = "#46B4AF", 'ZFTA-Cluster 4' = "#4682B4",'ST-YAP1' = "#B4AF46")
colors_Cycling <- c("Non-cycling" = "black", 'Cycling' = "red")
colors_Gender <- c("Female" = "#D95F02", "Male" = "#1B9E77", "Unknown" = "#999999")
colors_Sampling <- c("Primary" = "#1B9E77", "Recurrence" = "#D95F02", "Unknown" = "#999999")
colors_Fusion <- c("CCND1::RP11; KDM2A::POLA2" = "#4E79A7", "EPN-YAP" = "#59A14F", "SYVN1::MAJIN; NFIB::VPS51" = "#E15759", "ZFTA::MAML2" = "#F28E2B", "ZFTA::NCOA2" = "#B07AA1", "ZFTA::RELA" = "#EDC948", "ZFTA::SS18" = "#76B7B2")
colors_NewClassification <- c("No prediction" = "#BDBDBD", "No prediction/ZFTA-Cluster 4 (subclass E)" = "#878787", "ST-YAP1" = "#1F78B4", "ZFTA-Cluster 1 (subclass B)" = "#33A02C", "ZFTA-Cluster 2 (subclass C)" = "#FF7F00", "ZFTA-Cluster 3 (subclass D)" = "#E31A1C", "ZFTA-Cluster 4 (subclass E)" = "#6A3D9A", "ZFTA-RELA (subclass A)" = "#B15928")
colors_sc_snRNAseq <- c("scRNA-seq" = "#E69F00", "snRNA-seq" = "#0072B2")



Reductions_Shiny <- function(df, red, anno, text_size_scale = 1) {
  colors_to_use <- switch (anno,
                           'Metaprogram' = colors_Metaprograms,
                           'Metaprogram_noCC' = colors_Metaprograms,
                           'Sample' = colors_Sample,
                           'Subtype' = colors_Subtype, 
                           'Gender' = colors_Gender,
                           'Sampling' = colors_Sampling,
                           'Fusion' = colors_Fusion,
                           'Cycling_prop' = colors_Cycling,
                           'New_classification' = colors_NewClassification,
                           'sc_snRNAseq' = colors_sc_snRNAseq,
                           NA
  )
  
  pt <- ggplot(df, mapping = aes(x = !!sym(glue('{red}_1')), y = !!sym(glue('{red}_2')), colour = !!sym(anno))) + 
    geom_point(size = 0.5) + 
    scale_color_manual(values = colors_to_use) + 
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) + 
    theme_bw() +
    theme(panel.background = element_rect(colour = "black", linewidth = 0.1),
          plot.title = element_text(size = 14 * text_size_scale, hjust = 0.5, angle = 0, face = "bold", vjust = 1),
          axis.ticks.length=unit(.2, "cm"), 
          axis.text = element_text(size = 12 * text_size_scale, family = NULL, colour = "black"), ,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          text = element_text(size = 12 * text_size_scale, family = NULL, color = "black"), 
          axis.title = element_text(size = 13 * text_size_scale, family = NULL, colour = "black"), 
          strip.text = element_text(size = 12.5 * text_size_scale, family = NULL, colour = "black", hjust = 0.5, margin = margin(3, 3, 3, 3)), 
          legend.title = element_text(size = 12 * text_size_scale, family = NULL, colour = "black", hjust = 0),
          legend.text = element_text(size = 11 * text_size_scale, family = NULL, colour = "black"), 
          plot.subtitle = element_text(size = 13 * text_size_scale, family = NULL, hjust = 0, margin = margin(b = 3))) + 
    labs(x = toupper(glue('{red}-1')), y = toupper(glue('{red}-2')))
  
  return(pt)
}


FeaturePlot_Shiny <- function(counts, features, red) {
  normalized_values <- log1p((t(as.matrix(counts[features,,drop=F])) / nCount_RNA) * 10000) |> as.matrix()
  normalized_values <- cbind(dataset_reductions, normalized_values)
  
  pt_list <- lapply(features, function(x) {
    normalized_values |> arrange(!!sym(x)) |> 
      ggplot(aes(x = !!sym(glue('{red}_1')), y = !!sym(glue('{red}_2')), color = !!sym(x))) + 
      geom_point(size = 0.5) + 
      scale_color_gradient(low = "lightgrey", high = "darkblue") +
      labs(x = toupper(glue('{red}-1')), y = toupper(glue('{red}-2')), color = '', title = x) + 
      theme_bw() +
      theme(panel.background = element_rect(colour = "black", linewidth = 1),
            plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
            axis.ticks.length = unit(0, "cm"), axis.text = element_text(size = 0),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            plot.caption = element_text(hjust = 0.5, size = 16))
  })
  plot <- ggarrange(plotlist = pt_list)
  
  return(plot)
}


FastRowScale_Seurat <- function(mat, center = TRUE, scale = TRUE, scale_max = 10) 
{
  if (center) {
    rm <- matrixStats::rowMeans2(x = mat, na.rm = TRUE)
  }
  if (scale) {
    if (center) {
      rsd <- matrixStats::rowSds(mat, center = rm)
    }
    else {
      rsd <- sqrt(x = rowSums2(x = mat^2)/(ncol(x = mat) -1))
    }
  }
  if (center) {
    mat <- mat - rm
  }
  if (scale) {
    mat <- mat/rsd
  }
  if (scale_max != Inf) {
    mat[mat > scale_max] <- scale_max
  }
  return(mat)
}



HierarchyPlot_Shiny <- function(counts, features = NULL, includeOriginal = FALSE, subgroups = NULL, originalOnly = FALSE) {
  if (!is.null(subgroups)) {
    pt_list <- lapply(subgroups, function(x) {
      tmp <- hierarchy |> mutate(tmp = ifelse(subtype == x, x, 'other'))
      
      ggplot(tmp |> arrange(desc(tmp))) +
        geom_point(aes(x = set_x, y = set_y,color = tmp), size = 1) + 
        scale_colour_manual(values = c('grey90', unname(colors_Subtype[which(names(colors_Subtype) == x)])), name = "") + 
        labs(title = x) +
        geom_vline(xintercept = 0, linetype = "dotted") + 
        geom_hline(yintercept = 0, linetype = "dotted") + 
        labs(x = '', y = '') +
        theme_bw() +
        theme(panel.background = element_rect(colour = "black", linewidth = 0.1),
              plot.title = element_text(hjust = 0.5, angle = 0, size = 11, face = "bold", vjust = 1),
              axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              legend.position = "none")
    })
    
    combined_plot <- ggarrange(plotlist = pt_list, ncol = 3, nrow = 2, legend = "none")
    combined_plot <- annotate_figure(
      combined_plot,
      left = text_grob("Neuronal-like  <----->  Neuroepithelial-like", rot = 90, vjust = 0.5, hjust = 0.5, size = 14, face = 'bold'),
      right = text_grob("Embryonic-like  <----->  Ependymal-like", rot = 270, vjust = -1, hjust = 0.5, size = 14, face = 'bold'),
      bottom = text_grob("Neuronal-like  <----->  Ependymal-like", vjust = 0.5, hjust = 0.5, size = 14, face = 'bold'), 
      top = text_grob("Neuroepithelial-like  <----->  Embryonic-like", vjust = -0.5, hjust = 0.5, size = 14, face = 'bold')
    ) + theme(plot.margin = margin(10, 15, 5, 5))
    
    return(combined_plot)
    
  } else if (originalOnly) {
    pt <- ggplot(hierarchy, aes(x = set_x, y = set_y, color = Metaprogram)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors_Metaprograms) +
      geom_vline(xintercept = 0, linetype = "dotted") + 
      geom_hline(yintercept = 0, linetype = "dotted") + 
      labs(color = '') + 
      scale_y_continuous(name = "Neuronal-like  <----->  Neuroepithelial-like", sec.axis = sec_axis(~ ., name = "Embryonic-like  <----->  Ependymal-like")) +
      scale_x_continuous(name = "Neuronal-like  <----->  Ependymal-like", sec.axis = sec_axis(~ ., name = "Neuroepithelial-like  <----->  Embryonic-like")) + 
      guides(color = guide_legend(override.aes = list(size = 3))) +
      theme_bw() +
      theme(panel.background = element_rect(colour = "black", linewidth = 0.1),
            plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
            axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.position = "bottom", legend.key = element_blank())
    return(pt)
  } else {
    ## Normalize data first
    normalized_values <- t(log1p((t(as.matrix(counts[features,,drop=F])) / nCount_RNA) * 10000))
    
    ## Then: Scale data
    object <- normalized_values[features, , drop = FALSE] |> as.matrix()
    object.names <- dimnames(x = object)
    
    scaled.data <- matrix(data = NA_real_, nrow = nrow(x = object), ncol = ncol(x = object), dimnames = object.names)
    
    arg.list <- list(mat = object[features, colnames(object), drop = FALSE], scale = T, center = T, scale_max = 10, display_progress = FALSE)
    arg.list <- arg.list[intersect(x = names(x = arg.list), y = names(x = formals(fun = FastRowScale_Seurat)))]
    
    data.scale <- do.call(what = FastRowScale_Seurat, args = arg.list)
    dimnames(data.scale) <- dimnames(object[features, colnames(object)])
    scaled.data[features, colnames(object)] <- data.scale
    rm(data.scale)
    
    dimnames(x = scaled.data) <- object.names
    scaled.data[is.na(x = scaled.data)] <- 0
    
    mat <- cbind(hierarchy, t(scaled.data))
    
    ## Include original hierarchy
    if (includeOriginal) {
      pt_original <- ggplot(mat, aes(x = set_x, y = set_y, color = subtype)) +
        geom_point(size = 1) +
        scale_color_manual(values = colors_Subtype) +
        geom_vline(xintercept = 0, linetype = "dotted") + 
        geom_hline(yintercept = 0, linetype = "dotted") + 
        labs(color = '') + 
        scale_y_continuous(name = "Neuronal-like  <----->  Neuroepithelial-like", sec.axis = sec_axis(~ ., name = "Embryonic-like  <----->  Ependymal-like")) +
        scale_x_continuous(name = "Neuronal-like  <----->  Ependymal-like", sec.axis = sec_axis(~ ., name = "Neuroepithelial-like  <----->  Embryonic-like")) + 
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme_bw() +
        theme(panel.background = element_rect(colour = "black", linewidth = 0.1),
              plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
              axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              legend.position = "bottom", legend.key = element_blank())
    }
    
    ## Generate Hierarchy by gene expression individually
    if (length(features) > 1) {
      pt_gene_list <- lapply(features, function(x) {
        ggplot(mat |> arrange(!!sym(x)), aes(x = set_x, y = set_y)) +
          geom_point(aes(fill = .data[[x]]), color = "black", shape = 21, size = 2, stroke = 0) +
          scale_fill_gradient2(low = '#00008b', mid = 'oldlace', high = '#8b0000', midpoint = 0, guide = guide_colorbar(ticks.colour = "black", frame.colour = "black"), name = "Score") + 
          geom_vline(xintercept = 0, linetype = "dotted") +
          geom_hline(yintercept = 0, linetype = "dotted") + 
          theme_bw() +
          theme(panel.background = element_rect(colour = "black", linewidth = 0.1),
                plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
                axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
          labs(title = x, x = '', y = '')
        
      })
      
      if (includeOriginal) {
        pt_gene_list <- c(list(pt_original), pt_gene_list)
        pt_combined <- ggarrange(plotlist = pt_gene_list)
        return(pt_combined)
      } else {
        pt_combined <- ggarrange(plotlist = pt_gene_list)
        return(pt_combined)
      }
      
    } else {
      pt_gene <- ggplot(mat |> arrange(!!sym(features)), aes(x = set_x, y = set_y)) +
        geom_point(aes(fill = .data[[features]]), color = "black", shape = 21, size = 2, stroke = 0) +
        scale_fill_gradient2(low = '#00008b', mid = 'oldlace', high = '#8b0000', midpoint = 0, guide = guide_colorbar(ticks.colour = "black", frame.colour = "black"), name = "Score") + 
        geom_vline(xintercept = 0, linetype = "dotted") +
        geom_hline(yintercept = 0, linetype = "dotted") + 
        theme_bw() +
        theme(panel.background = element_rect(colour = "black", linewidth = 0.1),
              plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
              axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        labs(x = '', y = '', title = features)
      
      if (includeOriginal) {
        pt_combined <- ggarrange(pt_original, pt_gene, ncol = 2)
        return(pt_combined)
      } else {
        return(pt_gene)
      }
    }
  }
}



################################################################################
## Functions for ui.R to create sections
################################################################################
## Generate download section for each type of plot
generateDownloadUI <- function(tabCondition, inputCondition, outputId) {
  conditionalPanel(
    condition = paste0("input.tab == ", tabCondition, " && ", inputCondition),
    wellPanel(
      h4("Download:"),
      selectInput(inputId = "downloadPlotFileType", label = strong("Select download file type"),
                  choices = list("PDF" = "pdf", "TIFF" = "tiff", "PNG" = "png")
      ),
      br(),
      strong("Set download image dimensions"),
      helpText("(units are inches for PDF, pixels for all other formats)"),
      div(class="row",
          div(class="col-xs-6",
              numericInput(inputId = "downloadPlotWidth", label = "Width (in)", value = 7, min = 1, max = 100)
          ),
          div(class="col-xs-6",
              numericInput(inputId = "downloadPlotHeight", label = "Height (in)", value = 7, min = 1, max = 100)
          )
      ),
      conditionalPanel(
        condition = "input.downloadPlotFileType != 'pdf'",
        div(class="row",
            div(class="col-xs-6",
                numericInput(inputId = "downloadPlotRes", label = "Resolution (ppi)", 
                             value = 72, min = 72, max = 600)
            )
        )
      ),
      br(),
      downloadButton(outputId = outputId, label = "Download", class= "btn-primary")
    )
  )
}

