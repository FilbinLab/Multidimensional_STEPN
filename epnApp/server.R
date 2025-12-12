#' server.R

# Copyright (C) Carlos Biagi Jr
#
# This is a free software; you can redistribute it and/or modify it under the
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

shinyServer(
  function(input, output, session) {
    
    options(shiny.maxRequestSize=200*1024^2)
    
    
    ############################################################################
    ##### Setting download plot options
    ############################################################################
    #' Get the selected download file type.
    download_plot_file_type <- reactive(input$downloadPlotFileType)
    observe({
      plotUnit <- "in"  # Use inches for all file types
      plotUnitDef <- 7  # Default size in inches
      plotUnitMin <- 1  # Minimum size in inches
      plotUnitMax <- 12 # Maximum size in inches
      plotUnitStep <- 0.1 # Step size in inches
      
      updateNumericInput(session, inputId = "downloadPlotHeight", label = sprintf("Height (%s)", plotUnit), value = plotUnitDef, min = plotUnitMin, max = plotUnitMax, step = plotUnitStep)
      updateNumericInput(session, "downloadPlotWidth", label = sprintf("Width (%s)", plotUnit), value = plotUnitDef, min = plotUnitMin, max = plotUnitMax, step = plotUnitStep)
    })
    
    #' Get the download dimensions.
    download_plot_height <- reactive(input$downloadPlotHeight)
    download_plot_width <- reactive(input$downloadPlotWidth)
    download_plot_res <- reactive(input$downloadPlotRes)
    
    
    ############################################################################
    ##### Reductions
    ############################################################################
    observe({
      updateSelectizeInput(session, inputId = "reductionType", choices = c('UMAP', 't-SNE'), selected = 'UMAP')
    })
    
    reductionType <- reactive({
      input$reductionType
    })
    
    
    observe({
      updateSelectizeInput(session, inputId = "reductionAnnotation", choices = grep('umap|tsne', colnames(dataset_reductions), value = T, invert = T), selected = 'Subtype')
    })
    
    reductionAnnotation <- reactive({
      input$reductionAnnotation
    })
    
    observe({
      updateSliderInput(session, inputId = "text_size_scale_dot")
    })
    
    text_size_scale_dot <- reactive({
      input$text_size_scale_dot
    })
    
    
    dpi <- 96
    plot_width_reductions_px <- reactive(input$plot_width_reductions * dpi)
    plot_height_reductions_px <- reactive(input$plot_height_reductions * dpi)
    
    
    plot_Reductions <- reactive({
      req(input$reductionType, input$reductionAnnotation)
      
      pt <- Reductions_Shiny(df = dataset_reductions, 
                             red = ifelse(input$reductionType == 't-SNE', 'tsne', 'umap'), 
                             anno = input$reductionAnnotation, 
                             text_size_scale = input$text_size_scale_dot)
      
      return(pt)
    })
    
    output$reductions <- renderPlot({
      plot_Reductions()
    }, width = plot_width_reductions_px, height = plot_height_reductions_px)
    
    output$downloadPlot_Reductions <- downloadHandler(
      filename = function() {
        glue("{Sys.Date()}_Reductions.{download_plot_file_type()}")
      },
      
      # The argument content below takes filename as a function and returns what's printed to it.
      content = function(file) {
        # Gets the name of the function to use from the downloadFileType reactive element.
        plotFunction <- match.fun(download_plot_file_type())
        plotFunction(file, width = download_plot_width(), height = download_plot_height())
        if(download_plot_file_type() != "pdf"){
          plotFunction(file, width = download_plot_width(), height = download_plot_height(), res = download_plot_res())
        }
        print(plot_Reductions())
        dev.off()
      }
    )
    
    
    ############################################################################
    ##### FeaturePlots
    ############################################################################
    observe({
      updateSelectizeInput(session, inputId = "featureplotReduction", choices = c('UMAP', 't-SNE'), selected = 'UMAP')
    })
    
    featureplotReduction <- reactive({
      input$featureplotReduction
    })
    
    
    observe({
      updateSelectizeInput(session, inputId = "featureplotGenes", choices = rownames(counts), selected = NULL, server = T)
    })
    
    featureplotGenes <- reactive({
      input$featureplotGenes
    })
    
    
    dpi <- 96
    plot_width_featureplot_px <- reactive(input$plot_width_featureplot * dpi)
    plot_height_featureplot_px <- reactive(input$plot_height_featureplot * dpi)
    
    
    plot_FeaturePlot <- reactive({
      req(input$featureplotReduction, input$featureplotGenes)
      
      pt <- FeaturePlot_Shiny(counts = counts, 
                              features = input$featureplotGenes, 
                              red = ifelse(input$featureplotReduction == 't-SNE', 'tsne', 'umap'))
      
      return(pt)
    })
    
    output$featureplot <- renderPlot({
      plot_FeaturePlot()
    }, width = plot_width_featureplot_px, height = plot_height_featureplot_px)
    
    output$downloadPlot_FeaturePlots <- downloadHandler(
      filename = function() {
        glue("{Sys.Date()}_FeaturePlot.{download_plot_file_type()}")
      },
      # The argument content below takes filename as a function and returns what's printed to it.
      content = function(file) {
        # Gets the name of the function to use from the downloadFileType reactive element.
        plotFunction <- match.fun(download_plot_file_type())
        plotFunction(file, width = download_plot_width(), height = download_plot_height())
        if(download_plot_file_type() != "pdf"){
          plotFunction(file, width = download_plot_width(), height = download_plot_height(), res = download_plot_res())
        }
        print(plot_FeaturePlot())
        dev.off()
      }
    )
    
    
    ############################################################################
    ##### Cell State
    ############################################################################
    observe({
      updateSelectizeInput(session, inputId = "typeHierarchy", choices = c('Metaprogram', 'Subgroup', 'Gene Expression'), selected = 'Metaprogram')
    })
    
    typeHierarchy <- reactive({
      input$typeHierarchy
    })
    
    
    observeEvent(typeHierarchy(), {
      req(input$typeHierarchy)
      if (typeHierarchy() == 'Subgroup') {
        observe({
          updateSelectizeInput(session, inputId = "subgroupsHierarchy", choices = names(colors_Subtype), selected = names(colors_Subtype))
        })
      } else if (typeHierarchy() == 'Gene Expression') {
        updateSelectizeInput(session, "genesHierarchy", choices = rownames(counts), selected = '', server = T)
      }
    })
    
    observe({
      updateCheckboxInput(session, "includeOriginalHierarchy", value = FALSE)
    })
    
    includeOriginalHierarchy <- reactive({
      input$includeOriginalHierarchy
    })
    
    
    dpi <- 96
    plot_width_hierarchy_px <- reactive(input$plot_width_hierarchy * dpi)
    plot_height_hierarchy_px <- reactive(input$plot_height_hierarchy * dpi)
    
    
    plot_Hierarchy <- reactive({
      if (input$typeHierarchy == 'Subgroup') {
        req(input$typeHierarchy, input$subgroupsHierarchy)
        return(HierarchyPlot_Shiny(counts = counts, 
                                   subgroups = input$subgroupsHierarchy))
      } else if (input$typeHierarchy == 'Gene Expression') {
        req(input$typeHierarchy, input$genesHierarchy)
        return(HierarchyPlot_Shiny(counts = counts, 
                                   features = input$genesHierarchy, 
                                   includeOriginal = includeOriginalHierarchy()))
      } else if (input$typeHierarchy == 'Metaprogram') {
        req(input$typeHierarchy)
        return(HierarchyPlot_Shiny(counts = counts, 
                                   originalOnly = TRUE))
      }
    })
    
    output$hierarchy <- renderPlot({
      plot_Hierarchy()
    }, width = plot_width_hierarchy_px, height = plot_height_hierarchy_px)
    
    output$downloadPlot_CellState <- downloadHandler(
      filename = function() {
        glue("{Sys.Date()}_CellState.{download_plot_file_type()}")
      },
      # The argument content below takes filename as a function and returns what's printed to it.
      content = function(file) {
        # Gets the name of the function to use from the downloadFileType reactive element.
        plotFunction <- match.fun(download_plot_file_type())
        plotFunction(file, width = download_plot_width(), height = download_plot_height())
        if(download_plot_file_type() != "pdf"){
          plotFunction(file, width = download_plot_width(), height = download_plot_height(), res = download_plot_res())
        }
        print(plot_Hierarchy())
        dev.off()
      }
    )
})