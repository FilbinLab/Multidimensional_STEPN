# UI-elements for Explore tab
tabPanel(title = "Explore", icon = icon("magnifying-glass-chart"), id = "explore",
         
         # link to the style.css file.
         tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
         
         sidebarLayout(fluid = T,
                       sidebarPanel(
                         #img(src = "GlioSyn_logo.png", class="responsive-image"),
                         #br(),
                         wellPanel(     
                           
                           ## Side tab for Reductions
                           conditionalPanel(
                             condition = "input.tab == 1",
                             
                             selectizeInput(inputId = "reductionType", label = "Select a reduction:", choices = NULL, selected = NULL, multiple = F, options = list(placeholder = 'Select a reduction', plugins = list('remove_button', 'drag_drop'))),
                             br(),
                             
                             selectizeInput(inputId = "reductionAnnotation", label = "Select an annotation:", choices = NULL, selected = NULL, multiple = F, options = list(placeholder = 'Select an annotation', plugins = list('remove_button', 'drag_drop'))),
                             br(),
                             
                             wellPanel( 
                               h4("Plot options:"),
                               p(style = "background-image: linear-gradient(to right, #c3d2e6 0%, #d0dccc 100%); padding-left: 10px; border: 1px solid #E3E3E3;", 
                                 strong("Plot size")
                               ),
                               div(class="row",
                                   div(class="col-xs-6",
                                       numericInput("plot_width_reductions", label = "Width (in):", min = NA, max = 15, step = 0.5, value = 8)
                                   ),
                                   div(class="col-xs-6",
                                       numericInput("plot_height_reductions", label = "Height (in):", min = 3, max = 15, step = 0.5, value = 6)
                                   )
                               ),
                               
                             ),
                           ),
                           
                           
                           ## Side tab for FeaturePlot
                           conditionalPanel(
                             condition = "input.tab == 2",
                             
                             selectizeInput(inputId = "featureplotReduction", label = "Select a reduction:", choices = NULL, selected = NULL, multiple = F, options = list(placeholder = 'Select a reduction', plugins = list('remove_button', 'drag_drop'))),
                             br(),
                             
                             selectizeInput(inputId = "featureplotGenes", label = glue("Select gene(s)"), choices = NULL, selected = NULL, multiple = T, options = list(placeholder = glue('{nrow(counts)} in total'), plugins = list('remove_button', 'drag_drop'))),
                             br(),
                             
                             wellPanel( 
                               h4("Plot options:"),
                               p(style = "background-image: linear-gradient(to right, #c3d2e6 0%, #d0dccc 100%); padding-left: 10px; border: 1px solid #E3E3E3;", 
                                 strong("Plot size")
                               ),
                               div(class="row",
                                   div(class="col-xs-6",
                                       numericInput("plot_width_featureplot", label = "Width (in):", min = NA, max = 15, step = 0.5, value = 8)
                                   ),
                                   div(class="col-xs-6",
                                       numericInput("plot_height_featureplot", label = "Height (in):", min = 3, max = 15, step = 0.5, value = 6)
                                   )
                               ),
                             ),
                           ),
                           
                           
                           ## Side tab for Hierarchy
                           conditionalPanel(
                             condition = "input.tab == 3",
                             
                             # First selection: Hierarchy by gene expression or subgroup
                             selectizeInput(inputId = "typeHierarchy", label = "Select type:", choices = NULL, multiple = F, options = list(placeholder = 'Select type:', plugins = list('remove_button','drag_drop'))),
                             br(),
                             
                             # If typeHierarchy is Subgroup
                             conditionalPanel(
                               condition = "input.typeHierarchy === 'Subgroup'",
                               selectizeInput(inputId = "subgroupsHierarchy", label = "Subgroup:", choices = NULL, multiple = T, options = list(placeholder = 'Subgroup:', plugins = list('remove_button','drag_drop'))),
                               br(),
                             ),
                             
                             # If typeHierarchy is Gene Expression
                             conditionalPanel(
                               condition = "input.typeHierarchy === 'Gene Expression'",
                               selectizeInput(inputId = "genesHierarchy", label = glue("Select gene(s):"), choices = NULL, multiple = T, options = list(placeholder = glue('{nrow(counts)} in total'), plugins = list('remove_button','drag_drop'))),
                               br(),
                               checkboxInput(inputId = "includeOriginalHierarchy", label = "Include original cell state.", value = FALSE),
                               br()
                             ),
                             
                             wellPanel( 
                               h4("Plot options:"),
                               p(style = "background-image: linear-gradient(to right, #c3d2e6 0%, #d0dccc 100%); padding-left: 10px; border: 1px solid #E3E3E3;", 
                                 strong("Plot size")
                               ),
                               div(class="row",
                                   div(class="col-xs-6",
                                       numericInput("plot_width_hierarchy", label = "Width (in):", min = NA, max = 15, step = 0.5, value = 8)
                                   ),
                                   div(class="col-xs-6",
                                       numericInput("plot_height_hierarchy", label = "Height (in):", min = 3, max = 15, step = 0.5, value = 6)
                                   )
                               ),
                             ),
                           ),
                           
                           
                           ## Side tab to include export plot options
                           conditionalPanel(
                             condition = "input.tab == 1 | input.tab == 2 | input.tab == 3",
                             
                             # New conditionalPanel to show the download section for Reduction only when the last selection is done
                             generateDownloadUI(
                               tabCondition = 1,
                               inputCondition = "input.reductionAnnotation !== ''",
                               outputId = "downloadPlot_Reductions"
                             ),
                             
                             # New conditionalPanel to show the download section for FeaturePlot only when the last selection is done
                             generateDownloadUI(
                               tabCondition = 2,
                               inputCondition = "input.featureplotGenes.length > 0",
                               outputId = "downloadPlot_FeaturePlots"
                             ),
                             
                             # New conditionalPanel to show the download section for Hierarchy only when the last selection is done
                             #generateDownloadUI(
                            #   tabCondition = 3,
                            #   inputCondition = "input.genesHierarchy !== null | input.subgroupsHierarchy !== '' | input.typeHierarchy == 'Metaprogram",
                            #   outputId = "downloadPlot_CellState"
                            # ),
                            generateDownloadUI(
                              tabCondition = 3,
                              inputCondition = "
                              (input.typeHierarchy === 'Gene Expression' &&
                              input.genesHierarchy &&
                              input.genesHierarchy.length > 0) ||
                              (input.typeHierarchy === 'Subgroup' &&
                              input.subgroupsHierarchy &&
                              input.subgroupsHierarchy.length > 0) ||
                              (input.typeHierarchy === 'Metaprogram')
                              ",
                              outputId = 'downloadPlot_CellState'
                            ),
                           ),
                         ),
                       ),
                       
                       mainPanel(
                         tabsetPanel(id = "tab",
                                     
                                     tabPanel(title = "Reductions", icon = icon("puzzle-piece"), value = 1,
                                              tabPanel(title = "Reductions", value = "reductions",
                                                       plotOutput(outputId = "reductions", width = '1000px', height = '1200px')
                                              ),
                                     ), 
                                     
                                     tabPanel(title = "Feature Plot", icon = icon("circle-dot"), value = 2,
                                              tabPanel(title = "Feature Plot", value = "featureplot",
                                                       plotOutput(outputId = "featureplot", width = '1000px', height = '1200px')
                                              ),
                                     ),
                                     
                                     tabPanel(title = "Cell State Plot", icon = icon("sitemap"), value = 3,
                                              tabPanel(title = "Cell State Plot", value = "hierarchy",
                                                       plotOutput(outputId = "hierarchy", width = '1000px', height = '1200px')
                                              ),
                                     )
                         )
                       )
         )
)


