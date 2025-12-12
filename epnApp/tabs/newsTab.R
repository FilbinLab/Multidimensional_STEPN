# UI-elements for News tab
tabPanel(title = "News", icon = icon("newspaper"),
         div(id = "home",
             
             h4(class = "outer", "v1.0 (12-12-2025)"),
             h5(class = "outer", "Home:"),
             tags$ul(
               tags$li('Description of the App and information related to citation, number of samples, where to obtain processed data, and contact details.')
             ),
             h5(class = "outer", "Reductions:"),
             tags$ul(
               tags$li('Addition of UMAP and tSNE reductions;'),
               tags$li('Annotations available: Metaprogram, Metaprogram_noCC, Sample, Subtype, Gender, Sampling, Fusion, Cycling_prop, New_classification, sc_snRNAseq.')
             ),
             h5(class = "outer", "Feature Plot:"),
             tags$ul(
               tags$li('Addition of UMAP and tSNE reductions;'),
               tags$li('10k most variable genes available to plot.')
             ),
             h5(class = "outer", "Cell State Plot:"),
             tags$ul(
               tags$li('Gene Expression: 10k most variable genes available;'),
               tags$li('Subgroup: coloring by subgroup;'),
               tags$li('Metaprogram: coloring by metaprograms.')
             ),
             
             hr()
         )
)