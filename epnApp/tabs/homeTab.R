# UI-elements for Home tab
tabPanel(
  title = "Home", 
  icon = icon("home"),
  # Main container
  div(id = "home",
      br(),
      p(class = "lead", HTML("<strong><span style='color:#3a69ae;'>Single-cell multidimensional profiling of tumor cell heterogeneity in supratentorial ependymomas</span></strong>")),
      div(class="intro-divider"),
      
      p(class = "outer"," This Shiny application provides an interactive resource to explore single-cell Smart-seq2 data from supratentorial ependymomas. The app focuses on tumor cell heterogeneity and cell states across patients and platforms."),
      br(),
      h4(class = "outer", "Manuscript"),
      tags$ul(
        tags$li(strong("Title:"), 'Single-cell multidimensional profiling of tumor cell heterogeneity in supratentorial ependymomas'),
        tags$li(strong("Current version:"), a("Preprint available on bioRxiv", href="https://www.biorxiv.org/content/10.1101/2024.08.07.607066v1", target="_blank")),
        tags$li(strong("Final peer-reviewed article:"), 'A link to the journal publication will be added here as soon as it becomes available.'),
      ),
      br(),
      h4(class = "outer", "Datasets included"),
      tags$ul(
        tags$li(strong("42 tumor samples"), 'profiled by single-cell RNA sequencing (Smart-seq2).', HTML("<strong><span style='color:#FF0000;'>[This App!]</span></strong>")),
        tags$li('Integrated', strong("single-cell and single-nucleus RNA-seq"), 'datasets.'),
        tags$li('All processed data have been deposited in', strong("GEO:"), a("scRNA-seq (Smart-seq2 and 10x)", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300150", target="_blank"), 'and', a("Spatial transcriptomics (Xenium)", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300146", target="_blank")),
      ),
      br(),
      h4(class = "outer", "Main features of the App"),
      tags$ol(
        tags$li(strong('Reductions viewer:'), 'UMAP or t-SNE representations coloring by multiple annotations (Metaprogram, Sample, Subtype, Gender, Sampling, Fusion, Cycling proportion, New classification, sc vs snRNA-seq)'),
        tags$li(strong('Feature plot:'), 'UMAP or t-SNE colored by the expression of any gene present in the dataset.'),
        tags$li(strong('Cell state plot:'), '(1) Cell state–level visualization of gene expression; (2) Coloring by subgroup and by metaprogram annotation'),
      ),
      br(),
      h4(class = "outer", "Contact"),
      p(class = "outer", 'For questions about the data or the application, please contact:', 
        a("Daeun Jeong", href="mailto:djeong3@mgh.harvard.edu"), ',',
        a("Sara Danielli", href="mailto:sara.g.danielli@gmail.com"), ', or',
        a("Carlos Biagi Junior", href="mailto:carlosa_biagijunior@dfci.harvard.edu")),
  )
)