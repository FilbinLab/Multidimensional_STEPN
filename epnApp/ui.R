shinyUI(
  fluidPage(
    
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    
    #shinyjs::useShinyjs(),
    
    navbarPage(title = strong("STEPN"), 
               windowTitle = "Multidimensional STEPN", 
               fluid = TRUE, 
               id = "nav",
               ## Tabs
               source("tabs/homeTab.R", local = TRUE)$value,
               source("tabs/exploreTab.R", local = TRUE)$value,
               source("tabs/newsTab.R", local = TRUE)$value,
    )
  )
)
