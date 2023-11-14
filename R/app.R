# Load the required libraries
library(shiny)
library(dplyr)
library(Seurat)
library(colourpicker)  # For color picker input
library(paletteer)
library(ggplot2)
options(shiny.maxRequestSize = 3000*1024^2)  # Set the maximum upload size to 3gb

#' @export 
dim.pals <- function(seurat_object = NULL) {

# Define the UI for the Shiny app
ui <- fluidPage(
  titlePanel("Dim Pals"),
  sidebarLayout(
    sidebarPanel(
      tags$h3("Re-Color Individual Cluster"),
      selectInput("cluster_color", "Select Cluster for Custom Color", choices = NULL),
      colourInput("col", "Select Cluster Color", "purple"),
      actionButton("assign_cluster_color_button", "Assign Cluster Color"),
      tags$h3("Re-Color Groups of Clusters"),
      selectInput("grouping_column", "Select Grouping Column", choices = NULL),
      selectInput("group_value", "Select Group", choices = NULL),
      selectInput("pal_package", "Select Palette Package", choices = unique(palettes_d_names$package)),
      selectInput("pal_name", "Select Palette Name", choices = unique(palettes_d_names$palette)),
      actionButton("assign_group_color_button", "Assign Group Color"),
      actionButton("plot_button", "Generate UMAP Plot")
    ),
    mainPanel(
      plotOutput("umap_plot"),
      verbatimTextOutput("color_list")
    )
  )
)

# Define the server for the Shiny app
server <- function(input, output, session) {
  seurat_object <- reactiveVal(seurat_object)
  color_assignments <- reactiveVal(data.frame(Cluster = character(0), Color = character(0)))
  colors_use <- reactiveVal(c(""))
  
  
  observe({
    seurat <- seurat_object()
    if (!is.null(seurat)) {
      seurat@meta.data$all_clusters <- "All_Clusters"
      seurat@meta.data <- seurat@meta.data %>% select(all_clusters, everything()) %>% arrange(seurat_clusters)
      updateSelectInput(session, "cluster_color", choices = unique(seurat$seurat_clusters))
      updateSelectInput(session, "grouping_column", choices = colnames(seurat@meta.data))
      seurat_object(seurat)
    }
  })
  
  # Observe event for updating palette choices based on the selected palette package
  observe({
    pal_package <- input$pal_package
    if (!is.null(pal_package)) {
      # Filter palettes based on the selected package
      palette_choices <- palettes_d_names %>%
        filter(package == pal_package) %>%
        pull(palette)
      
      updateSelectInput(session, "pal_name", choices = palette_choices)
    }
  })
  
  observeEvent(input$assign_cluster_color_button, {
    cluster <- input$cluster_color
    color <- input$col
    color_assignments_data <- data.frame(Cluster = cluster, Color = color)
    color_assignments_data <- rbind(color_assignments(), color_assignments_data)
    color_assignments(color_assignments_data)
  })
  
  observeEvent(input$assign_group_color_button, {
    seurat <- seurat_object()
    grouping_column <- input$grouping_column
    group_value <- input$group_value
    package <- input$pal_package
    palette_name <- input$pal_name
    
    if (!is.null(grouping_column) && !is.null(group_value)) {
      color_assignments_data <- seurat@meta.data %>% 
        select(seurat_clusters, grouping_column) %>% 
        unique() %>% 
        filter(.data[[grouping_column]] == group_value) %>% 
        select(seurat_clusters) 
        color_assignments_data <-color_assignments_data %>%
        mutate(Color = paletteer::paletteer_d(paste0(package,"::",palette_name), 
                                              nrow(color_assignments_data), 
                                              type = "continuous")) %>% 
        rename("Cluster" = "seurat_clusters")
      
      color_assignments_data <- rbind(color_assignments(), color_assignments_data)
      color_assignments(color_assignments_data)
    }
  })
  
  observeEvent(input$plot_button, {
    seurat <- seurat_object()
    color_assignments_data <- color_assignments()
    
    
    output$umap_plot <- renderPlot({
      if (!is.null(seurat) && nrow(color_assignments_data) > 0) {
        cols <- rep("#000000", length(unique(seurat$seurat_clusters)))
        for (i in 1:nrow(color_assignments_data)) {
          cluster_name <- color_assignments_data$Cluster[i]
          color <- color_assignments_data$Color[i]
          cols[cluster_name] <- color
        }
        colors_use(cols)
        DimPlot(seurat, group.by = "seurat_clusters", cols = cols) + NoAxes() + ggtitle("") 
      } else {
        DimPlot(seurat, group.by = "seurat_clusters") + NoAxes() + ggtitle("")
      }
    })
    
    nr_cols <- reactive({
      isolate(colors_use())
    })
    
    output$color_list <- renderPrint({
      cat("palette <- ","c(", paste(shQuote(unlist(nr_cols())), collapse = ", "), ")\n")
    })
  })
  
  observe({
    seurat <- seurat_object()
    grouping_column <- input$grouping_column
    if (!is.null(seurat) && !is.null(grouping_column)) {
      unique_values <- unique(seurat@meta.data[[grouping_column]])
      updateSelectInput(session, "group_value", choices = unique_values)
    }
  })
}

# Run the Shiny app
shinyApp(ui, server)
}
