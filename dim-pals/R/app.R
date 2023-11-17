library(shiny)
library(dplyr)
library(Seurat)
library(colourpicker)  # For color picker input
library(paletteer) # Package containing color palettes
library(ggplot2)

dim.pals <- function(seurat_object = NULL) {
  ui <- fluidPage(
    titlePanel("dim-pals"),
    sidebarLayout(
      sidebarPanel(
        tags$h3("re-color individual clusters"),
        selectInput("cluster_to_color","Select Cluster for Custom Color", choices = NULL),
        colourInput("col", "Select Cluster Color", "purple"),
        actionButton("assign_cluster_color_button", "Assign Cluster Color"),
        tags$h3("Apply Palette to Groups of Clusters"),
        tags$h4("Select Clusters"),
        selectInput("grouping_column", "Select Grouping Column", choices = NULL),
        selectInput("group_value", "Select Group", choices = NULL),
        tags$h4("Palette"),
        selectInput("pal_types", "Select Palette Type", choices = c("Discrete","Continuous","Dynamic")),
        selectInput("pal_package", "Select Palette Package", choices = NULL),
        selectInput("pal_name", "Select Palette Name", choices = NULL),
        actionButton("assign_group_color_button", "Assign Palette"),
        tags$h3("Refresh Plot"),
        actionButton("plot_button","Generate Plot"),
        checkboxInput("show_labels", "Show Cluster Labels", value = FALSE)
      ),
      mainPanel(
        plotOutput("umap_plot"),
        verbatimTextOutput("color_list")
      )
    )
  )
  
  server <- function(input, output, session){
    ### INITIALIZE PALETTE
    initialized_palette <- rep("grey", length(unique(seurat_object$seurat_clusters)))
    names(initialized_palette) <- c(0:(max(as.numeric(data$seurat_clusters)) -1))
    
    ##DEFINE REACTIVE VARIABLES
    seurat_object <- reactiveVal(seurat_object)
    palette <- reactiveVal(initialized_palette)
    plotted <- reactiveVal(NULL)
    
    ##UPDATE SELECTION BOX INPUT OPTIONS
    observe({
      seurat <- seurat_object()
      if (!is.null(seurat)) {
        seurat@meta.data <- seurat@meta.data %>% 
          mutate(all_clusters = "Select All") %>% 
          select(all_clusters, everything()) %>% 
          arrange(seurat_clusters)
        updateSelectInput(session, "cluster_color", choices = unique(seurat$seurat_clusters))
        updateSelectInput(session, "grouping_column", choices = colnames(seurat@meta.data))
        updateSelectInput(session, "cluster_to_color", choices = unique(seurat$seurat_clusters))
        seurat_object(seurat)
      }
    })
    
    ##UPDATE THE GROUP SELECTION MENU
    observe({
      seurat <- seurat_object()
      grouping_column <- input$grouping_column
      if (!is.null(seurat) && !is.null(grouping_column)) {
        unique_values <- unique(seurat@meta.data[[grouping_column]])
        updateSelectInput(session, "group_value", choices = unique_values)
        seurat_object(seurat)
      }
    })
    
    # Observe event for updating palette choices based on the selected palette package
    observe({
      
      if (input$pal_types == "Discrete"){
           pkg_choices <- palettes_d_names %>% pull(package)
           updateSelectInput(session, "pal_package", choices = pkg_choices)
      }
      
      if (input$pal_types == "Continuous"){
        pkg_choices <- palettes_c_names %>% pull(package)
        updateSelectInput(session, "pal_package", choices = pkg_choices)
      }
      
      if (input$pal_types == "Dynamic"){
        pkg_choices <- palettes_dynamic_names %>% pull(package)
        updateSelectInput(session, "pal_package", choices = pkg_choices)
      }
    })
    
    ##update pal type selection
    observe({
      pkg_select <- input$pal_package
      
      if (input$pal_types == "Discrete" && !is.null(pkg_select)){
        pal_choices <- palettes_d_names %>% filter(package == pkg_select) %>% pull(palette)
        updateSelectInput(session, "pal_name", choices = pal_choices)
      }
      
      if (input$pal_types == "Continuous" && !is.null(pkg_select)){
        pal_choices <- palettes_c_names %>% filter(package == pkg_select) %>% pull(palette)
        updateSelectInput(session, "pal_name", choices = pal_choices)
      }
      
      if (input$pal_types == "Dynamic" && !is.null(pkg_select)){
        pal_choices <- palettes_dynamic_names %>% filter(package == pkg_select) %>% pull(palette)
        updateSelectInput(session, "pal_name", choices = pal_choices)
      }
    })
      
    
    ##SELECT CLUSTER RECOLOR BUTTON
    observeEvent(input$assign_cluster_color_button, {
      cluster <- input$cluster_to_color
      color <- input$col
      newpal <- palette()
      newpal[[paste(cluster)]] <- color
      palette(newpal)
    })
    
    ##GENERATE GROUP PALETTE BUTTON
    observeEvent(input$assign_group_color_button, {
      seurat <- seurat_object()
      grouping_col  <- input$grouping_column
      group_val  <- input$group_value
      package <- input$pal_package
      palette_name <- input$pal_name
      newpal <- isolate(palette())
      
      if(!is.null(seurat) && !is.null(grouping_col) && !is.null(group_val) && !is.null(package)){
      
      group_clusters <- seurat@meta.data %>% 
        select(seurat_clusters, grouping_col) %>% 
        unique() %>% 
        filter(.data[[grouping_col]] == group_val) %>% 
        pull(seurat_clusters)

      
      pal_name <- paste0(package,"::",palette_name)
      
      tryCatch({
        if (input$pal_types == "Discrete" && !is.null(package)) {
          pal_vect <- paletteer::paletteer_d(pal_name, length(group_clusters), type = "continuous")
        }
        
        if (input$pal_types == "Continuous" && !is.null(package)) {
          pal_vect <- paletteer::paletteer_c(pal_name, length(group_clusters))
        }
        
        if (input$pal_types == "Dynamic" && !is.null(package)) {
          pal_vect <- paletteer::paletteer_dynamic(pal_name, length(group_clusters))
        }
        
        names(pal_vect) <- as.character(group_clusters)
        match_indices <- match(names(pal_vect), names(newpal))
        newpal[match_indices[!is.na(match_indices)]] <- pal_vect[!is.na(match_indices)]
        palette(newpal)
        
        grouping_col <- NULL
        group_val <- NULL
        package <- NULL
        palette_name <- NULL
        newpal <- NULL
      }, error = function(e) {
        # Display a warning without crashing the app
        showNotification("Selected palette has fewer colors than needed. Please select a different palette.", type = "warning")
      })
    }
    })
    
    ##DEFAULT UMAP PLOT GENERATION ON APP LAUNCH
    observe({
      pal <- palette()
      seurat <- seurat_object()
      if (is.null(plotted())){
        output$umap_plot <- renderPlot(
        DimPlot(seurat, group.by = "seurat_clusters", cols = pal) + NoAxes() + ggtitle("")
        )
      }
    })
    
    ##BUTTON: GENERATE UMAP
    observeEvent(input$plot_button,{
      plotted(TRUE)
      seurat <- seurat_object()
      pal <- palette()
      output$umap_plot <- renderPlot(
        if (input$show_labels) {
          DimPlot(seurat, group.by = "seurat_clusters", cols = pal, label = TRUE, label.size = 7, label.box = TRUE) + NoAxes() + ggtitle("")
        } else {
          DimPlot(seurat, group.by = "seurat_clusters", cols = pal) + NoAxes() + ggtitle("")
        }
      )
    })
    
    ## VECTOR TEXT OUTPUT
    output$color_list <- renderPrint({
      cat("Copy and paste this palette for use in the seurat 'cols' argument", "\n","\n")
      cat("palette <- ","c(", paste(shQuote(unlist(palette())), collapse = ", "), ")\n")
    })
  }
  
  shinyApp(ui,server)
}

