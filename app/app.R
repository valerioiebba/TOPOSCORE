library(DT)
library(ggplot2)
library(dplyr)
# Add at the beginning with other library imports
library(readr)
library(ggnewscale)
# Add to library imports at the beginning
library(ggrepel)
# Add after get_data_path function

# Source the calculation function
source("R/toposcore_calc.R")
source("R/kde_plot.R")
ui <- fluidPage(
  titlePanel("Toposcore Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      # File upload and calculate button group
      fluidRow(
        column(8,
               fileInput("file", "Upload microbiome tab-separated data",
                         accept = c("text/tsv", "text/tab-separated-values,text/plain", ".txt/.tsv"))
        ),
        selectInput("Level", "Select taxonomy level and database", choices = c("Species - MetaPhlAn 4 vJan21", "SGB - MetaPhlAn 4 vJan21", "Species - GTDB r207", "Species - GTDB r220")),
        column(4,
               br(), # Add some spacing to align with file input
               actionButton("calculate", "Calculate Toposcore", 
                            class = "btn-primary")
        )
      ),
      
      # Input data template section
      h4("Input Data Requirements:"),
      HTML("<p>The input table should be a multi-level taxa (rows) x samples (columns) data frame.</p>
      <p><strong>Required column:</strong></p>
           <ul>
             <li>clade_name: Full level taxonomy for each taxon at different taxonomic levels</li>
           </ul>
      <p><strong>For MetaPhlAn 4 options, the required input is the output of the <a href='https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4#merging-tables'>merge_metaphlan_tables.py</a> command</strong></p>"),
      downloadButton("download_template", "Download MetaPhlAn 4 Template File"),
      HTML("<p><strong>For GTDB options, the required input can be output of running the <a href='https://github.com/bluenote-1577/sylph/wiki/Taxonomic-profiling-with-the-GTDB‐R214-database'>Sylph wiki</a> and merging profiles </strong></p>"),
      downloadButton("download_template_gtdb", "Download Sylph GTDB Template File"),
      
      # Display summary stats
      hr(),
      h4("Summary Statistics:"),
      verbatimTextOutput("summary"),
      
      hr(),
      # Display species lists
      h4("SIG1 Species:"),
      verbatimTextOutput("sig1_species"),
      
      h4("SIG2 Species:"),
      verbatimTextOutput("sig2_species"),
      
      # Download buttons
      downloadButton("download", "Download Results"),
      
    ),
    
    # In the UI section, modify the plotOutput height
    mainPanel(
      tabsetPanel(
        tabPanel("Results", 
                 DTOutput("results_table"),
                 downloadButton("download_results", "Download Results Table"),
                 hr(),
                 # Sample selection controls
                 fluidRow(
                   column(4,
                          selectizeInput("sample_select", "Select Sample(s):",
                                         choices = NULL,
                                         multiple = TRUE,
                                         options = list(
                                           placeholder = 'Click or type to select samples'
                                         ))
                   ),
                   column(4,
                          checkboxInput("select_all", 
                                        "Select All Samples", 
                                        value = FALSE)
                   )
                 ),
                 # Score distribution plot
                 plotOutput("score_dist"),
                 downloadButton("download_plot", "Download Plot"),
                 hr(),
                 plotOutput("species_heatmap", height = "800px"))  # Increased from 600px to 800px
      )
    )
  )
)

# Add helper function after libraries
get_data_path <- function(filename) {
  file.path("data", filename)
}

server <- function(input, output, session) {
  
  # Update color_values function with OS12 mapping
  color_values <- function(type, data) {
    if (type == "SIG_class") {
      return(c("SIG1" = "red", "Gray" = "gray", "SIG2" = "darkgreen"))
    } else if (type == "OS12") {
      return(c("R" = "green", "NR" = "red"))  # R for Responders, NR for Non-Responders
    }
    # Add other color mappings if needed
  }
  # Handle "Select All" checkbox
  observeEvent(input$select_all, {
    req(results())
    if (input$select_all) {
      updateSelectizeInput(session, "sample_select",
                           choices = results()$Sample_id,
                           selected = results()$Sample_id)
    } else {
      updateSelectizeInput(session, "sample_select",
                           choices = results()$Sample_id,
                           selected = NULL)
    }
  })
  
  # Update sample selection when results change
  observe({
    req(results())
    updateSelectizeInput(session, "sample_select",
                         choices = results()$Sample_id,
                         selected = if(input$select_all) results()$Sample_id else NULL)
  })
  
  # Add download handler for results table
  output$download_results <- downloadHandler(
    filename = function() {
      paste("toposcore_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(results(), file, row.names = FALSE)
    }
  )
  # Modified heatmap output
  output$species_heatmap <- renderPlot({
    req(results(), input$sample_select)
    selected_samples <- input$sample_select
    
    tryCatch({
      # Read original data and signature species
      original_data <- read_csv(input$file$datapath)
      sig1_species <- read.table(get_data_path("sig1.txt"), stringsAsFactors = FALSE, sep = "\t", skip = 1)[,1]
      sig2_species <- read.table(get_data_path("sig2.txt"), stringsAsFactors = FALSE, sep = "\t", skip = 1)[,1]
      
      # Add Akkermansia as a separate group
      akk_species <- "Akkermansia_muciniphila"
      
      # Create data frame for heatmap with three groups
      heatmap_data <- data.frame(
        Species = rep(c(sig1_species, akk_species, sig2_species), length(selected_samples)),
        Group = rep(c(rep("SIG1", length(sig1_species)), 
                      rep("AKK", 1),
                      rep("SIG2", length(sig2_species))), length(selected_samples)),
        Sample = rep(selected_samples, each = length(c(sig1_species, akk_species, sig2_species))),
        Abundance = 0  # Default value
      )
      
      # Set factor levels to control order
      heatmap_data$Group <- factor(heatmap_data$Group, levels = c("SIG1", "AKK", "SIG2"))
      
      # Update abundance values where available
      for(sample_id in selected_samples) {
        sample_data <- original_data[original_data$Sample_id == sample_id, ]
        for(species in c(sig1_species, akk_species, sig2_species)) {
          if(species %in% colnames(sample_data)) {
            heatmap_data$Abundance[heatmap_data$Species == species & 
                                     heatmap_data$Sample == sample_id] <- 
              as.numeric(sample_data[[species]])
          }
        }
      }
      
      # Calculate global maximum abundance across all species
      max_global <- max(heatmap_data$Abundance)
      
      # Create enhanced heatmap with different color scales for AKK
      ggplot(heatmap_data, aes(x = Sample, y = Species)) +
        # Add tiles with different color scales based on group
        geom_tile(data = subset(heatmap_data, Group != "AKK"),
                  aes(fill = Abundance)) +
        # Different color scales for SIG species
        scale_fill_gradient2(low = "white", high = "red", 
                             mid = "pink", midpoint = max_global/2,
                             limits = c(0, max_global),
                             name = "Relative\nAbundance") +
        new_scale_fill() +
        geom_tile(data = subset(heatmap_data, Group == "AKK"),
                  aes(fill = Abundance)) +
        scale_fill_gradient(low = "white", high = "blue",
                            limits = c(0, max_global),
                            name = "Akkermansia\nAbundance") +
        facet_grid(Group ~ ., scales = "free", space = "free") +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          strip.text = element_text(size = 12, face = "bold"),
          strip.background = element_rect(
            fill = c("lightblue", "lightgray", "lightgreen"),
            colour = "black"
          )
        ) +
        labs(title = "Species Abundance Heatmap",
             subtitle = paste("Selected Sample(s):", 
                              paste(selected_samples, collapse = ", ")))
      
    }, error = function(e) {
      # Return an error plot
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = paste("Error generating heatmap:", e$message)) +
        theme_void()
    })
  })
  
  # Add template download handler
  output$download_template <- downloadHandler(
    filename = function() {
      "template_mpa4_vJan21_test_data.tsv"
    },
    content = function(file) {
      file.copy("data/mpa4_vJan21_test_data.tsv", file)
    }
  )
  
  output$download_template_gtdb <- downloadHandler(
    filename = function() {
      "template_sylph_GTDBr220_test_data.tsv"
    },
    content = function(file) {
      file.copy("data/sylph_GTDBr220_test_data.tsv", file)
    }
  )
  
  # Add plot download handler
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("toposcore_distribution_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = {
        # Create base reference KDE plot - exactly the same as UI plot
        base_plot <- plot_toposcoreb01_density(scores_disc, clin_disc, lims = c(0.5351, 0.7911)) +
          # Expand the plot limits to make room for labels
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)), limits = c(-0.2, NA))
        
        # Add vertical lines for selected samples if any
        if (!is.null(input$sample_select) && length(input$sample_select) > 0) {
          selected_data <- results()[results()$Sample_id %in% input$sample_select, ]
          
          # Add vertical lines and labels
          base_plot <- base_plot +
            # Add vertical lines
            geom_vline(data = selected_data,
                       aes(xintercept = S_score),
                       color = "red",
                       linetype = "dashed",
                       linewidth = 1) +
            # Add text labels
            annotate("text",
                     x = selected_data$S_score,
                     y = rep(3.1, nrow(selected_data)),
                     label = selected_data$Sample_id,
                     color = "red",
                     angle = 45,
                     hjust = 1,
                     size = 4)
        }
        
        base_plot
        
      }, width = 10, height = 6)
    }
  )
  
  # Reactive value to store results
  results <- reactiveVal(NULL)
  
  # Calculate results when button is clicked
  observeEvent(input$calculate, {
    req(input$file)
    
    # Add error handling
    tryCatch({
      convert_options <- c("species_Jan21" = "Species - MetaPhlAn 4 vJan21" , 
                           "SGB_Jan21" = "SGB - MetaPhlAn 4 vJan21", 
                           "GTDB_r207"= "Species - GTDB r207", 
                           "GTDB_r220" = "Species - GTDB r220")
      res <- calculate_toposcore(input$file$datapath, names(convert_options)[convert_options == input$Level])
      results(res)
      showNotification("Calculation completed successfully!", type = "message")
    }, error = function(e) {
      showNotification(
        paste("Error in calculation:", e$message), 
        type = "error",
        duration = NULL
      )
    })
  })
  
  # Display species lists
  output$sig1_species <- renderPrint({
    sig1_species <- read.table("data/sig1.txt", stringsAsFactors = FALSE, sep = "\t", skip = 1)[,1]
    cat(paste(sig1_species, collapse = "\n"))
  })
  
  output$sig2_species <- renderPrint({
    sig2_species <- read.table("data/sig2.txt", stringsAsFactors = FALSE, sep = "\t", skip = 1)[,1]
    cat(paste(sig2_species, collapse = "\n"))
  })
  
  # Make the results table interactive with row selection
  output$results_table <- renderDT({
    req(results())
    tryCatch({
      datatable(results(), 
                options = list(pageLength = 10),
                selection = 'single',  # Allow single row selection
                rownames = FALSE)
    }, error = function(e) {
      showNotification(
        paste("Error displaying results:", e$message),
        type = "error"
      )
      return(NULL)
    })
  })
  
  
  
  # Display summary statistics with error handling
  output$summary <- renderPrint({
    req(results())
    tryCatch({
      # First line: S-Score classifications
      cat("S-Score Classifications:\n")
      sig_counts <- table(results()$SIG_class)
      cat(sprintf("SIG1(%d) | Gray(%d) | SIG2(%d)\n\n",
                  sig_counts["SIG1"], sig_counts["Gray"], sig_counts["SIG2"]))
      
      # Second line: Toposcore classifications
      cat("Toposcore Classifications:\n")
      topo_counts <- table(results()$Toposcore)
      cat(sprintf("SIG1+(%d) | SIG2+(%d)",
                  topo_counts["SIG1+"], topo_counts["SIG2+"]))
    }, error = function(e) {
      cat("Error generating summary statistics:", e$message)
    })
  })
  
  # Load reference data
  scores_disc <- read_csv("data/scores_disc.csv")
  clin_disc <- read_csv("data/clin_disc.csv")
  
  # Modify score_dist to use only sample_select
  output$score_dist <- renderPlot({
    # Create base reference KDE plot
    base_plot <- plot_toposcoreb01_density(scores_disc, clin_disc, lims = c(0.5351, 0.7911)) +
      # Expand the plot limits to make room for labels
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)), limits = c(-0.5, NA))
    
    # Add vertical lines for selected samples if any
    if (!is.null(input$sample_select) && length(input$sample_select) > 0) {
      selected_data <- results()[results()$Sample_id %in% input$sample_select, ]
      
      # Add vertical lines and labels
      base_plot <- base_plot +
        # Add vertical lines
        geom_vline(data = selected_data,
                   aes(xintercept = S_score),
                   color = "red",
                   linetype = "dashed",
                   linewidth = 1) +
        # Add text labels
        annotate("text",
                 x = selected_data$S_score,
                 y = rep(3.1, nrow(selected_data)),
                 label = selected_data$Sample_id,
                 color = "red",
                 angle = 45,
                 hjust = 1,
                 size = 4)
    }
    
    base_plot
  }, height = 400)  # Set a fixed height for the plot
  
  
  # Download handler with error handling
  output$download <- downloadHandler(
    filename = function() {
      paste("toposcore_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      tryCatch({
        write.csv(results(), file, row.names = FALSE)
      }, error = function(e) {
        showNotification(
          paste("Error downloading results:", e$message),
          type = "error"
        )
      })
    }
  )
}


# Run the app
shinyApp(ui = ui, server = server)