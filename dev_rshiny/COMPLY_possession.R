
#######################

#    C O M P L Y v2.0

#   Version:02/04/26
#   macdonald.megan@epa.gov

#######################





# This is a Shiny web application. You can run the application by clicking the 'Run App' button above.



setwd("/home/kbaker/shiny/comply")


library(dplyr)
library(shiny)
library(leaflet)
library(tidygeocoder)
library(DT)
library(reshape2)
library(ggplot2)
library(tidyr)
library(plotly)
library(bslib)
library(tibble)

library(leaflet.extras)

library(sf)




rads <- read.csv(file="comply171_original_possession_data.csv",skip=0,header=T)
nrow(rads)
long_rads <- melt(rads, 
                  id.vars = "Nuclide", 
                  variable.name = "Form", 
                  value.name = "Allowance")
long_rads$Radionuclide <- gsub(" ", "", long_rads$Nuclide)
long_rads$Type <- gsub(" ", "", long_rads$Form)
long_rads <- long_rads[,c(4,5,3)]






# Prepare blank dataframes


rads <- data.frame(
  `Radionuclide` = as.character(NA),
  `Type` =  as.character(NA),
  `Emission` = round(as.numeric(NA),1),
  `Allowance` = round(as.numeric(NA),1),
  `Fraction` = round(as.numeric(NA),1)
)

xsum <- data.frame(
  `Fraction` = round(as.numeric(NA),1)
)




#### APP CODE STARTS HERE:


#### UI

ui <- page_sidebar(
  title = "COMPLY v2.0 Possession Tables",
  # Define the sidebar panel
  sidebar = sidebar(
    h6("Enter Facility Details:"),
    textInput("CompanyName", "Company Name:", ""),
    textInput("CompanyContact", "Company Contact Name:", ""),
    textInput("CompanyContactNum", "Company Contact Number:", ""),
    textInput("FacilityName", "Facility Name:", ""),
    h6("Enter Coordinates or Address:"),
    fluidRow(
      column(6,
             textInput("Latitude", "Latitude:", "")
      ),
      column(6,
             textInput("Longitude", "Longitude:", "")
      )
    ),
    textInput("address_input", "Facility Address:", "109 T.W. Alexander Drive, Durham, NC, 27709"),
    actionButton("geocode_button", "Plot Location"),
    #actionButton("add_heatmap_btn", "HEATMAP"),
    h6("Coordinates:"),
    textOutput("coords")
  ),
  
  # Main panel content 
  # Main row: Tabs within a card using navset_card_tab
  card(
    full_screen = TRUE,# Optional: allows the card to go full screen
    height = "1000px", # Sets a fixed height
    navset_card_tab(
      
      nav_panel("Rads",
                actionButton("add_rad_button", "Add Radionuclide", class = "btn-sm", width = '20%'),
                tags$head(
                  tags$style(HTML("
      /* Add a 1px solid black border to all table cells and headers */
      .dataTable td, .dataTable th {
        border: 1px solid black;
      }
        /* Change background and text color of all DT buttons */
      .dt-button {
        background-color: #f0f0f0 !important; /* White background */
        color: black !important;             /* Black text */
        border: 1px solid #000500 !important; /* Black border */
      }
    "))
                ),
                DTOutput("rad_table")  #(dynamic table)
      ),
      
      nav_panel("Compliance Summary",
                #DTOutput("sums"),
                textOutput("fraction_display"),
                textOutput("fraction_display_i")
      ),
      
      nav_panel("Export Data",
                downloadButton("downloadCSVData", "Download .CSV file", class = "btn-sm", width = '20%'),
                downloadButton("downloadReport", "Generate HTML Report", class = "btn-sm", width = '20%')
      )
    )
  ), #end condition for main row
  
  # Ancillary information: Two cards using layout_columns for a grid
  layout_columns(
    # Card 1: Value box or plot
    card(
      full_screen = TRUE, # Optional: allows the card to go full screen
      height = "1000px", # Sets a fixed height
      card_header("Map"),
      leafletOutput("map")
    ),
    # Card 2: Value box or plot
    
    card(
      full_screen = TRUE, # Optional: allows the card to go full screen
      height = "1000px", # Sets a fixed height
      card_header("Plot"),
      plotlyOutput(outputId = "graph")
      #      full_screen = TRUE,# Optional: allows the card to go full screen
#      height = "1000px", # Sets a fixed height
#     card_header("Compliance Summary"),
#      "This is the text box content inside the new Shiny card."
    ),
    
    # Use col_widths to specify that both cards should be side-by-side in one row
    col_widths = c(6, 6) # Divides the row into two equal columns (6+6=12)
  )
  
)

server <- function(input, output, session) {
  


  
  
  # RADS ------------------------------------------------------------------
  
  # Reactive value to store the data frame
  rad_data <- reactiveVal(rads)
  rad_sum <- reactiveVal(xsum)
  i_sum <- reactiveVal(xsum)
  
  # Render the data table
  output$rad_table <- renderDT({
   tab <-  datatable(rad_data(), 
                     class = 'cell-border stripe', # Add the 'cell-borders' class
                     extensions = 'Buttons',
                     options = list( 
                       ordering = FALSE,  # Disables sorting functionality on all columns
                       searching = FALSE, # Removes the global search box
                       info = FALSE,       # Removes the "Showing X to Y of Z entries" display
                       lengthChange = FALSE, # This removes the 'Show entries' dropdown
                       dom = 'Blfti', # Excludes 'p' for pagination
                       columnDefs = list(list(className = 'dt-center', targets = "_all")),
                       rownames = FALSE,
                       pageLength = -1,      # Shows all rows (-1)
                       paging = FALSE,
                       filter = "top",
                       buttons = c('copy', 'csv', 'excel', 'pdf'), # Specify the export buttons
                       # Optional: customize export options, e.g.
                       exportOptions = list(
                         columns = ':visible' # Only export visible columns
                       )))
   tab
  })

  
  # Observe the "Add Row" button click to show the modal
  observeEvent(input$add_rad_button, {
    showModal(
      modalDialog(
        title = "Enter New Radionuclide Data",
        #textInput("input_ID", "Stack ID:", ""),
        selectInput(
          inputId = "input_RN",
          label = "Radionuclide:",
          choices = unique(long_rads$Radionuclide),
          selected = NULL
        ),
        selectInput(
          inputId = "input_TYPE",
          label = "Type:",
          choices = unique(long_rads$Type),
          selected = NULL
        ),
        numericInput("input_EMISSION", "Emission:", value = NULL, min = 0),
        footer = tagList(
          actionButton("save_rad_btn", "Save"),
          modalButton("Cancel")
        )
      )
    )
  })
  
  
  
  # Observe the "Save" button click within the modal
  observeEvent(input$save_rad_btn, {
    # Get the current data
    current_data <- rad_data()
  
    
    # Create the new row
    new_row <- data.frame(
      `Radionuclide` = input$input_RN,
      `Type` = input$input_TYPE,
      `Emission` = input$input_EMISSION,
      `Allowance` = subset(long_rads, Type == input$input_TYPE & Radionuclide == input$input_RN)[1,3],
      `Fraction` = round((input$input_EMISSION/(subset(long_rads, Type == input$input_TYPE & Radionuclide == input$input_RN)[1,3]) * 1),5),

      stringsAsFactors = FALSE
    )
    
    # Append the new row using bind_rows for safety with types
    updated_data <- bind_rows(current_data, new_row)
    cleaned_data <- updated_data[rowSums(is.na(updated_data)) < ncol(updated_data), ]
    print(cleaned_data)
    temp <- cleaned_data
    temp$element <- substr(temp$Radionuclide,1,1)
    iodine_data <- subset(temp,element=="I")
    print(iodine_data)
    
    new_row2 <- data.frame(
      `Fraction` = sum(cleaned_data$Fraction),
      stringsAsFactors = FALSE
    )    
    
    new_row3 <- data.frame(
      `Fraction` = sum(iodine_data$Fraction),
      stringsAsFactors = FALSE
    )    
  
    # Update the reactive value, which automatically re-renders the table
    rad_data(cleaned_data)
    rad_sum(new_row2)
    i_sum(new_row3)
    # Update the reactive value, which automatically re-renders the table
    print(rad_data())
    
    # Close the modal
    removeModal()
  })
  
  

  
  
# BUILD COMPLETE DATA FRAME   ----------------------------------------------

    output$sums <- renderDT({
    datatable(rad_sum(), options = list(dom = 't')) # 'dom = t' removes search/pagination
  })
    
    combined_data <- reactive({
      req(rad_data())
      # Use bind_rows to stack dataframes vertically (assumes same column names)
      combined_df <- rad_data

      # Remove rows with any NA values
      cleaned <- as.data.frame(combined_df[rowSums(is.na(combined_df[, -1])) < ncol(combined_df[, -1]), ])
      print(cleaned)
      
      rad_distances_doses <- sum(cleaned$Fraction)
      print("sum")
      print(rad_distances_doses)
      return(rad_distances_doses)
    })
 

    output$fraction_display <- renderText({
      selected_row <- rad_sum()
      selected_age <- selected_row$Fraction
      # Create the text string using paste0() or paste()
      paste0("Total possession fraction is ", selected_age, ".")
    })        
    
    output$fraction_display_i <- renderText({
      selected_row <- i_sum()
      selected_age <- selected_row$Fraction
      # Create the text string using paste0() or paste()
      paste0("Total iodine possession fraction is ", selected_age, ".")
    })        
    
        
    # Render the interactive plot using renderPlotly and ggplotly
    output$graph <- renderPlotly({
      # Use req() to ensure data is available before plotting
      req(rad_data())
      # Call the reactive expression like a function: combined_data()
      p <- ggplot(rad_data(), aes(x = 1, y = Fraction, fill = Radionuclide)) +
        geom_col(alpha = 0.6) +
        labs(title = "Possession Fraction Plot", color = "Source Data") +
        theme_classic()+ 
        xlab(" ")+
        ylab("Fraction")+
        theme(
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white", color = NA)
        ) +
        labs(color = "Radionuclides")
      
      ggplotly(p) # Transform the ggplot object into an interactive plotly object
    })
  
    

# Data Downloads ----------------------------------------------------------
    
    # Handler for the CSV download button
    output$downloadCSVData <- downloadHandler(
      filename = function() {
        # Dynamic filename generation
        paste0("COMPLY_", Sys.Date(), ".csv") 
      },
      content = function(file) {
        # Write the content to a temporary file Shiny provides
        write.csv(rad_data, file, row.names = FALSE) 
      }
    )    
  
    ##### HTML Report
    # Reactive data for both app display and report
    report_data <-  reactive({rad_data()})
    report_plot <- reactive({
     p <-  ggplot(rad_data(), aes(x = 1, y = Fraction, fill = Radionuclide)) +
        geom_col(alpha = 0.6) +
        labs(title = "Possession Fraction Plot", color = "Source Data") +
        theme_classic()+
        xlab(" ")+
        ylab("Fraction")+
        theme(
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white", color = NA)
        ) +
        labs(color = "Radionuclides")

      ggplotly(p) # Transform the ggplot object into an interactive plotly object
    })
#     
#     

    output$HTML_table <- renderDT({ report_data() })
    output$HTML_plot <- renderPlotly({ report_plot() })


    output$downloadReport <- downloadHandler(
      filename = function() { paste0("report_", Sys.Date(), ".html") },
      content = function(file) {
        # Use a temporary directory for the .Rmd to ensure write permissions
        temp_rmd <- file.path(tempdir(), "COMPLY_Output.Rmd")
        file.copy("COMPLY_Output.Rmd", temp_rmd, overwrite = TRUE)

        # Pass reactive objects to the R Markdown template
        rmarkdown::render(temp_rmd, output_file = file,
                          params = list(report_data = report_data(),report_plot = report_plot()),
                          envir = new.env(parent = globalenv())
        )
      })


    
    
    

# MAPPING -----------------------------------------------------------------

    # Reactive expression to perform geocoding
    geocoded_address <- eventReactive(input$geocode_button, {
      if (input$Latitude == "") {
        
        full_address <- input$address_input
        some_addresses <- tibble::tribble(
          ~name,                  ~addr,
          "company",          full_address
        )
        print(some_addresses)
        
        # Use tidygeocoder to geocode the address string
        geocode_result <- some_addresses %>%
          geocode(addr, method = 'osm', lat = latitude , long = longitude)
        print(geocode_result)
        
        # Check if geocoding was successful
        if (nrow(geocode_result) > 0 && !is.na(geocode_result$latitude) && !is.na(geocode_result$longitude)) {
          geocode_result
        } else {
          showModal(modalDialog(
            title = "Geocoding Error",
            "Could not find coordinates for the entered address. Please check the address and try again.",
            easyClose = TRUE
          ))
          NULL
        }
      } else {
        latitude <- input$Latitude
        longitude <- input$Longitude
        name <- "company"
        addr <- ""
        geocode_result <- as.tibble(cbind(name, addr, as.numeric(latitude), as.numeric(longitude)))
        names(geocode_result) <- c("name","addr","latitude","longitude")
        print(geocode_result)
      }
    })
    
    
    # Render the Leaflet map
    output$map <- renderLeaflet({
      # Initialize an empty Leaflet map
      leaflet() %>%
        addProviderTiles(providers$Esri.WorldTopoMap) %>%
        setView(lng = -98.5833, lat = 39.8333, zoom = 5) # Default view (center of the US)
    })
    
    # Update the map when a new address is geocoded
    observe({
      result <- geocoded_address()
      result$latitude <- as.numeric(result$latitude)
      result$longitude <- as.numeric(result$longitude)
      print(result)
      if (!is.null(result)) {
        # Clear existing markers
        leafletProxy("map") %>% clearMarkers()
        
        # Add a marker at the new coordinates and center the map
        leafletProxy("map") %>%
          addMarkers(
            lng = result$longitude,
            lat = result$latitude,
            popup = result$addr
          ) %>%
          setView(
            lng = result$longitude,
            lat = result$latitude,
            zoom = 10
          )
      }
    })
    
    # Render the coordinates text
    output$coords <- renderText({
      result <- geocoded_address()
      result$latitude <- as.numeric(result$latitude)
      result$longitude <- as.numeric(result$longitude)
      if (!is.null(result)) {
        paste("Latitude:", round(result$latitude, 6), ", Longitude:", round(result$longitude, 6))
      } else {
        "Enter an address and click 'Plot Location'."
      }
    })    



}


# Define the viewer with specific width and height (in pixels)
vwr <- dialogViewer('COMPLY v2.0', width = 3200, height = 7400)
# 
# # Run the app using runGadget and the defined viewer
#
#runGadget(shinyApp(ui = ui, server = server), viewer = vwr)
shinyApp(ui,server)