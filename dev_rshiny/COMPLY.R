
#######################

#    C O M P L Y

#   Version:02/04/26
#   macdonald.megan@epa.gov

#######################





# This is a Shiny web application. You can run the application by clicking the 'Run App' button above.



setwd("C:/Users/MMACDO01/OneDrive - Environmental Protection Agency (EPA)/MM Projects & Code/COMPLY Shiny App/COMPLY")


library(dplyr)
library(shiny)
library(leaflet)
library(tidygeocoder)
library(DT)
library(reshape2)
library(ggplot2)
library(tidyr)
library(plotly)

library(leaflet.extras)

library(sf)

library(dplyr)
library(shiny)
library(leaflet)
library(tidygeocoder)
library(DT)
library(reshape2)
library(ggplot2)
library(tidyr)
library(plotly)

library(leaflet.extras)



# define these lookup tables for use in the app ... need to make tehse available with app code in folder
doses_info <- read.csv(file="C:/Users/MMACDO01/OneDrive - Environmental Protection Agency (EPA)/MM Projects & Code/COMPLY Shiny App/COMPLY/comply171_original_input_data.csv",skip=0,header=T)
doses_info$Radionuclide <- trimws(doses_info$Nuclide)


rads <- read.csv(file="comply171_original_possession_tables.csv",skip=0,header=T)
nrow(rads)
long_rads <- melt(rads, 
                  id.vars = "Nuclide", 
                  variable.name = "Form", 
                  value.name = "Allowance")
long_rads$Radionuclide <- gsub(" ", "", long_rads$Nuclide)
long_rads$Type <- gsub(" ", "", long_rads$Form)
long_rads <- long_rads[,c(4,5,3)]








# Prepare blank dataframes

stacks <- data.frame(
  `Stack ID` = as.character(NA),
  `Release Height` = round(as.numeric(NA),1),
  `Building Height` = round(as.numeric(NA),1),
  `Building Width` = round(as.numeric(NA),1),
  `Source Distance`= round(as.numeric(NA),1),
  `Veg Distance`= round(as.numeric(NA),1),
  `Cow Distance`= round(as.numeric(NA),1),
  `Mean WS`= round(as.numeric(NA),1)
)

rads <- data.frame(
  `Stack ID`= as.character(NA),
  `Radionuclide` = as.character(NA),
  `Type` =  as.character(NA),
  `Emission` = round(as.numeric(NA),1),
  `Allowance` = round(as.numeric(NA),1),
  `Percent` = round(as.numeric(NA),1)
  
)

doses <- data.frame(
  `Stack.ID` = as.character(NA),
  `Radionuclide` = as.character(NA),
  `Type` =  as.character(NA),
  `Source.Distance` =  round(as.numeric(NA),1),
  `dose`= round(as.numeric(NA),1),
  `vdose` = round(as.numeric(NA),1),
  `cdose` =  round(as.numeric(NA),1),
  `total` = round(as.numeric(NA),1)

)




#### APP CODE STARTS HERE:


#### UI



ui <- page_sidebar(
  title = "COMPLY",
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
    actionButton("geocode_button", "Geocode Address"),
    actionButton("add_heatmap_btn", "HEATMAP"),
    h6("Coordinates:"),
    textOutput("coords")
  ),
  
  # Main panel content 
  
  # First row: Two cards using layout_columns for a grid
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
      
    ),
    # Use col_widths to specify that both cards should be side-by-side in one row
    col_widths = c(6, 6) # Divides the row into two equal columns (6+6=12)
  ),
  
  # Second row: Tabs within a card using navset_card_tab
  card(
    full_screen = TRUE,# Optional: allows the card to go full screen
    height = "1000px", # Sets a fixed height
    navset_card_tab(
      nav_panel("Stacks",
                actionButton("add_data_button", "Add Stack", class = "btn-sm", width = '20%'),
                DTOutput("stack_table")  #(dynamic table)
                
      ),
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
      nav_panel("Doses",
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
                
                DTOutput("dose_table")  #(dynamic table)
      ),
      nav_panel("Export Data",
                downloadButton("downloadCSVData", "Download .CSV file", class = "btn-sm", width = '20%'),
                downloadButton("downloadReport", "Generate HTML Report", class = "btn-sm", width = '20%')
      )
    )
  )
)

server <- function(input, output, session) {
  

# STACKS ------------------------------------------------------------------
  
  # Reactive value to store the data frame
  stack_data <- reactiveVal(stacks)
  
  # Render the data table
  output$stack_table <- renderDT({
    datatable(stack_data(), 
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
  })
  
  # Observe the "Add Row" button click to show the modal
  observeEvent(input$add_data_button, {
    showModal(
      modalDialog(
        title = "Enter New Row Data",
        textInput("input_ID", "Stack ID:", ""),
        numericInput("input_RH", "Release Height:", value = NULL, min = 0),
        numericInput("input_BH", "Building Height:", value = NULL, min = 0),
        numericInput("input_BW", "Building Width:", value = NULL, min = 0),
        numericInput("input_SD", "Source Distance:", value = NULL, min = 0),
        numericInput("input_VD", "Veg Distance:", value = NULL, min = 0),
        numericInput("input_CD", "Cow Distance:", value = NULL, min = 0),
        numericInput("input_WS", "Mean WS:", value = NULL, min = 0),
        footer = tagList(
          actionButton("save_btn", "Save"),
          modalButton("Cancel")
        )
      )
    )
  })
  
  
  # Observe the "Save" button click within the modal
  observeEvent(input$save_btn, {
    # Get the current data
    current_data <- stack_data()
    
    # Create the new row
    new_row <- data.frame(
      `Stack ID` = input$input_ID,
      `Release Height` = input$input_RH,
      `Building Height` = input$input_BH,
      `Building Width` = input$input_BW,
      `Source Distance` = input$input_SD,
      `Veg Distance` = input$input_VD,
      `Cow Distance` = input$input_CD,
      `Mean WS`= input$input_WS,

      stringsAsFactors = FALSE
    )
    # Append the new row using bind_rows for safety with types
    updated_data <- bind_rows(current_data, new_row)
    cleaned_data <- updated_data[rowSums(is.na(updated_data)) < ncol(updated_data), ]

    
    # Update the reactive value, which automatically re-renders the table
    stack_data(cleaned_data)
    print(stack_data())
    
    
    # Close the modal
    removeModal()
  })

 
  
  
  # RADS ------------------------------------------------------------------
  
  # Reactive value to store the data frame
  rad_data <- reactiveVal(rads)
  
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
        textInput("input_ID", "Stack ID:", ""),
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
      `Stack ID` = input$input_ID,
      `Radionuclide` = input$input_RN,
      `Type` = input$input_TYPE,
      `Emission` = input$input_EMISSION,
      `Allowance` = subset(long_rads, Type == input$input_TYPE & Radionuclide == input$input_RN)[1,3],
      `Percent` = round((input$input_EMISSION/(subset(long_rads, Type == input$input_TYPE & Radionuclide == input$input_RN)[1,3]) * 100),2),

      stringsAsFactors = FALSE
    )
    print("new_row!")
    print(new_row)
    # Append the new row using bind_rows for safety with types
    updated_data <- bind_rows(current_data, new_row)
    cleaned_data <- updated_data[rowSums(is.na(updated_data)) < ncol(updated_data), ]
    
    
    # Update the reactive value, which automatically re-renders the table
    rad_data(cleaned_data)
    # Update the reactive value, which automatically re-renders the table
    print(rad_data())
    # Close the modal
    removeModal()
  })
  
  

# DOSES -------------------------------------------------------------------
  
  # Render the data table using the result of the eventReactive expression
  output$dose_table <- renderDT({
    req(dose_data())
    datatable(dose_data(), 
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
                )
              ))
  })
  
  
  
  
  
# BUILD COMPLETE DATA FRAME   ----------------------------------------------

 
    
    combined_data <- reactive({
      req(rad_data())
      # Use bind_rows to stack dataframes vertically (assumes same column names)
      combined_df <-  merge(stack_data(), rad_data(), by= "Stack.ID")

      # Remove rows with any NA values
      cleaned <- as.data.frame(combined_df[rowSums(is.na(combined_df[, -1])) < ncol(combined_df[, -1]), ])
      print(cleaned)
      
      
      #distances
      rep_length <- nrow(cleaned)
      receptors <- seq(5,2000,5)
      
      rads_distances <- cleaned %>%
        mutate(Source.Distance = NA) %>%
        mutate(row_id = row_number()) %>%
        uncount(400) %>% # Duplicates each row 400 times
        group_by(row_id) %>%
        #mutate(across(-Radionuclide, ~ ifelse(row_number() > 1, NA, .))) %>%
        ungroup() %>%
        mutate(Source.Distance = rep(seq(5,2000,5), times = rep_length))%>%
        select(-row_id) %>% mutate(
        area = Building.Width * Building.Height,
        xcrit = 2.5 * sqrt(area),
        xarea = pmin(Building.Width,Building.Height),
        sigmaz  = (0.06*Source.Distance)/sqrt(1.0+(0.0015*Source.Distance)),
        capsig  =  sqrt((sigmaz^2) + (area/3.14)),
        pp  =   round(2.032 / (Source.Distance*capsig),2), 
        conc = round(ifelse(Source.Distance < xcrit, (0.25*as.numeric(Emission))/(3.14*Mean.WS*1.0*xarea),
                            (0.25*as.numeric(Emission)*pp) / Mean.WS),2)
        #note:conc will be INF if WS = 0...
      )

      print(head(rads_distances))
      rad_distances_doses <- left_join(rads_distances, doses_info, by = "Radionuclide")

      # For each distance, calculate a total, inhalation, vegetable, and a cowmilk dose
      rad_distances_doses$dose <- rad_distances_doses$WBdose.y * rad_distances_doses$conc
      rad_distances_doses$vdose <- rad_distances_doses$VegDose * rad_distances_doses$conc
      rad_distances_doses$cdose <- rad_distances_doses$CowMilkDose * rad_distances_doses$conc
      rad_distances_doses$total <-   rad_distances_doses$dose +  rad_distances_doses$vdose + rad_distances_doses$cdose

      return(rad_distances_doses)
    
    
    })
  
  
  
  
  # make cut-down dose-specific version
  dose_data <- reactive({
    # Use bind_rows to stack dataframes vertically (assumes same column names)
    full_data <-  combined_data()
    df_subset <- full_data %>% select(`Stack.ID`, `Radionuclide`, `Type`, `Source.Distance`, `dose`, `vdose`, `cdose`, `total`)
    
    print(df_subset)
    return(df_subset)
  })
  
  
    
    # Render the interactive plot using renderPlotly and ggplotly
    output$graph <- renderPlotly({
      # Use req() to ensure data is available before plotting
      req(combined_data())
      # Call the reactive expression like a function: combined_data()
      p <- ggplot(combined_data(), aes(x = Source.Distance, y = total, color = Radionuclide)) +
        geom_point(alpha = 0.6) +
        labs(title = "Combined Reactive Data Plot", color = "Source Data") +
        theme_classic()+ 
        xlab("Distance from Source (m)")+
        ylab("Total Dose (mrem/y)")+
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
        write.csv(combined_data(), file, row.names = FALSE) 
      }
    )    
  
    ##### HTML Report
    # Reactive data for both app display and report
    report_data <-  reactive({combined_data()})
    report_plot <- reactive({
     p <-  ggplot(combined_data(), aes(x = Source.Distance, y = total, color = Radionuclide)) +
        geom_point(alpha = 0.6) +
        labs(title = "Combined Reactive Data Plot", color = "Source Data") +
        theme_classic()+
        xlab("Distance from Source (m)")+
        ylab("Total Dose (mrem/y)")+
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
    output$HTML_plot <- renderPlotly({ report_plot() })
    output$HTML_table <- renderDT({ report_data() })



    output$downloadReport <- downloadHandler(
      filename = function() { paste0("report_", Sys.Date(), ".html") },
      content = function(file) {
        # Use a temporary directory for the .Rmd to ensure write permissions
        temp_rmd <- file.path(tempdir(), "COMPLY_Output.Rmd")
        file.copy("COMPLY_Output.Rmd", temp_rmd, overwrite = TRUE)

        # Pass reactive objects to the R Markdown template
        rmarkdown::render(temp_rmd, output_file = file,
                          params = list(report_plot = report_plot(), report_data = report_data()),
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
        addProviderTiles(providers$Esri.WorldImagery) %>%
        setView(lng = -98.5833, lat = 39.8333, zoom = 4) # Default view (center of the US)
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
            zoom = 15
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
        "Enter an address and click 'Geocode'."
      }
    })    

    # Add heatmap 
    observeEvent(input$add_heatmap_btn, {
      req(combined_data())
      req(geocoded_address())
      
      result <- geocoded_address()
      result$latitude <- as.numeric(result$latitude)
      result$longitude <- as.numeric(result$longitude)
      print("heatmap button pressed!")
      print(result)
      
      
      rad_distances_doses <- combined_data()
      
      rad_distances_doses$lat <- result$latitude
      rad_distances_doses$long <- result$longitude
      
      R <- 6378137 # Earth radius in meters
      lat_rad <- result$latitude * pi / 180
      rad_distances_doses$new_lat <- rad_distances_doses$lat + (rad_distances_doses$Source.Distance / R) * (180 / pi)
      rad_distances_doses$new_lon <- rad_distances_doses$long + (rad_distances_doses$Source.Distance / R) * (180 / pi) / cos(lat_rad)
      print(rad_distances_doses)
      
      leafletProxy("map") %>%
        addHeatmap(
          data = rad_distances_doses,
          lng = ~rad_distances_doses$new_lon,
          lat = ~rad_distances_doses$new_lat,
          intensity = ~rad_distances_doses$total , # Optional: weighting by magnitude
          radius = 20,       # Radius of influence in pixels
          blur = 15,         # Blur effect
          max = 0.05
        )
      
    })


}


# Define the viewer with specific width and height (in pixels)
vwr <- dialogViewer('COMPLY', width = 3200, height = 7400)
# 
# # Run the app using runGadget and the defined viewer
#
runGadget(shinyApp(ui = ui, server = server), viewer = vwr)