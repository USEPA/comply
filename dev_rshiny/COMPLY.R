
### nuclide names // rates csv read in to drop down 
### option to also upload lat // long 
### 
setwd("C:/Users/MMACDO01/OneDrive - Environmental Protection Agency (EPA)/MM Projects & Code/COMPLY Shiny App/COMPLY")


library(shiny)
library(leaflet)
library(tidygeocoder)
library(DT)
library(reshape2)

# define these lookup tables for use in the app 
rads <- read.csv(file="comply171_original_possession_tables.csv",skip=0,header=T)
long_rads <- melt(rads, 
                  id.vars = "Nuclide", 
                  variable.name = "Form", 
                  value.name = "Value")
long_rads$Nuclide <- gsub(" ", "", long_rads$Nuclide)
long_rads$Form <- gsub(" ", "", long_rads$Form)



# Define the user interface (UI)
ui <- fluidPage(
  titlePanel("COMPLY Sample"),
  sidebarLayout(
    sidebarPanel(
      h4("Enter Facility Details:"),
      textInput("CompanyName", "Company Name:", ""),
      textInput("CompanyContact", "Company Contact Name:", ""),
      textInput("CompanyContactNum", "Company Contact Number:", ""),
      textInput("FacilityName", "Facility Name:", ""),
      h4("Enter Coordinates or Address:"),
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
      hr(),
      h4("Add Stacks:"),
      numericInput("releaseheight", "Release Height (m):", value = 0),
      numericInput("buildingheight", "Building Height (m):", value = 0),
      numericInput("buildingwidth", "Building Width (m):", value = 0),
      numericInput("distancesource", "Distance to Source (m):", value = 0),
      numericInput("distanceveg", "Distance to Vegetation (m):", value = 0),
      numericInput("distancecow", "Distance to Cows (m):", value = 0),
      numericInput("meanws", "Mean Wind Speed (m/s):", value = 0),
      actionButton("add_row_btn", "Add Stack"),
      fluidRow(
        column(6,
      selectInput(
        inputId = "radionuclide",
        label = "Radionuclide:",
        choices = unique(long_rads$Nuclide)
      )),
      column(6,
             selectInput(
               inputId = "Type",
               label = "Type:",
               choices = c("Solid", "Liquid", "Gas")
             ))
      ),
      numericInput("emission", "Emission(Ci/y):", value = 0),
      numericInput("stacknumber", "Stack #:", value = 0),
      
      actionButton("add_row_btn_rn", "Add Radionuclide")
      
      ### option for more radionuclides
      ),
    mainPanel(
      # Output: A Leaflet map showing the geocoded location
      leafletOutput("map"),
      # Output: Display coordinates
      h4("Coordinates:"),
      textOutput("coords"),
      h4("Stacks:"),
      DTOutput("data_table"), # Use DTOutput for interactive table
      h4("Radionuclides:"),
      DTOutput("rn_table"), # Use DTOutput for interactive table
      textOutput("possesion")
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  # Reactive expression to perform geocoding
  geocoded_address <- eventReactive(input$geocode_button, {
    # Combine the separate address components into a single string
    # full_address <- paste(
    #   input$street,
    #   input$city,
    #   input$state,
    #   input$zip,
    #   sep = ", "
    # )
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
  

  # Initialize reactive values for the table data: STACKS
  table_data <- reactiveValues(df = data.frame(Name = character(0)))
  
  # Observe the "Add Row" button click
  observeEvent(input$add_row_btn, {
    counter_value <- nrow(table_data$df) + 1
    
    # Create a new row from input values
    new_row <- data.frame(Stack.No = counter_value,
                          Release.Height = input$releaseheight, 
                          Building.Height = input$buildingheight,
                          Building.Width = input$buildingwidth,
                          Distance.Source = input$distancesource,
                          Distance.Veg = input$distanceveg,
                          Distance.Cow = input$distancecow,
                          Mean.WS = input$meanws)
    
    # Append the new row to the existing data frame
    table_data$df <- rbind(table_data$df, new_row)
    
    # Clear input fields after adding (optional)
    updateNumericInput(session, "releaseheight", value = 0)
    updateNumericInput(session, "buildingheight", value = 0)
    updateNumericInput(session, "buildingwidth", value = 0)
    updateNumericInput(session, "distancesource", value = 0)
    updateNumericInput(session, "distanceveg", value = 0)
    updateNumericInput(session, "distancecow", value = 0)
    updateNumericInput(session, "meanws", value = 0)
    
    counter_value <- nrow(table_data$df) + 1
  })
  
  # Render the table
  output$data_table <- renderDT({ # Use renderDT for interactive table
    datatable(
      table_data$df, rownames = FALSE,
      options = list(
        searching = FALSE, # Removes the global search box
        info = FALSE,       # Removes the "Showing X to Y of Z entries" display
        lengthChange = FALSE, # This removes the 'Show entries' dropdown
        dom = 'lfti' # Excludes 'p' for pagination
      )
    )

  })
  
  
 
  # Initialize reactive values for the table data: RADIONUCLIDES
  table_rn <- reactiveValues(df = data.frame(Name = character(0), Age = numeric(0)))
  
  # Observe the "Add Row" button click
  observeEvent(input$add_row_btn_rn, {
    # Create a new row from input values
    type <- input$Type
    print(type)
    nuclide <-input$radionuclide
    print(nuclide)
    allowance_val <- subset(long_rads, Form == type & Nuclide == nuclide)[1,3]
    allowance_val <- as.numeric(allowance_val)
    print(allowance_val)
    new_row_rn <- data.frame(stack.Number = input$stacknumber,
                          Radionuclide = input$radionuclide, 
                          Phase = type,
                          Emission = input$emission,
                          Allowance = allowance_val,
                          Percent = round((input$emission/allowance_val * 100),2))
    names(new_row_rn) <- c("Stack", "Nuclide","Phase","Possession Amount (Ci/y)","Legal Limit (Ci/y)","Percent of Limit")
    print(new_row_rn)
    
    
   
    
    
        
    # Append the new row to the existing data frame
    table_rn$df <- rbind(table_rn$df, new_row_rn)
    
    # Clear input fields after adding (optional)
    updateNumericInput(session, "decayrate", value = 0)
    updateNumericInput(session, "stacknumber", value = 0)
    updateTextInput(session, "radionuclide", value = " ")
    
  })
  
  # Render the table
  output$rn_table <- renderDT({ # Use renderDT for interactive table
    # table_data$df
    datatable(
      table_rn$df,
      options = list(
        searching = FALSE, # Removes the global search box
        info = FALSE,       # Removes the "Showing X to Y of Z entries" display
        lengthChange = FALSE, # This removes the 'Show entries' dropdown
        dom = 'lfti' # Excludes 'p' for pagination
      )
    )
  })
  
  
  
  # Observe changes in the radius input and update the circle
  observeEvent(input$distancesource, {
    
    
    result <- geocoded_address()
  
    
    
    
    leafletProxy("map") %>%
     # clearShapes() %>% # Clear existing circles
      addCircles(
        lng = result$longitude,
        lat = result$latitude,
        radius = input$distancesource,
        color = "blue",
        fillColor = "lightblue",
        fillOpacity = 0.5,
        popup = paste("Distance:", input$distancesource, "meters")
      )
  })
  
  output$possesion <- renderText({
    type <- input$Type
    nuclide <-input$radionuclide
    allowance_val <- subset(long_rads, Form == type & Nuclide == nuclide)[1,3]
    allowance_val <- as.numeric(allowance_val)
    Emission = input$emission
    Percent = input$emission/allowance_val * 100
    
    if (Percent > 100) {
      paste("You possess", round(input$emission/allowance_val, 1), " times the allowable amount given in the possession table.")
    } else {
      "You are within possession limits :D ."
    }
  })
  # 
  # 
  
  
  
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

