# =========================================================
# COMPLY v2.0
# Author: MacDonald, Megan | RTP, NC
# Requires (in working directory):
#   - comply171_original_possession_data.csv
#   - comply171_original_input_data.csv
#   - COMPLY_Output.Rmd
# =========================================================

library(shiny)
library(bslib)
library(dplyr)
library(tibble)
library(DT)
library(plotly)
library(reshape2)
library(leaflet)
library(tidygeocoder)
library(openair)
library(rmarkdown)
library(sf)
library(raster)
library(akima)


# Load Files --------------------------------------------------------------
POSSESSION_CSV <- "comply171_original_possession_data.csv"  # possession table
DOSES_CSV      <- "comply171_original_input_data.csv"        # dose coefficients
REPORT_RMD     <- "COMPLY_Output.Rmd"                        # report template


# DT Styling --------------------------------------------------------------
dt_css <- "
  .dataTable td, .dataTable th { border: 1px solid black; }
  .dt-button {
    background-color: #f0f0f0 !important;
    color: black !important;
    border: 1px solid #000500 !important;
  }
"


# Prep Data Frames --------------------------------------------------------
rads_src <- read.csv(POSSESSION_CSV,
                     header = TRUE,
                     stringsAsFactors = FALSE)
long_rads <- melt(
  rads_src,
  id.vars = "Nuclide",
  variable.name = "Form",
  value.name = "Allowance"
) %>%
  mutate(Radionuclide = gsub(" ", "", Nuclide),
         Type         = gsub(" ", "", Form)) %>%
  dplyr::select(Radionuclide, Type, Allowance)

doses_info <- read.csv(DOSES_CSV, header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(Radionuclide = trimws(Nuclide))


rads <- data.frame(
  Stack.ID = character(),
  Radionuclide = character(),
  Type = character(),
  Emission = numeric(),
  Allowance = numeric(),
  Fraction = numeric(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stacks <- data.frame(
  Stack.ID = character(),
  Release.Height = numeric(),
  Building.Height = numeric(),
  Building.Width = numeric(),
  Source.Distance = numeric(),
  Veg.Distance = numeric(),
  Cow.Distance = numeric(),
  Mean.WS = numeric(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

xsum <- data.frame(Fraction = as.numeric(NA))




# ---- UI ----
ui <- page_sidebar(
  title = "COMPLY v2.0",
  sidebar = sidebar(
    h6("Enter Facility Details:"),
    textInput("CompanyName", "Company Name:", ""),
    textInput("CompanyContact", "Company Contact Name:", ""),
    textInput("CompanyContactNum", "Company Contact Email:", ""),
    textInput("FacilityName", "Facility Name:", ""),
    h6("Enter Coordinates or Address:"),
    fluidRow(column(6, textInput(
      "Latitude", "Latitude:", ""
    )), column(
      6, textInput("Longitude", "Longitude:", "")
    )),
    textInput(
      "address_input",
      "Facility Address:",
      "109 T.W. Alexander Drive, Durham, NC, 27709"
    ),
    actionButton("geocode_button", "Plot Location"),
    h6("Coordinates:"),
    textOutput("coords")
  ),
  
  ## Level 1 ---------------------------------------------------------
  navset_tab(
    id = "level_tabs",
    selected = "level 1",
    
    nav_panel(
      "level 1",
      card(
        full_screen = TRUE,
        height = "500px",
        navset_card_tab(
          nav_panel(
            "Rads",
            tags$div(
              class = "d-flex gap-2 mb-2",
              actionButton("add_rad_button", "Add Radionuclide", class = "btn-sm"),
              actionButton("delete_rad_button", "Delete Selected Row", class = "btn-sm")
            ),
            tags$head(tags$style(HTML(dt_css))),
            DTOutput("rad_table")
          ),
          nav_panel(
            "Compliance Summary",
            textOutput("fraction_display"),
            textOutput("fraction_display_i")
          ),
          nav_panel(
            "Export Data",
            downloadButton(
              "downloadReport",
              "Generate HTML Report",
              class = "btn-sm",
              width = "20%"
            )
          )
        )
      ),
      layout_columns(
        card(
          full_screen = TRUE,
          height = "500px",
          card_header("Map"),
          leafletOutput("map")
        ),
        card(
          full_screen = TRUE,
          height = "500px",
          card_header("Plot"),
          plotlyOutput(outputId = "graph")
        ),
        col_widths = c(6, 6)
      )
    ),
    
    ## Level 2 / 3 ---------------------------------------------------------
    nav_panel(
      "level 2/3",
      layout_columns(
        card(
          full_screen = TRUE,
          height = "500px",
          card_header("Map"),
          leafletOutput("map_l2")
        ),
        card(
          full_screen = TRUE,
          height = "500px",
          card_header("Plot"),
          plotlyOutput(outputId = "graph_l2")
        ),
        col_widths = c(6, 6)
      ),
      card(
        full_screen = TRUE,
        height = "500px",
        navset_card_tab(
          nav_panel(
            "Stacks",
            tags$div(
              class = "d-flex gap-2 mb-1",
              actionButton("add_data_button", "Add Stack", class = "btn-sm"),
              actionButton("delete_stack_button", "Delete Selected Row", class = "btn-sm")
            ),
            tags$small(
              class = "text-muted",
              "Note: For Level 2, use same distance values for Source Distance, Veg Distance, and Cow Distance."
            ),
            DTOutput("stack_table")
          ),
          nav_panel(
            "Emissions",
            tags$div(
              class = "d-flex gap-2 mb-2",
              actionButton("add_rad_button_l2", "Add Radionuclide", class = "btn-sm"),
              actionButton("delete_rad_button_l2", "Delete Selected Row", class = "btn-sm")
            ),
            tags$head(tags$style(HTML(dt_css))),
            DTOutput("rad_table_l2")
          ),
          nav_panel(
            "Diagnostic Results",
            tags$head(tags$style(HTML(dt_css))),
            DTOutput("dose_table")
          ),
          nav_panel(
            "Compliance Summary",
            textOutput("fraction_display_l2"),
            textOutput("fraction_display_i_l2")
          ),
          nav_panel(
            "Export Data",
            downloadButton(
              "downloadCSVData",
              "Download .CSV file",
              class = "btn-sm",
              width = "20%"
            ),
            downloadButton(
              "downloadReport",
              "Generate HTML Report",
              class = "btn-sm",
              width = "20%"
            )
          )
        )
      )
    ),
    
    ## Level 4 / Radon ---------------------------------------------------------
    nav_panel(
      "level 4/Radon",
      card(
        full_screen = TRUE,
        height = "160px",
        card_header("Radon Data Inputs & Parameters"),
        fluidRow(
          column(
            width = 4,
            fileInput(
              "level4_input_file",
              "Wind Data File (.txt, .csv)",
              accept = c(".txt", ".csv")
            )
          ),
          column(
            width = 4,
            tags$div(
              class = "d-flex align-items-center",
              actionButton("radon_params_btn", "Set Radon Parameters", class = "btn-sm")
            )
          ),
          column(
            width = 4,
            tags$div(
              class = "d-flex align-items-center",
              actionButton("calculate_radon", "Calculate Radon Sectors", class = "btn-sm")
            )
          )
        )
      ),
      
      navset_card_tab(
        nav_panel("Radon Data & Parameters", DTOutput("radon_allout_table")),
        nav_panel("Meterological Data", DTOutput("wind_data_table")),
        
        nav_panel(
          "Compliance Summary",
          # add placeholder or actual outputs
          textOutput("radon_summary")
        ),
        
        nav_panel(
          "Export Data"
          # download buttons commented out
          # downloadButton("downloadRadonCSV", "Download .CSV", class="btn-sm")
        )
        ),
                  layout_columns(
                    card(
                      full_screen = TRUE,
                      height = "450px",
                      card_header("Wind Rose"),
                      plotOutput("radon_wind_rose", height = "400px")
                    ),
                    card(
                      full_screen = TRUE,
                      height = "450px",
                      card_header("Radon Dose Map"),
                      leafletOutput("radon_dose_map", height = "400px")
                    ),
                    col_widths = c(6, 6)
                  )
        )
      )
    )
  

    
    
    # ---- Server ----
    server <- function(input, output, session) {
      # Store Empty Vars
      stack_diameter    <- reactiveVal(NA_real_)
      exit_velocity     <- reactiveVal(NA_real_)
      release_height    <- reactiveVal(NA_real_)
      stack_temp        <- reactiveVal(NA_real_)
      air_temp          <- reactiveVal(NA_real_)
      vertical_momentum <- reactiveVal(NA_character_)
      emission          <- reactiveVal(NA_real_)
      percent_calm      <- reactiveVal(NA_real_)
      
      
      # LEVEL 1 -----------------------------------------------------------------
      
      
      ## ===== STACKS (Level 1/2) =====
      stack_data <- reactiveVal(stacks)
      
      ### Stacks Table ------------------------------------------------------------
      #Stacks tabs in Level 1 and Level 2/3
      output$stack_table <- renderDT({
        datatable(
          stack_data(),
          class = "cell-border stripe",
          extensions = "Buttons",
          selection = "single",
          options = list(
            ordering = FALSE,
            searching = FALSE,
            info = FALSE,
            lengthChange = FALSE,
            dom = "Blfti",
            columnDefs = list(list(
              className = "dt-center", targets = "_all"
            )),
            rownames = FALSE,
            pageLength = -1,
            paging = FALSE,
            filter = "top",
            buttons = c("copy", "csv", "excel", "pdf"),
            exportOptions = list(columns = ":visible")
          )
        )
      })
      
      ### Stacks Modal ------------------------------------------------------------
      observeEvent(input$add_data_button, {
        showModal(
          modalDialog(
            title = "Enter New Stack Data",
            textInput("input_ID", "Stack ID:", ""),
            numericInput(
              "input_RH",
              "Release Height:",
              value = NULL,
              min = 0
            ),
            numericInput(
              "input_BH",
              "Building Height:",
              value = NULL,
              min = 0
            ),
            numericInput(
              "input_BW",
              "Building Width:",
              value = NULL,
              min = 0
            ),
            numericInput(
              "input_SD",
              "Source Distance:",
              value = NULL,
              min = 0
            ),
            numericInput(
              "input_VD",
              "Veg Distance:",
              value = NULL,
              min = 0
            ),
            numericInput(
              "input_CD",
              "Cow Distance:",
              value = NULL,
              min = 0
            ),
            numericInput("input_WS", "Mean WS:", value = NULL, min = 0),
            tags$small(
              class = "text-muted",
              "Note: For Level 2, use same distance values for Source Distance, Veg Distance, and Cow Distance."
            ),
            footer = tagList(actionButton("save_btn", "Save"), modalButton("Cancel")),
            easyClose = TRUE
          )
        )
      })
      observeEvent(input$save_btn, {
        current <- stack_data()
        new_row <- data.frame(
          Stack.ID = input$input_ID,
          Release.Height = input$input_RH,
          Building.Height = input$input_BH,
          Building.Width = input$input_BW,
          Source.Distance = input$input_SD,
          Veg.Distance = input$input_VD,
          Cow.Distance = input$input_CD,
          Mean.WS = input$input_WS,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        updated <- bind_rows(current, new_row)
        cleaned <- updated[rowSums(is.na(updated)) < ncol(updated), ]
        stack_data(cleaned)
        removeModal()
      })
      observeEvent(input$delete_stack_button, {
        sel <- input$stack_table_rows_selected
        if (is.null(sel) || !length(sel)) {
          showModal(
            modalDialog(
              title = "No row selected",
              "Click a row in the Stacks table to select it, then press 'Delete Selected Row'.",
              easyClose = TRUE
            )
          )
        } else {
          current <- stack_data()
          updated <- current[-sel, , drop = FALSE]
          cleaned <- updated[rowSums(is.na(updated)) < ncol(updated), , drop = FALSE]
          stack_data(cleaned)
        }
      })
      
      
      ## ===== RADS (POSSESSION) =====
      rad_data           <- reactiveVal(rads)
      rad_sum_possession <- reactiveVal(xsum)
      i_sum_possession   <- reactiveVal(xsum)
      
      ### Rads Table --------------------------------------------------------------
      output$rad_table <- renderDT({
        datatable(
          rad_data(),
          class = "cell-border stripe",
          extensions = "Buttons",
          selection = "single",
          options = list(
            ordering = FALSE,
            searching = FALSE,
            info = FALSE,
            lengthChange = FALSE,
            dom = "Blfti",
            columnDefs = list(list(
              className = "dt-center", targets = "_all"
            )),
            rownames = FALSE,
            pageLength = -1,
            paging = FALSE,
            filter = "top",
            buttons = c("copy", "csv", "excel", "pdf"),
            exportOptions = list(columns = ":visible")
          )
        )
      })
      
      
      ### Rads Modal --------------------------------------------------------------
      observeEvent(input$add_rad_button, {
        showModal(
          modalDialog(
            title = "Enter New Radionuclide Data",
            textInput("input_STACK", "Stack ID:", ""),
            selectInput(
              "input_RN",
              "Radionuclide:",
              choices = unique(long_rads$Radionuclide)
            ),
            selectInput("input_TYPE", "Type:", choices = unique(long_rads$Type)),
            numericInput(
              "input_EMISSION",
              "Emission:",
              value = NULL,
              min = 0
            ),
            footer = tagList(
              actionButton("save_rad_btn", "Save"),
              modalButton("Cancel")
            )
          )
        )
      })
      observeEvent(input$delete_rad_button, {
        sel <- input$rad_table_rows_selected
        if (is.null(sel) || !length(sel)) {
          showModal(
            modalDialog(
              title = "No row selected",
              "Click a row in the table to select it, then press 'Delete Selected Row'.",
              easyClose = TRUE
            )
          )
        } else {
          current <- rad_data()
          updated <- current[-sel, , drop = FALSE]
          cleaned <- updated[rowSums(is.na(updated)) < ncol(updated), , drop = FALSE]
          rad_data(cleaned)
          if (!nrow(cleaned)) {
            rad_sum_possession(data.frame(Fraction = 0))
            i_sum_possession(data.frame(Fraction = 0))
          } else {
            temp <- cleaned
            temp$element <- substr(temp$Radionuclide, 1, 1)
            iodine <- subset(temp, element == "I")
            rad_sum_possession(data.frame(Fraction = sum(temp$Fraction, na.rm = TRUE)))
            i_sum_possession(data.frame(Fraction = sum(iodine$Fraction, na.rm = TRUE)))
          }
        }
      })
      observeEvent(input$save_rad_btn, {
        current <- rad_data()
        
        allow_row <- subset(long_rads,
                            Type == input$input_TYPE &
                              Radionuclide == input$input_RN)
        if (!nrow(allow_row) || is.na(allow_row$Allowance[1])) {
          showModal(
            modalDialog(
              title = "No matching allowance",
              "No allowance found for that Radionuclide and Type. Please check your selection.",
              easyClose = TRUE
            )
          )
          return(invisible(NULL))
        }
        allow_val <- allow_row$Allowance[1]
        
        new_row <- data.frame(
          Stack.ID = input$input_STACK,
          Radionuclide = input$input_RN,
          Type = input$input_TYPE,
          Emission = input$input_EMISSION,
          Allowance = allow_val,
          Fraction  = round((input$input_EMISSION / allow_val), 5),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        
        updated <- bind_rows(current, new_row)
        cleaned <- updated[rowSums(is.na(updated)) < ncol(updated), ]
        temp <- cleaned
        temp$element <- substr(temp$Radionuclide, 1, 1)
        iodine <- subset(temp, element == "I")
        
        rad_data(cleaned)
        rad_sum_possession(data.frame(Fraction = sum(cleaned$Fraction, na.rm = TRUE)))
        i_sum_possession(data.frame(Fraction = sum(iodine$Fraction, na.rm = TRUE)))
        removeModal()
      })
      
      
      
      
      ## L1 Compliance Summary  -------------------------------------------------
      output$fraction_display <- renderText({
        #l1
        selected <- rad_sum_possession()
        val <- round(selected$Fraction, 2)
        paste0(
          "Stack possession fraction is ",
          val,
          " times the allowable amount given in the possession table."
        )
      })
      output$fraction_display_i <- renderText({
        selected <- i_sum_possession()
        val <- round(selected$Fraction, 2)
        paste0("Total iodine possession fraction is ", val, ".")
      })
      
      
      ## L1 Bar Graph ------------------------------------------------------------
      output$graph <- renderPlotly({
        df <- as.data.frame(rad_data()) %>%
          dplyr::filter(
            !is.na(Fraction),
            is.finite(Fraction),
            !is.na(Radionuclide),
            nzchar(Radionuclide)
          )
        if (!nrow(df)) {
          return(
            plotly::plotly_empty(type = "scatter", mode = "markers") |>
              layout(title = "No valid fraction data to plot yet.")
          )
        }
        agg <- df %>%
          group_by(Radionuclide) %>%
          summarise(Fraction = sum(Fraction, na.rm = TRUE),
                    .groups = "drop")
        plotly::plot_ly(
          data = agg,
          x = ~ Radionuclide,
          y = ~ Fraction,
          type = "bar",
          hoverinfo = "x+y"
        ) |>
          layout(
            title = "Possession Fraction by Radionuclide",
            xaxis = list(title = "Radionuclide"),
            yaxis = list(title = "Fraction"),
            showlegend = FALSE
          )
      })
      
      
      ## L1 Export ---------------------------------------------------------------
      report_meta <- reactive({
        res   <- tryCatch({
          # geocode if lat blank
          if (identical(input$Latitude, "")) {
            some_addresses <- tibble::tibble(name = "company", addr = input$address_input)
            tidygeocoder::geocode(
              some_addresses,
              address = addr,
              method = "osm",
              lat = latitude,
              long = longitude
            )
          } else {
            tibble::tibble(
              name = "company",
              addr = input$address_input,
              latitude = as.numeric(input$Latitude),
              longitude = as.numeric(input$Longitude)
            )
          }
        }, error = function(e)
          NULL)
        lat_in <- suppressWarnings(as.numeric(input$Latitude))
        lon_in <- suppressWarnings(as.numeric(input$Longitude))
        lat <- if (!is.null(res) &&
                   "latitude"  %in% names(res))
          as.numeric(res$latitude)
        else
          lat_in
        lon <- if (!is.null(res) &&
                   "longitude" %in% names(res))
          as.numeric(res$longitude)
        else
          lon_in
        list(
          company_name = input$CompanyName,
          company_contact = input$CompanyContact,
          company_contact_email = input$CompanyContactNum,
          facility_name = input$FacilityName,
          address = input$address_input,
          latitude = lat,
          longitude = lon
        )
      })
      report_data <- reactive({
        rad_data()
      })
      report_plot <- reactive({
        df <- as.data.frame(rad_data()) %>%
          filter(
            !is.na(Fraction),
            is.finite(Fraction),
            !is.na(Radionuclide),
            nzchar(Radionuclide)
          )
        if (!nrow(df)) {
          return(
            plotly::plotly_empty(type = "scatter", mode = "markers") |>
              layout(title = "Possession Fraction Plot")
          )
        }
        agg <- df %>%
          group_by(Radionuclide) %>%
          summarise(Fraction = sum(Fraction, na.rm = TRUE),
                    .groups = "drop")
        plotly::plot_ly(
          data = agg,
          x = ~ Radionuclide,
          y = ~ Fraction,
          type = "bar",
          hoverinfo = "x+y"
        ) |>
          layout(
            title = "Possession Fraction by Radionuclide",
            xaxis = list(title = "Radionuclide"),
            yaxis = list(title = "Fraction"),
            showlegend = FALSE
          )
      })
      output$downloadReport <- downloadHandler(
        filename = function() {
          paste0("report_", Sys.Date(), ".html")
        },
        content  = function(file) {
          if (!file.exists(REPORT_RMD))
            stop("Missing file: ", REPORT_RMD)
          temp_rmd <- file.path(tempdir(), "COMPLY_Output.Rmd")
          file.copy(REPORT_RMD, temp_rmd, overwrite = TRUE)
          rmarkdown::render(
            temp_rmd,
            output_file = file,
            params = list(
              report_data = report_data(),
              report_plot = report_plot(),
              meta = report_meta()
            ),
            envir = new.env(parent = globalenv())
          )
        }
      )     
      
      
      
      
      ## LEAFLET: L1/L2 ----------------------------------------------------------
      output$map <- renderLeaflet({
        leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap) %>%
          setView(lng = -98.5833,
                  lat = 39.8333,
                  zoom = 5)
      })
      output$map_l2 <- renderLeaflet({
        leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap) %>%
          setView(lng = -98.5833,
                  lat = 39.8333,
                  zoom = 5)
      })
      
      observeEvent(input$geocode_button, {
        # Attempt address geocode unless lat/lon provided
        if (identical(input$Latitude, "")) {
          some_addresses <- tibble::tibble(name = "company", addr = input$address_input)
          geocode_result <- tidygeocoder::geocode(
            some_addresses,
            address = addr,
            method = "osm",
            lat = latitude,
            long = longitude
          )
          if (nrow(geocode_result) &&
              !is.na(geocode_result$latitude) &&
              !is.na(geocode_result$longitude)) {
            lat <- as.numeric(geocode_result$latitude)
            lon <- as.numeric(geocode_result$longitude)
            for (id in c("map", "map_l2")) {
              leafletProxy(id) %>%
                clearMarkers() %>%
                addMarkers(lng = lon,
                           lat = lat,
                           popup = geocode_result$addr) %>%
                setView(lng = lon,
                        lat = lat,
                        zoom = 10)
            }
            output$coords <- renderText({
              paste("Latitude:",
                    round(lat, 6),
                    ", Longitude:",
                    round(lon, 6))
            })
          } else {
            showModal(
              modalDialog(
                title = "Geocoding Error",
                "Could not find coordinates for the entered address. Please check the address and try again.",
                easyClose = TRUE
              )
            )
          }
        } else {
          lat <- as.numeric(input$Latitude)
          lon <- as.numeric(input$Longitude)
          for (id in c("map", "map_l2")) {
            leafletProxy(id) %>%
              clearMarkers() %>%
              addMarkers(lng = lon,
                         lat = lat,
                         popup = input$address_input) %>%
              setView(lng = lon,
                      lat = lat,
                      zoom = 10)
          }
          output$coords <- renderText({
            paste("Latitude:",
                  round(lat, 6),
                  ", Longitude:",
                  round(lon, 6))
          })
        }
      })
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      # LEVEL 2/3 ---------------------------------------------------------------
      
      
      ## Rad Emissions Table ---------------------------------------------------------
      # Emissions tab
      output$rad_table_l2 <- renderDT({
        datatable(
          rad_data(),
          class = "cell-border stripe",
          extensions = "Buttons",
          selection = "single",
          options = list(
            ordering = FALSE,
            searching = FALSE,
            info = FALSE,
            lengthChange = FALSE,
            dom = "Blfti",
            columnDefs = list(list(
              className = "dt-center", targets = "_all"
            )),
            rownames = FALSE,
            pageLength = -1,
            paging = FALSE,
            filter = "top",
            buttons = c("copy", "csv", "excel", "pdf"),
            exportOptions = list(columns = ":visible")
          )
        )
      })
      
      
      ## Rad Modal (L2) ----------------------------------------------------------
      observeEvent(input$add_rad_button_l2, {
        sd <- isolate(stack_data())
        stack_choices <- character(0)
        if (!is.null(sd) && "Stack.ID" %in% names(sd)) {
          ids <- sd$`Stack.ID`
          stack_choices <- unique(ids[!is.na(ids) && nzchar(ids)])
        }
        if (!length(stack_choices)) {
          showModal(
            modalDialog(
              title = "No stacks found",
              "Please add at least one Stack in the Stacks tab before adding emissions.",
              easyClose = TRUE
            )
          )
          return(invisible(NULL))
        }
        showModal(
          modalDialog(
            title = "Enter New Radionuclide Data",
            selectInput("input_STACK", "Stack ID:", choices = stack_choices),
            selectInput(
              "input_RN",
              "Radionuclide:",
              choices = unique(long_rads$Radionuclide)
            ),
            selectInput("input_TYPE", "Type:", choices = unique(long_rads$Type)),
            numericInput(
              "input_EMISSION",
              "Emission:",
              value = NULL,
              min = 0
            ),
            footer = tagList(
              actionButton("save_rad_btn", "Save"),
              modalButton("Cancel")
            )
          )
        )
      })
      observeEvent(input$delete_rad_button_l2, {
        sel <- input$rad_table_l2_rows_selected
        if (is.null(sel) || !length(sel)) {
          showModal(
            modalDialog(
              title = "No row selected",
              "Click a row in the Emissions table to select it, then press 'Delete Selected Row'.",
              easyClose = TRUE
            )
          )
        } else {
          current <- rad_data()
          updated <- current[-sel, , drop = FALSE]
          cleaned <- updated[rowSums(is.na(updated)) < ncol(updated), , drop = FALSE]
          rad_data(cleaned)
          
          if (!nrow(cleaned)) {
            rad_sum_possession(data.frame(Fraction = 0))
            i_sum_possession(data.frame(Fraction = 0))
          } else {
            temp <- cleaned
            temp$element <- substr(temp$Radionuclide, 1, 1)
            iodine <- subset(temp, element == "I")
            rad_sum_possession(data.frame(Fraction = sum(temp$Fraction, na.rm = TRUE)))
            i_sum_possession(data.frame(Fraction = sum(iodine$Fraction, na.rm = TRUE)))
          }
        }
      })
      
      
      
      ## L2 Dose-by-distance grid/table ---------------------------------------
      rad_sum_dose <- reactiveVal(data.frame(total = 0))
      i_sum_dose   <- reactiveVal(data.frame(total = 0))
      
      combined_data_l2 <- reactive({
        s <- stack_data()
        r <- rad_data()
        req(!is.null(s), !is.null(r))
        req(nrow(s) > 0, nrow(r) > 0)
        names(s) <- make.names(names(s))
        names(r) <- make.names(names(r))
        s <- s[rowSums(is.na(s)) < ncol(s), , drop = FALSE]
        r <- r[rowSums(is.na(r)) < ncol(r), , drop = FALSE]
        req("Stack.ID" %in% names(s), "Stack.ID" %in% names(r))
        merged <- merge(s, r, by = "Stack.ID")
        req(nrow(merged) > 0)
        
        n_rows <- nrow(merged)
        df <- merged %>%
          dplyr::slice(rep(dplyr::row_number(), each = 400)) %>%
          dplyr::mutate(Source.Distance = rep(seq(5, 2000, 5), times = n_rows)) %>%
          dplyr::mutate(
            area  = Building.Width * Building.Height,
            xcrit = 2.5 * sqrt(area),
            xarea = pmin(Building.Width, Building.Height),
            sigmaz = (0.06 * Source.Distance) / sqrt(1.0 + (0.0015 * Source.Distance)),
            capsig = sqrt((sigmaz ^ 2) + (area / 3.14)),
            pp     = 2.032 / (Source.Distance * capsig),
            conc   = ifelse(
              Source.Distance < xcrit,
              (0.25 * as.numeric(Emission)) / (3.14 * Mean.WS * 1.0 * xarea),
              (0.25 * as.numeric(Emission) * pp) / Mean.WS
            )
          ) %>%
          dplyr::left_join(doses_info, by = "Radionuclide") %>%
          dplyr::mutate(
            airconc = conc,
            dose    = WBdose.y    * conc,
            vdose   = VegDose     * conc,
            cdose   = CowMilkDose * conc,
            total   = dose + vdose + cdose
          )
        df
      })
      #Dignostic Results tab
      output$dose_table <- renderDT({
        req(combined_data_l2())
        datatable(
          combined_data_l2(),
          class = "cell-border stripe",
          extensions = "Buttons",
          selection = "none",
          options = list(
            ordering = FALSE,
            searching = FALSE,
            info = FALSE,
            lengthChange = FALSE,
            dom = "Blfti",
            columnDefs = list(list(
              className = "dt-center", targets = "_all"
            )),
            rownames = FALSE,
            pageLength = -1,
            paging = FALSE,
            filter = "top",
            buttons = c("copy", "csv", "excel", "pdf"),
            exportOptions = list(columns = ":visible")
          )
        )
      })
      
      dose_data <- reactive({
        combined_data_l2() %>%
          dplyr::select(Stack.ID,
                        Radionuclide,
                        Source.Distance,
                        airconc,
                        dose,
                        vdose,
                        cdose,
                        total)
      })
      
      
      ## Level 2 CSV Export ------------------------------------------------------
      output$downloadCSVData <- downloadHandler(
        filename = function() {
          paste0("comply_l2_dose_", Sys.Date(), ".csv")
        },
        content  = function(file) {
          write.csv(dose_data(), file, row.names = FALSE)
        }
      )
      
      distance_data <- reactive({
        s <- stack_data()
        names(s) <- make.names(names(s))
        r <- rad_data()
        names(r) <- make.names(names(r))
        s <- s[rowSums(is.na(s)) < ncol(s), , drop = FALSE]
        r <- r[rowSums(is.na(r)) < ncol(r), , drop = FALSE]
        req("Stack.ID" %in% names(s), "Stack.ID" %in% names(r))
        xrads <- merge(s, r, by = "Stack.ID")
        req(nrow(xrads) > 0)
        xrads$area  <- xrads$Building.Width * xrads$Building.Height
        xrads$xcrit <- 2.5 * sqrt(xrads$area)
        xrads$xarea <- pmin(xrads$Building.Width, xrads$Building.Height)
        calc_conc <- function(dist,
                              area,
                              xcrit,
                              xarea,
                              Emission,
                              Mean.WS) {
          sigmaz <- (0.06 * dist) / sqrt(1.0 + (0.0015 * dist))
          capsig <- sqrt((sigmaz ^ 2) + (area / 3.14))
          pp     <- 2.032 / (dist * capsig)
          ifelse(
            dist < xcrit,
            (0.25 * as.numeric(Emission)) / (3.14 * Mean.WS * 1.0 * xarea),
            (0.25 * as.numeric(Emission) * pp) / Mean.WS
          )
        }
        xrads$conc.mei <- calc_conc(
          xrads$Source.Distance,
          xrads$area,
          xrads$xcrit,
          xrads$xarea,
          xrads$Emission,
          xrads$Mean.WS
        )
        xrads$conc.veg <- calc_conc(
          xrads$Veg.Distance,
          xrads$area,
          xrads$xcrit,
          xrads$xarea,
          xrads$Emission,
          xrads$Mean.WS
        )
        xrads$conc.cow <- calc_conc(
          xrads$Cow.Distance,
          xrads$area,
          xrads$xcrit,
          xrads$xarea,
          xrads$Emission,
          xrads$Mean.WS
        )
        xrad_doses <- left_join(xrads, doses_info, by = "Radionuclide")
        xrad_doses$airconc <- xrad_doses$conc.mei
        xrad_doses$dose    <- xrad_doses$WBdose.y    * xrad_doses$conc.mei
        xrad_doses$vdose   <- xrad_doses$VegDose     * xrad_doses$conc.veg
        xrad_doses$cdose   <- xrad_doses$CowMilkDose * xrad_doses$conc.cow
        xrad_doses$total   <- xrad_doses$dose + xrad_doses$vdose + xrad_doses$cdose
        total_sum  <- sum(xrad_doses$total, na.rm = TRUE)
        temp <- xrad_doses
        temp$element <- substr(temp$Radionuclide, 1, 1)
        iodine <- subset(temp, element == "I")
        iodine_sum <- sum(iodine$total, na.rm = TRUE)
        rad_sum_dose(data.frame(total = total_sum))
        i_sum_dose(data.frame(total = iodine_sum))
        xrad_doses
      })
      
      
      ## Compliance Summary (L2) -------------------------------------------------
      output$fraction_display_l2 <- renderText({
        s <- rad_sum_dose()
        val <- if (!is.null(s$total))
          s$total
        else
          NA_real_
        paste0("Sum of total dose across modeled receptors is ",
               signif(val, 4),
               " mrem/y.")
      })
      output$fraction_display_i_l2 <- renderText({
        s <- i_sum_dose()
        val <- if (!is.null(s$total))
          s$total
        else
          NA_real_
        paste0("Sum of iodine-only dose across modeled receptors is ",
               signif(val, 4),
               " mrem/y.")
      })
      
      
      
      
      
      ## L2 Line Graph -----------------------------------------------------------
      output$graph_l2 <- renderPlotly({
        df <- as.data.frame(combined_data_l2())
        if (is.null(df)) {
          return(
            plotly::plotly_empty(type = "scatter", mode = "markers") |>
              layout(title = "No data to plot yet.")
          )
        }
        df <- df %>%
          mutate(Source.Distance = suppressWarnings(as.numeric(Source.Distance)),
                 total           = suppressWarnings(as.numeric(total))) %>%
          filter(
            !is.na(Source.Distance),
            is.finite(Source.Distance),
            !is.na(total),
            is.finite(total),
            !is.na(Radionuclide),
            nzchar(Radionuclide)
          ) %>%
          arrange(Radionuclide, Source.Distance)
        if (!nrow(df)) {
          return(
            plotly::plotly_empty(type = "scatter", mode = "markers") |>
              layout(title = "No valid dose data to plot yet.")
          )
        }
        plotly::plot_ly(
          data  = df,
          x     = ~ Source.Distance,
          y     = ~ total,
          color = ~ Radionuclide,
          type  = "scatter",
          mode  = "markers"
        ) |>
          layout(
            title = "Dose by Distance",
            xaxis = list(title = "Distance from Source (m)"),
            yaxis = list(title = "Total Dose (mrem/y)")
          )
      })
      
      
      
      
      
      
      
      # Level 4 / RADON ---------------------------------------------------------
      
      ## Radon params Modal ------------------------------------------------------
      observeEvent(input$radon_params_btn, {
        # UNUSED: directions <- c("N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW")
        showModal(
          modalDialog(
            title = "Radon Stack Parameters",
            numericInput(
              "radon_diam",
              "Diameter (m):",
              value = stack_diameter(),
              min = 0.01
            ),
            numericInput(
              "radon_vel",
              "Exit Velocity (m³/s):",
              value = exit_velocity(),
              min = 0
            ),
            numericInput(
              "radon_ht",
              "Release Height (m):",
              value = release_height(),
              min = 0
            ),
            numericInput(
              "radon_tstack",
              "Stack Temperature (°C):",
              value = stack_temp(),
              min = -100
            ),
            numericInput(
              "radon_tair",
              "Air Temperature (°C):",
              value = air_temp(),
              min = -100
            ),
            selectInput(
              "radon_vertmom",
              "Vertical Momentum?",
              choices = c("YES", "NO"),
              selected = vertical_momentum()
            ),
            numericInput(
              "radon_emit",
              "Radon Emissions (Ci/year):",
              value = emission(),
              min = 0
            ),
            numericInput(
              "radon_calm",
              "Percent Calms (%):",
              value = percent_calm(),
              min = 0,
              max = 100
            ),
            tags$hr(),
            footer = tagList(
              actionButton("radon_save", "Save", class = "btn-primary"),
              modalButton("Cancel")
            ),
            easyClose = TRUE
          )
        )
      })
      
      
      ## Met params (input file) -------------------------------------------------
      # Reactive: Parse wind data from input.level4r.txt (lines 21–36)
      wind_data_table <- reactive({
        req(input$level4_input_file)
        lines <- readLines(input$level4_input_file$datapath, warn = FALSE)
        # Wind data is lines 21–36 (R is 1-based)
        wind_lines <- lines[21:36]
        # Remove the "Dir,Freq,Speed(m/s)|" prefix if present
        wind_lines <- gsub("^Dir,Freq,Speed\\(m/s\\)\\|", "", wind_lines)
        # Split each line by "|" and then by ","
        wind_df <- do.call(rbind, lapply(wind_lines, function(line) {
          # Remove any leading/trailing whitespace and split by "|"
          parts <- strsplit(line, "\\|")[[1]]
          # The actual data is in parts[1]
          vals <- strsplit(parts[1], ",")[[1]]
          # Drop the last column (fifth value, sector-to)
          vals <- vals[1:4]
          # Clean up whitespace
          vals <- trimws(vals)
          return(vals)
        }))
        wind_df <- as.data.frame(wind_df, stringsAsFactors = FALSE)
        colnames(wind_df) <- c("Direction", "Frequency", "Speed", "Distance")
        # Convert numeric columns
        wind_df$Frequency <- as.numeric(wind_df$Frequency)
        wind_df$Speed <- as.numeric(wind_df$Speed)
        wind_df$Distance <- as.numeric(wind_df$Distance)
        wind_df
      })
      # Meterlogical Table
      output$wind_data_table <- renderDT({
        wind_data_table()
      })
      
      ## Wind Rose ---------------------------------------------------------------
      output$radon_wind_rose <- renderPlot({
        df <- wind_data_table()
        
        df <- df[complete.cases(df[, c("Frequency", "Speed", "Direction")]), ]
        
        req(!is.null(df), nrow(df) > 0)
        
        # Map directions to angles for openair
        dir_map <- c(
          "N" = 0,
          "NNE" = 22.5,
          "NE" = 45,
          "ENE" = 67.5,
          "E" = 90,
          "ESE" = 112.5,
          "SE" = 135,
          "SSE" = 157.5,
          "S" = 180,
          "SSW" = 202.5,
          "SW" = 225,
          "WSW" = 247.5,
          "W" = 270,
          "WNW" = 292.5,
          "NW" = 315,
          "NNW" = 337.5
        )
        # Add angle column
        df$angle <- dir_map[toupper(trimws(df$Direction))]
        # Expand rows by frequency (for openair windRose)
        freq_counts <- pmax(round(df$Frequency * 1000), 1)
        expanded <- df[rep(seq_len(nrow(df)), freq_counts), ]
        # openair expects columns: ws (wind speed), wd (wind direction in degrees)
        windrose_df <- data.frame(ws = expanded$Speed, wd = expanded$angle)
        windrose_df <- windrose_df[!is.na(windrose_df$ws) &
                                     !is.na(windrose_df$wd), ]
        
        openair::windRose(
          windrose_df,
          ws = "ws",
          wd = "wd",
          paddle = FALSE,
          key.title = "Wind Speed (m/s)",
          main = "Wind Rose",
          cols = c("yellow", "chartreuse3", "cadetblue", "deepskyblue")
        )
      })
      
      
      
      ## ---- Dose Calculation Function ----
      calculate_radon_allout <- function(wind_df,
                                         releaseht,
                                         exitvelocity,
                                         diameter,
                                         stacktemp,
                                         airtemp,
                                         verticalmomentum,
                                         rademis,
                                         pctcalm,
                                         nwindbins) {
        # Initialize output
        allout <- data.frame()
        maxdose <- 0.0
        maxsector <- 0
        conc <- numeric(nwindbins)
        ihdose <- numeric(nwindbins)
        
        for (nnuc in 1:nwindbins) {
          windbinname <- wind_df$Direction[nnuc]
          windbinfreq <- wind_df$Frequency[nnuc]
          windbinspeed <- wind_df$Speed[nnuc]
          windbinsrdist <- wind_df$Distance[nnuc]
          sectorname <- windbinname
          
          # --- Calculations ---
          HEFF <- releaseht
          if (releaseht > 0.0) {
            RHOAIR <- 0.0735 * (80.0 + 460.0) / (airtemp + 460.0)
            RHOAIR <- RHOAIR * 454.0 / (0.3048 ^ 3)
            RHOSTK <- RHOAIR * (airtemp + 460.0) / (stacktemp + 460.0)
            AREA <- pi * (diameter ^ 2.0) / 4.0
            STKVEL <- exitvelocity / AREA
            R <- STKVEL / windbinspeed
            FM <- (RHOSTK / RHOAIR) * STKVEL * exitvelocity / pi
            TERM1 <- 75.0 * FM / (((1.0 + 3.0 / R) * windbinspeed) ^ 2.0)
            TERM2 <- releaseht + 7.79 * sqrt(FM) / ((1.0 + 3.0 / R) * windbinspeed)
            DELHMO <- if (toupper(verticalmomentum) == "YES")
              0.93 * (TERM1 ^ (3 / 7)) * (TERM2 ^ (1 / 7))
            else
              0.0
            FB <- 9.81 * (RHOAIR - RHOSTK) * exitvelocity / (pi * RHOAIR)
            DELHBU <- 23.67 * ((FB / (windbinspeed ^ 3.0)) ^ 0.6) * (releaseht + 195.1 * FB / (windbinspeed ^
                                                                                                 3.0)) ^ 0.4
            HEFF <- releaseht + max(DELHMO, DELHBU, na.rm = TRUE)
          }
          SIGMAZ <- (0.06 * windbinsrdist) / sqrt(1.0 + 0.0015 * windbinsrdist)
          CAPSIG <- sqrt(SIGMAZ ^ 2 + 0.0 / pi)
          zz <- 0.5 * ((HEFF / CAPSIG) ^ 2)
          pp <- if (zz > 0.001)
            2.032 * exp(-zz) / (windbinsrdist * CAPSIG)
          else
            2.032 / (windbinsrdist * CAPSIG)
          denom <- (1.0 - pctcalm)
          denom <- ifelse(denom <= 0, 1e-6, denom)
          conc[nnuc] <- ((windbinfreq / denom) * rademis * pp) / windbinspeed
          if (!is.na(conc[nnuc]) &&
              conc[nnuc] < 1e-32)
            conc[nnuc] <- 1e-32
          
          # Inhalation dose
          INHFAC <- 0.0092 * 100. * 1000.
          TINIT <- 11.8
          T <- TINIT + windbinsrdist / (windbinspeed * 60.)
          EQILOUT <- 1 - 0.0479 * exp(-T / 4.39) - 2.1963 * exp(-T / 38.6) + 1.2442 * exp(-T / 28.4)
          EQILOUT <- min(EQILOUT, 1.0)
          EQILIN <- 0.35 * (1. + EQILOUT)
          EQILIN <- min(EQILIN, 0.7)
          INHDOS <- conc[nnuc] * 1E12 * (1 / 1000.) * (1 / 100.) * (8760 / 170.) * INHFAC
          FINDOOR <- 0.75
          ihdose[nnuc] <- INHDOS * (FINDOOR * EQILIN + (1.0 - FINDOOR) * EQILOUT)
          
          if (ihdose[nnuc] > maxdose) {
            maxsector <- nnuc
            maxdose <- ihdose[nnuc]
          }
          
          # Store results
          allout <- rbind(
            allout,
            data.frame(
              Stack = 1,
              Emissions_Ci_y = rademis,
              Sector = nnuc,
              Max_Sector = maxsector,
              Sector_Name = sectorname,
              Distance_m = windbinsrdist,
              Total_Dose_mrem_y = ihdose[nnuc],
              Conc_Ci_m3 = conc[nnuc],
              Wind_Frequency = windbinfreq,
              Wind_Speed_m_s = windbinspeed
            )
          )
        }
        return(allout)
      }
      
      
      radon_allout <- reactiveVal(NULL)
      
      observeEvent(input$calculate_radon, {
        req(input$level4_input_file)
        req(
          input$radon_ht,
          input$radon_vel,
          input$radon_diam,
          input$radon_tstack,
          input$radon_tair,
          input$radon_vertmom,
          input$radon_emit,
          input$radon_calm
        )
        
        wind_df <- wind_data_table()
        nwindbins <- nrow(wind_df)
        
        allout <- calculate_radon_allout(
          wind_df = wind_df,
          releaseht = input$radon_ht,
          exitvelocity = input$radon_vel,
          diameter = input$radon_diam,
          stacktemp = input$radon_tstack,
          airtemp = input$radon_tair,
          verticalmomentum = input$radon_vertmom,
          rademis = input$radon_emit,
          pctcalm = input$radon_calm,
          nwindbins = nwindbins
        )
        
        radon_allout(allout)  # Store for use in dose_results
        
        
        ## Radon Dose Output Table -------------------------------------------------
        
        
        output$radon_allout_table <- renderDT({
          datatable(allout,
                    extensions = "Buttons",
                    options = list(
                      dom = "Blfti",
                      buttons = c("copy", "csv", "excel", "pdf")
                    ))
        })
      })
      
      
      observeEvent(input$radon_save, {
        stack_diameter(input$radon_diam)
        exit_velocity(input$radon_vel)
        release_height(input$radon_ht)
        stack_temp(input$radon_tstack)
        air_temp(input$radon_tair)
        vertical_momentum(input$radon_vertmom)
        emission(input$radon_emit)
        percent_calm(input$radon_calm)
        removeModal()
      })
      
      
      ### Radon Dose Map Output ---------------------------------------------------
      # GET lat/long
      site_coords <- reactive({
        # Try user-entered lat/lon first
        lat_in <- suppressWarnings(as.numeric(input$Latitude))
        lon_in <- suppressWarnings(as.numeric(input$Longitude))
        if (!is.na(lat_in) &&
            !is.na(lon_in) &&
            nzchar(input$Latitude) && nzchar(input$Longitude)) {
          return(list(lat = lat_in, lon = lon_in))
        }
        # If not available, try geocoded address
        some_addresses <- tibble::tibble(name = "company", addr = input$address_input)
        geocode_result <- tryCatch(
          tidygeocoder::geocode(
            some_addresses,
            address = addr,
            method = "osm",
            lat = latitude,
            long = longitude
          ),
          error = function(e)
            NULL
        )
        if (!is.null(geocode_result) &&
            !is.na(geocode_result$latitude) &&
            !is.na(geocode_result$longitude)) {
          return(list(
            lat = as.numeric(geocode_result$latitude),
            lon = as.numeric(geocode_result$longitude)
          ))
        }
        # If neither is available, return NA
        return(list(lat = NA_real_, lon = NA_real_))
      })
      
      # get data
      dose_results <- reactive({
        allout <- radon_allout()
        req(allout)
        # Use the wind data for directions/distances (should match allout order)
        wind_df <- wind_data_table()
        dirs <- wind_df$Direction
        dist <- wind_df$Distance
        dose <- allout$Total_Dose_mrem_y
        data.frame(dirs = dirs,
                   dist = dist,
                   dose = dose)
      })
      
      #make map
      output$radon_dose_map <- renderLeaflet({
        req(dose_results())
        coords <- site_coords()
        req(!is.na(coords$lat), !is.na(coords$lon))
        
        dir_map <- c(
          "N" = 0,
          "NNE" = 22.5,
          "NE" = 45,
          "ENE" = 67.5,
          "E" = 90,
          "ESE" = 112.5,
          "SE" = 135,
          "SSE" = 157.5,
          "S" = 180,
          "SSW" = 202.5,
          "SW" = 225,
          "WSW" = 247.5,
          "W" = 270,
          "WNW" = 292.5,
          "NW" = 315,
          "NNW" = 337.5
        )
        
        dr <- dose_results() %>%
          mutate(
            angle = dir_map[dirs],
            x = dist * sin(angle * pi / 180),
            y = dist * cos(angle * pi / 180)
          )
        
        g <- st_as_sf(
          dr,
          coords = c("x", "y"),
          crs = paste0(
            "+proj=lcc +lat_1=",
            coords$lat,
            " +lat_2=",
            coords$lat,
            " +lat_0=",
            coords$lat,
            " +lon_0=",
            coords$lon,
            " +datum=NAD83"
          )
        )
        p <- st_transform(g, "+proj=longlat +datum=NAD83")
        dr <- cbind(dr, st_coordinates(p))
        names(dr)[names(dr) == "X"] <- "lon"
        names(dr)[names(dr) == "Y"] <- "lat"
        
        x_out <- seq(-1000, 1000, 5)
        y_out <- x_out
        interp <- with(dr,
                       akima::interp(
                         x,
                         y,
                         dose,
                         xo = x_out,
                         yo = y_out,
                         linear = TRUE,
                         extrap = TRUE
                       ))
        interp$z[interp$z > 10] <- 9.99
        
        grid_pts <- expand.grid(x = interp$x, y = interp$y)
        g_grid <- st_as_sf(
          grid_pts,
          coords = c("x", "y"),
          crs = paste0(
            "+proj=lcc +lat_1=",
            coords$lat,
            " +lat_2=",
            coords$lat,
            " +lat_0=",
            coords$lat,
            " +lon_0=",
            coords$lon,
            " +datum=NAD83"
          )
        )
        p_grid <- st_transform(g_grid, "+proj=longlat +datum=NAD83")
        xy <- st_coordinates(p_grid)
        
        r <- raster(
          ncols = length(interp$x),
          nrows = length(interp$y),
          xmn = min(xy[, 1]),
          xmx = max(xy[, 1]),
          ymn = min(xy[, 2]),
          ymx = max(xy[, 2])
        )
        r1 <- rasterize(cbind(xy[, 1], xy[, 2]), r)
        values(r1) <- as.vector(interp$z)
        
        obins <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
        my.colors <- colorRampPalette(
          c(
            grey(.8),
            "mediumpurple",
            "darkorchid4",
            "#002FFF",
            "lightblue",
            "green",
            "yellow",
            "orange",
            "red",
            "brown"
          )
        )
        pal1 <- colorBin(my.colors(length(obins)), bins = obins, na.color = "transparent")
        
        leaflet() %>%
          addProviderTiles("Esri.WorldStreetMap", group = "Street") %>%
          addProviderTiles("Esri.WorldImagery", group = "Imagery") %>%
          addProviderTiles("OpenTopoMap", group = "Topo") %>%
          addCircleMarkers(
            lng = dr$lon,
            lat = dr$lat,
            color = "black",
            fillColor = "yellow",
            radius = 7,
            fillOpacity = 0.7,
            group = "MEI"
          ) %>%
          addRasterImage(
            flip(r1, "y"),
            colors = pal1,
            opacity = 0.6,
            group = "Interpolated"
          ) %>%
          addLayersControl(
            baseGroups = c("Street", "Imagery", "Topo"),
            overlayGroups = c("MEI", "Interpolated"),
            options = layersControlOptions(collapsed = FALSE)
          )
      })
      
    }
    
    # ---- Viewer / run ----
    vwr <- dialogViewer("COMPLY v2.0", width = 3200, height = 7400)
    runGadget(shinyApp(ui = ui, server = server), viewer = vwr)
    # shinyApp(ui, server)  # alternative