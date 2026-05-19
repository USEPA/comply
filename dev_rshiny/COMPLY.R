
setwd("C:/Users/MMACDO01/OneDrive - Environmental Protection Agency (EPA)/Github/comply/dev_rshiny")



#######################

#    C O M P L Y

#   Version:05/19/26
#   macdonald.megan@epa.gov

#######################

# This is a Shiny web application. You can run the application by clicking the 'Run App' button above.

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

# Helper to ensure column names are syntactic
make_syntactic <- function(df) { names(df) <- make.names(names(df)); df }

# Used in Level 1: load possession table
rads_src <- read.csv(file="comply171_original_possession_data.csv", skip=0, header=TRUE)
long_rads <- melt(rads_src, id.vars = "Nuclide", variable.name = "Form", value.name = "Allowance")
long_rads$Radionuclide <- gsub(" ", "", long_rads$Nuclide)
long_rads$Type <- gsub(" ", "", long_rads$Form)
long_rads <- long_rads[, c(4, 5, 3)]
# long_rads columns: Radionuclide, Type, Allowance

# Used in Level 2: dose coefficients/lookup
doses_info <- read.csv(file="comply171_original_input_data.csv", skip=0, header=TRUE)
doses_info$Radionuclide <- trimws(doses_info$Nuclide)
doses_info <- make_syntactic(doses_info)  # ensure consistent names for joins/mutate

# align expected dose coefficient column names if needed
cn <- names(doses_info)
if ("WB.dose.y" %in% cn && !"WBdose.y" %in% cn) {
  names(doses_info)[names(doses_info) == "WB.dose.y"] <- "WBdose.y"
}
if ("Veg.dose" %in% cn && !"VegDose" %in% cn) {
  names(doses_info)[names(doses_info) == "Veg.dose"] <- "VegDose"
}
if ("CowMilk.dose" %in% cn && !"CowMilkDose" %in% cn) {
  names(doses_info)[names(doses_info) == "CowMilk.dose"] <- "CowMilkDose"
}

# Blank templates as zero-row data frames
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
doses <- data.frame(
  Stack.ID = as.character(NA),
  Radionuclide = as.character(NA),
  Type = as.character(NA),
  Source.Distance = as.numeric(NA),
  airconc = as.numeric(NA),
  dose = as.numeric(NA),
  vdose = as.numeric(NA),
  cdose = as.numeric(NA),
  total = as.numeric(NA),
  check.names = FALSE
)

radon <- data.frame(
  Stack.ID = as.character(NA),
  Emissions = as.numeric(NA),
  Sector = as.character(NA),
  Max.Sector = as.numeric(NA),
  Sector.Name = as.numeric(NA),
  Distance = as.numeric(NA),
  Total.Dose = as.numeric(NA),
  Conc = as.numeric(NA),
  Wind.Frequency = as.numeric(NA),
  Wind.Speed = as.numeric(NA),
  
  check.names = FALSE
)


xsum <- data.frame(Fraction = as.numeric(NA))

#### UI
ui <- page_sidebar(
  title = "COMPLY v2.0",
  sidebar = sidebar(
    h6("Enter Facility Details:"),
    textInput("CompanyName", "Company Name:", ""),
    textInput("CompanyContact", "Company Contact Name:", ""),
    textInput("CompanyContactNum", "Company Contact Email:", ""),
    textInput("FacilityName", "Facility Name:", ""),
    h6("Enter Coordinates or Address:"),
    fluidRow(
      column(6, textInput("Latitude", "Latitude:", "")),
      column(6, textInput("Longitude", "Longitude:", ""))
    ),
    textInput("address_input", "Facility Address:", "109 T.W. Alexander Drive, Durham, NC, 27709"),
    actionButton("geocode_button", "Plot Location"),
    h6("Coordinates:"),
    textOutput("coords")
  ),
  
  navset_tab(
    id = "level_tabs",
    selected = "level 1",
    
    nav_panel("level 1",
              card(
                full_screen = TRUE,
                height = "500px",
                navset_card_tab(
                  nav_panel("Rads",
                            tags$div(class = "d-flex gap-2 mb-2",
                                     actionButton("add_rad_button", "Add Radionuclide", class = "btn-sm"),
                                     actionButton("delete_rad_button", "Delete Selected Row", class = "btn-sm")
                            ),
                            tags$head(tags$style(HTML("
              .dataTable td, .dataTable th { border: 1px solid black; }
              .dt-button {
                background-color: #f0f0f0 !important;
                color: black !important;
                border: 1px solid #000500 !important;
              }"))),
                            DTOutput("rad_table")
                  ),
                  nav_panel("Compliance Summary",
                            textOutput("fraction_display"),
                            textOutput("fraction_display_i")
                  ),
                  nav_panel("Export Data",
                            downloadButton("downloadReport", "Generate HTML Report", class = "btn-sm", width = '20%')
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
    ###### LEVEL 2
    nav_panel("level 2/3",
              # Map and Plot
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
                  nav_panel("Stacks",
                            tags$div(
                              class = "d-flex gap-2 mb-1",
                              actionButton("add_data_button", "Add Stack", class = "btn-sm"),
                              actionButton("delete_stack_button", "Delete Selected Row", class = "btn-sm")
                            ),
                            # Footnote under the buttons
                            tags$small(
                              class = "text-muted",
                              "Note: For Level 2, use same distance values for Source Distance, Veg Distance, and Cow Distance."
                            ),
                            DTOutput("stack_table")
                  ),
                  nav_panel("Emissions",
                            # Buttons for Level 2 emissions table
                            tags$div(class = "d-flex gap-2 mb-2",
                                     actionButton("add_rad_button_l2", "Add Radionuclide", class = "btn-sm"),
                                     actionButton("delete_rad_button_l2", "Delete Selected Row", class = "btn-sm")
                            ),
                            tags$head(
                              tags$style(HTML("
      .dataTable td, .dataTable th { border: 1px solid black; }
      .dt-button {
        background-color: #f0f0f0 !important;
        color: black !important;
        border: 1px solid #000500 !important;
      }
    "))
                            ),
                            DTOutput("rad_table_l2")
                  ),
                  nav_panel("Diagnostic Results",
                            tags$head(tags$style(HTML("
              .dataTable td, .dataTable th { border: 1px solid black; }
              .dt-button {
                background-color: #f0f0f0 !important;
                color: black !important;
                border: 1px solid #000500 !important;
              }"))),
                            DTOutput("dose_table")
                  ),
                  nav_panel("Compliance Summary",
                            textOutput("fraction_display_l2"),
                            textOutput("fraction_display_i_l2")
                  ),
                  nav_panel("Export Data",
                            downloadButton("downloadCSVData", "Download .CSV file", class = "btn-sm", width = '20%'),
                            downloadButton("downloadReport", "Generate HTML Report", class = "btn-sm", width = '20%')
                  )
                )
              )
    ),
    
    # ---- RADON TAB: Button + summary output ----
    nav_panel("Radon",
              h4("Radon content goes here."),
              tags$div(class = "d-flex gap-2 mb-2",
                       actionButton("radon_params_btn", "Set Radon Parameters", class = "btn-sm")
              ),
              tags$small(class = "text-muted",
                         "Enter and save Radon stack parameters here for use in downstream calculations."
              ),
              textOutput("radon_params_summary")
    )
  )
)

server <- function(input, output, session) {
  
  # ========== STACKS ==========
  stack_data <- reactiveVal(stacks)
  
  output$stack_table <- renderDT({
    datatable(
      stack_data(), 
      class = 'cell-border stripe',
      extensions = 'Buttons',
      selection = 'single',            
      options = list( 
        ordering = FALSE,
        searching = FALSE,
        info = FALSE,
        lengthChange = FALSE,
        dom = 'Blfti',
        columnDefs = list(list(className = 'dt-center', targets = "_all")),
        rownames = FALSE,
        pageLength = -1,
        paging = FALSE,
        filter = "top",
        buttons = c('copy', 'csv', 'excel', 'pdf'),
        exportOptions = list(columns = ':visible')
      )
    )
  })
  
  observeEvent(input$add_data_button, {
    showModal(
      modalDialog(
        title = "Enter New Stack Data",
        textInput("input_ID", "Stack ID:", ""),
        numericInput("input_RH", "Release Height:", value = NULL, min = 0),
        numericInput("input_BH", "Building Height:", value = NULL, min = 0),
        numericInput("input_BW", "Building Width:", value = NULL, min = 0),
        numericInput("input_SD", "Source Distance:", value = NULL, min = 0),
        numericInput("input_VD", "Veg Distance:", value = NULL, min = 0),
        numericInput("input_CD", "Cow Distance:", value = NULL, min = 0),
        numericInput("input_WS", "Mean WS:", value = NULL, min = 0),
        
        # Footnote inside the modal
        tags$small(
          class = "text-muted",
          "Note: For Level 2, use same distance values for Source Distance, Veg Distance, and Cow Distance."
        ),
        
        footer = tagList(
          actionButton("save_btn", "Save"),
          modalButton("Cancel")
        ),
        easyClose = TRUE
      )
    )
  })
  
  observeEvent(input$save_btn, {
    current_data <- stack_data()
    
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
    
    updated_data <- bind_rows(current_data, new_row)
    cleaned_data <- updated_data[rowSums(is.na(updated_data)) < ncol(updated_data), ]
    stack_data(cleaned_data)
    removeModal()
  })
  
  observeEvent(input$delete_stack_button, {
    sel <- input$stack_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      showModal(modalDialog(
        title = "No row selected",
        "Click a row in the Stacks table to select it, then press 'Delete Selected Row'.",
        easyClose = TRUE
      ))
    } else {
      current <- stack_data()
      # Remove the selected row
      updated <- current[-sel, , drop = FALSE]
      # Drop any all-NA template rows
      cleaned <- updated[rowSums(is.na(updated)) < ncol(updated), , drop = FALSE]
      stack_data(cleaned)
    }
  })
  
  
  # ========== RADS (POSSESSION) ==========
  rad_data <- reactiveVal(rads)
  rad_sum_possession <- reactiveVal(xsum)
  i_sum_possession <- reactiveVal(xsum)
  
  output$rad_table <- renderDT({
    datatable(
      rad_data(),
      class = 'cell-border stripe',
      extensions = 'Buttons',
      selection = 'single',
      options = list(
        ordering = FALSE, searching = FALSE, info = FALSE, lengthChange = FALSE,
        dom = 'Blfti',
        columnDefs = list(list(className = 'dt-center', targets = "_all")),
        rownames = FALSE,
        pageLength = -1, paging = FALSE, filter = "top",
        buttons = c('copy', 'csv', 'excel', 'pdf'),
        exportOptions = list(columns = ':visible')
      )
    )
  })
  
  observeEvent(input$add_rad_button, {
    showModal(
      modalDialog(
        title = "Enter New Radionuclide Data",
        textInput("input_STACK", "Stack ID:", ""),
        selectInput("input_RN", "Radionuclide:", choices = unique(long_rads$Radionuclide), selected = NULL),
        selectInput("input_TYPE", "Type:", choices = unique(long_rads$Type), selected = NULL),
        numericInput("input_EMISSION", "Emission:", value = NULL, min = 0),
        footer = tagList(
          actionButton("save_rad_btn", "Save"),
          modalButton("Cancel")
        )
      )
    )
  })
  
  observeEvent(input$delete_rad_button, {
    sel <- input$rad_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      showModal(modalDialog(
        title = "No row selected",
        "Click a row in the table to select it, then press 'Delete Selected Row'.",
        easyClose = TRUE
      ))
    } else {
      current <- rad_data()
      updated <- current[-sel, , drop = FALSE]
      cleaned <- updated[rowSums(is.na(updated)) < ncol(updated), , drop = FALSE]
      rad_data(cleaned)
      
      if (nrow(cleaned) == 0) {
        rad_sum_possession(data.frame(Fraction = 0))
        i_sum_possession(data.frame(Fraction = 0))
      } else {
        temp <- cleaned
        temp$element <- substr(temp$Radionuclide, 1, 1)
        iodine_data <- subset(temp, element == "I")
        rad_sum_possession(data.frame(Fraction = sum(temp$Fraction)))
        i_sum_possession(data.frame(Fraction = sum(iodine_data$Fraction)))
      }
    }
  })
  
  # CHANGED: safer allowance lookup and fraction calc
  observeEvent(input$save_rad_btn, {
    current_data <- rad_data()
    
    allow_row <- subset(long_rads, Type == input$input_TYPE & Radionuclide == input$input_RN)
    if (nrow(allow_row) == 0 || is.na(allow_row$Allowance[1])) {
      showModal(modalDialog(
        title = "No matching allowance",
        "No allowance found for that Radionuclide and Type. Please check your selection.",
        easyClose = TRUE
      ))
      return(invisible(NULL))
    }
    allow_val <- allow_row$Allowance[1]
    
    new_row <- data.frame(
      Stack.ID = input$input_STACK,
      Radionuclide = input$input_RN,
      Type = input$input_TYPE,
      Emission = input$input_EMISSION,
      Allowance = allow_val,
      Fraction = round((input$input_EMISSION / allow_val), 5),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    updated_data <- bind_rows(current_data, new_row)
    cleaned_data <- updated_data[rowSums(is.na(updated_data)) < ncol(updated_data), ]
    temp <- cleaned_data
    temp$element <- substr(temp$Radionuclide, 1, 1)
    iodine_data <- subset(temp, element == "I")
    
    rad_data(cleaned_data)
    rad_sum_possession(data.frame(Fraction = sum(cleaned_data$Fraction, na.rm = TRUE)))
    i_sum_possession(data.frame(Fraction = sum(iodine_data$Fraction, na.rm = TRUE)))
    removeModal()
  })
  
  # CHANGED: Level 2 "Add" button with safe choices and guard for no stacks
  observeEvent(input$add_rad_button_l2, {
    sd <- isolate(stack_data())
    stack_choices <- character(0)
    if (!is.null(sd) && "Stack.ID" %in% names(sd)) {
      ids <- sd$`Stack.ID`
      stack_choices <- unique(ids[!is.na(ids) && nzchar(ids)])
    }
    
    if (length(stack_choices) == 0) {
      showModal(modalDialog(
        title = "No stacks found",
        "Please add at least one Stack in the Stacks tab before adding emissions.",
        easyClose = TRUE
      ))
      return(invisible(NULL))
    }
    
    showModal(
      modalDialog(
        title = "Enter New Radionuclide Data",
        selectInput("input_STACK", "Stack ID:", choices = stack_choices, selected = NULL),
        selectInput("input_RN", "Radionuclide:", choices = unique(long_rads$Radionuclide), selected = NULL),
        selectInput("input_TYPE", "Type:", choices = unique(long_rads$Type), selected = NULL),
        numericInput("input_EMISSION", "Emission:", value = NULL, min = 0),
        footer = tagList(
          actionButton("save_rad_btn", "Save"),
          modalButton("Cancel")
        )
      )
    )
  })
  
  observeEvent(input$delete_rad_button_l2, {
    sel <- input$rad_table_l2_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      showModal(modalDialog(
        title = "No row selected",
        "Click a row in the Emissions table to select it, then press 'Delete Selected Row'.",
        easyClose = TRUE
      ))
    } else {
      current <- rad_data()
      updated <- current[-sel, , drop = FALSE]
      cleaned <- updated[rowSums(is.na(updated)) < ncol(updated), , drop = FALSE]
      rad_data(cleaned)
      
      # Optional: update possession summaries if you show them in Level 1
      if (exists("rad_sum_possession") && exists("i_sum_possession")) {
        if (nrow(cleaned) == 0) {
          rad_sum_possession(data.frame(Fraction = 0))
          i_sum_possession(data.frame(Fraction = 0))
        } else {
          temp <- cleaned
          temp$element <- substr(temp$Radionuclide, 1, 1)
          iodine_data <- subset(temp, element == "I")
          rad_sum_possession(data.frame(Fraction = sum(temp$Fraction, na.rm = TRUE)))
          i_sum_possession(data.frame(Fraction = sum(iodine_data$Fraction, na.rm = TRUE)))
        }
      }
    }
  })
  
  
  # Level 2 emissions table 
  output$rad_table_l2 <- renderDT({
    datatable(
      rad_data(),
      class = 'cell-border stripe',
      extensions = 'Buttons',
      selection = 'single',
      options = list(
        ordering = FALSE,
        searching = FALSE,
        info = FALSE,
        lengthChange = FALSE,
        dom = 'Blfti',
        columnDefs = list(list(className = 'dt-center', targets = "_all")),
        rownames = FALSE,
        pageLength = -1,
        paging = FALSE,
        filter = "top",
        buttons = c('copy', 'csv', 'excel', 'pdf'),
        exportOptions = list(columns = ':visible')
      )
    )
  })
  
  # ========== BUILD COMPLETE DATA FRAME (Level 1 summary) ==========
  output$sums <- renderDT({
    datatable(rad_sum_possession(), options = list(dom = 't'))
  })
  
  combined_data <- reactive({
    df <- rad_data()
    req(nrow(df) > 0)
    cleaned <- df[rowSums(is.na(df[, -1, drop = FALSE])) < ncol(df[, -1, drop = FALSE]), ]
    sum(cleaned$Fraction)
  })
  
  # ========== LEVEL 2: Dose-by-distance calculations ==========
  rad_sum_dose <- reactiveVal(data.frame(total = 0))
  i_sum_dose   <- reactiveVal(data.frame(total = 0))
  
  combined_data_l2 <- reactive({
    s <- stack_data() %>% make_syntactic()
    r <- rad_data()   %>% make_syntactic()
    
    # Remove fully NA rows
    s <- s[rowSums(is.na(s)) < ncol(s), , drop = FALSE]
    r <- r[rowSums(is.na(r)) < ncol(r), , drop = FALSE]
    
    # Require keys to exist and non-empty data
    req("Stack.ID" %in% names(s), "Stack.ID" %in% names(r))
    merged <- merge(s, r, by = "Stack.ID")
    req(nrow(merged) > 0)
    
    # Expand to distances (5..2000 by 5) per row without uncount()
    n_rows <- nrow(merged)
    df <- merged %>%
      dplyr::slice(rep(dplyr::row_number(), each = 400)) %>%
      dplyr::mutate(Source.Distance = rep(seq(5, 2000, 5), times = n_rows)) %>%
      dplyr::mutate(
        area  = Building.Width * Building.Height,
        xcrit = 2.5 * sqrt(area),
        xarea = pmin(Building.Width, Building.Height),
        sigmaz = (0.06 * Source.Distance) / sqrt(1.0 + (0.0015 * Source.Distance)),
        capsig = sqrt((sigmaz^2) + (area / 3.14)),
        pp     = 2.032 / (Source.Distance * capsig),
        conc   = ifelse(
          Source.Distance < xcrit,
          (0.25 * as.numeric(Emission)) / (3.14 * Mean.WS * 1.0 * xarea),
          (0.25 * as.numeric(Emission) * pp) / Mean.WS
        )
      )
    
    # Join dose coefficients
    df <- df %>% dplyr::left_join(doses_info, by = "Radionuclide")
    
    # Calculate doses
    df <- df %>%
      dplyr::mutate(
        airconc = conc,
        dose  = WBdose.y    * conc,   # columns must exist in doses_info (aligned earlier)
        vdose = VegDose     * conc,
        cdose = CowMilkDose * conc,
        total = dose + vdose + cdose
      )
    
    df
  })
  
  # Dose table (full)
  output$dose_table <- renderDT({
    req(combined_data_l2())
    datatable(combined_data_l2(),
              class = 'cell-border stripe',
              extensions = 'Buttons',
              options = list(
                ordering = FALSE, searching = FALSE, info = FALSE, lengthChange = FALSE,
                dom = 'Blfti',
                columnDefs = list(list(className = 'dt-center', targets = "_all")),
                rownames = FALSE,
                pageLength = -1, paging = FALSE, filter = "top",
                buttons = c('copy', 'csv', 'excel', 'pdf'),
                exportOptions = list(columns = ':visible')
              )
    )
  })
  
  # Cut-down dose-specific version for display/export
  dose_data <- reactive({
    full_data_l2 <- combined_data_l2()
    full_data_l2 %>%
      select(Stack.ID, Radionuclide, Source.Distance, airconc, dose, vdose, cdose, total)
  })
  
  # Specific receptors (MEI / veg / cow distances) for compliance
  distance_data <- reactive({
    s <- stack_data() %>% make_syntactic()
    r <- rad_data() %>% make_syntactic()
    
    s <- s[rowSums(is.na(s)) < ncol(s), , drop = FALSE]
    r <- r[rowSums(is.na(r)) < ncol(r), , drop = FALSE]
    
    req("Stack.ID" %in% names(s), "Stack.ID" %in% names(r))
    xrads <- merge(s, r, by = "Stack.ID")
    req(nrow(xrads) > 0)
    
    xrads$area  <- xrads$Building.Width * xrads$Building.Height
    xrads$xcrit <- 2.5 * sqrt(xrads$area)
    xrads$xarea <- pmin(xrads$Building.Width, xrads$Building.Height)
    
    # MEI distance
    xrads$sigmaz  <- (0.06 * xrads$Source.Distance)/sqrt(1.0 + (0.0015 * xrads$Source.Distance))
    xrads$capsig  <- sqrt((xrads$sigmaz^2) + (xrads$area/3.14))
    xrads$pp      <- 2.032 / (xrads$Source.Distance * xrads$capsig)
    xrads$conc.mei <- ifelse(
      xrads$Source.Distance < xrads$xcrit,
      (0.25 * as.numeric(xrads$Emission))/(3.14 * xrads$Mean.WS * 1.0 * xrads$xarea),
      (0.25 * as.numeric(xrads$Emission) * xrads$pp) / xrads$Mean.WS
    )
    
    # Veg receptor
    xrads$sigmaz  <- (0.06 * xrads$Veg.Distance)/sqrt(1.0 + (0.0015 * xrads$Veg.Distance))
    xrads$capsig  <- sqrt((xrads$sigmaz^2) + (xrads$area/3.14))
    xrads$pp      <- 2.032 / (xrads$Veg.Distance * xrads$capsig)
    xrads$conc.veg <- ifelse(
      xrads$Veg.Distance < xrads$xcrit,
      (0.25 * as.numeric(xrads$Emission))/(3.14 * xrads$Mean.WS * 1.0 * xrads$xarea),
      (0.25 * as.numeric(xrads$Emission) * xrads$pp) / xrads$Mean.WS
    )
    
    # Cow receptor
    xrads$sigmaz  <- (0.06 * xrads$Cow.Distance)/sqrt(1.0 + (0.0015 * xrads$Cow.Distance))
    xrads$capsig  <- sqrt((xrads$sigmaz^2) + (xrads$area/3.14))
    xrads$pp      <- 2.032 / (xrads$Cow.Distance * xrads$capsig)
    xrads$conc.cow <- ifelse(
      xrads$Cow.Distance < xrads$xcrit,
      (0.25 * as.numeric(xrads$Emission))/(3.14 * xrads$Mean.WS * 1.0 * xrads$xarea),
      (0.25 * as.numeric(xrads$Emission) * xrads$pp) / xrads$Mean.WS
    )
    
    xrad_doses <- left_join(xrads, doses_info, by = "Radionuclide")
    
    xrad_doses$airconc <- xrad_doses$conc.mei
    xrad_doses$dose    <- xrad_doses$WBdose.y   * xrad_doses$conc.mei
    xrad_doses$vdose   <- xrad_doses$VegDose    * xrad_doses$conc.veg
    xrad_doses$cdose   <- xrad_doses$CowMilkDose * xrad_doses$conc.cow
    xrad_doses$total   <- xrad_doses$dose + xrad_doses$vdose + xrad_doses$cdose
    
    # Update Level 2 summary reactives
    total_sum <- sum(xrad_doses$total, na.rm = TRUE)
    temp <- xrad_doses
    temp$element <- substr(temp$Radionuclide, 1, 1)
    iodine_data <- subset(temp, element == "I")
    iodine_sum <- sum(iodine_data$total, na.rm = TRUE)
    
    rad_sum_dose(data.frame(total = total_sum))
    i_sum_dose(data.frame(total = iodine_sum))
    
    xrad_doses
  })
  
  compliance_data <- reactive({
    full_data <- distance_data()
    full_data %>% select(Stack.ID, Radionuclide, Source.Distance, airconc, dose, vdose, cdose, total)
  })
  
  # ========== COMPLIANCE SUMMARIES ==========
  output$fraction_display <- renderText({
    selected_row <- rad_sum_possession()
    selected_age <- round(selected_row$Fraction, 2)
    paste0("Stack possession fraction is ", selected_age,
           " times the allowable amount given in the possession table.")
  })
  
  output$fraction_display_i <- renderText({
    selected_row <- i_sum_possession()
    selected_age <- round(selected_row$Fraction, 2)
    paste0("Total iodine possession fraction is ", selected_age, ".")
  })
  
  output$fraction_display_l2 <- renderText({
    s <- rad_sum_dose()
    val <- if (!is.null(s$total)) s$total else NA_real_
    paste0("Sum of total dose across modeled receptors is ", signif(val, 4), " mrem/y.")
  })
  
  output$fraction_display_i_l2 <- renderText({
    s <- i_sum_dose()
    val <- if (!is.null(s$total)) s$total else NA_real_
    paste0("Sum of iodine-only dose across modeled receptors is ", signif(val, 4), " mrem/y.")
  })
  
  # ========== PLOTS ==========
  output$graph <- renderPlotly({
    df <- as.data.frame(rad_data())
    
    # keep only valid rows
    df <- df %>%
      dplyr::filter(
        !is.na(Fraction),
        is.finite(Fraction),
        !is.na(Radionuclide),
        nzchar(Radionuclide)
      )
    
    if (nrow(df) == 0) {
      return(plotly::plotly_empty(type = "scatter", mode = "markers") |>
               layout(title = "No valid fraction data to plot yet."))
    }
    
    # aggregate by radionuclide (if multiple rows per RN exist)
    agg <- df %>%
      dplyr::group_by(Radionuclide) %>%
      dplyr::summarise(Fraction = sum(Fraction, na.rm = TRUE), .groups = "drop")
    
    plotly::plot_ly(
      data = agg,
      x = ~Radionuclide,
      y = ~Fraction,
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
  
  # In the report section
  report_plot <- reactive({
    df <- as.data.frame(rad_data())
    
    df <- df %>%
      dplyr::filter(
        !is.na(Fraction),
        is.finite(Fraction),
        !is.na(Radionuclide),
        nzchar(Radionuclide)
      )
    
    if (nrow(df) == 0) {
      return(plotly::plotly_empty(type = "scatter", mode = "markers") |>
               layout(title = "Possession Fraction Plot"))
    }
    
    agg <- df %>%
      dplyr::group_by(Radionuclide) %>%
      dplyr::summarise(Fraction = sum(Fraction, na.rm = TRUE), .groups = "drop")
    
    plotly::plot_ly(
      data = agg,
      x = ~Radionuclide,
      y = ~Fraction,
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
  
  output$graph_l2 <- renderPlotly({
    # Get the computed L2 data frame
    df <- as.data.frame(combined_data_l2())
    
    # If the reactive returns NULL early, bail out
    if (is.null(df)) {
      return(plotly::plotly_empty(type = "scatter", mode = "markers") |>
               layout(title = "No data to plot yet."))
    }
    
    # Coerce types (defensive) and keep only valid rows
    df <- df %>%
      dplyr::mutate(
        Source.Distance = suppressWarnings(as.numeric(Source.Distance)),
        total           = suppressWarnings(as.numeric(total))
      ) %>%
      dplyr::filter(
        !is.na(Source.Distance), is.finite(Source.Distance),
        !is.na(total),           is.finite(total),
        !is.na(Radionuclide),    nzchar(Radionuclide)
      ) %>%
      dplyr::arrange(Radionuclide, Source.Distance)
    
    if (nrow(df) == 0) {
      return(plotly::plotly_empty(type = "scatter", mode = "markers") |>
               layout(title = "No valid dose data to plot yet."))
    }
    
    plotly::plot_ly(
      data  = df,
      x     = ~Source.Distance,
      y     = ~total,
      color = ~Radionuclide,
      type  = "scatter",
      mode  = "markers"      # switch to "lines+markers" if you prefer lines
    ) |>
      layout(
        title = "Dose by Distance",
        xaxis = list(title = "Distance from Source (m)"),
        yaxis = list(title = "Total Dose (mrem/y)")
      )
  })
  
  # ========== REPORT DATA ==========
  report_meta <- reactive({
    res <- tryCatch(geocoded_address(), error = function(e) NULL)
    lat_in <- suppressWarnings(as.numeric(input$Latitude))
    lon_in <- suppressWarnings(as.numeric(input$Longitude))
    lat <- if (!is.null(res) && !is.null(res$latitude)) as.numeric(res$latitude) else lat_in
    lon <- if (!is.null(res) && !is.null(res$longitude)) as.numeric(res$longitude) else lon_in
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
  
  report_data <- reactive({ rad_data() })
  report_plot <- reactive({
    p <- ggplot(rad_data(), aes(x = 1, y = Fraction, fill = Radionuclide)) +
      geom_col(alpha = 0.6) +
      labs(title = "Possession Fraction Plot", color = "Source Data") +
      theme_classic() + xlab(" ") + ylab("Fraction") +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(fill = "white", color = NA)
      ) +
      labs(color = "Radionuclides")
    ggplotly(p)
  })
  
  output$HTML_table <- renderDT({ report_data() })
  output$HTML_plot <- renderPlotly({ report_plot() })
  
  output$downloadReport <- downloadHandler(
    filename = function() { paste0("report_", Sys.Date(), ".html") },
    content = function(file) {
      temp_rmd <- file.path(tempdir(), "COMPLY_Output.Rmd")
      file.copy("COMPLY_Output.Rmd", temp_rmd, overwrite = TRUE)
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
  
  # ========== MAPPING ==========
  geocoded_address <- eventReactive(input$geocode_button, {
    if (input$Latitude == "") {
      full_address <- input$address_input
      some_addresses <- tibble::tribble(~name, ~addr, "company", full_address)
      geocode_result <- some_addresses %>%
        geocode(addr, method = 'osm', lat = latitude, long = longitude)
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
      name <- "company"; addr <- ""
      geocode_result <- tibble::as_tibble(cbind(name, addr, as.numeric(latitude), as.numeric(longitude)))
      names(geocode_result) <- c("name","addr","latitude","longitude")
      geocode_result
    }
  })
  
  output$map <- renderLeaflet({
    leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap) %>%
      setView(lng = -98.5833, lat = 39.8333, zoom = 5)
  })
  output$map_l2 <- renderLeaflet({
    leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap) %>%
      setView(lng = -98.5833, lat = 39.8333, zoom = 5)
  })
  
  # Update both maps when a new address is geocoded
  observe({
    result <- geocoded_address()
    if (is.null(result)) return()
    result$latitude  <- as.numeric(result$latitude)
    result$longitude <- as.numeric(result$longitude)
    for (id in c("map", "map_l2")) {
      leafletProxy(id) %>%
        clearMarkers() %>%
        addMarkers(lng = result$longitude, lat = result$latitude, popup = result$addr) %>%
        setView(lng = result$longitude, lat = result$latitude, zoom = 10)
    }
  })
  
  output$coords <- renderText({
    result <- geocoded_address()
    if (is.null(result)) return("Enter an address and click 'Plot Location'.")
    result$latitude <- as.numeric(result$latitude)
    result$longitude <- as.numeric(result$longitude)
    paste("Latitude:", round(result$latitude, 6), ", Longitude:", round(result$longitude, 6))
  })
  
  # ---- RADON PARAMETERS (stored as numeric reactive values) ----
  # These can be called later as radon_exit_velocity() and radon_diameter()
  radon_exit_velocity <- reactiveVal(NA_real_)
  radon_diameter      <- reactiveVal(NA_real_)
  
  # Open the modal when user clicks the button on the Radon tab
  observeEvent(input$radon_params_btn, {
    showModal(
      modalDialog(
        title = "Radon Stack Parameters",
        # Use the current saved values (if any) as defaults
        numericInput(
          "radon_exit_velocity_input",
          "Exit velocity (m³/s):",
          value = radon_exit_velocity(),
          min = 0,
          step = 0.01
        ),
        numericInput(
          "radon_diameter_input",
          "Diameter (m):",
          value = radon_diameter(),
          min = 0.0001,   # positive diameter
          step = 0.01
        ),
        footer = tagList(
          actionButton("save_radon_btn", "Save"),
          modalButton("Cancel")
        ),
        easyClose = TRUE
      )
    )
  })
  
  # Validate and save numeric values
  observeEvent(input$save_radon_btn, {
    ev <- suppressWarnings(as.numeric(input$radon_exit_velocity_input))
    d  <- suppressWarnings(as.numeric(input$radon_diameter_input))
    
    # basic validation
    if (is.na(ev) || is.na(d) || ev < 0 || d <= 0) {
      showModal(modalDialog(
        title = "Invalid values",
        "Please enter a non-negative exit velocity (m³/s) and a positive diameter (m).",
        easyClose = TRUE
      ))
      return(invisible(NULL))
    }
    
    radon_exit_velocity(ev)
    radon_diameter(d)
    removeModal()
  })
  
  # Small summary shown on the Radon tab
  output$radon_params_summary <- renderText({
    ev <- radon_exit_velocity()
    d  <- radon_diameter()
    if (is.na(ev) || is.na(d)) {
      "No radon parameters saved yet."
    } else {
      paste0(
        "Saved Radon parameters — Exit velocity: ",
        signif(ev, 6), " m³/s; Diameter: ",
        signif(d, 6), " m."
      )
    }
  })
}

# Define the viewer with specific width and height (in pixels)
vwr <- dialogViewer('COMPLY v2.0', width = 3200, height = 7400)
runGadget(shinyApp(ui = ui, server = server), viewer = vwr)
#shinyApp(ui, server)

