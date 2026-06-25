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
library(openair)
library(stringr)

# Helper to ensure column names are syntactic
make_syntactic <- function(df) {
  names(df) <- make.names(names(df))
  df
}

# Used in Level 1: load possession table
rads_src <- read.csv(file = "comply171_original_possession_data.csv",
                     skip = 0,
                     header = TRUE)
long_rads <- melt(
  rads_src,
  id.vars = "Nuclide",
  variable.name = "Form",
  value.name = "Allowance"
)
long_rads$Radionuclide <- gsub(" ", "", long_rads$Nuclide)
long_rads$Type <- gsub(" ", "", long_rads$Form)
long_rads <- long_rads[, c(4, 5, 3)]
# long_rads columns: Radionuclide, Type, Allowance

# Used in Level 2: dose coefficients/lookup
doses_info <- read.csv(file = "comply171_original_input_data.csv",
                       skip = 0,
                       header = TRUE)
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

#
# # Wind Rose Funcs
# convert_dir_to_deg <- function(x) {
#   lookup <- c(
#     "N"=0, "NNE"=22.5, "NE"=45, "ENE"=67.5,
#     "E"=90, "ESE"=112.5, "SE"=135, "SSE"=157.5,
#     "S"=180, "SSW"=202.5, "SW"=225, "WSW"=247.5,
#     "W"=270, "WNW"=292.5, "NW"=315, "NNW"=337.5
#   )
#   lookup[toupper(trimws(x))]
# }
# parse_wind_file <- function(path) {
#   raw <- readLines(path)
#   # Extract part after "|"
#   clean <- str_split_fixed(raw, "\\|", 2)[, 2] |> str_trim()
#   df <- str_split_fixed(clean, ",", 5) |>
#     as.data.frame(stringsAsFactors = FALSE)
#   names(df) <- c("Dir", "Freq", "Speed", "Dist", "Toward")
#   df$Freq  <- as.numeric(df$Freq)
#   df$Speed <- as.numeric(df$Speed)
#   # Expand rows by frequency to approximate observations
#   mult <- 1000
#   # Convert Freq to numeric safely
#   df$Freq <- suppressWarnings(as.numeric(df$Freq))
#   # Replace NA or invalid frequencies with zero
#   df$Freq[is.na(df$Freq) | df$Freq < 0] <- 0
#   # Compute times while guaranteeing valid integers
#   times <- pmax(round(df$Freq * 1000), 1)
#   times[is.na(times) | times < 1] <- 1   # absolute safety net
#   expanded <- df[rep(seq_len(nrow(df)), times = times), ]
#   data.frame(
#     ws = expanded$Speed,
#     wd = convert_dir_to_deg(expanded$Dir)
#   ) |>
#     dplyr::filter(!is.na(wd))   # remove impossible directions
#
#
# }



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
            tags$head(tags$style(
              HTML(
                "
              .dataTable td, .dataTable th { border: 1px solid black; }
              .dt-button {
                background-color: #f0f0f0 !important;
                color: black !important;
                border: 1px solid #000500 !important;
              }"
              )
            )),
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
              width = '20%'
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
    ###### LEVEL 2
    nav_panel(
      "level 2/3",
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
          nav_panel(
            "Stacks",
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
          nav_panel(
            "Emissions",
            # Buttons for Level 2 emissions table
            tags$div(
              class = "d-flex gap-2 mb-2",
              actionButton("add_rad_button_l2", "Add Radionuclide", class = "btn-sm"),
              actionButton("delete_rad_button_l2", "Delete Selected Row", class = "btn-sm")
            ),
            tags$head(tags$style(
              HTML(
                "
      .dataTable td, .dataTable th { border: 1px solid black; }
      .dt-button {
        background-color: #f0f0f0 !important;
        color: black !important;
        border: 1px solid #000500 !important;
      }
    "
              )
            )),
            DTOutput("rad_table_l2")
          ),
          nav_panel(
            "Diagnostic Results",
            tags$head(tags$style(
              HTML(
                "
              .dataTable td, .dataTable th { border: 1px solid black; }
              .dt-button {
                background-color: #f0f0f0 !important;
                color: black !important;
                border: 1px solid #000500 !important;
              }"
              )
            )),
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
              width = '20%'
            ),
            downloadButton(
              "downloadReport",
              "Generate HTML Report",
              class = "btn-sm",
              width = '20%'
            )
          )
        )
      )
    ),
    
    
    nav_panel(
      "Radon",
      
      h4("Radon Modeling Tools"),
      
      layout_columns(
        # ---- Wind Rose card ----
        card(
          full_screen = TRUE,
          height = "450px",
          card_header("Wind Rose"),
          plotOutput("radon_wind_rose", height = "380px")
        ),
        
        # ---- Additional Plot card ----
        card(
          full_screen = TRUE,
          height = "450px",
          # card_header("Radon Plot"),
          # plotOutput("radon_plot", height = "380px")
        ),
        
        col_widths = c(6, 6)
      ),
      
      # ---- Data & Parameters (bottom box) ----
      card(
        full_screen = TRUE,
        height = "620px",
        
        card_header("Radon Data & Parameters"),
        
        
        # # Put the params button inside the data box (per your request)
        # tags$div(class = "d-flex gap-2 mb-2",
        #          actionButton("radon_params_btn", "Set Radon Parameters", class = "btn-sm")
        # ),
        # 
        # tags$small(class = "text-muted",
        #            "Upload wind data and configure radon stack parameters here."
        # ),
        # 
        # # Wind file upload (sectors)
        # fileInput("radon_wind_file", "Upload Wind Data File",
        #           accept = c(".txt", ".csv")),
        # 
        
        tags$div(
          class = "d-flex align-items-end gap-3 mb-2",
          tags$div(
            class = "flex-fill",
            # d-grid + w-100 makes the button fill its half of the row
            tags$div(class = "d-grid",
                     actionButton("radon_params_btn", "Set Radon Parameters", class = "btn-sm btn-primary")
            )
          ),
          tags$div(
            class = "flex-fill",
            fileInput("radon_wind_file", "Upload Wind Data File",
                      accept = c(".txt", ".csv"),
                      width = "100%")
          )
        ),
        # Parsed/expanded wind data
        DTOutput("radon_allout_table"),
        
        tags$hr(),
        
        # ---- Final dose/conc results (allout) ----
        h5("Final Sector Results"),
        br(),
        textOutput("radon_total_dose")
      )
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
          columnDefs = list(list(
            className = 'dt-center', targets = "_all"
          )),
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
          
          # Footnote inside the modal
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
        showModal(
          modalDialog(
            title = "No row selected",
            "Click a row in the Stacks table to select it, then press 'Delete Selected Row'.",
            easyClose = TRUE
          )
        )
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
          ordering = FALSE,
          searching = FALSE,
          info = FALSE,
          lengthChange = FALSE,
          dom = 'Blfti',
          columnDefs = list(list(
            className = 'dt-center', targets = "_all"
          )),
          rownames = FALSE,
          pageLength = -1,
          paging = FALSE,
          filter = "top",
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
          selectInput(
            "input_RN",
            "Radionuclide:",
            choices = unique(long_rads$Radionuclide),
            selected = NULL
          ),
          selectInput(
            "input_TYPE",
            "Type:",
            choices = unique(long_rads$Type),
            selected = NULL
          ),
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
      if (is.null(sel) || length(sel) == 0) {
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
      
      allow_row <- subset(long_rads,
                          Type == input$input_TYPE & Radionuclide == input$input_RN)
      if (nrow(allow_row) == 0 || is.na(allow_row$Allowance[1])) {
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
          selectInput(
            "input_STACK",
            "Stack ID:",
            choices = stack_choices,
            selected = NULL
          ),
          selectInput(
            "input_RN",
            "Radionuclide:",
            choices = unique(long_rads$Radionuclide),
            selected = NULL
          ),
          selectInput(
            "input_TYPE",
            "Type:",
            choices = unique(long_rads$Type),
            selected = NULL
          ),
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
      if (is.null(sel) || length(sel) == 0) {
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
        
        # Optional: update possession summaries if you show them in Level 1
        if (exists("rad_sum_possession") &&
            exists("i_sum_possession")) {
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
          columnDefs = list(list(
            className = 'dt-center', targets = "_all"
          )),
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
          capsig = sqrt((sigmaz ^ 2) + (area / 3.14)),
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
          dose  = WBdose.y    * conc,
          # columns must exist in doses_info (aligned earlier)
          vdose = VegDose     * conc,
          cdose = CowMilkDose * conc,
          total = dose + vdose + cdose
        )
      
      df
    })
    
    # Dose table (full)
    output$dose_table <- renderDT({
      req(combined_data_l2())
      datatable(
        combined_data_l2(),
        class = 'cell-border stripe',
        extensions = 'Buttons',
        options = list(
          ordering = FALSE,
          searching = FALSE,
          info = FALSE,
          lengthChange = FALSE,
          dom = 'Blfti',
          columnDefs = list(list(
            className = 'dt-center', targets = "_all"
          )),
          rownames = FALSE,
          pageLength = -1,
          paging = FALSE,
          filter = "top",
          buttons = c('copy', 'csv', 'excel', 'pdf'),
          exportOptions = list(columns = ':visible')
        )
      )
    })
    
    # Cut-down dose-specific version for display/export
    dose_data <- reactive({
      full_data_l2 <- combined_data_l2()
      full_data_l2 %>%
        select(Stack.ID,
               Radionuclide,
               Source.Distance,
               airconc,
               dose,
               vdose,
               cdose,
               total)
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
      xrads$sigmaz  <- (0.06 * xrads$Source.Distance) / sqrt(1.0 + (0.0015 * xrads$Source.Distance))
      xrads$capsig  <- sqrt((xrads$sigmaz ^ 2) + (xrads$area / 3.14))
      xrads$pp      <- 2.032 / (xrads$Source.Distance * xrads$capsig)
      xrads$conc.mei <- ifelse(
        xrads$Source.Distance < xrads$xcrit,
        (0.25 * as.numeric(xrads$Emission)) / (3.14 * xrads$Mean.WS * 1.0 * xrads$xarea),
        (0.25 * as.numeric(xrads$Emission) * xrads$pp) / xrads$Mean.WS
      )
      
      # Veg receptor
      xrads$sigmaz  <- (0.06 * xrads$Veg.Distance) / sqrt(1.0 + (0.0015 * xrads$Veg.Distance))
      xrads$capsig  <- sqrt((xrads$sigmaz ^ 2) + (xrads$area / 3.14))
      xrads$pp      <- 2.032 / (xrads$Veg.Distance * xrads$capsig)
      xrads$conc.veg <- ifelse(
        xrads$Veg.Distance < xrads$xcrit,
        (0.25 * as.numeric(xrads$Emission)) / (3.14 * xrads$Mean.WS * 1.0 * xrads$xarea),
        (0.25 * as.numeric(xrads$Emission) * xrads$pp) / xrads$Mean.WS
      )
      
      # Cow receptor
      xrads$sigmaz  <- (0.06 * xrads$Cow.Distance) / sqrt(1.0 + (0.0015 * xrads$Cow.Distance))
      xrads$capsig  <- sqrt((xrads$sigmaz ^ 2) + (xrads$area / 3.14))
      xrads$pp      <- 2.032 / (xrads$Cow.Distance * xrads$capsig)
      xrads$conc.cow <- ifelse(
        xrads$Cow.Distance < xrads$xcrit,
        (0.25 * as.numeric(xrads$Emission)) / (3.14 * xrads$Mean.WS * 1.0 * xrads$xarea),
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
      full_data %>% select(Stack.ID,
                           Radionuclide,
                           Source.Distance,
                           airconc,
                           dose,
                           vdose,
                           cdose,
                           total)
    })
    
    # ========== COMPLIANCE SUMMARIES ==========
    output$fraction_display <- renderText({
      selected_row <- rad_sum_possession()
      selected_age <- round(selected_row$Fraction, 2)
      paste0(
        "Stack possession fraction is ",
        selected_age,
        " times the allowable amount given in the possession table."
      )
    })
    
    output$fraction_display_i <- renderText({
      selected_row <- i_sum_possession()
      selected_age <- round(selected_row$Fraction, 2)
      paste0("Total iodine possession fraction is ", selected_age, ".")
    })
    
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
        return(
          plotly::plotly_empty(type = "scatter", mode = "markers") |>
            layout(title = "No valid fraction data to plot yet.")
        )
      }
      
      # aggregate by radionuclide (if multiple rows per RN exist)
      agg <- df %>%
        dplyr::group_by(Radionuclide) %>%
        dplyr::summarise(Fraction = sum(Fraction, na.rm = TRUE),
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
        return(
          plotly::plotly_empty(type = "scatter", mode = "markers") |>
            layout(title = "Possession Fraction Plot")
        )
      }
      
      agg <- df %>%
        dplyr::group_by(Radionuclide) %>%
        dplyr::summarise(Fraction = sum(Fraction, na.rm = TRUE),
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
    
    output$graph_l2 <- renderPlotly({
      # Get the computed L2 data frame
      df <- as.data.frame(combined_data_l2())
      
      # If the reactive returns NULL early, bail out
      if (is.null(df)) {
        return(
          plotly::plotly_empty(type = "scatter", mode = "markers") |>
            layout(title = "No data to plot yet.")
        )
      }
      
      # Coerce types (defensive) and keep only valid rows
      df <- df %>%
        dplyr::mutate(Source.Distance = suppressWarnings(as.numeric(Source.Distance)),
                      total           = suppressWarnings(as.numeric(total))) %>%
        dplyr::filter(
          !is.na(Source.Distance),
          is.finite(Source.Distance),
          !is.na(total),
          is.finite(total),
          !is.na(Radionuclide),
          nzchar(Radionuclide)
        ) %>%
        dplyr::arrange(Radionuclide, Source.Distance)
      
      if (nrow(df) == 0) {
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
      res <- tryCatch(
        geocoded_address(),
        error = function(e)
          NULL
      )
      lat_in <- suppressWarnings(as.numeric(input$Latitude))
      lon_in <- suppressWarnings(as.numeric(input$Longitude))
      lat <- if (!is.null(res) &&
                 !is.null(res$latitude))
        as.numeric(res$latitude)
      else
        lat_in
      lon <- if (!is.null(res) &&
                 !is.null(res$longitude))
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
      p <- ggplot(rad_data(), aes(x = 1, y = Fraction, fill = Radionuclide)) +
        geom_col(alpha = 0.6) +
        labs(title = "Possession Fraction Plot", color = "Source Data") +
        theme_classic() + xlab(" ") + ylab("Fraction") +
        theme(
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(
            color = "black",
            fill = NA,
            size = 1
          ),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white", color = NA)
        ) +
        labs(color = "Radionuclides")
      ggplotly(p)
    })
    
    output$HTML_table <- renderDT({
      report_data()
    })
    output$HTML_plot <- renderPlotly({
      report_plot()
    })
    
    output$downloadReport <- downloadHandler(
      filename = function() {
        paste0("report_", Sys.Date(), ".html")
      },
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
        some_addresses <- tibble::tribble( ~ name, ~ addr, "company", full_address)
        geocode_result <- some_addresses %>%
          geocode(addr,
                  method = 'osm',
                  lat = latitude,
                  long = longitude)
        if (nrow(geocode_result) > 0 &&
            !is.na(geocode_result$latitude) &&
            !is.na(geocode_result$longitude)) {
          geocode_result
        } else {
          showModal(
            modalDialog(
              title = "Geocoding Error",
              "Could not find coordinates for the entered address. Please check the address and try again.",
              easyClose = TRUE
            )
          )
          NULL
        }
      } else {
        latitude <- input$Latitude
        longitude <- input$Longitude
        name <- "company"
        addr <- ""
        geocode_result <- tibble::as_tibble(cbind(name, addr, as.numeric(latitude), as.numeric(longitude)))
        names(geocode_result) <- c("name", "addr", "latitude", "longitude")
        geocode_result
      }
    })
    
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
    
    # Update both maps when a new address is geocoded
    observe({
      result <- geocoded_address()
      if (is.null(result))
        return()
      result$latitude  <- as.numeric(result$latitude)
      result$longitude <- as.numeric(result$longitude)
      for (id in c("map", "map_l2")) {
        leafletProxy(id) %>%
          clearMarkers() %>%
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
    
    output$coords <- renderText({
      result <- geocoded_address()
      if (is.null(result))
        return("Enter an address and click 'Plot Location'.")
      result$latitude <- as.numeric(result$latitude)
      result$longitude <- as.numeric(result$longitude)
      paste(
        "Latitude:",
        round(result$latitude, 6),
        ", Longitude:",
        round(result$longitude, 6)
      )
    })
    
    #  ---- RADON PARAMETERS (stored as numeric reactive values) ----
    #  These can be called later as radon_exit_velocity() and radon_diameter()
    exit_velocity <- reactiveVal(NA_real_)
    stack_diameter      <- reactiveVal(NA_real_)
    stack_temp      <- reactiveVal(NA_real_)
    air_temp      <- reactiveVal(NA_real_)
    vertical_momentum      <- reactiveVal(NA_real_)
    emission      <- reactiveVal(NA_real_)
    percent_calm      <- reactiveVal(NA_real_)
    release_height <- reactiveVal(NA_real_)
    
    
    # Open the modal when user clicks the button on the Radon tab
    observeEvent(input$radon_params_btn, {
      showModal(
        modalDialog(
          title = "Radon Stack Parameters",
          #textInput("stack_id", "Stack ID:", value = ""),
          
          numericInput("radon_diam", "Diameter (m):", value = stack_diameter(), min = 0.01),
          numericInput("radon_vel",  "Exit Velocity (m³/s):", value = exit_velocity(), min = 0),
          
          numericInput("radon_ht",   "Release Height (m):", value = release_height(), min = 0),
          numericInput("radon_tstack", "Stack Temperature (°C):", value = stack_temp(), min = -100),
          numericInput("radon_tair",   "Air Temperature (°C):", value = air_temp(), min = -100),
          
          selectInput("radon_vertmom", "Vertical Momentum?", choices = c("YES","NO"),
                      selected = vertical_momentum()),
          
          numericInput("radon_emit", "Radon Emissions (Ci/year):", value = emission(), min = 0),
          numericInput("radon_calm", "Percent Calms (%):", value = percent_calm(), min = 0, max = 100),
          
          footer = tagList(
            actionButton("radon_save", "Save", class = "btn-primary"),
            modalButton("Cancel")
          ),
          easyClose = TRUE
        )
      )
    })
    
    observeEvent(input$radon_save, {
      #stack_id(input$radon_id)
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
    
    # ---- Convert compass (N, NNE, NE, ...) to degrees ----
    convert_dir_to_deg <- function(x) {
      lookup <- c(
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
      # Trim/uppercase to ensure matching
      out <- lookup[toupper(trimws(x))]
      # Return NA for any unknown tokens
      as.numeric(out)
    }
    
    # ---- Safe expansion by frequency (prevents invalid rep(times)) ----
    safe_expand_by_freq <- function(df,
                                    freq_col = "Freq",
                                    mult = 1000L) {
      f <- suppressWarnings(as.numeric(df[[freq_col]]))
      f[is.na(f) | f < 0] <- 0
      times <- pmax(round(f * mult), 1L)
      times[is.na(times) | times < 1L] <- 1L
      df[rep.int(seq_len(nrow(df)), times = times), , drop = FALSE]
    }
    
    # ---- Parse uploaded wind file (handles your exact line format) ----
    parse_wind_file <- function(path) {
      raw <- readLines(path, warn = FALSE)
      rhs <- vapply(strsplit(raw, "\\|"), function(z)
        if (length(z) >= 2)
          z[2]
        else
          NA_character_, "")
      rhs <- trimws(rhs)
      rhs <- rhs[!is.na(rhs) & nzchar(rhs)]
      m <- do.call(rbind, strsplit(rhs, ","))
      df <- as.data.frame(m, stringsAsFactors = FALSE)
      if (ncol(df) < 5) {
        stop("Wind file lines must have 5 fields: Dir, Freq, Speed, Dist, Toward.")
      }
      names(df) <- c("Dir", "Freq", "Speed", "Dist", "Toward")
      df$Freq  <- suppressWarnings(as.numeric(df$Freq))
      df$Speed <- suppressWarnings(as.numeric(df$Speed))
      expanded <- safe_expand_by_freq(df, freq_col = "Freq", mult = 1000L)
      wd_raw <- trimws(expanded$Dir)
      wd_num <- suppressWarnings(as.numeric(wd_raw))
      if (all(!is.na(wd_num))) {
        wd_deg <- wd_num
      } else {
        wd_deg <- convert_dir_to_deg(wd_raw)
      }
      keep <- !is.na(wd_deg) &
        is.finite(wd_deg) & wd_deg >= 0 & wd_deg <= 360
      expanded <- expanded[keep, , drop = FALSE]
      wd_deg   <- wd_deg[keep]
      out <- data.frame(ws = expanded$Speed,
                        wd = wd_deg,
                        stringsAsFactors = FALSE)
      out <- out[!is.na(out$ws) &
                   is.finite(out$ws) & out$ws >= 0, , drop = FALSE]
      if (nrow(out) == 0) {
        stop("No valid wind rows after parsing/cleaning (check Freq/Speed/Dir fields).")
      }
      
      out
    }
    
    
    radon_wind_data <- reactive({
      req(input$radon_wind_file)
      validate(need(
        file.exists(input$radon_wind_file$datapath),
        "Uploaded file not found. Please re-upload."
      ))
      parse_wind_file(input$radon_wind_file$datapath)
    })
    
    
    output$radon_data_table <- renderDT({
      req(radon_wind_data())
      DT::datatable(
        radon_wind_data(),
        class = 'cell-border stripe',
        options = list(dom = 'Blfti', pageLength = 10)
      )
    })
    
    
    output$radon_wind_rose <- renderPlot({
      req(radon_wind_data())
      openair::windRose(
        radon_wind_data(),
        ws = "ws",
        wd = "wd",
        paddle = FALSE,
        key.header = "Wind Speed (m/s)",
        key.footer = " ",
        cols = "hue"
      )
    })
    
    # --------- helpers ----------
    # Map compass to degrees (for the wind rose only)
    convert_dir_to_deg <- function(x) {
      lookup <- c(
        "N"=0, "NNE"=22.5, "NE"=45, "ENE"=67.5,
        "E"=90, "ESE"=112.5, "SE"=135, "SSE"=157.5,
        "S"=180, "SSW"=202.5, "SW"=225, "WSW"=247.5,
        "W"=270, "WNW"=292.5, "NW"=315, "NNW"=337.5
      )
      out <- lookup[toupper(trimws(x))]
      as.numeric(out)
    }
    
    # Safe expansion by frequency to create enough points for windRose()
    safe_expand_by_freq <- function(df, mult = 1000L) {
      f <- suppressWarnings(as.numeric(df$Freq))
      f[is.na(f) | f < 0] <- 0
      times <- pmax(round(f * mult), 1L)
      times[is.na(times) | times < 1L] <- 1L
      df[rep.int(seq_len(nrow(df)), times = times), , drop = FALSE]
    }
    
    # Parse your wind-sector file (RHS after the '|')
    parse_wind_bins <- function(path) {
      raw <- readLines(path, warn = FALSE)
      rhs <- vapply(strsplit(raw, "\\|"), function(z) if (length(z) >= 2) z[2] else NA_character_, "")
      rhs <- trimws(rhs)
      rhs <- rhs[!is.na(rhs) & nzchar(rhs)]
      
      m <- do.call(rbind, strsplit(rhs, ","))
      df <- as.data.frame(m, stringsAsFactors = FALSE)
      if (ncol(df) < 5) stop("Each line must be: Dir,Freq,Speed,Distance,Sector")
      names(df) <- c("DirFrom", "Freq", "Speed", "Dist", "SectorTo")
      
      df$Freq  <- suppressWarnings(as.numeric(df$Freq))
      df$Speed <- suppressWarnings(as.numeric(df$Speed))
      df$Dist  <- suppressWarnings(as.numeric(df$Dist))
      
      # keep only valid rows
      df <- df[!is.na(df$Freq) & !is.na(df$Speed) & !is.na(df$Dist), , drop = FALSE]
      if (nrow(df) == 0) stop("No valid wind rows found after parsing.")
      df
    }
    
    # For the wind rose
    to_windrose_df <- function(bins_df) {
      expanded <- safe_expand_by_freq(
        data.frame(Dir = bins_df$DirFrom, Freq = bins_df$Freq, Speed = bins_df$Speed),
        mult = 1000L
      )
      wd <- convert_dir_to_deg(expanded$Dir)
      out <- data.frame(ws = expanded$Speed, wd = wd)
      out <- out[!is.na(out$wd) & is.finite(out$wd) & out$wd >= 0 & out$wd <= 360 &
                   !is.na(out$ws) & is.finite(out$ws) & out$ws >= 0, , drop = FALSE]
      out
    }
    
    # --------- reactives ----------
    # Read user wind file (sector bins)
    radon_wind_bins <- reactive({
      req(input$radon_wind_file)
      parse_wind_bins(input$radon_wind_file$datapath)
    })
    
    # Show parsed wind bins
    output$radon_data_table <- DT::renderDT({
      req(radon_wind_bins())
      DT::datatable(radon_wind_bins(),
                    class = 'cell-border stripe',
                    options = list(dom = 'Blfti', pageLength = 16))
    })
    
    # Wind rose using expanded freq samples
    output$radon_wind_rose <- renderPlot({
      req(radon_wind_bins())
      df <- to_windrose_df(radon_wind_bins())
      openair::windRose(df, ws = "ws", wd = "wd",
                        paddle = TRUE, key.header = "Wind Speed (m/s)",
                        main = "Wind Rose")
    })
    
    # Simple scatter demo (you can customize later)
    output$radon_plot <- renderPlot({
      req(radon_wind_bins())
      df <- radon_wind_bins()
      plot(df$Dist, df$Speed,
           xlab = "Receptor Distance (m)", ylab = "Wind Speed (m/s)",
           main = "Wind Speed vs Distance", pch = 19, col = "#2C7FB8")
    })
    
    # --------- RADON SECTOR CALC (your formulas) ----------
    radon_allout <- reactive({
      bins <- radon_wind_bins()
      req(nrow(bins) > 0)
      
      # Pull user-entered params (reactiveVals you already have)
      releaseht        <- suppressWarnings(as.numeric(release_height()))   
      exitvelocity     <- suppressWarnings(as.numeric(exit_velocity()))
      diameter         <- suppressWarnings(as.numeric(stack_diameter()))
      stacktemp        <- suppressWarnings(as.numeric(stack_temp()))
      airtemp          <- suppressWarnings(as.numeric(air_temp()))
      verticalmomentum <- toupper(trimws(vertical_momentum()))
      rademis          <- suppressWarnings(as.numeric(emission()))
      pctcalm_input    <- suppressWarnings(as.numeric(percent_calm()))
      
      # Correct: release height should not be stack_temp; use your actual release height source
      # If you store Release Height elsewhere, replace with that reactive.
      # For now, defensively set to 0 if missing:
      releaseht <- if (!is.na(releaseht)) releaseht else 0
      
      # Safety defaults
      if (is.na(exitvelocity) || exitvelocity <= 0) exitvelocity <- 0.3
      if (is.na(diameter)     || diameter <= 0)     diameter     <- 0.1
      if (is.na(stacktemp)) stacktemp <- 20
      if (is.na(airtemp))   airtemp   <- 20
      if (is.na(rademis) || rademis < 0) rademis <- 0
      
      calm_frac <- pctcalm_input
      # If user supplies % calm, convert to fraction:
      if (!is.na(calm_frac) && calm_frac > 1) calm_frac <- calm_frac / 100
      if (is.na(calm_frac)) calm_frac <- 0
      calm_frac <- pmin(pmax(calm_frac, 0), 0.99)  # clamp to [0, 0.99]
      
      # constants
      pi_ <- 3.14159265
      g   <- 9.81
      
      # accumulators
      maxdose   <- 0.0
      maxsector <- 0L
      maxdist   <- max(bins$Dist, na.rm = TRUE)
      totaldose <- 0.0
      
      # vectors to fill
      nwindbins <- nrow(bins)
      conc    <- numeric(nwindbins)
      ihdose  <- numeric(nwindbins)
      
      # loop over sectors
      for (i in seq_len(nwindbins)) {
        wind_freq  <- bins$Freq[i]
        wind_speed <- bins$Speed[i]
        rdist      <- bins$Dist[i]
        
        # avoid divide by zero
        wind_speed <- ifelse(is.na(wind_speed) || wind_speed <= 0, 0.3, wind_speed)
        
        HEFF <- releaseht
        
        if (releaseht > 0) {
          # Air & stack density
          RHOAIR <- 0.0735 * (80.0 + 460.0) / (airtemp + 460.0)
          RHOAIR <- RHOAIR * 454.0 / (0.3048^3)
          RHOSTK <- RHOAIR * (airtemp + 460.0) / (stacktemp + 460.0)
          
          # Geometry & momentum terms
          AREA   <- pi_ * (diameter^2.0) / 4.0
          STKVEL <- exitvelocity / AREA
          R      <- STKVEL / wind_speed
          FM     <- (RHOSTK / RHOAIR) * STKVEL * exitvelocity / pi_
          
          TERM1  <- 75.0 * FM / (((1.0 + 3.0 / R) * wind_speed)^2.0)
          TERM2  <- releaseht + 7.79 * sqrt(FM) / ((1.0 + 3.0 / R) * wind_speed)
          
          if (!is.na(TERM1) && !is.na(TERM2) && verticalmomentum == "YES") {
            DELHMO <- 0.93 * (TERM1^(3/7)) * (TERM2^(1/7))
          } else {
            DELHMO <- 0.0
          }
          
          FB     <- g * (RHOAIR - RHOSTK) * exitvelocity / (pi_ * RHOAIR)
          DELHBU <- 23.67 * ((FB / (wind_speed^3.0))^0.6) *
            (releaseht + 195.1 * FB / (wind_speed^3.0))^0.4
          
          HEFF <- releaseht + max(DELHMO, DELHBU, na.rm = TRUE)
          if (is.infinite(HEFF) || is.na(HEFF)) HEFF <- releaseht
        }
        
        SIGMAZ <- (0.06 * rdist) / sqrt(1.0 + 0.0015 * rdist)
        CAPSIG <- sqrt(SIGMAZ^2 + 0.0 / pi_)
        
        zz <- 0.5 * ((HEFF / CAPSIG)^2)
        
        if (!is.na(zz) && zz > 0.001) {
          pp <- 2.032 * exp(-zz) / (rdist * CAPSIG)
        } else {
          pp <- 2.032 / (rdist * CAPSIG)
        }
        
        # concentration
        denom <- (1.0 - calm_frac)
        denom <- ifelse(denom <= 0, 1e-6, denom)
        conc_i <- ( (wind_freq / denom) * rademis * pp ) / wind_speed
        if (!is.na(conc_i) && conc_i < 1e-32) conc_i <- 1e-32
        conc[i] <- conc_i
        
        # inhalation dose for sector i
        INHFAC <- 0.0092 * 100.0 * 1000.0
        TINIT  <- 11.8
        T      <- TINIT + rdist / (wind_speed * 60.0)
        
        EQILOUT <- 1 - 0.0479 * exp(-T / 4.39) - 2.1963 * exp(-T / 38.6) + 1.2442 * exp(-T / 28.4)
        if (!is.na(EQILOUT)) EQILOUT <- min(EQILOUT, 1.0)
        
        EQILIN <- 0.35 * (1 + EQILOUT)
        if (!is.na(EQILIN)) EQILIN <- min(EQILIN, 0.7)
        
        INHDOS  <- conc_i * 1e12 * (1/1000) * (1/100) * (8760/170) * INHFAC
        FINDOOR <- 0.75
        
        ihdose_i <- INHDOS * (FINDOOR * EQILIN + (1 - FINDOOR) * EQILOUT)
        ihdose[i] <- ihdose_i
        
        if (!is.na(ihdose_i) && ihdose_i > maxdose) {
          maxdose   <- ihdose_i
          maxsector <- i
        }
      } # end sector loop
      
      totaldose <- totaldose + ifelse(is.na(maxdose), 0, maxdose)
      
      # build allout rows
      out <- data.frame(
        Stack              = 1L,
        `Emissions (Ci/y)` = rademis,
        Sector             = seq_len(nwindbins),
        `Max Sector`       = maxsector,
        `Sector Name`      = bins$SectorTo,
        `Distance (m)`     = bins$Dist,
        `Total Dose (mrem/y)` = ihdose,
        `Conc (Ci/m3)`        = conc,
        `Wind Frequency`      = bins$Freq,
        `Wind Speed (m/s)`    = bins$Speed,
        check.names = FALSE
      )
      
      attr(out, "totaldose") <- totaldose
      out
    })
    
    # Final table in the bottom Radon box
    output$radon_allout_table <- DT::renderDT({
      req(radon_allout())
      DT::datatable(radon_allout(),
                    class = 'cell-border stripe',
                    extensions = 'Buttons',
                    options = list(
                      dom = 'Blfti',
                      buttons = c('copy','csv','excel','pdf'),
                      pageLength = 12
                    ))
    })
    
    # Show total dose summary under the table
    output$radon_total_dose <- renderText({
      df <- radon_allout()
      td <- attr(df, "totaldose")
      paste0("Sum of sector MEI inhalation dose (selects max sector per stack): ",
             signif(td, 6), " mrem/y")
    })
    
  }
  
  # Define the viewer with specific width and height (in pixels)
  vwr <- dialogViewer('COMPLY v2.0', width = 3200, height = 7400)
  runGadget(shinyApp(ui = ui, server = server), viewer = vwr)
  #shinyApp(ui, server)
  