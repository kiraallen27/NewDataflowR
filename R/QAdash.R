#'@name QAdash
#'@title Creates interactive dashboard for data QA
#'@description Creates interactive dashboard for data QA
#'@param yearmon numeric designation of survey date formatted as yyyymm
#'@param original.data specify whether to read in original streamcleaned data (TRUE) or already QAd data (FALSE)
#'@param fdir character file path to local data directory
#'@export
#'@import shiny
#'@import shinyWidgets
#'@import DT
#'@import plotly
#'@import dplyr
#'@import tibble
#'@import lubridate
#'@import readr
#'@import rlang
#'@examples \dontrun{
#'QAdash(yearmon=202510, original.data=TRUE)
#'}

QAdash <- function(yearmon, original.data=FALSE, fdir = getOption("fdir")) {
  trip_lookup_raw <- read.csv(file.path(fdir, "trip_lookup.csv"))


  # ==========================================================
  # STEP 1: Start of GUI Function for filtering on instrument
  # ==========================================================

  prepare_lookup_tables <- function(instrument,
                                    trip_lookup_raw) {   #{ indicates start of function

    # SAFETY: validate instrument (can be one OR two instruments)
    valid_instruments <- c("YSI", "C6", "Manta", "DF")

    # check all provided instruments are valid
    if (!all(instrument %in% valid_instruments)) {
      stop("Instrument must be one or more of: YSI, C6, Manta, DF")
    }


    #=============================
    #STEP 2: Filter by instrument
    #=============================
    # Clean column names - may need to change if get updated in lookup table; added cleaning/normalizing parameter synonyms
    #here just renaming
    trip_lookup <- trip_lookup_raw %>%
      dplyr::filter(.data$Instrument_General %in% instrument) %>%    #added this in to filter on Instrument_General column in lookup table-is only: YSI or Manta or DF or +C6
      dplyr::rename(
        #Instrument*updatedname            = Instrument, #don't need but these are placeholders in case gets changed somehow
        #Standard_Parameter*updatedname    = Standard_Parameter, #placeholder same as "Instrument" note
        parameter_names_raw   = Parameter_column_name
        # Units*updatedname                 = Units, #placeholder for column
        #Probe_ID_Final*updated name        = Probe_ID_Final #placeholder column so is in here in case need to update/rename
      )

    #======================
    #STEP 3: Clean synonyms
    #======================
    #synonym cleaning (synonyms = alternate names for columns of one param; eg. sal, salinity, salinity.pss)
    trip_lookup <- trip_lookup %>%
      dplyr::mutate(
        # Clean synonyms split by comma, trim whitespace, tolower() for matching, remove empty strings, remove stray punctuation (., .., etc.); turns
        # Parameter_column_name in trip_lookup_raw to new column called "parameter_synonyms" that reads as a list (e.g., c("salinity", "salinity.pss"))
        parameter_synonyms = parameter_names_raw %>%
          gsub("[\r\n\t]", " ", .) %>%                # kill control chars and replaces w/ spaces
          strsplit(",") %>%                           # splits multiple synonyms separated by commas into individual (eg. sal, salinity, salinity.pss))
          lapply(function(x) {                        # for each vector of synonyms, apply the following cleaning steps/rules:
            x <- trimws(x)                            # removes leading/trailing whitespace
            x <- tolower(x)                           # makes everything lowercase
            x <- gsub("\\.+$", "", x)                 # remove trailing periods or sequences of periods ".."
            x <- gsub("[^a-z0-9_.]", "", x)           # remove weird chars that aren't lowercase letters, underscores, periods (weird characters/stray punctuation)
            x <- x[x != ""]                           # removes empty strings
            unique(x)                                 # removes duplicates
          })
      )

    names(trip_lookup)
    glimpse(trip_lookup$parameter_synonyms)

    #==============================
    #STEP 4: Group synonyms by parameter + instrument (Instrument key here b/c of possible redundancy of synonyms betweein instruments)
    #==============================
    ##PARAMETER LOOKUP TABLES (each column from trip_lookup becomes list), synonym checks, duplicate checks, checks in actual dataset too; columns=
    #Standard_Parameter, Instrument serial#, columns (list for different synonyms), units, probes (list for probe #s)

    trip_lookup <- trip_lookup %>%
      dplyr::group_by(Standard_Parameter, Instrument_General) %>%  #Instrument_General is YSI, Manta, DF, or C6-no differentiation by serial #
      dplyr::mutate(parameter_synonyms = list(unique(unlist(parameter_synonyms)))) %>%
      dplyr::ungroup()

    #==========================
    #STEP 5: Build param_lookup
    #============================
    param_lookup <- trip_lookup %>%
      dplyr::group_by(Standard_Parameter, Instrument_General) %>%   #Instrument_General is YSI, Manta, DF, or C6-no differentiation by serial #
      dplyr::summarise(
        columns = list(unique(unlist(parameter_synonyms))),
        Units   = unique(units)[1],   #"[1]" indicates use first one, also lowercase units
        probes  = list(unique(Probe_ID_Final)),
        .groups = "drop"
      )
    #param_lookup tab
    #trip_lookup$parameter_synonyms
    str(param_lookup)
    head(param_lookup)
    print(param_lookup, n=73) #73 should cover all entries but may need to adjust # if want to update

    #==========================
    #STEP 6: Build column map
    #============================
    #column map is list of Instruments with synonyms and indication of # of synonyms (e.g., character[3] then 'do' 'dots' 'do.ts')
    column_map <- param_lookup %>%
      dplyr::group_by(Standard_Parameter) %>%
      dplyr::summarise(
        all_possible_columns = list(unique(unlist(columns))),
        .groups = "drop"
      ) %>%
      tibble::deframe()

    str(column_map) #shows column map list

    #param_lookup %>% select(Standard_Parameter, Units)
    #names(param_lookup)
    #names(trip_lookup)
    #glimpse(param_lookup)

    #===============================================================
    #STEP 7: build unit_map, check empty synonym sets, check duplicate synonyms
    #====================================================

    #7a. unit map is list of Standard_Parameter with units e.g., "RFU", "per_sat"
    unit_map <- param_lookup %>%
      dplyr::group_by(Standard_Parameter) %>%
      dplyr::summarise(unit = unique(Units)[1], .groups = "drop") %>% #note shift from "units" lowercase to "Units" uppercase here
      tibble::deframe()

    print(unit_map) #check lookup system is not reading something as "NULL" or wonky

    #7b. check empty synonym sets; ie.e., every Standard_Parameter has >= 1 recognized synonym string (i.e., is internally consistent)
    which_empty <- sapply(column_map, length) == 0
    column_map[which_empty] #should not show anything in list

    str(param_lookup) #make sure there are no weird designations not reading as character <chr>

    #troubleshoot character so makes function instead of string
    #param_lookup %>% select(Standard_Parameter, Units)%>%
    # print(n=73)-change # based on str results

    #7c. check if duplicate synonym sets contain overlapping/ambiguous synonyms (i.e., same raw column name under multiple params)
    dup_synonyms <- tibble(
      synonym = unlist(column_map),
      parameter = rep(names(column_map), lengths(column_map))
    ) %>%
      dplyr::count(synonym, sort = TRUE) %>%
      dplyr::filter(n > 1)

    dup_synonyms #****should be SOME in trip_lookup table, e.g., parameter_synonyms has n=3 b/c "temp" for YSI OR Manta OR DF; this is OK b/c only one instrument
    #active/trip but would be an issue otherwise


    #==========================
    #STEP 8: Return results
    #============================
    return(list(
      trip_lookup  = trip_lookup,
      param_lookup = param_lookup,
      column_map   = column_map,
      unit_map     = unit_map,
      dup_synonyms = dup_synonyms,
      empty_sets   = which_empty

    ))

    #names(trip_lookup)
    #glimpse(param_lookup)

  }

  #2- flag constant

  #Global flag constants-will be defined per parameter interactively in the Shiny GUI
  FLAG_GOOD    <- 1
  FLAG_SUSPECT <- 2
  FLAG_BAD     <- 3
  #not handling NA in here; do so later and change to "Missing" to avoid trouble with handling/reading NA vs <NA>

  #define labels
  flag_labels <- c("Good",  "Suspect",  "Bad")

  # Global color map for ALL parameters - be sure to refer to this "flag_colors" in plots
  flag_colors <- c(
    "Good"    = "#2ECC71",
    "Suspect" = "orange",
    "Bad"     = "#E74C3C",
    "Missing" = "darkgray"
  )

  # Optional map of numeric->label (convenient for logs) **YES DO THIS
  flag_map <- c(
    `1` = "Good",
    `2` = "Suspect",
    `3` = "Bad"
  )

  #3 "Define" dynamic threshold function (set up format to be applied in GUI) #
  flag_by_dynamic_thresholds <- function(x,
                                         bad_min     = NA,
                                         suspect_min = NA,
                                         suspect_max = NA,
                                         bad_max     = NA)
  {
    out <- rep("Good", length(x))

    # BAD (outer)
    if (!is.na(bad_min)) out[x < bad_min] <- "Bad"
    if (!is.na(bad_max)) out[x > bad_max] <- "Bad"

    # SUSPECT (inner)
    if (!is.na(suspect_min)) out[out == "Good" & x < suspect_min] <- "Suspect"
    if (!is.na(suspect_max)) out[out == "Good" & x > suspect_max] <- "Suspect"

    # Missing
    out[is.na(x)] <- "Missing"

    out
  }

  #4 Define universal .csv loader function - this is to define "load_dataflow_csv", determine date/time column names and create "timestamp" and unique ID#
  #(vs do it by hand each time); also R is having issues reading Date Time column from a lot of combined streamcleaned files so best to create new column

  #define load_dataflow_csv function
  load_dataflow_csv <- function(path) {

    df <- read.csv(path, check.names = FALSE)

    # Required: "date" column
    if (!"date" %in% names(df)) {
      stop("ERROR: The dataset does not contain a 'date' column.") #checks there is a date column for new datetime column
    }

    # Detect time column name
    time_col <- NULL
    if ("time" %in% names(df)) time_col <- "time"
    if ("time..hh.mm.ss." %in% names(df)) time_col <- "time..hh.mm.ss."

    if (is.null(time_col)) {
      stop("ERROR: Could not find a time column ('time' or 'time..hh.mm.ss.')") #checks there is a time column to use for new datetime column
    }


    # Build new timestamp safely (MDY HMS, UTC)
    df <- df %>%
      dplyr::mutate(
        timestamp = mdy_hms(
          paste(date, .data[[time_col]]),
          tz = "UTC"
        )
      )

    # Create ID for QA tracking if not present - is a unique #
    if (!"id" %in% names(df)) {
      df <- df %>% dplyr::arrange(timestamp) %>% dplyr::mutate(id = row_number())
    }

    return(df)
  }

  #5 load dataflow data
  if(original.data == FALSE && file.exists(file.path(fdir, "DF_FullDataSets", "QA datasets", paste(yearmon, "_qa.csv", sep = "")))){
    df <- load_dataflow_csv(file.path(fdir, "DF_FullDataSets", "QA datasets", paste(yearmon, "_qa.csv", sep = "")))
  } else {
    #Read original streamcleaned file if no QA yet
    fdir_fd <- file.path(fdir, "DF_FullDataSets")
    flist <- list.files(fdir_fd, include.dirs = T, full.names = T)
    flist <- flist[substring(basename(flist),1,6) == yearmon]
    df <- load_dataflow_csv(flist)
  }

  #####################################
  #CREATE QA LOG to track corrections made to raw data
  #####################################
  #####################################
  # LOAD EXISTING QA LOG OR CREATE NEW ONE
  #####################################

  qa_log_path <- file.path(
    fdir,
    "DF_FullDataSets",
    "QA datasets",
    paste0(yearmon, "_QA_log.csv")
  )

  #-----------------------------------
  # Define empty QA log structure
  #-----------------------------------

  empty_qa_log <- tibble::tibble(
    id = integer(),
    parameter = character(),
    original_value = numeric(),
    action = character(),
    flag = character(),
    corrected_value = numeric(),
    reason = character(),
    context = character(),
    reviewer = character(),
    date = as.Date(character())
  )

  #-----------------------------------
  # Load existing QA log if available
  #-----------------------------------

  if (file.exists(qa_log_path)) {

    qa_log <- readr::read_csv(
      qa_log_path,
      show_col_types = FALSE
    )

    # Make sure the expected columns exist
    # in case an older log is missing something.
    missing_cols <- setdiff(
      names(empty_qa_log),
      names(qa_log)
    )

    if (length(missing_cols) > 0) {
      for (col in missing_cols) {
        qa_log[[col]] <- empty_qa_log[[col]]
      }
    }

    # Put columns into the expected order
    qa_log <- qa_log %>%
      dplyr::select(
        dplyr::all_of(names(empty_qa_log))
      )

  } else {

    # No existing log: start a new one
    qa_log <- empty_qa_log
  }

  # Helper: append one QA log entry to qa_log and also make sure new row matches column types in log
  append_qa_log <- function(
    qa_log, id, parameter, original_value,
    action, flag, corrected_value, reason,
    context, reviewer) {
    new_entry <- tibble::tibble(
      id              = as.integer(id),
      parameter       = as.character(parameter),
      original_value  = if(is.null(original_value)) NA_real_ else
        as.numeric(original_value),
      action          = as.character(action),          # "flag", "correct"| "set_na"
      flag            = if (is.null(flag)) NA_character_ else as.character(flag),
      corrected_value = if (is.null(corrected_value)) NA_real_ else
        as.numeric(corrected_value),
      reason          = as.character(reason),  #user-entered text only
      context         = as.character(context), #column for diagnostic levels
      reviewer        = as.character(reviewer), #fills in with reviewer name
      date            = as.Date(Sys.Date())     #fills in with date of review
    )
    dplyr::bind_rows(qa_log, new_entry)
  }

  ###Get username to populate "reviewer" column of QA Log-should work in both Shiny and batch
  #--- reviewer resolver helpers (done outside shiny server) ---
  resolve_friendly <- function(raw) {
    raw <- trimws(if (is.null(raw)) "" else raw)
    if (!nzchar(raw)) return("unknown")
    raw
  }

  get_shiny_reviewer <- function(session, prefer_full_name = TRUE) {
    # Optional override for testing/service accounts
    override <- getOption("qa.reviewer_override", NULL)
    if (!is.null(override) && nzchar(override)) return(override)

    # Shiny-hosted user if available (RStudio Connect/Shiny Server Pro can set session$user)
    if (!is.null(session$user) && nzchar(session$user)) {
      return(if (prefer_full_name) resolve_friendly(session$user) else session$user)
    }

    # Env fallbacks (batch/non-Shiny) - might not need if end up not using batch
    user <- Sys.getenv("USERNAME")
    if (!nzchar(user)) user <- Sys.getenv("USER")
    if (nzchar(user)) return(if (prefer_full_name) resolve_friendly(user) else user)

    "unknown"
  }

  get_reviewer_unified <- function(session = NULL, prefer_full_name = TRUE) {
    if (!is.null(session)) {
      return(list(reviewer = get_shiny_reviewer(session, prefer_full_name), source = "shiny"))
    }
    # batch/local
    user <- Sys.getenv("USERNAME"); if (!nzchar(user)) user <- Sys.getenv("USER")
    list(reviewer = if (nzchar(user)) resolve_friendly(user) else "unknown", source = "env")
  }

  #if want to check has reviewer name:
  id <- get_reviewer_unified(NULL, prefer_full_name = TRUE)

  ##########################################
  ################## SHINY #################
  ##########################################

  #=====================#
  #===========UI========#
  #=====================#

  ui <- fluidPage( #define buttons and inputs; fluidPage=responsive layout that automatically scales with browser window size

    titlePanel("Dataflow QA/QC"), #title-tricky to update dynamically so good w/ this for now

    # ---- Parameter Selector ----

    selectInput(
      "parameter",
      "Choose Parameter:",
      choices  = NULL        # UI starts empty, gets populated later server-side w/ updateSelectInput() once data and available params are known
    ),

    # ----Plotly outputs: RAW/CORRECTED plot vs MAP plots side by side----
    fluidRow(                                                # ensures UI elements are side-by-side inside a 12-coolumn grid
      column(6, plotlyOutput("time_plot", height="400px")),  # first number (column) refers to width; 6 is half the page width. If screen is narrow, stack
      column(6, plotlyOutput("map_plot", height="400px"))    # first number (column) refers to width - same as time_plot but for map figure
    ),
    br(),                                                    # line break
    textOutput("selection_count"),
    # ----Panels showing selected points & debug flag counts----
    DT::dataTableOutput("selected_info"),  # shows selected points, raw text formatting; helpful for debugging
    DT::dataTableOutput("debug_flags"),    # to deal with NAs; diagnostic output-"verbatim" prevents Shiny trying to format/escape

    #=======================================================
    # Instrument selectors
    # =====================================================
    #---Instrument selector buttons------
    radioButtons(
      inputId = "primary_instrument",
      label = "Select Main Instrument",
      choices = c("YSI", "Manta", "DF"),
      selected = "YSI",
      inline = TRUE                       #keeps buttons on one horizontal line instead of stacked vertically
    ),

    checkboxInput(
      inputId = "use_C6",
      label = "Include C6?",
      value = FALSE            #TRUE/FALSE toggle used later in filtering/instrument selection logic
    ),

    # =========================================================
    # RANGE-BASED (DYNAMIC) THRESHOLDING CONTROLS
    # =========================================================

    # ---- Range-based/dynamic thresholds: Good / Suspect / Bad ----
    h3("Range-based Thresholds (Good / Suspect / Bad)"),


    # Dynamic thresholds for the currently selected parameter & trip
    numericInput("range_bad_min",
                 "BAD minimum (value < this):",
                 value = NA_real_),

    numericInput("range_suspect_min",
                 "SUSPECT minimum (value < this, but >= bad_min):",
                 value = NA_real_),

    numericInput("range_suspect_max",
                 "SUSPECT maximum (value > this, but <= bad_max):",
                 value = NA_real_),

    numericInput("range_bad_max",
                 "BAD maximum (value > this):",
                 value = NA_real_),

    actionButton("apply_range_thresholds",
                 "Apply Threshold Flags"),    #triggers server logic that applies parameter-specific range rules

    #=========================================================
    # STEP-CHANGE DIAGNOSTICS CONTROLS (difference b/n pts-lets set threshold, visually mark pts, and select flagged ones w/ step change w/ a click)
    # =========================================================

    h3("Step-change Diagnostics"),

    checkboxInput(
      "enable_step",
      "Enable step-change diagnostics",
      value = FALSE
    ),

    conditionalPanel(
      condition = "input.enable_step == true",

      radioButtons(
        "step_method",
        "Step-change method:",
        choices = c(
          "Raw value change" = "raw",
          "Percent change" = "percent"
        ),
        selected = "raw",
        inline = TRUE
      ),

      numericInput(
        "step_threshold",
        "Step-change threshold:",
        value = 0.3,
        min = 0,
        step = 0.1
      ),

      checkboxInput(
        "highlight_step",
        "Highlight step-change points",
        value = TRUE
      ),

      actionButton(
        "select_step_flagged",
        "Select flagged + step-change points"
      )
    ),

    #========================
    #---QA Action Controls---
    #========================

    h3("Select Data"),

    checkboxInput(
      "enable_flag_select",
      "Enable QA flag selection",
      value = FALSE
    ),

    selectInput(
      "qa_flag_select",
      "Select points by QA flag:",
      choices = c(
        "Good",
        "Suspect",
        "Bad",
        "Missing"
      ),
      selected = "Suspect"
    ),

    actionButton(
      "select_qa_flag",
      "Select points with this QA flag"
    ),

    #Slider bar
    checkboxInput(
      "enable_value_select",
      "Enable value-range selection",
      value = FALSE
    ),

    sliderInput(
      "value_select",
      "Select points by value:",
      min = 0,
      max = 100,
      value = c(0,100),
      step = 0.1
    ),

    # Date/time selection
    checkboxInput(
      "enable_datetime_select",
      "Enable date/time selection",
      value = FALSE
    ),

    shinyWidgets::airDatepickerInput(
      "datetime_start",
      "Start date/time:",
      value = NULL,
      timepicker = TRUE,
      dateFormat = "yyyy-MM-dd",
      autoClose = TRUE
    ),

    shinyWidgets::airDatepickerInput(
      "datetime_end",
      "End date/time:",
      value = NULL,
      timepicker = TRUE,
      dateFormat = "yyyy-MM-dd",
      autoClose = TRUE
    ),

    actionButton(
      "select_datetime",
      "Select date/time range"
    ),

    # ---Action selector - these are actions to perform on selected points: flag, correct, set NA
    h3("QA Actions"),

    radioButtons(
      "action", "Action:",
      choices = c(
        "Flag only"             = "flag",
        "Correct value"         = "correct",
        "Set NA (Missing)"      = "set_na"
      ),
      selected = "flag",
      inline   = TRUE          #naming format: label=internal_value; server logic reacts to input$action using these internal values
    ),

    # ---Show QC flag selector/choice ONLY for "flag" and "correct" (Good, Suspect, Bad only)---
    conditionalPanel(
      condition = "input.action == 'flag' || input.action == 'correct'",           # =='correct' means condition evaluated in JavaScript, not R.So use input.action, not input$action
      selectInput("flag", "Set Flag:", choices = c("Good", "Suspect", "Bad"))
    ),

    #--- Corrected numeric value — shows only for "correct"---
    conditionalPanel(
      condition = "input.action == 'correct'",
      numericInput("new_value", "Correct Value:", value = NA_real_)        #conditional so field only appears when user selects "Correct Value" action
    ),

    #----Reason field — always shown; required for NA & Correct----
    tagList(
      tags$label("Reason:"),
      textInput("reason", NULL, ""),
      tags$small(style = "color:#666;",
                 "Required for 'Set NA' and 'Correct value'.")        #tagList groups elements w/o adding extra layout containers
    ),

    # ----Button creation for applying actions, clearing selection, undo, etc.----
    actionButton("apply", "Apply QA Change"),
    actionButton("clear_sel", "Clear Selection"),
    actionButton("undo", "Undo Last Action", icon = icon("undo")),  # add undo last action button and uses FontAwesome icons
    #checkboxInput("accumulate","Add to existing selection", TRUE),  # TRUE means new selections added to old ones


    #=========================================================
    # SAVE / DOWNLOAD CONTROLS **server portion needs to define downloadHandler() for each control for these to work
    # =========================================================

    actionButton("save_all", "Save QA Data"),
    downloadButton("download_corrected", "Download Corrected Data"),
    downloadButton("download_log",       "Download QA Log"),


    # =========================================================
    # QA OUTPUT LOG & COUNT DISPLAYS
    # =========================================================

    DT::dataTableOutput("qa_log_table"),      # interactive, sortable table of logged QA actions
    textOutput("qa_log_count")                  # shows number of QA actions applied
  )

  #========================
  #=========SERVER=========
  #========================

  #server logic:
  # 1) Reactive values (rv) - is where working memory is, 2) Parameter selector reactive (active_param), 3) Data‑processing reactives (step_cols, etc.)
  # 4) Observers (these respond to UI events:; step observer, apply observer, Undo observer, Clear selection observer, step-change "select flagged" button,
  # 5) Outputs (tables, plots, etc.)

  #start/initialize server ("rv" only active in this server portion!! If try to do something outside of bottom bracket, won't work and will lead to your demise)
  server <- function(input, output, session) {


    ###### A. INITIALIZATION & CORE REACTIVES#############

    #===============================================================================================
    #A.1 Reactive values storage (live working dataset in the app)-i.e., creates a centralized state object for the entire Shiny session
    #================================================================
    #stores the live working dataset and user selections i.e.,  stores the app’s "live state”; need to have so is update-able, etc.
    #is not a dataframe, or a list, each session gets its own, saved when "Save All" writes files
    #
    # rv$data     - main working dataset (raw + corrected + flags); this is what plots and selection tools visualize
    # rv$qa_log   - audit trail table, row by row QA actions appended over time (eventually exported) - pairs with UI qa_log_table
    # rv$selected_ids - current Plotly-selected point IDs    #for both time_plot and map_plot, these are the IDs selected with Plotly
    # rv$use_range_colors - plotting logic: TRUE → use range-based threshold flags; FALSE → use QC flag_plot values; pairs with "Apply Range Thresholds" in UI
    # --------------------------------------------------------

    rv <- reactiveValues(
      data             = df %>% dplyr::mutate(across(everything(), ~ .)),
      qa_log           = qa_log,
      qa_log_saved_n = nrow(qa_log),
      selected_ids     = integer(),
      reviewer         = NULL,
      reviewer_source  = NULL,
      undo_snapshot    = NULL,
      use_range_colors = FALSE
    )

    #testing in here-will also remove
    observe({
      print("APPLY VALUE:")
      print(input$apply)
    })


    #================================
    #A.2 Dynamic column scaffolding
    #================================

    #-------------------------------
    #A.2a ensure_dynamic_cols()
    #-------------------------------
    # inside server b/c observer calls it
    # initializes dynamic columns if don't exist yet for selected parameter
    # need in here b/c if don't, will glitch error b/c columns won't exist and will get "Error in : Variable 'CDOM_rfu_c6_corrected' not found in data"

    observeEvent(input$parameter, {          #pairs with UI selectInput, numericInput, plot outputs, and QC action radioButtons
      AP <- active_param()
      ensure_dynamic_cols(AP)
    })

    ensure_dynamic_cols <- function(AP) {

      # Make corrected column
      if (!AP$corrected_col %in% names(rv$data)) {
        rv$data[[AP$corrected_col]] <- suppressWarnings(
          as.numeric(rv$data[[AP$col]])
        )
      }

      # Internal QC flag column
      if (!AP$flag_label %in% names(rv$data)) {
        rv$data[[AP$flag_label]] <- NA_character_
      }

      # Plot-safe/export-safe QC flag column
      if (!AP$flag_plot %in% names(rv$data)) {
        rv$data[[AP$flag_plot]] <- factor(
          "Missing",
          levels = names(flag_colors)
        )
      }
      # Manual QA override column
      # NA means there is no manual override, so the automatic
      # range-based flag should be used.
      if (!AP$manual_flag_col %in% names(rv$data)) {

        rv$data[[AP$manual_flag_col]] <- NA_character_

      }

      # Range diagnostics
      if (!AP$range_value_col %in% names(rv$data)) {
        rv$data[[AP$range_value_col]] <- NA_real_
      }
      if (!AP$range_flag_col %in% names(rv$data)) {
        rv$data[[AP$range_flag_col]] <- NA_character_
      }
      if (!AP$range_label_col %in% names(rv$data)) {
        rv$data[[AP$range_label_col]] <- factor(
          "Missing",
          levels = names(flag_colors)
        )
      }

    }

    #----------------------------
    #A.2b - observer that calls 2a
    #----------------------------
    #when a new parameter is selected in UI

    observeEvent(input$apply_range_thresholds, {   #this makes sure "Apply Range Thresholds" functions and colors get updated then
      rv$use_range_colors <- TRUE
      rv$data <- rv$data   # triggers re-render
    })


    #===========================
    #A.3 new helper functions - inside server b/c depends on reactives (lookup_results()) (similar was previously outside server)
    #===========================

    #------------------
    #A.3a - find param
    #------------------
    #given a Standard_Parmeter (e.g., DO_mgL), find corresponding raw data column in rv$data by matching synonyms;
    #indirectly tied to selectInput ("parameter") or plot/QC action that works on "active parameter" in the UI
    #pairs with lookup_results() reactive, which uses instrument selector UI portion (:radioButtons: and :checkboxInput:)

    find_param_column <- function(df, standard_param) {

      # Get synonyms for this parameter FROM instrument-aware lookup
      synonyms <- lookup_results()$column_map[[standard_param]]

      # No synonyms found for this parameter
      if (is.null(synonyms)) {
        return(NA_character_)
      }

      # Normalize df column names for matching
      df_cols_lower <- tolower(names(df))

      # Find matches in the dataset columns
      hits <- names(df)[df_cols_lower %in% synonyms]

      # No match = parameter not present in this trip
      if (length(hits) == 0) {
        return(NA_character_)
      }

      # Return the first match (consistent with your original logic)
      hits[1]
    }


    #----------------------------------------------------------------------------
    #A.3b create final dataset to be exported when user clicks "Download Corrected Data"
    #-------------------------------------------------------------------------------
    # raw values+corrected values+plot-safe flags+core metadata columns (lat, long, timestamp, etc)
    # NO Diagnostic or internal flags

    build_corrected_subset <- function(df) {

      # Standard parameters known from lookup
      params <- unique(lookup_results()$param_lookup$Standard_Parameter)

      # Dynamic raw columns: depends on which exist in dataset
      raw_cols <- sapply(params, function(p) find_param_column(df, p))
      raw_cols <- raw_cols[!is.na(raw_cols)]

      # Corrected columns (if present)
      corrected_cols <- paste0(params, "_corrected")
      corrected_cols <- corrected_cols[corrected_cols %in% names(df)]

      # Exported, plot-safe flags
      flag_plot_cols <- paste0(params, "_flag_plot")
      flag_plot_cols <- flag_plot_cols[flag_plot_cols %in% names(df)]

      # Core metadata to keep (you may adjust)
      core_cols <- c(
        "id", "timestamp", "date",
        "time", "time..hh.mm.ss.",
        "name", "lon_dd", "lat_dd"
      )

      keep_cols <- unique(c(
        core_cols,
        raw_cols,
        corrected_cols,
        flag_plot_cols
      ))

      # Final subset with timestamp in local time
      df %>%
        dplyr::mutate(timestamp = with_tz(timestamp, "America/New_York")) %>%
        dplyr::select(dplyr::any_of(keep_cols))
    }


    ######### B. USER IDENTITY AND INSTRUMENT SELECTION #############

    #=========================================
    #B.1.  OBSERVER "R" Resolve reviewer identity
    #=========================================
    # does 1x per session, for documentation purposes
    observe({
      id <- get_reviewer_unified(session, prefer_full_name = TRUE)
      rv$reviewer        <- id$reviewer
      rv$reviewer_source <- id$source

      message(sprintf(
        "QA Reviewer (Shiny): %s [source=%s]",
        rv$reviewer, rv$reviewer_source
      ))
    })


    #===========================================================
    # B.2 Instrument selection and lookup chain
    #===========================================================

    #-------------------------------------------------------------
    # B.2a Resolve instrument selected + lookup_results reactives
    #--------------------------------------------------------------
    # builds parameter lookup tables based on instrument selection
    # pairs with radioButtons, checkboxInput from UI + indirectly with selectInput
    # b/c instrument choice affects parameters, synonym mapping, and default thresholds
    # C, below, depends on this so need it set up before

    # 1. Create the instrument vector from GUI controls
    instrument_selected <- reactive({
      primary <- input$primary_instrument

      if (isTRUE(input$use_C6)) {
        c(primary, "C6")
      } else {
        primary
      }
    })

    #---------------------------------
    # B.2b Assign active_instrument
    #---------------------------------
    active_instrument <- instrument_selected


    #--------------------------------
    # B.2c Run lookup table builder
    #--------------------------------
    lookup_results <- reactive({
      out <- prepare_lookup_tables(instrument_selected(), trip_lookup_raw)

      threshold_defaults <- setNames(
        lapply(unique(out$param_lookup$Standard_Parameter), function(p) {
          list(
            bad_min     = NA_real_,
            suspect_min = NA_real_,
            suspect_max = NA_real_,
            bad_max     = NA_real_
          )
        }),
        unique(out$param_lookup$Standard_Parameter)
      )

      list(
        trip_lookup        = out$trip_lookup,
        param_lookup       = out$param_lookup,
        column_map         = out$column_map,
        unit_map           = out$unit_map,
        dup_synonyms       = out$dup_synonyms,
        empty_sets         = out$empty_sets,
        threshold_defaults = threshold_defaults   # <-- NEW
      )
    })


    #-----------------------
    # B2.d Debug observers
    #-----------------------
    # Optional debug output - for during development, only
    observe({
      message("Selected instrument(s): ", paste(instrument_selected(), collapse = ", "))
      message("Parameters found: ", paste(lookup_results()$param_lookup$Standard_Parameter, collapse = ", "))
    })

    #some checks
    observe({
      cat("\n=== DEBUG lookup_results()$param_lookup ===\n")
      print(lookup_results()$param_lookup)
    })
    #end checks


    ######### C. OBSERVER AUTO-SELECT FIRST AVAILABLE STANDARD_PARAMETER IN THE DATASET #############
    #then auto-updates; pairs with selectInput("parameter") from UI

    #==========================
    # C.1 DEBUG PARAMETER LIST
    #==========================
    #ties to UI radioButtons and checkboxInput
    #ONCE APP IS STABLE, CAN BLANK THIS OUT
    observeEvent(active_instrument(), {
      req(lookup_results())

      cat("\n=== DEBUG param list ===\n")
      print(lookup_results()$param_lookup$Standard_Parameter)
    })

    #==========================
    # C.2 MAIN PARAMETER SELECTOR OBSERVER
    #==========================
    #ties to selectInput("parameter"...) in UI

    observeEvent(active_instrument(), {
      req(lookup_results())

      all_params <- sort(unique(lookup_results()$param_lookup$Standard_Parameter))

      available <- sapply(all_params, function(p) {
        col <- find_param_column(rv$data, p)
        !is.na(col)
      })

      if (any(available)) {
        updateSelectInput(
          session,
          "parameter",
          choices  = all_params[available],
          selected = all_params[which(available)[1]]
        )
      }
    })

    #================================
    # C.3 Active Parameter reactive
    #================================
    #HEART OF APP HERE - this is the central hub used by almost all parts of server
    #----Create/Resolve all relevant column names for the active parameter----
    # This function used throughout the server for naming raw, corrected, flag, flag_plot columns + internal diagnostic columns (for Plotly internal use)
    # NOTE: active_param() is used extensively inside server logic.
    # It "returns a list of canonical names used throughout QC".
    # ----------------------------------------------------------

    #define active param
    active_param <- reactive({
      req(input$parameter)

      std_param <- input$parameter
      raw_col   <- find_param_column(rv$data, std_param) # raw data column name in rv$data
      req(!is.na(raw_col))

      #troubleshooting time_plot and map_plot not showing up
      print(paste("std_param:", std_param, "| raw_col:", raw_col))

      # Build AP as an object
      AP <- list(
        std_param       = std_param,
        col             = raw_col,

        # corrected and flag columns (parameter-specific)
        corrected_col   = paste0(std_param, "_corrected"),
        flag_label      = paste0(std_param, "_flag"),
        flag_plot       = paste0(std_param, "_flag_plot"),

        # Manual QA override:
        # NA = no manual override; use automatic range flag
        manual_flag_col = paste0(std_param, "_manual_flag"),

        # step diagnostic columns-INTERNAL(not exported) columns
        step_raw_col    = paste0(std_param, "_raw_value"),
        step_delta_col  = paste0(std_param, "_delta_value"),
        step_flag_col   = paste0(std_param, "_step_flag"),
        step_label_col  = paste0(std_param, "_step_flag_label"),

        # Range diagnostics-INTERNAL(not exported) columns
        range_value_col = paste0(std_param, "_range_value"),
        range_flag_col  = paste0(std_param, "_range_flag"),
        range_label_col = paste0(std_param, "_range_flag_label")
      )

      ensure_dynamic_cols(AP)   #this syntax added as a fix to previous issues with maps not showing up

      return(AP)
    })


    #=======================================================
    # C.3a Reset selection + reset value-range selector
    #=======================================================
    observeEvent(input$parameter, {

      # A parameter change always starts with NO selected points
      rv$selected_ids <- integer()

      # Value-range selection must be explicitly re-enabled
      updateCheckboxInput(
        session,
        "enable_value_select",
        value = FALSE
      )
      # Date/time selection must also be explicitly re-enabled
      updateCheckboxInput(
        session,
        "enable_datetime_select",
        value = FALSE
      )

      AP <- active_param()

      rawvals <- suppressWarnings(
        as.numeric(rv$data[[AP$col]])
      )

      if (any(is.finite(rawvals))) {

        slider_min <- min(rawvals, na.rm = TRUE)
        slider_max <- max(rawvals, na.rm = TRUE)

        # Prevent this programmatic slider update from being
        # interpreted as a user selection
        freezeReactiveValue(input, "value_select")

        updateSliderInput(
          session,
          "value_select",
          min = slider_min,
          max = slider_max,
          value = c(slider_min, slider_max)
        )
      }
    })

    #=======================================================
    # C.3b Select points by value range using slider
    #=======================================================
    observeEvent(
      input$value_select,
      ignoreInit = TRUE,
      {

        # Do nothing unless value-range selection is enabled
        if (!isTRUE(input$enable_value_select)) {
          return(NULL)
        }

        AP <- active_param()

        rawvals <- suppressWarnings(
          as.numeric(rv$data[[AP$col]])
        )

        sel_ids <- rv$data$id[
          !is.na(rawvals) &
            rawvals >= input$value_select[1] &
            rawvals <= input$value_select[2]
        ]

        # A new slider selection REPLACES the old selection
        rv$selected_ids <- sel_ids

        showNotification(
          paste(
            "Selected",
            length(sel_ids),
            "points by value range."
          ),
          type = "message"
        )
      }
    )

    #=======================================================
    # C.3c Enabling/disabling value-range selection
    #=======================================================
    observeEvent(
      input$enable_value_select,
      {

        # If the user turns value selection OFF,
        # immediately clear any slider-based selection.
        if (!isTRUE(input$enable_value_select)) {

          rv$selected_ids <- integer()

        }
      },
      ignoreInit = TRUE
    )
    #=======================================================
    # C.3d Select points by date/time range
    #=======================================================
    observeEvent(
      input$select_datetime,
      {

        # Do nothing unless date/time selection is enabled
        if (!isTRUE(input$enable_datetime_select)) {
          showNotification(
            "Enable date/time selection first.",
            type = "warning"
          )
          return(NULL)
        }

        # Make sure both date/time values exist
        req(input$datetime_start, input$datetime_end)

        # Convert inputs to POSIXct
        start_time <- as.POSIXct(
          input$datetime_start,
          tz = "UTC"
        )
        #Issue with time not matching up, selected times seems to be 4 hours off so sutract four hours to make up difference
        start_time <- start_time - 14400

        end_time <- as.POSIXct(
          input$datetime_end,
          tz = "UTC"
        )
        end_time <- end_time - 14400

        # Make sure start comes before end
        if (is.na(start_time) || is.na(end_time)) {
          showNotification(
            "Please provide both a valid start and end date/time.",
            type = "error"
          )
          return(NULL)
        }

        if (start_time > end_time) {
          showNotification(
            "Start date/time must be before end date/time.",
            type = "error"
          )
          return(NULL)
        }

        # Use the timestamp column directly
        timestamps <- rv$data$timestamp

        # Select IDs within the requested date/time range
        sel_ids <- rv$data$id[
          !is.na(timestamps) &
            timestamps >= start_time &
            timestamps <= end_time
        ]

        # A new date/time selection REPLACES the old selection
        rv$selected_ids <- sel_ids

        showNotification(
          paste0(
            "Selected ",
            length(sel_ids),
            " point",
            ifelse(length(sel_ids) == 1, "", "s"),
            " from ",
            format(start_time, "%Y-%m-%d %H:%M"),
            " to ",
            format(end_time, "%Y-%m-%d %H:%M"),
            "."
          ),
          type = "message"
        )
      }
    )
    #=======================================================
    # C.3e Enable/disable date/time selection
    #=======================================================
    observeEvent(
      input$enable_datetime_select,
      {

        # Turning date/time selection OFF clears the selection
        if (!isTRUE(input$enable_datetime_select)) {
          rv$selected_ids <- integer()
        }
      },
      ignoreInit = TRUE
    )
    #=======================================================
    # C.3f Select points by final QA flag
    #=======================================================

    observeEvent(
      input$select_qa_flag,
      {

        # Do nothing unless explicitly enabled
        if (!isTRUE(input$enable_flag_select)) {
          showNotification(
            "Enable QA flag selection first.",
            type = "warning"
          )
          return(NULL)
        }

        AP <- active_param()

        flag_col <- AP$flag_label

        req(flag_col %in% names(rv$data))

        selected_flag <- input$qa_flag_select

        # Use the FINAL QA flag:
        # manual override if present, otherwise automatic range flag.
        flag_values <- as.character(
          rv$data[[flag_col]]
        )

        sel_ids <- rv$data$id[
          !is.na(flag_values) &
            flag_values == selected_flag
        ]

        # Replace current selection
        rv$selected_ids <- sel_ids

        showNotification(
          paste0(
            "Selected ",
            length(sel_ids),
            " point",
            ifelse(length(sel_ids) == 1, "", "s"),
            " with QA flag: ",
            selected_flag,
            "."
          ),
          type = "message"
        )
      }
    )
    observeEvent(
      input$enable_flag_select,
      {

        if (!isTRUE(input$enable_flag_select)) {
          rv$selected_ids <- integer()
        }

      },
      ignoreInit = TRUE
    )
    #=======================================================
    # C.4 reset threshold UI when param or instrument changes
    #======================================================
    #matches with "numericInput("Range_bad_min",...)" etc. from UI

    observeEvent(
      list(input$parameter, active_instrument()),
      {
        req(lookup_results())
        AP <- isolate(active_param())

        # Find default values for this parameter
        thr <- lookup_results()$threshold_defaults[[AP$std_param]]
        if (is.null(thr)) return()

        # Update all four UI fields
        updateNumericInput(session, "range_bad_min",     value = thr$bad_min)
        updateNumericInput(session, "range_suspect_min", value = thr$suspect_min)
        updateNumericInput(session, "range_suspect_max", value = thr$suspect_max)
        updateNumericInput(session, "range_bad_max",     value = thr$bad_max)
      },
      ignoreInit = TRUE
    )

    #===========================================
    # C.5 Reset range coloring when param changes:
    #===========================================
    observeEvent(input$parameter, {
      rv$use_range_colors <- FALSE
    })


    ########## D. STEP-CHANGE DIAGNOSTICS #######################

    #=======================
    # D.0 define Step-change plot overlay color settings
    #=======================
    step_color <- "#A100FF"   # vivid purple
    #step_outline <- "#FFFFFF" # white outline for vis.


    #========================================================
    # D.1 Compute Step-change diagnostics for active parameter
    #=========================================================
    #creates columns for raw_value, delta_value, step_flag, step_flag_label

    step_cols <- reactive({

      # ---------------------------------------------------------
      # Only calculate diagnostics when explicitly enabled
      # ---------------------------------------------------------

      req(isTRUE(input$enable_step))

      AP <- active_param()
      req(AP$col)

      req(!is.null(input$step_threshold))
      req(input$step_threshold >= 0)

      # ---------------------------------------------------------
      # Sort data
      # ---------------------------------------------------------

      d <- rv$data %>%
        dplyr::arrange(id, timestamp)

      # ---------------------------------------------------------
      # Extract raw values as numeric
      # ---------------------------------------------------------

      raw_num <- suppressWarnings(
        as.numeric(d[[AP$col]])
      )

      # ---------------------------------------------------------
      # Calculate raw change
      # ---------------------------------------------------------

      previous_value <- dplyr::lag(raw_num)

      delta <- raw_num - previous_value

      # ---------------------------------------------------------
      # Calculate percent change
      # ---------------------------------------------------------

      percent_change <- ifelse(
        !is.na(previous_value) &
          previous_value != 0,

        (delta / previous_value) * 100,

        NA_real_
      )

      # ---------------------------------------------------------
      # Choose calculation method
      # ---------------------------------------------------------

      if (identical(input$step_method, "percent")) {

        change_value <- percent_change

      } else {

        change_value <- delta

      }

      # ---------------------------------------------------------
      # Flag step changes
      # ---------------------------------------------------------

      step_flag <- !is.na(change_value) &
        abs(change_value) > input$step_threshold

      # ---------------------------------------------------------
      # Build output
      # ---------------------------------------------------------

      out <- d %>%
        dplyr::mutate(

          !!AP$step_raw_col := raw_num,

          !!AP$step_delta_col := change_value,

          !!AP$step_flag_col := step_flag,

          !!AP$step_label_col := factor(
            dplyr::case_when(

              is.na(.data[[AP$step_raw_col]]) ~
                NA_character_,

              .data[[AP$step_flag_col]] ~
                "Step change > threshold",

              TRUE ~
                "Good"
            ),

            levels = c(
              "Good",
              "Step change > threshold"
            )
          )
        ) %>%

        dplyr::select(
          id,
          !!AP$step_raw_col,
          !!AP$step_delta_col,
          !!AP$step_flag_col,
          !!AP$step_label_col
        )

      out
    })

    #================================================
    # D.2 Merge step-change diagnostics into rv$data
    #=================================================
    # updates main working dataset (rv$data) each time user changes parameter, threshold, underlying data updates, when step_cols() re-evalutated
    # Plotly and QC logic need rv$data to contain these columns, not just in a separate dataframe
    observe({

      req(isTRUE(input$enable_step))

      AP <- active_param()
      sc <- step_cols()

      sc[[AP$step_raw_col]] <-
        as.numeric(sc[[AP$step_raw_col]])

      sc[[AP$step_delta_col]] <-
        as.numeric(sc[[AP$step_delta_col]])

      sc[[AP$step_flag_col]] <-
        as.logical(sc[[AP$step_flag_col]])

      rv$data <- rv$data %>%

        dplyr::select(
          -dplyr::any_of(c(
            AP$step_raw_col,
            AP$step_delta_col,
            AP$step_flag_col,
            AP$step_label_col
          ))
        ) %>%

        dplyr::left_join(
          sc,
          by = "id"
        )
    })

    #======================================
    # D.3 Select step-change flagged points
    #=====================================
    #interacts with Plotly selection
    #corresponding UI element is "actionButton("select_step_flagged", "Select flagged + step-change points")"

    #NEW CODE:
    #---------------------------------------------------
    # D.3a  Retreive active parameter + diagnostic column
    #----------------------------------------------------

    observeEvent(input$select_step_flagged, {

      AP <- active_param()
      flag_col <- AP$step_flag_col

      #--------------------------------------------------------------------
      # D.3b  Ensure the diagnostic column exists (step-change computed for this parameter)
      #---------------------------------------------------------------------------

      if (!flag_col %in% names(rv$data)) {
        showNotification("No step-change diagnostics available for this parameter.",
                         type = "warning")
        return(NULL)
      }

      #----------------------------------------
      # D.3c  Extract flagged IDs
      #----------------------------------------
      step_ids <- rv$data %>%
        dplyr::filter(.data[[flag_col]] == TRUE) %>%
        dplyr::pull(id)

      #----------------------------------------
      # D.3d  Update rv$selected_ids (accumulate vs replace)
      #----------------------------------------

      rv$selected_ids <- step_ids

      #--------------------
      # D.3e  Notify user
      #--------------------

      showNotification(
        paste0("Selected ", length(step_ids), " step-change flagged point",
               ifelse(length(step_ids) == 1, "", "s"), "."),
        type = "message"
      )
    })

    ########### E. RANGE-BASED THRESHOLD DIAGNOSTICS ##########

    #=====================================================
    # E.1 Subset data to selected parameter & trip = "slice"
    #=====================================================
    range_dat <- reactive({
      req(rv$data, input$parameter)
      rv$data %>%
        dplyr::mutate(Value = suppressWarnings(as.numeric(.data[[active_param()$col]])),
                      Date   = timestamp) %>%
        # no "trip" filtering – whole dataset is one trip so no trip column or anything like that
        dplyr::select(id, Date, Value)
    })

    #========================================================
    # E.2 Collect threshold values from UI (namespaced inputs)
    #=========================================================
    # has "range" in here b/c is absolute values, different from "step"
    # retrieves 4 numeric threshold values from UI (bad, suspect, min and max)

    range_thresholds <- reactive({                  #these values are defined/reset up in C.4
      list(
        bad_min     = input$range_bad_min,
        suspect_min = input$range_suspect_min,
        suspect_max = input$range_suspect_max,
        bad_max     = input$range_bad_max
      )
    })


    #===============================================================
    # E.3 Apply dynamic threshold function to the sliced dataset
    #==============================================================
    #classifies each data point into Good, Suspect, Bad based on E.1 range and E.2 threshold list

    range_flagged <- reactive({
      req(range_dat())
      dat <- range_dat()

      th <- range_thresholds()

      dat$range_flag <- flag_by_dynamic_thresholds(
        x           = dat$Value,
        bad_min     = th$bad_min,
        suspect_min = th$suspect_min,
        suspect_max = th$suspect_max,
        bad_max     = th$bad_max
      )

      dat
    })

    #===========================================
    # E.4 Output visualizations and Table
    #============================================
    #colors defined in beginning outside of ShinyApp before data even read in

    #-------------------------------------------
    # E.4a Time-series output with flagged points
    #--------------------------------------------
    #this is not actually time_plot but is auxiliary diagnostic

    output$range_ts <- renderPlot({
      req(range_flagged())
      dat <- range_flagged()

      plot(dat$Date, dat$Value,
           pch = 16, col = "black",
           xlab = "Date", ylab = "Value",
           main = paste("Dynamic Range Thresholds —", input$parameter))

      # Highlight "Bad"
      bad_pts <- dat %>% dplyr::filter(range_flag == "Bad")
      points(bad_pts$Date, bad_pts$Value,
             col = flag_colors["Bad"], pch = 16, cex = 1.3)

      # Highlight "Suspect"
      sus_pts <- dat %>% dplyr::filter(range_flag == "Suspect")
      points(sus_pts$Date, sus_pts$Value,
             col = flag_colors["Suspect"], pch = 16, cex = 1.3)
    })

    #--------------------------------------
    # E.4b Histogram with threshold lines
    #--------------------------------------
    #displays dist'n of values and overlays thresholds to be tuned interactively

    output$range_hist <- renderPlot({
      req(range_flagged(), range_thresholds())
      dat  <- range_flagged()
      th   <- range_thresholds()

      hist(dat$Value,
           breaks = 40, col = "grey85", border = "white",
           main = paste("Value Distribution:", input$parameter),
           xlab = "Value")

      # BAD limits
      if (!is.na(th$bad_min))
        abline(v = th$bad_min, col = flag_colors["Bad"], lwd = 2)
      if (!is.na(th$bad_max))
        abline(v = th$bad_max, col = flag_colors["Bad"], lwd = 2)

      # SUSPECT limits
      if (!is.na(th$suspect_min))
        abline(v = th$suspect_min, col = flag_colors["Suspect"], lwd = 2, lty = 2)
      if (!is.na(th$suspect_max))
        abline(v = th$suspect_max, col = flag_colors["Suspect"], lwd = 2, lty = 2)
    })

    #------------------------------------------------------
    # E.4c Table of flagged values ("Suspect" or "Bad" only)
    #-----------------------------------------------------
    output$range_flag_table <- DT::renderDT({
      req(range_flagged())
      range_flagged() %>%
        dplyr::filter(range_flag != "Good") %>%
        DT::datatable(
          options = list(pageLength = 12, scrollX = TRUE),
          rownames = FALSE
        )
    })

    #-----------------------------
    # E.4d Download flagged values
    #------------------------
    #exports only values that failed thresholds

    output$download_range_flags <- downloadHandler(
      filename = function() {
        paste0("range_flags_", input$parameter, ".csv")
      },
      content = function(file) {
        write.csv(
          range_flagged() %>% dplyr::filter(range_flag != "Good"),
          file, row.names = FALSE
        )
      }
    )

    #=================================================
    # E.5 Merge range diagnostics + automatic QA flags
    #=================================================

    observeEvent(
      list(
        input$parameter,
        input$range_bad_min,
        input$range_suspect_min,
        input$range_suspect_max,
        input$range_bad_max
      ),
      {
        AP <- active_param()
        req(AP$col)

        rd <- range_flagged()
        req(nrow(rd) > 0)

        # Dynamic column names
        range_value_col  <- AP$range_value_col
        range_flag_col   <- AP$range_flag_col
        range_label_col  <- AP$range_label_col
        manual_flag_col  <- AP$manual_flag_col
        flag_label_col   <- AP$flag_label
        flag_plot_col    <- AP$flag_plot

        # --------------------------------------------------
        # Build range diagnostic columns
        # --------------------------------------------------

        new_cols <- rd %>%
          dplyr::mutate(
            !!range_value_col := as.numeric(Value),

            !!range_flag_col := as.character(range_flag),

            !!range_label_col := factor(
              range_flag,
              levels = names(flag_colors)
            )
          ) %>%
          dplyr::select(
            id,
            !!rlang::sym(range_value_col),
            !!rlang::sym(range_flag_col),
            !!rlang::sym(range_label_col)
          )

        # --------------------------------------------------
        # Merge range diagnostics into working dataset
        # --------------------------------------------------

        rv$data <- rv$data %>%
          dplyr::select(
            -dplyr::any_of(c(
              range_value_col,
              range_flag_col,
              range_label_col
            ))
          ) %>%
          dplyr::left_join(new_cols, by = "id")

        # --------------------------------------------------
        # Apply automatic range flag ONLY when there is
        # no manual QA override.
        # --------------------------------------------------

        automatic_flag <- as.character(
          rv$data[[range_flag_col]]
        )

        manual_flag <- as.character(
          rv$data[[manual_flag_col]]
        )

        # Manual override takes precedence.
        final_flag <- ifelse(
          !is.na(manual_flag) & nzchar(manual_flag),
          manual_flag,
          automatic_flag
        )

        # Write final QA designation to internal flag
        rv$data[[flag_label_col]] <- final_flag

        # Write final QA designation to plot/export flag
        rv$data[[flag_plot_col]] <- factor(
          final_flag,
          levels = names(flag_colors)
        )

        # Force reactive update
        rv$data <- rv$data
      }
    )


    ######## F. QC Actions################

    #======================================================
    # F.1 Apply QA Change (multiple substeps to complete logic)
    #========================================================

    #--------------------------
    # F.1a Validate user input
    #--------------------------
    #prevents incomplete/accidental modifications and enforces QA documentation "discipline"

    observeEvent(input$apply, {

      AP <- active_param()
      flag_plot_col <- AP$flag_plot_col
      selected_ids <- rv$selected_ids

      print("========== APPLY QA DEBUG ==========")
      print(paste("Selected IDs:", length(selected_ids)))
      print(selected_ids)
      print(paste("QA log BEFORE:", nrow(rv$qa_log)))

      # If nothing is selected, stop
      if (length(selected_ids) == 0) {
        showNotification("No points selected.", type = "error")
        return(NULL)
      }

      # Force re-render
      rv$data <- rv$data

      # Reason field is required for set NA and Correct
      if (input$action %in% c("correct", "set_na") &&
          !nzchar(input$reason)) {
        showNotification(
          "A reason is required for 'Correct' or 'Set NA'.",
          type = "error"
        )
        return(NULL)
      }

      #----------------------------
      # F.1b Creates undo snapshot
      #----------------------------
      # "structured" snapshot only for specific IDs and params-not snapshot of entire dataset; done before saves

      rv$undo_snapshot <- list(
        ids            = selected_ids,
        ap             = AP,

        prev_corrected =
          rv$data[[AP$corrected_col]][
            rv$data$id %in% selected_ids
          ],

        prev_flags =
          rv$data[[AP$flag_label]][
            rv$data$id %in% selected_ids
          ],

        prev_manual_flags =
          rv$data[[AP$manual_flag_col]][
            rv$data$id %in% selected_ids
          ],

        log_n_before = nrow(rv$qa_log)
      )

      #------------------------------------
      # F.1c  Apply action row-by-row; THIS IS A LOOP COMPOSED OF STEPS i-v and iv has a-e subdivisions)
      #---------------------------------

      # F.1c (i)-Prep references and data object
      #Shortcut column references - these are necessary for Apply logic; are not part of snapshot; goes after snapshot b/c changes rv$data

      df <- rv$data

      flag_col         <- AP$flag_label
      plot_col         <- AP$flag_plot
      corr_col         <- AP$corrected_col
      manual_flag_col  <- AP$manual_flag_col
      param            <- AP$std_param


      # F.1c (ii)  Loop through all selected IDs and apply action
      #  #ID row being edited; iterates over selected datapts (id'd by unique IDs); skips IDs not found

      for (id_val in selected_ids) {

        row_idx <- which(df$id == id_val)
        if (length(row_idx) == 0) next

        # F.1c (iii) Initialize editable fields, extract raw + existing values for use in log

        raw_value      <- df[[AP$col]][row_idx]          # raw numeric value
        old_corrected  <- df[[corr_col]][row_idx]        # previous corrected

        #initialize editable fields (may be modified)
        new_corrected  <- old_corrected
        final_flag     <- df[[flag_col]][row_idx]
        final_plot     <- df[[plot_col]][row_idx]


        # F.1c (iv) Apply chosen action

        # F.1c (iv.a) FLAG ONLY (does not change corrected)-------------
        if (input$action == "flag") {

          user_flag <- input$flag

          final_flag <- user_flag
          final_plot <- user_flag

          # Record explicit manual QA decision
          df[[manual_flag_col]][row_idx] <- user_flag

          # corrected value unchanged
        }

        #F.1 c (iv.b) CORRECT VALUE----------------------
        if (input$action == "correct") {

          req(!is.na(input$new_value))

          user_flag     <- input$flag
          new_corrected <- input$new_value

          final_flag <- user_flag
          final_plot <- user_flag

          # Correct + flag is also an explicit manual QA decision
          df[[manual_flag_col]][row_idx] <- user_flag
        }

        #F.1 c (iv.c) Correct Value--------------------------
        if (input$action == "set_na") {

          new_corrected <- NA_real_

          final_flag <- "Missing"
          final_plot <- "Missing"

          # Explicit manual designation
          df[[manual_flag_col]][row_idx] <- "Missing"
        }


        # F.1 c (v) Apply updates to dataframe------
        df[[corr_col]][row_idx]  <- new_corrected
        df[[flag_col]][row_idx]  <- final_flag
        df[[plot_col]][row_idx]  <- final_plot


        #---------------------------------
        #F.1d Capture diagnostic context
        #---------------------------------
        #records metadata; does not modify the data; is used only for log creation in F.1e
        #key is that it captures technical conditions when the action was made (e.g., bad and suspect ranges)
        #--------------------------
        # Extract diagnostic context safely
        #--------------------------

        # Step-change diagnostics
        if (!is.null(AP$step_flag_col) &&
            AP$step_flag_col %in% names(df) &&
            length(df[[AP$step_flag_col]]) >= row_idx) {
          step_flag_val <- df[[AP$step_flag_col]][row_idx]
        } else {
          step_flag_val <- NA_character_
        }

        if (!is.null(AP$step_delta_col) &&
            AP$step_delta_col %in% names(df) &&
            length(df[[AP$step_delta_col]]) >= row_idx) {
          step_delta_val <- df[[AP$step_delta_col]][row_idx]
        } else {
          step_delta_val <- NA_real_
        }

        # Range diagnostics
        if (!is.null(AP$range_flag_col) &&
            AP$range_flag_col %in% names(df) &&
            length(df[[AP$range_flag_col]]) >= row_idx) {
          range_flag_val <- df[[AP$range_flag_col]][row_idx]
        } else {
          range_flag_val <- NA_character_
        }

        if (!is.null(AP$range_value_col) &&
            AP$range_value_col %in% names(df) &&
            length(df[[AP$range_value_col]]) >= row_idx) {
          range_value_val <- df[[AP$range_value_col]][row_idx]
        } else {
          range_value_val <- NA_real_
        }

        diag_context <- sprintf(
          " | step_flag=%s step_delta=%s | range_flag=%s range_value=%s",
          step_flag_val,
          ifelse(is.na(step_delta_val), "NA", sprintf("%.4f", step_delta_val)),
          range_flag_val,
          ifelse(is.na(range_value_val), "NA", sprintf("%.4f", range_value_val))
        )

        # If user supplied a reason, append diagnostics.
        # If blank (e.g., flag-only), diagnostics still get stored.

        #final_reason <- paste0(input$reason, diag_context) #this puts reason typed in as well as thresholds into same cell; separates reason from context

        final_reason <- input$reason           # ONLY what user typed
        final_context <- diag_context          # NEW variable storing diagnostics under which decision was made


        #--------------------------
        # F.1e Append QA LOG entry - THIS IS THE AUDIT TRAIL
        #--------------------------
        print("========== QA VALUE LENGTH DEBUG ==========")

        print(paste("id_val length:", length(id_val)))
        print(paste("param length:", length(param)))
        print(paste("raw_value length:", length(raw_value)))
        print(paste("input$action length:", length(input$action)))
        print(paste("final_flag length:", length(final_flag)))
        print(paste("new_corrected length:", length(new_corrected)))
        print(paste("final_reason length:", length(final_reason)))
        print(paste("final_context length:", length(final_context)))
        print(paste("reviewer length:", length(rv$reviewer)))

        print("Actual structures:")
        str(id_val)
        str(param)
        str(raw_value)
        str(input$action)
        str(final_flag)
        str(new_corrected)
        str(final_reason)
        str(final_context)
        str(rv$reviewer)

        print("============================================")


        rv$qa_log <- append_qa_log(
          rv$qa_log,
          id              = id_val,
          parameter       = param,
          original_value  = raw_value,
          action          = input$action,
          flag            = final_flag,
          corrected_value = new_corrected,
          reason          = final_reason,
          context         = final_context,
          reviewer        = rv$reviewer
        )

      } # end for loop; runs once for each selected ID


      #-----------------------------
      #F.1f Finalize updates
      #-----------------------------

      # F.1f (i)  Actual Update reactive dataset-----------
      rv$data <- df  #writes all changes in F.1c into main wordking dataset
      rv$data <- rv$data   #force reactive invalidation so map_plot fully redraws b/c scattermapbox caches traces (i.e., when a point is deleted, it is removed from map)


      # F.1f (ii)  Clear selection after every applied change-------------
      rv$selected_ids <- integer()

      # F.1 f (iii) Notify user of completion
      #critical b/c confirms acion succeeded and shows how many pts edited
      showNotification(
        sprintf("Applied '%s' to %d points.", input$action, length(selected_ids)),
        type = "message"
      )

      # F.1 f (iv) Debug prints (optional)
      print(rv$selected_ids) #troubleshoot data format issue/data column name
      print("INSIDE APPLY") #troubleshoot data format issue/data column name
    }) #these brackets indicate end of observeEvent


    #=======================
    # F.2 Undo last action
    #========================

    observeEvent(input$undo, {

      #----------------------------
      #F.2a Require a stored snapshot - created in F.1b
      #-----------------------------
      if (is.null(rv$undo_snapshot)) {
        showNotification("Nothing to undo.", type = "warning")
        return(NULL)
      }

      #----------------------------
      #F.2b Extract snapshot fields
      #-----------------------------
      snap <- rv$undo_snapshot

      ids            <- snap$ids               # list of IDs modified
      AP             <- snap$ap                # active parameter structure (AP)
      prev_corrected <- snap$prev_corrected    # original corrected values
      prev_flags     <- snap$prev_flags        # original internal QC flags
      prev_manual_flags <- snap$prev_manual_flags
      log_n_before   <- snap$log_n_before      # number of QA log rows before this action

      df <- rv$data


      #--------------------------------------
      #F.2c Restore corrected values and QC flags
      #---------------------------------------
      # Apply restoration for each affected ID edited in F1.c
      for (i in seq_along(ids)) {
        id_val <- ids[i]
        idx    <- which(df$id == id_val)
        if (length(idx) == 0) next

        df[[AP$corrected_col]][idx] <- prev_corrected[i]
        df[[AP$flag_label]][idx]    <- prev_flags[i]

        # Restore manual override state
        df[[AP$manual_flag_col]][idx] <- prev_manual_flags[i]

        # Restore plot-safe flag based on internal flag
        if (!is.null(prev_flags[i]) && !is.na(prev_flags[i])) {
          df[[AP$flag_plot]][idx] <- factor(
            prev_flags[i],
            levels = names(flag_colors)
          )
        } else {
          df[[AP$flag_plot]][idx] <- factor(
            "Missing",
            levels = names(flag_colors)
          )
        }
      }


      #-------------------------------
      #F.2d Restore truncated QA log
      #-------------------------------
      if (nrow(rv$qa_log) > log_n_before) {
        rv$qa_log <- rv$qa_log[seq_len(log_n_before), , drop = FALSE]      #any new entries added in last Apply are removed
      }


      #-------------------------------
      #F.2e Update dataset rv$data and clear snapshot
      #-------------------------------
      # Update dataset
      rv$data <- df

      # Clear undo snapshot
      rv$undo_snapshot <- NULL  #cannot undo an undo


      #------------------
      #F.2f Notify user
      #----------------------
      showNotification("Undo: last change reverted.", type = "message")
    })


    ######## G. SELECTION LOGIC - these are the observers################
    #Plotly click/lasso + enables switching params

    ###TIME PLOT OBSERVER##

    #================================
    # G.2 Time plot Selection Logic
    #=================================

    #----------------------------------------------------------
    # G.2a Time_plot: Click-to-select (single point)
    #----------------------------------------------------------
    #simple click collects a single point
    #"accumulate" checked adds to selection; otherwise replaces

    observeEvent(
      event_data("plotly_click", source = "time_plot", priority = "event"),
      {
        click <- event_data(
          "plotly_click",
          source = "time_plot",
          priority = "event"
        )

        req(click)

        clicked_id <- click$customdata

        req(!is.null(clicked_id))

        clicked_id <- as.integer(clicked_id)

        if (isTRUE(input$accumulate)) {

          rv$selected_ids <- unique(
            c(rv$selected_ids, clicked_id)
          )

        } else {

          rv$selected_ids <- clicked_id

        }
      }
    )

    #-------------------------------------------------
    # G.2b Time_plot: Lasso / box multipt selection
    #-------------------------------------------------
    observeEvent(
      event_data("plotly_selected", source = "time_plot", priority = "event"),
      {
        sel <- event_data(
          "plotly_selected",
          source = "time_plot",
          priority = "event"
        )

        req(sel)

        req("customdata" %in% names(sel))

        ids <- as.integer(sel$customdata)

        ids <- ids[!is.na(ids)]

        if (isTRUE(input$accumulate)) {

          rv$selected_ids <- unique(
            c(rv$selected_ids, ids)
          )

        } else {

          rv$selected_ids <- ids

        }
      }
    )


    #----------------------------------------------------------
    #  # G.2c Time_plot: Deselect event (click in empty space)
    #----------------------------------------------------------

    observeEvent(
      event_data("plotly_deselect", source = "time_plot", priority = "event"),
      {
        rv$selected_ids <- integer()
      }
    )


    #======================================
    # G.3 Map plot Selection Logic Observer
    #======================================

    #--------------------------------------------------
    #  # G.3a Map_plot: Click-to-select (single point)
    #--------------------------------------------------
    #same logic as time_plot

    observeEvent(
      event_data("plotly_click", source = "map_plot", priority = "event"),
      {
        click <- event_data(
          "plotly_click",
          source = "map_plot",
          priority = "event"
        )

        req(click)

        clicked_id <- click$customdata

        req(!is.null(clicked_id))

        clicked_id <- as.integer(clicked_id)

        if (isTRUE(input$accumulate)) {

          rv$selected_ids <- unique(
            c(rv$selected_ids, clicked_id)
          )

        } else {

          rv$selected_ids <- clicked_id

        }
      }
    )


    #----------------------------------------------------------
    #  G.3b Map_plot: Lasso / box multipt selection
    #----------------------------------------------------------
    observeEvent(
      event_data("plotly_selected", source = "map_plot", priority = "event"),
      {
        sel <- event_data(
          "plotly_selected",
          source = "map_plot",
          priority = "event"
        )

        req(sel)

        req("customdata" %in% names(sel))

        ids <- as.integer(sel$customdata)

        ids <- ids[!is.na(ids)]

        if (isTRUE(input$accumulate)) {

          rv$selected_ids <- unique(
            c(rv$selected_ids, ids)
          )

        } else {

          rv$selected_ids <- ids

        }
      }
    )


    #-------------------------------------------------------
    # G.3c Map_plot: Deselect event (click in empty space)
    #-------------------------------------------------------

    observeEvent(input$clear_sel, {

      rv$selected_ids <- integer()

      showNotification(
        "Selection cleared.",
        type = "message"
      )
    })


    #============================================
    # G.4 Clear selection - wipes out Part 5 selection
    #============================================
    #pairs with actionButton("clear_sel", "Clear_Selection") in UI

    observeEvent(input$clear_sel, {

      # Reset the selected_id list to empty
      rv$selected_ids <- integer()

      # Notify user
      showNotification("Selection cleared.", type = "message")
    })


    ######################## H - PLOT CONSTRUCTION #################################

    #====================================================
    # H.1 TIME SERIES PLOT (Corrected values + QC colors)
    #====================================================
    output$time_plot <- renderPlotly({

      #-----------------------------------------
      # H.1a Retrieve active parameter metadata
      #----------------------------------------

      # H.1a (i) Call active_param()
      AP <- active_param()    #retrieves list

      # H.1a (ii) validate required dynamic columns
      req(AP$col, AP$corrected_col, AP$flag_plot)  #ensures column names exist before plotting

      # H.1a (iii) Extract units safely
      unit_label <- lookup_results()$param_lookup %>%          # gets param lookup frame
        dplyr::filter(Standard_Parameter == AP$std_param) %>%  # filters one row=selected param
        dplyr::pull(Units) %>%                                 # extracts Units
        as.character() %>%                                     # converts pulls into chr string
        .[1]                                                   # first value ensures 1 unit string

      # H.1a (iv) Fallback missing units
      if (is.null(unit_label) || is.na(unit_label)) unit_label <- ""

      # H.1a (v) Assign parameter-specific column names
      d <- rv$data
      raw_col  <- AP$col
      corr_col <- AP$corrected_col

      # Step-change diagnostic columns
      step_flag_col  <- AP$step_flag_col
      step_delta_col <- AP$step_delta_col
      step_label_col <- AP$step_label_col

      # Make sure step diagnostic columns exist
      if (!(step_flag_col %in% names(d))) {
        d[[step_flag_col]] <- FALSE
      }

      if (!(step_delta_col %in% names(d))) {
        d[[step_delta_col]] <- NA_real_
      }


      # ============================================================
      # BEFORE range thresholds are applied:
      # all points are neutral/light gray
      # ============================================================

      if (!isTRUE(rv$use_range_colors)) {

        p <- plot_ly(
          data = d,
          x = ~timestamp,
          y = ~.data[[corr_col]],
          type = "scatter",
          mode = "markers",
          source = "time_plot",
          customdata = ~id,

          marker = list(
            size = 5,
            opacity = 0.7,
            color = "lightgray"
          ),

          text = ~paste0(
            "Basin: ", name,
            "<br>", AP$std_param, " raw: ",
            ifelse(
              is.na(.data[[raw_col]]),
              "NA",
              sprintf(
                "%.3f %s",
                as.numeric(.data[[raw_col]]),
                unit_label
              )
            ),
            "<br>", AP$std_param, " corrected: ",
            ifelse(
              is.na(.data[[corr_col]]),
              "NA",
              sprintf(
                "%.3f %s",
                as.numeric(.data[[corr_col]]),
                unit_label
              )
            ),
            "<br>ΔValue: ",
            ifelse(
              is.na(.data[[step_delta_col]]),
              "NA",
              sprintf(
                "%.3f %s",
                as.numeric(.data[[step_delta_col]]),
                unit_label
              )
            ),
            "<br>Time: ",
            format(timestamp, "%Y-%m-%d %H:%M %Z")
          ),

          hovertemplate = "%{text}<extra></extra>"
        )


        # ============================================================
        # AFTER range thresholds are applied:
        # Good / Suspect / Bad + Missing / Deleted
        # ============================================================

      } else {

        d <- d %>%
          dplyr::mutate(
            plot_flag = dplyr::case_when(
              as.character(.data[[AP$flag_plot]]) %in% c("Good", "Suspect", "Bad") ~
                as.character(.data[[AP$flag_plot]]),
              is.na(.data[[AP$corrected_col]]) ~ "Missing",
              TRUE ~ as.character(.data[[AP$range_label_col]])
            )
          )

        d$plot_flag <- factor(
          d$plot_flag,
          levels = names(flag_colors)
        )

        p <- plot_ly(
          data = d,
          x = ~timestamp,
          y = ~.data[[corr_col]],
          type = "scatter",
          mode = "markers",
          source = "time_plot",
          customdata = ~id,

          color = ~plot_flag,
          colors = flag_colors,

          marker = list(
            size = 5,
            opacity = 0.7
          ),

          text = ~paste0(
            "Basin: ", name,
            "<br>", AP$std_param, " raw: ",
            ifelse(
              is.na(.data[[raw_col]]),
              "NA",
              sprintf(
                "%.3f %s",
                as.numeric(.data[[raw_col]]),
                unit_label
              )
            ),
            "<br>", AP$std_param, " corrected: ",
            ifelse(
              is.na(.data[[corr_col]]),
              "NA",
              sprintf(
                "%.3f %s",
                as.numeric(.data[[corr_col]]),
                unit_label
              )
            ),
            "<br>ΔValue: ",
            ifelse(
              is.na(.data[[step_delta_col]]),
              "NA",
              sprintf(
                "%.3f %s",
                as.numeric(.data[[step_delta_col]]),
                unit_label
              )
            ),
            "<br>Time: ",
            format(timestamp, "%Y-%m-%d %H:%M %Z"),
            "<br>Flag: ",
            as.character(.data[[AP$range_label_col]])
          ),

          hovertemplate = "%{text}<extra></extra>"
        )
      }


      # ============================================================
      # DEBUG
      # ============================================================

      print("=== TIME PLOT DEBUG ===")
      print("use_range_colors:")
      print(rv$use_range_colors)

      if (isTRUE(rv$use_range_colors)) {
        print("Plot flag distribution:")
        print(table(d$plot_flag, useNA = "always"))
      }

      print("Step flag distribution:")
      print(table(d[[step_flag_col]], useNA = "always"))

      print("======================")

      #----------------------------------------------------------
      # H.1 d Step-change overlay (conditional-i.e., if enabled)
      #---------------------------------------------------------


      # ---- STEP-CHANGE OVERLAY (only one trace, shown in legend) ----
      if (isTRUE(input$enable_step) && isTRUE(input$highlight_step)) {

        d_step <- d %>% dplyr::filter(.data[[step_flag_col]] == TRUE)

        if (nrow(d_step) > 0) {

          p <- p %>% add_markers(
            data = d_step,
            x = ~timestamp,
            y = ~.data[[corr_col]],
            customdata = ~id,
            marker = list(
              size = 8,
              color = "#A100FF"
            ),
            name = "Step-change > threshold",     # ONE clear legend label
            legendgroup = "step_change",          # group for this trace only
            showlegend = TRUE,                    # <-- IMPORTANT
            inherit = FALSE,
            text = ~paste0(
              "Basin: ", name,
              "<br>", AP$std_param, " raw: ", ifelse(is.na(.data[[raw_col]]), "NA",
                                                     sprintf("%.3f %s", as.numeric(.data[[raw_col]]), unit_label)),
              "<br>ΔValue: ", ifelse(is.na(.data[[step_delta_col]]), "NA",
                                     sprintf("%.3f", as.numeric(.data[[step_delta_col]]))),
              "<br>QC: Step-change > ", sprintf("%.2f", input$step_threshold)
            ),
            hovertemplate = "%{text}<extra></extra>"
          )
        }
      }


      #-------------------------------
      # H.1e Layout and highlight behavior
      #-------------------------------

      # H.1e (i) Configure plot (layout, axes, drag)
      p <- p %>%
        layout(
          dragmode = "lasso",
          yaxis = list(
            title = paste0(AP$std_param, " (corrected, ", unit_label, ")"),
            hoverformat = ".3f"
          ),
          xaxis = list(
            type = "date",
            rangebreaks = list(list(pattern = "hour", bounds = c(18, 10))),
            tickformat = "%Y-%m-%d\n%H:%M",
            hoverformat = "%Y-%m-%d %H:%M %Z"
          )
        ) %>%

        # H.1e(ii) Define interactive highlighting rules
        highlight(
          on = "plotly_selected",
          off = "plotly_deselect",
          persistent = TRUE,
          opacityDim = 0.25,
          selected = attrs_selected(marker = list(size = 12, opacity = 1))
        )

      #-------------------------------
      # H.1f Register plotly events - TIME PLOT
      #-------------------------------
      # Register events AFTER all traces and layout
      #NOTE: you will  get a few "'time_plot' is not registered" messages early in the console. This surprisingly doesn't mean that the event
      #is not registered. We know it is b/c we can click on pts, it's from an earlier part where it is not registering while cycling and registers at the end!!!!
      #-------------------------------
      p <- event_register(p, "plotly_click")
      p <- event_register(p, "plotly_selected")
      p <- event_register(p, "plotly_deselect")

      # H.1g Return final plot
      p
    })


    #=========================================================
    # H.2 MAP PLOT (Corrected values + QC colors + Step overlay)
    #==========================================================

    #-------------------------------------
    # H.2a Retrieve active parameter data
    #-------------------------------------
    output$map_plot <- renderPlotly({

      # H.2a (i) Get active parameter + required columns---
      AP <- active_param()
      req(AP$col, AP$corrected_col, AP$flag_plot)
      req(rv$data)   # <-- line added to prevent point from still being plotted when it has been deleted; means map will re-render after delete

      # H.2a (ii) Safe unit extraction-get units for this parameter to include in hovertext
      #doing b/c stays consisent w/ timeplot and can throw errors b/c unit becomes "NULL"
      unit_label <- lookup_results()$param_lookup %>%
        dplyr::filter(Standard_Parameter == AP$std_param) %>%
        dplyr::pull(Units) %>%
        as.character() %>%
        .[1]
      if (is.null(unit_label) || is.na(unit_label)) unit_label <- ""  #fallback for missing units

      # H.2a (iii)  Assign dynamic column names (metadata)
      raw_col  <- AP$col
      corr_col <- AP$corrected_col

      # ALWAYS use the actual QA/QC flag for map colors.
      # Range thresholds are diagnostic only and should not overwrite
      # the user's Good/Suspect/Bad/Missing/Deleted status.
      flag_col <- AP$flag_plot


      #-----------------
      # H2.b Data prep
      #-----------------

      # H.2b (i) Create working dataset + numeric lat/lon
      d <- rv$data %>%
        dplyr::mutate(
          lon_num = suppressWarnings(as.numeric(lon_dd)),
          lat_num = suppressWarnings(as.numeric(lat_dd))
        )


      # H.2b (ii) Create map-specific color category
      #
      # Normal points:
      #   use the existing range/QC flag (Good/Suspect/Bad)
      #
      # Manually QA'd points:
      #   Missing -> Missing color
      #   Deleted -> Deleted color

      if (isTRUE(rv$use_range_colors)) {

        # Once thresholds have been applied:
        # range diagnostic colors, with Missing/Deleted overrides.

        d <- d %>%
          dplyr::mutate(
            map_flag = dplyr::case_when(
              as.character(.data[[AP$flag_plot]]) %in% c("Good", "Suspect", "Bad") ~
                as.character(.data[[AP$flag_plot]]),
              is.na(.data[[AP$corrected_col]]) ~ "Missing",
              TRUE ~ as.character(.data[[AP$range_label_col]])
            )
          )

        d$map_flag <- factor(
          d$map_flag,
          levels = names(flag_colors)
        )

      } else {

        # Before thresholds are applied:
        # no categorical coloring.

        d$map_flag <- "Unclassified"
      }

      # H.2b (iii) Debug print
      # ===== DEBUG MAP COLORS =====
      print("========== MAP COLOR DEBUG ==========")

      print("Range label column:")
      print(AP$range_label_col)

      print("Range diagnostic distribution:")
      print(table(
        as.character(d[[AP$range_label_col]]),
        useNA = "always"
      ))

      print("Manual flag distribution:")
      print(table(
        as.character(d[[AP$flag_plot]]),
        useNA = "always"
      ))

      print("Final map_flag distribution:")
      print(table(
        as.character(d$map_flag),
        useNA = "always"
      ))

      print("Flag colors:")
      print(flag_colors)

      print("=====================================")
      # ===== END DEBUG =====

      # H.2b(iv) Ensure/define step-change diagnostic columns exist (fallback if missing)
      step_flag_col  <- AP$step_flag_col
      step_delta_col <- AP$step_delta_col
      step_label_col <- AP$step_label_col

      #Ensure missing diagnostics don't break the plot
      if (!(step_flag_col %in% names(d))) d[[step_flag_col]] <- FALSE
      if (!(step_delta_col %in% names(d))) d[[step_delta_col]] <- NA_real_

      #troubleshooting if plots aren't appearing
      #print(attributes(d$timestamp))
      #print(paste("Rows in d:", nrow(d)))
      #print(AP)


      #-----------------------------------------
      # H.2c Base map scatter plot
      #-----------------------------------------
      #corrected values w/ QC plot flags + setup for step-change overlay
      if (isTRUE(rv$use_range_colors)) {
        p <- plot_ly(
          #p <- p %>% plot_ly(    #calls in previous time trace to match that behavior
          data = d,
          type = "scattermapbox",
          mode = "markers",
          lon = ~lon_num,
          lat = ~lat_num,
          source = "map_plot",
          customdata = ~id,

          #mapping
          color = ~map_flag,
          colors = flag_colors,       #same pallete in time_plot
          #inherit=FALSE,                #forces plotly to rebuild trace cleanly-makes it fully contained so traces can't be "inherited" from prior state when re-renders
          marker = list(size = 6, opacity = 0.8),
          text = ~paste0(
            "Basin: ", name,
            "<br>", AP$std_param, " raw: ",
            ifelse(is.na(.data[[raw_col]]), "NA",
                   sprintf("%.3f %s", as.numeric(.data[[raw_col]]), unit_label)),
            "<br>", AP$std_param, " corr: ",
            ifelse(is.na(.data[[corr_col]]), "NA",
                   sprintf("%.3f %s", as.numeric(.data[[corr_col]]), unit_label)),
            "<br>ΔValue: ",
            ifelse(is.na(.data[[step_delta_col]]), "NA",
                   sprintf("%.3f %s", as.numeric(.data[[step_delta_col]]), unit_label)),
            "<br>Time: ", format(timestamp, "%Y-%m-%d %H:%M"),
            "<br>Flag: ", as.character(.data[[AP$flag_plot]])
          ),
          hovertemplate = "%{text}<extra></extra>"
        )
      } else {
        p <- plot_ly(
          #p <- p %>% plot_ly(    #calls in previous time trace to match that behavior
          data = d,
          type = "scattermapbox",
          mode = "markers",
          lon = ~lon_num,
          lat = ~lat_num,
          source = "map_plot",
          customdata = ~id,
          marker = list(size = 6, opacity = 0.8, color="lightgray"),
          text = ~paste0(
            "Basin: ", name,
            "<br>", AP$std_param, " raw: ",
            ifelse(is.na(.data[[raw_col]]), "NA",
                   sprintf("%.3f %s", as.numeric(.data[[raw_col]]), unit_label)),
            "<br>", AP$std_param, " corr: ",
            ifelse(is.na(.data[[corr_col]]), "NA",
                   sprintf("%.3f %s", as.numeric(.data[[corr_col]]), unit_label)),
            "<br>ΔValue: ",
            ifelse(is.na(.data[[step_delta_col]]), "NA",
                   sprintf("%.3f %s", as.numeric(.data[[step_delta_col]]), unit_label)),
            "<br>Time: ", format(timestamp, "%Y-%m-%d %H:%M"),
            "<br>Flag: ", as.character(.data[[AP$flag_plot]])
          ),
          hovertemplate = "%{text}<extra></extra>"
        )
      }


      #-----------------------------------------
      # H.2d Step-change overlay (conditional)
      #-----------------------------------------
      if (isTRUE(input$highlight_step)) {

        d_step <- d %>%
          dplyr::filter(
            .data[[step_flag_col]] == TRUE,
            is.finite(lon_num),
            is.finite(lat_num)
          )

        if (nrow(d_step) > 0) {

          p <- p %>% add_trace(
            data = d_step,
            type = "scattermapbox",
            mode = "markers",
            lon = ~lon_num,
            lat = ~lat_num,
            customdata = ~id,
            marker = list(size = 8, color = "#A100FF"),
            name = "Step-change > threshold",
            legendgroup="step",
            showlegend = TRUE,

            text = ~paste0(
              "Basin: ", name,
              "<br>", AP$std_param, " raw: ",
              ifelse(is.na(.data[[raw_col]]), "NA",
                     sprintf("%.3f %s", as.numeric(.data[[raw_col]]), unit_label)),
              "<br>ΔValue: ",
              ifelse(is.na(.data[[step_delta_col]]), "NA",
                     sprintf("%.3f %s", as.numeric(.data[[step_delta_col]]), unit_label)),
              "<br>Time: ", format(timestamp, "%Y-%m-%d %H:%M"),
              "<br>QC: Step-change>", sprintf("%.2f", input$step_threshold)  #fun story need "Step-change>" not "Step-change >" to match
            ),
            hovertemplate = "%{text}<extra></extra>",
            inherit = FALSE  #prevents overlay caching
          )
        }
      }


      #-----------------------------------------
      # H.2e Map layout & interaction formatting
      #-----------------------------------------
      p <- p %>%
        layout(
          mapbox = list(
            style = "carto-positron",
            zoom = 9,
            center = list(
              lon = mean(d$lon_num, na.rm = TRUE),
              lat = mean(d$lat_num, na.rm = TRUE)
            )
          ),
          dragmode = "select"
        ) %>%
        highlight(
          on = "plotly_selected",
          off = "plotly_deselect",
          persistent = TRUE,
          opacityDim = 0.25,
          selected = attrs_selected(marker = list(size = 14, opacity = 1))
        )


      #-----------------------------------------
      # H.2f Register plotly events - MAP PLOT
      #-----------------------------------------

      p <- event_register(p, "plotly_click")
      p <- event_register(p, "plotly_selected")
      p <- event_register(p, "plotly_deselect")

      #-----------------------------------------
      # H.2g Return final plot
      #-----------------------------------------
      p
    })

    ########### I. Outputs ################

    #===========================
    # I.1 SELECTED INFO PANEL
    #============================
    #detailed info for selected points

    output$selected_info <- DT::renderDataTable({
      ids <- rv$selected_ids

      # If no selection, show nothing
      if (length(ids) == 0) {
        return(NULL)
      }

      AP <- active_param()
      req(AP$col, AP$corrected_col, AP$flag_label, AP$flag_plot)

      sel <- rv$data %>%
        dplyr::filter(id %in% ids) %>%
        dplyr::select(
          id,
          timestamp,
          name,
          lon_dd,
          lat_dd,
          !!rlang::sym(AP$col),
          !!rlang::sym(AP$corrected_col),
          !!rlang::sym(AP$range_flag_col),
          !!rlang::sym(AP$manual_flag_col),
          !!rlang::sym(AP$flag_label),
          !!rlang::sym(AP$flag_plot)
        ) %>%
        dplyr::rename(
          Raw = !!rlang::sym(AP$col),
          Corrected = !!rlang::sym(AP$corrected_col),
          Range_Flag = !!rlang::sym(AP$range_flag_col),
          Manual_Flag = !!rlang::sym(AP$manual_flag_col),
          QA_Flag = !!rlang::sym(AP$flag_label),
          QA_Flag_Plot = !!rlang::sym(AP$flag_plot)
        )

      DT::datatable(
        sel,
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollY = "250px",
          scrollX = TRUE,
          paging = TRUE
        )
      )
    })

    output$selection_count <- renderText({

      paste(
        "Selected points:",
        length(rv$selected_ids)
      )

    })


    #==================
    # I.2 QA Log table
    #====================
    # Displays full QA log with pagination

    output$qa_log_table <- DT::renderDataTable({
      req(nrow(rv$qa_log) >= 0)
      rv$qa_log
    }, options = list(pageLength = 6, scrollX = TRUE))

    #======================
    # I.3 QA LOG ROW COUNT
    #======================
    output$qa_log_count <- renderText({
      paste("QA log rows:", nrow(rv$qa_log))
    })

    ######### J. Helper functions (for K and L) ######

    #=============================
    # J.1 Build corrected dataset
    #============================
    # used for download (K) and Save-All (L)

    build_corrected_subset <- function(data) {

      AP <- active_param()

      data %>%
        dplyr::mutate(
          timestamp = with_tz(timestamp, "America/New_York"),
          flag      = .data[[AP$flag_plot]],
          corrected = .data[[AP$corrected_col]]
        ) %>%
        dplyr::select(
          id, timestamp, date, dplyr::any_of(c("time", "time..hh.mm.ss.")),
          name, lon_dd, lat_dd,
          !!AP$col,
          corrected,
          flag
        )
    }
    #================================================
    # J.2 Build full QA dataset
    #================================================
    # Creates a copy of the original-style dataset and
    # replaces only values that received QA actions.
    #
    # All un-QA'd values remain unchanged.
    # QA'd values are replaced by their final corrected value.
    #
    # This is different from build_corrected_subset(), which
    # creates a reduced output for the active parameter.

    build_full_qa_dataset <- function(data) {

      # Start with the complete working dataset
      qa_data <- data

      # If there are no QA actions, return the data unchanged
      if (nrow(rv$qa_log) == 0) {
        return(qa_data)
      }

      # Use the QA log as the authoritative record of which
      # ID + parameter combinations were actually QA'd.
      #
      # If a point was QA'd more than once, keep the most
      # recent action because that represents its final state.
      latest_qa <- rv$qa_log %>%
        dplyr::mutate(
          qa_order = dplyr::row_number()
        ) %>%
        dplyr::group_by(id, parameter) %>%
        dplyr::slice_max(
          order_by = qa_order,
          n = 1,
          with_ties = FALSE
        ) %>%
        dplyr::ungroup()

      # Process each parameter independently
      for (param in unique(latest_qa$parameter)) {

        # Find the actual raw-data column corresponding to
        # this Standard_Parameter.
        raw_col <- find_param_column(qa_data, param)

        # Skip if the parameter cannot be found
        if (is.na(raw_col) || !raw_col %in% names(qa_data)) {
          warning(
            paste(
              "Could not find raw data column for parameter:",
              param
            )
          )
          next
        }

        # Get QA records for this parameter
        param_qa <- latest_qa %>%
          dplyr::filter(parameter == param)

        # Match QA records to rows in the full dataset using ID
        match_idx <- match(param_qa$id, qa_data$id)

        # Only keep successfully matched IDs
        valid <- !is.na(match_idx)

        if (!any(valid)) next

        match_idx <- match_idx[valid]
        corrected_values <- param_qa$corrected_value[valid]

        # Replace the original/raw value with the final QA value
        qa_data[[raw_col]][match_idx] <- corrected_values
      }

      # Remove internal QA columns before saving.
      #
      # These columns are created by the dashboard and are not
      # part of the original Dataflow dataset.
      internal_cols <- grep(
        "(_corrected|_flag|_flag_plot|_manual_flag|_raw_value|_delta_value|_step_flag|_step_flag_label|_range_value|_range_flag|_range_flag_label)$",
        names(qa_data),
        value = TRUE
      )

      qa_data <- qa_data %>%
        dplyr::select(-dplyr::any_of(internal_cols))

      return(qa_data)
    }


    ######## K. Downloads (Client Side) ##############

    #===============================
    # K.1 — Download Corrected Subset
    #================================

    output$download_corrected <- downloadHandler(

      filename = function() {
        AP <- active_param()
        paste0("dataflow_corrected_", AP$std_param, "_",
               format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },

      content = function(file) {

        # Build the complete QA'd dataset using the QA log
        corrected_data <- build_full_qa_dataset(rv$data)

        readr::write_csv(
          corrected_data,
          file
        )
      }
    )

    #===============================
    # K.2 — Download QA Log
    #===============================

    output$download_log <- downloadHandler(

      filename = function() {
        paste0(
          "dataflow_QA_log_",
          yearmon,
          ".csv"
        )
      },

      content = function(file) {
        readr::write_csv(
          rv$qa_log,
          file
        )
      }
    )

    ######### L. Save All (Network Output) ############

    #================================================
    # L.1 — Save full QA dataset + QA log to network
    #================================================

    observeEvent(input$save_all, {

      print("========== SAVE DEBUG ==========")
      print(paste("QA log rows:", nrow(rv$qa_log)))
      print(paste("Previously saved rows:", rv$qa_log_saved_n))
      print("Last QA log record:")
      print(tail(rv$qa_log, 1))
      print("================================")

      # QA output folder
      base_dir <- file.path(
        fdir,
        "DF_FullDataSets",
        "QA datasets"
      )

      # Check that folder exists
      if (!dir.exists(base_dir)) {
        showNotification(
          paste("Output folder does not exist:", base_dir),
          type = "error",
          duration = 6
        )
        return(invisible(NULL))
      }

      #----------------------------------------------
      # Build exact filenames
      #----------------------------------------------

      qa_data_path <- file.path(
        base_dir,
        paste0(yearmon, "_qa.csv")
      )

      qa_log_path <- file.path(
        base_dir,
        paste0(yearmon, "_QA_log.csv")
      )

      #----------------------------------------------
      # Build and save QA dataset
      #----------------------------------------------

      tryCatch({

        # Create full dataset with only QA'd values replaced
        full_qa_data <- build_full_qa_dataset(rv$data)

        # Write full QA dataset
        readr::write_csv(
          full_qa_data,
          qa_data_path
        )

        #--------------------------------------------
        # Save QA LOG
        #--------------------------------------------

        if (file.exists(qa_log_path)) {

          # Existing QA log: append only NEW records
          if (nrow(rv$qa_log) > rv$qa_log_saved_n) {

            new_log_records <- rv$qa_log[
              seq.int(
                from = rv$qa_log_saved_n + 1,
                to   = nrow(rv$qa_log)
              ),
              ,
              drop = FALSE
            ]

            readr::write_csv(
              new_log_records,
              qa_log_path,
              append = TRUE,
              col_names = FALSE
            )

            rv$qa_log_saved_n <- nrow(rv$qa_log)

            log_message <- paste(
              "Appended",
              nrow(new_log_records),
              "new QA log record(s)."
            )

          } else {

            log_message <- "No new QA log records to append."

          }

        } else {

          # No QA log exists yet: create a new one
          readr::write_csv(
            rv$qa_log,
            qa_log_path
          )

          rv$qa_log_saved_n <- nrow(rv$qa_log)

          log_message <- paste(
            "Created new QA log with",
            nrow(rv$qa_log),
            "record(s)."
          )
        }
        #--------------------------------------------
        # Confirmation
        #--------------------------------------------

        showNotification(
          paste0(
            "QA data saved successfully:\n",
            qa_data_path,
            "\n\nQA log saved to:\n",
            qa_log_path
          ),
          type = "message",
          duration = 10
        )

      }, error = function(e) {

        showNotification(
          paste(
            "Save failed:",
            conditionMessage(e)
          ),
          type = "error",
          duration = 10
        )

      })
    })


  }  # END OF SERVER

  #cat("\n\n***** STARTING SHINY APP NOW *****\n\n") #IF TROUBLESHOOTING, indicates start of where "loop" is so can catch it

  #========================
  #=========Launch App======
  #========================
  shinyApp(ui, server)

}
