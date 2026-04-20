#'@name QAflags
#'@title Create plots that flag data based on values or change
#'@description Create plots that flag data based on values or change
#'@param yearmon numeric designation of survey date formatted as yyyymm
#'@param param name of parameter to be examined
#'@param bad_min values less than this are flagged as bad
#'@param bad_max values more than this are flagged as bad
#'@param sus_min values less than this are flagged as suspect
#'@param sus_max values more than this are flagged as suspect
#'@param step_threshold values flagged if change between susequent measurements are greater than this threshold
#'@param min_time data after this time (24 hour clock) are not plotted
#'@param max_time data before this time (24 hour clock) are not plotted
#'@param plots user specifies which plots to create (flag points, flag map, step flag points, step flag map), can specify 2 at a time
#'@param original.data specify whether to read in original streamcleaned data (TRUE) or already QAd data (FALSE)
#'@param fdir character file path to local data directory
#'@export
#'@import dplyr
#'@import zoo
#'@import rlang
#'@import lubridate
#'@import ggplot2
#'@import scales
#'@import purrr
#'@import tibble
#'@import plotly
#'@import htmltools
#'@examples \dontrun{
#'d <- QAflags(yearmon=202510, param="tal.pc.rfu", bad_min=0, bad_max=10, sus_min=0.01, sus_max=4, step_threshold=2,plots=c("flag points", "step flag points"), original.data=F)
#'}

QAflags <- function (yearmon, param, bad_min, bad_max, sus_min, sus_max, step_threshold,
                     min_time=17, max_time=9, plots, original.data=FALSE,
                     fdir = getOption("fdir")) {
  #Load data
  if(original.data == FALSE && file.exists(file.path(fdir, "DF_FullDataSets", "QA datasets", paste(yearmon, "j_qa.csv", sep = "")))){
    df <- read.csv(file.path(fdir, "DF_FullDataSets", "QA datasets", paste(yearmon, "j_qa.csv", sep = "")))
  } else {
    #Read original streamcleaned file if no QA yet
    fdir_fd <- file.path(fdir, "DF_FullDataSets")
    flist <- list.files(fdir_fd, include.dirs = T, full.names = T)
    flist <- flist[substring(basename(flist),1,6) == yearmon]
    df <- read.csv(flist, stringsAsFactors = FALSE)
  }

  #Define flag ranges
  FLAG_GOOD    <- "Good"
  FLAG_SUSPECT <- "Suspect"
  FLAG_BAD     <- "Bad"

  #Function to define ranges to be flagged
  range_rule <- function(min_value, max_value, flag_value) {
    list(
      flag = flag_value,
      condition = function(x) {
        is.na(x) | x < min_value | x > max_value
      }
    )
  }

  rules <- list(
    # SUSPECT: unusually low or high but possible
    range_rule(min_value = sus_min, max_value = sus_max, flag_value = FLAG_SUSPECT),  #CHANGE values based on time of year
    # BAD: physically impossible or sensor failure or unlikely high
    range_rule(min_value = bad_min, max_value = bad_max, flag_value = FLAG_BAD) #CHANGE values based on time of year
  )

  #Function that applies rules
  flag_by_rules <- function(x, rules) {
    flag <- rep(FLAG_GOOD, length(x))

    priority <- c(Good=1, Suspect=2, Bad=3)

    for (r in rules) {
      idx <- r$condition(x)
      # Only update if new flag has higher priority
      flag[idx] <- ifelse(
        priority[r$flag] > priority[flag[idx]],
        r$flag,
        flag[idx]
      )
    }
    flag
  }

  #Combine full date and time into one column
  #if keeping timestamp column in QA data could add check to see if already exists when reading in data set
  df <- df %>%
    mutate(timestamp=mdy_hms(paste(date, time..hh.mm.ss.), tz="UTC"))  #Doing this instead of using datetime column because datetime doesn't include seconds
  df$timestamp <- as.POSIXct(df$timestamp)

  #Apply rules
  df$flag <- flag_by_rules(
    x =df[,param],
    rules = rules
  )

  #For step based change
  df <- df %>%
    arrange(timestamp) %>%
    ungroup()
  x <- df[[param]]
  df$delta_value <- x - dplyr::lag(x)
  df$step_flag <- !is.na(df$delta_value) & abs(df$delta_value) > step_threshold

  qc_plot <- function(df, time_col, value_col, flag_label_col, y_label, min_time, max_time) {

    # Guardrails
    stopifnot(
      time_col %in% names(df),
      value_col %in% names(df),
      flag_label_col %in% names(df)
    )

    # Ensure factor levels for colors
    if (length(unique(df[[flag_label_col]]))==3) {
      df[[flag_label_col]] <- factor(
        df[[flag_label_col]],
        levels = c("Good", "Suspect", "Bad")
      )
    } else if (length(unique(df[[flag_label_col]]))==2) {
      df[[flag_label_col]] <- factor(
        df[[flag_label_col]],
        levels = c(TRUE, FALSE)
      )
    }


    plot_ly(
      data = df,
      x = df[[time_col]],
      y = df[[value_col]],
      type = "scatter",
      mode = "markers",
      color = df[[flag_label_col]],
      colors = flag_colors,
      marker = list(size = 6),
      #can't use hovertemplate only b/c won't show info when hover correctly; need to create text then show text in hovertemplate:

      text = ~paste(
        "Time:", timestamp,
        "<br>",value_col, ":", round(df[[value_col]], 2),
        "<br>Name:", name,
        "<br>QC:", flag_label_col
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        title = paste(y_label, "with QC Flags"),
        xaxis = list(
          # title = "Time",
          type = "date",
          rangebreaks = list(
            list(bounds = c(min_time, max_time), pattern = "hour") #bounds here, adjust based on data trip times
          )
        ),
        yaxis = list(title = y_label)
      )
  }

  flag_points_plot <- NA
  flag_map_plot <- NA
  step_flag_points_plot <- NA
  step_flag_map_plot <- NA

  if (any(plots == "flag points")) {
    flag_colors <- c(
      "Good"    = "#2ECC71",
      "Suspect" = "orange",
      "Bad"     = "#E74C3C"
    )

    flag_points_plot <- qc_plot(df, time_col="timestamp", value_col=param, flag_label_col="flag", y_label=param, min_time=min_time, max_time=max_time)

  }
  if (any(plots == "flag map")) {
    map_center <- list(
      lon = mean(df$lon_dd, na.rm = TRUE),
      lat = mean(df$lat_dd, na.rm = TRUE)
    )
    flag_map_plot<-plot_ly(
      data = df,
      type="scattermapbox",   #can be issues using "scattermapbox" type when combining map figure with another plot that is a different type
      mode = "markers",
      lon = ~lon_dd,
      lat = ~lat_dd,
      height = "100%", #this is needed in here to ensure when stacked, shows completely
      color = ~flag,
      colors = c(
        "Good" = "#2ECC71",
        "Suspect" = "orange",
        "Bad" = "#E74C3C"
      ),
      marker = list(size = 4, opacity = 0.8),
      text = ~paste(
        "Time:", timestamp,
        "<br>",param, ":", round(df[[param]], 1),
        "<br>QC:", flag,
        "<br>Name:", name
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        mapbox = list(
          style = "carto-positron",
          zoom = 9.5,
          center = map_center    #need to have done above to reference here
        ),
        legend = list(title = list(text = "QC Flag")),
        margin = list(l = 0, r = 0, t = 0, b = 0)
      )
  }
  if (any(plots == "step flag points")) {
    flag_colors <- c(
      "TRUE"    = "red",
      "FALSE" = "blue"
    )
    step_flag_points_plot <- qc_plot(df, time_col="timestamp", value_col=param, flag_label_col="step_flag", y_label=paste0(param, "change from previous value"), min_time=min_time, max_time=max_time)
  }
  if (any(plots == "step flag map")) {
    map_center <- list(
      lon = mean(df$lon_dd, na.rm = TRUE),
      lat = mean(df$lat_dd, na.rm = TRUE)
    )
    step_flag_map_plot <- plot_ly(
      data = df,
      type = "scattermapbox",
      mode = "markers",
      lon = ~lon_dd,
      lat = ~lat_dd,
      color = ~factor(step_flag),
      colors = c(
        "FALSE" = "blue",
        "TRUE"  = "#E74C3C"
      ),
      marker = list(
        size = 4,
        opacity = 0.8,
        line = list(
          width = ifelse(df$step_flag, 2, 0),
          color = "black"
        )
      ),
      text = ~paste(
        "Time:", timestamp,
        "<br>", param, ":", round(df[[param]], 2),
        "<br>ΔValue:", round(delta_value, 2),
        "<br>Name:", name,
        "<br>Step flag:", step_flag
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        mapbox = list(
          style = "carto-positron",
          zoom = 9.5,
          center = list(
            lon = mean(df$lon, na.rm = TRUE),
            lat = mean(df$lat, na.rm = TRUE)
          )
        ),
        legend = list(title = list(text = "Step-change Flag")),
        margin = list(l = 0, r = 0, t = 0, b = 0)
      )
  }
  all_plots <- list(flag_points_plot, flag_map_plot, step_flag_points_plot, step_flag_map_plot)
  use_plots <-  which(!is.na(all_plots))
  if (length(use_plots) > 2) {
    stop("Error: only two plots can be displayed at once")
  }
  if (length(use_plots) == 1) {
    plot <- all_plots[use_plots]
  } else {
    plot <- htmltools::browsable(
      tags$div(
        style = "
      display: flex;
      flex-direction: column;
      height: 100vh;
      width: 100%;
    ",
        tags$div(
          style = "flex: 0 0 40vh; width: 100%; position: relative;",
          all_plots[use_plots[1]]
        ),
        tags$div(
          style = "flex: 0 0 60vh; width: 100%; position: relative;",
          all_plots[use_plots[2]]
        )
      )
    )
  }
  return(list(df, plot))
}



#'@name replace_data
#'@title Replace bad data values
#'@description Replace bad data values based on user-specified criteria
#'@param df dataframe (if specifying values based on QA flags, the dataframe must contain the columns with those flags)
#'@param yearmon numeric designation of survey date formatted as yyyymm
#'@param value_var name of parameter of interest
#'@param replace_value what to replace selected values with (default is to set selected values to NA)
#'@param min_value replace values less than this
#'@param max_value replace values more than this
#'@param start_date replace values occuring after this datetime (must be used with end_date)
#'@param end_date replace values before this datetime (must be used with start_date)
#'@param date_var the date/time column that start_date and end_date reference (default is timestamp column created in QAflags function)
#'@param flag selects values with a given flag
#'@param flag_var column name that contains flags
#'@param step_flag selects values with given step flag
#'@param step_flag_var column name that contains step flags
#'@param logic whether to follow "and" or "or"conditions when selecting values
#'@param save whether to save output to QA datasets folder
#'@param fdir character file path to local data directory
#'@export
#'@import dplyr
#'@import lubridate
#'@examples \dontrun{
#'df <- replace_data(df=data, yearmon=202510, value_var="tal.pc.rfu", flag="Bad", step_flag=TRUE, logic="and", save=T)
#'df <- replace_data(df=data, yearmon=202510, value_var="tal.pc.rfu", min_value=0, replace_val=0, save=T)
#'dat <- replace_data(df=data, yearmon=202510, value_var="tal.pc.rfu", start_date="2025-10-21 12:57:42", end_date="2025-10-21 13:13:13", replace_val=20, save=F)
#'}

replace_data <- function(df=df, yearmon,
                         value_var,
                         replace_val =NA,
                         min_value = NULL, #replace values less than this
                         max_value = NULL, #replace values more than this
                         start_date = NULL, #replace values after this date
                         end_date = NULL, #replace values before this date
                         date_var = "timestamp",
                         flag = NULL, #replace values based on flag label
                         flag_var = "flag",
                         step_flag = NULL, #replace values based on step flag label
                         step_flag_var = "step_flag",
                         logic = "or",
                         save = TRUE,
                         fdir = getOption("fdir")) {

  #logic <- match.arg(logic)

  # If no criteria, return unchanged data
  if (all(sapply(list(min_value, max_value, start_date, end_date, flag, step_flag), is.null))) {
    return(df)
  }

  condition <- if (logic == "and") rep(TRUE, nrow(df)) else rep(FALSE, nrow(df))

  combine <- function(a, b) {
    if (logic == "and") a & b else a | b
  }

  # Numeric range
  if (!is.null(min_value) && !is.null(value_var)) {
    condition <- combine(condition, df[[value_var]] < min_value)
  }

  if (!is.null(max_value) && !is.null(value_var)) {
    condition <- combine(condition, df[[value_var]] > max_value)
  }

  # Date range
  if (!is.null(start_date) && !is.null(end_date) && !is.null(date_var)) {
    date_range <- seq.POSIXt(as.POSIXct(start_date), as.POSIXct(end_date), by="sec")
    date_range <- force_tz(date_range, tzone="UTC")
    condition <- combine(condition, df[[date_var]] %in% date_range)
  }

  # Flag conditions
  if (!is.null(flag) && !is.null(flag_var)) {
    condition <- combine(condition, df[[flag_var]] %in% flag)
  }

  if (!is.null(step_flag) && !is.null(step_flag_var)) {
    condition <- combine(condition, df[[step_flag_var]] %in% step_flag)
  }

  # Replace values with NA where condition is TRUE
  df <- df %>%
    mutate(
      across(all_of(value_var),
             ~ replace(., condition, replace_val))
    )

  if (save == TRUE) {
    col_rmv <- which(colnames(df) %in% c("flag", "flag_label", "delta_value", "step_flag", "step_flag_label"))
    df <- df[,-col_rmv]
    write.csv(df, file.path(fdir, "DF_FullDataSets", "QA datasets", paste(yearmon, "j_qa.csv", sep = "")))
  }

  return(df)

}





