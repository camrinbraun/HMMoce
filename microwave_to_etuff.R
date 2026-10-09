# ==============================================================================
# microwave_to_etuff.R
# Robust converter from Microwave Telemetry XLS/XLSX files to standardized eTUFF
# ==============================================================================

library(readxl)
library(dplyr)
library(reshape2)
library(data.table)
library(lubridate)
library(tags2etuff)

#' Convert Microwave Telemetry XLS/XLSX to eTUFF file
#'
#' @param excel_file Path to Microwave .xls or .xlsx file
#' @param meta_row Data frame containing tag metadata (instrument_name, platform, etc.)
#' @param out_file Destination path for the .txt eTUFF file
#' @param obs_types_url Optional custom URL/path to eTUFF-ObservationTypes.csv
#' @return An object of class \code{etuff}
#' @export
microwave_to_etuff <- function(excel_file, meta_row, out_file = NULL, obs_types_url = NULL) {
  
  if (!file.exists(excel_file)) {
    stop("Excel file not found: ", excel_file)
  }
  
  if (is.null(out_file)) {
    out_file <- file.path(dirname(excel_file), paste0(meta_row$instrument_name[1], "_eTUFF.txt"))
  }
  
  # 1. Fetch obsTypes mapping --------------------------------------------------
  if (is.null(obs_types_url)) {
    obs_types_url <- "https://raw.githubusercontent.com/camrinbraun/tagbase/master/eTUFF-ObservationTypes.csv"
  }
  
  obsTypes <- tryCatch(
    read.csv(obs_types_url, stringsAsFactors = FALSE),
    error = function(e) {
      stop("Unable to download obsTypes from: ", obs_types_url, "\nDetails: ", e$message)
    }
  )
  if (names(obsTypes)[1] != "VariableID") names(obsTypes)[1] <- "VariableID"

  # Standardize / populate required metadata attributes
  meta_row <- as.data.frame(meta_row, stringsAsFactors = FALSE)
  if (!"uid_no" %in% names(meta_row)) meta_row$uid_no <- 1
  if (!"instrument_type" %in% names(meta_row)) meta_row$instrument_type <- "PSAT"
  if (!"model" %in% names(meta_row)) meta_row$model <- "Microwave PSAT"
  if (!"owner_contact" %in% names(meta_row)) meta_row$owner_contact <- "contact@example.com"
  if (!"serial_number" %in% names(meta_row)) meta_row$serial_number <- as.character(meta_row$instrument_name[1])
  if (!"taxonomic_serial_number" %in% names(meta_row)) meta_row$taxonomic_serial_number <- "168566" # Rachycentron canadum
  if (!"end_details" %in% names(meta_row)) meta_row$end_details <- "Pop-off"
  if (!"end_type" %in% names(meta_row)) meta_row$end_type <- "Pop-off"
  if (!"waypoints_source" %in% names(meta_row)) meta_row$waypoints_source <- "GPS/Argos"
  if (!"found_problem" %in% names(meta_row)) meta_row$found_problem <- "no"
  if (!"person_qc" %in% names(meta_row)) meta_row$person_qc <- meta_row$person_owner[1]
  
  sheet_names <- readxl::excel_sheets(excel_file)
  all_data <- list()
  
  # Helper to standardize rows against obsTypes
  standardize_obs <- function(df, val_vars) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    m <- reshape2::melt(df, id.vars = "DateTime", measure.vars = intersect(val_vars, names(df)))
    m$VariableName <- as.character(m$variable)
    m$VariableValue <- as.numeric(m$value)
    m <- merge(m, obsTypes[, c("VariableID", "VariableName", "VariableUnits")], by = "VariableName", all.x = TRUE)
    m <- m[!is.na(m$VariableValue), c("DateTime", "VariableID", "VariableValue", "VariableName", "VariableUnits")]
    return(m)
  }

  # 2. Extract Time Series (Depth & Temp) ---------------------------------------
  if ("Press Data" %in% sheet_names && "Temp Data" %in% sheet_names) {
    message("-> Processing Depth and Temperature time series...")
    
    p_raw <- readxl::read_excel(excel_file, sheet = "Press Data", skip = 1)
    p_df <- data.frame(
      DateTime = as.POSIXct(p_raw[[1]], tz = "UTC"),
      pressure = as.numeric(p_raw[[2]]),
      depth = abs(as.numeric(p_raw[[3]]))
    )
    p_df <- p_df[!is.na(p_df$DateTime), ]
    
    t_raw <- readxl::read_excel(excel_file, sheet = "Temp Data", skip = 1)
    t_df <- data.frame(
      DateTime = as.POSIXct(t_raw[[1]], tz = "UTC"),
      temperature = as.numeric(t_raw[[3]])
    )
    t_df <- t_df[!is.na(t_df$DateTime), ]
    
    mti <- merge(p_df, t_df, by = "DateTime", all = TRUE)
    all_data[[length(all_data) + 1]] <- standardize_obs(mti, c("depth", "pressure", "temperature"))
  }

  # 3. Extract Min/Max Summaries -----------------------------------------------
  if ("Press Data (MinMax)" %in% sheet_names && "Temp Data (MinMax)" %in% sheet_names) {
    message("-> Processing Min/Max Depth and Temp summaries...")
    
    pm_raw <- readxl::read_excel(excel_file, sheet = "Press Data (MinMax)", skip = 1)
    pm_df <- data.frame(
      DateTime = as.POSIXct(pm_raw[[1]], tz = "UTC"),
      depthMin = abs(as.numeric(pm_raw[[4]])),
      depthMax = abs(as.numeric(pm_raw[[5]]))
    )
    pm_df <- pm_df[!is.na(pm_df$DateTime), ]
    
    tm_raw <- readxl::read_excel(excel_file, sheet = "Temp Data (MinMax)", skip = 1)
    tm_df <- data.frame(
      DateTime = as.POSIXct(tm_raw[[1]], tz = "UTC"),
      tempMin = as.numeric(tm_raw[[4]]),
      tempMax = as.numeric(tm_raw[[5]])
    )
    tm_df <- tm_df[!is.na(tm_df$DateTime), ]
    
    minmax <- merge(pm_df, tm_df, by = "DateTime", all = TRUE)
    all_data[[length(all_data) + 1]] <- standardize_obs(minmax, c("depthMin", "depthMax", "tempMin", "tempMax"))
  }

  # 4. Extract Sunrise / Sunset Times ------------------------------------------
  if ("Sunrise and Sunset Times" %in% sheet_names) {
    message("-> Processing Sunrise / Sunset observations...")
    srss_raw <- tryCatch(
      readxl::read_excel(excel_file, sheet = "Sunrise and Sunset Times", skip = 1),
      error = function(e) NULL
    )
    if (!is.null(srss_raw) && ncol(srss_raw) >= 4) {
      sr_dates <- as.POSIXct(srss_raw[[1]], tz = "UTC")
      # Extract depth at sunrise/sunset if available
      sr_depth <- if (ncol(srss_raw) >= 5) abs(as.numeric(srss_raw[[5]])) else rep(NA, length(sr_dates))
      ss_depth <- if (ncol(srss_raw) >= 6) abs(as.numeric(srss_raw[[6]])) else rep(NA, length(sr_dates))
      
      sr_df <- data.frame(DateTime = sr_dates, depthSunrise = sr_depth)
      ss_df <- data.frame(DateTime = sr_dates, depthSunset = ss_depth)
      
      all_data[[length(all_data) + 1]] <- standardize_obs(sr_df, "depthSunrise")
      all_data[[length(all_data) + 1]] <- standardize_obs(ss_df, "depthSunset")
    }
  }

  # 5. Extract Lat / Long Track Estimates ---------------------------------------
  if ("Lat&Long" %in% sheet_names) {
    message("-> Processing Light-based Geolocation estimates...")
    ll_raw <- tryCatch(
      readxl::read_excel(excel_file, sheet = "Lat&Long", skip = 1),
      error = function(e) NULL
    )
    if (!is.null(ll_raw) && ncol(ll_raw) >= 3) {
      # First column: date / excel serial; col 2: lat; col 3: lon
      dates_val <- suppressWarnings(as.numeric(ll_raw[[1]]))
      dt <- as.POSIXct(dates_val * 86400, origin = "1899-12-30", tz = "UTC")
      lat <- as.numeric(ll_raw[[2]])
      lon <- as.numeric(ll_raw[[3]])
      
      # Microwave west longitudes are typically reported as positive degrees
      if (all(lon >= 0, na.rm = TRUE) && any(lon > 50, na.rm = TRUE)) {
        lon <- lon * -1
      }
      
      loc_df <- data.frame(DateTime = dt, latitude = lat, longitude = lon)
      loc_df <- loc_df[!is.na(loc_df$DateTime), ]
      all_data[[length(all_data) + 1]] <- standardize_obs(loc_df, c("latitude", "longitude"))
    }
  }

  # 6. Combine and Format Observations -----------------------------------------
  returnData <- do.call(rbind, all_data)
  if (is.null(returnData) || nrow(returnData) == 0) {
    stop("No valid observations could be extracted from: ", excel_file)
  }
  
  returnData <- dplyr::distinct(returnData, DateTime, VariableName, .keep_all = TRUE)
  returnData <- returnData[order(returnData$DateTime, returnData$VariableID), ]
  returnData$DateTime <- format(returnData$DateTime, "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  # 7. Write Header and Data to eTUFF File -------------------------------------
  message("-> Writing eTUFF header to: ", out_file)
  tags2etuff::build_meta_head(meta_row = meta_row, filename = out_file, write_hdr = TRUE)
  
  message("-> Appending ", nrow(returnData), " records to: ", out_file)
  data.table::fwrite(
    returnData,
    file = out_file,
    sep = ",",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    append = TRUE
  )

  # 8. Return canonical etuff R object via read_etuff -------------------------
  message("-> Reading and verifying generated eTUFF object...")
  res <- tags2etuff::read_etuff(out_file)
  message("-> Conversion successful!")
  return(res)
}
