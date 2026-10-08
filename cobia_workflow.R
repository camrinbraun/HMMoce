# ==============================================================================
# Cobia PSAT Tag 36419 Workflow Script
# Direct processing from Microwave Telemetry XLS export to HMMoce
# ==============================================================================

library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
# devtools::load_all("c:/Users/benjamin.galuardi/Documents/GitHub/HMMoce")
invisible(lapply(list.files("c:/Users/benjamin.galuardi/Documents/GitHub/HMMoce/R", 
                            pattern = "\\.R$", full.names = TRUE), source))

# 1. Dynamic Metadata and Paths ------------------------------------------------
base_dir <- "G:/.shortcut-targets-by-id/1jN4NVXfNcxVqIPxVvl3k3e3e4W3J5HgU/Cobia Shared Files"

get_cobia_meta <- function(tag_id, base_dir) {
  tag_id <- as.character(tag_id)
  tag_dir <- file.path(base_dir, "Updated PSAT Reports", paste0(tag_id, "_updated(2013)"))
  excel_file <- file.path(tag_dir, paste0(tag_id, "_updated(2013).xls"))
  
  if (!file.exists(excel_file)) {
    stop(paste("Excel file not found for tag:", tag_id, "at", excel_file))
  }
  
  tracks_file <- file.path(base_dir, "Data files", "Cobia tracks.csv")
  tracks <- read.csv(tracks_file)
  tag_track <- tracks[tracks[, 1] == as.numeric(tag_id), ]
  
  if (nrow(tag_track) > 0) {
    tag_track$Date <- sprintf("%04d-%02d-%02d", tag_track$Year, tag_track$Month, tag_track$Day)
    tag_track <- tag_track[order(tag_track$Date), ]
    
    t_start <- as.POSIXct(paste(tag_track$Date[1], "00:00:00"), tz = "UTC")
    t_end   <- as.POSIXct(paste(tag_track$Date[nrow(tag_track)], "00:00:00"), tz = "UTC")
    lat_start <- tag_track$Lat_N[1]
    lon_start <- tag_track$Lon_E[1]
    lat_end   <- tag_track$Lat_N[nrow(tag_track)]
    lon_end   <- tag_track$Lon_E[nrow(tag_track)]
  } else {
    t_start <- NA; t_end <- NA; lat_start <- NA; lon_start <- NA; lat_end <- NA; lon_end <- NA
  }
  
  meta <- data.frame(
    instrument_name = tag_id,
    platform = "Rachycentron canadum",
    person_owner = "Jim Franks",
    manufacturer = "Microwave",
    time_coverage_start = t_start,
    time_coverage_end = t_end,
    geospatial_lat_start = lat_start,
    geospatial_lon_start = lon_start,
    geospatial_lat_end = lat_end,
    geospatial_lon_end = lon_end,
    excel_file = excel_file,
    tag_dir = tag_dir,
    stringsAsFactors = FALSE
  )
  return(meta)
}

# Select target tag ID (e.g. "36419", "36420", "36421", "36422", "40302", "49157", "55550", "55552")
instrument_id <- "36420"
tag_meta <- get_cobia_meta(instrument_id, base_dir)
excel_file <- tag_meta$excel_file

# 2. Extract Depth & Temperature Time Series from XLS ---------------------------
message("Reading Depth and Temperature time series from: ", excel_file)

# Depth (Pressure)
p_raw <- readxl::read_excel(excel_file, sheet = "Press Data", skip = 1)
depth_ts <- p_raw %>%
  dplyr::select(
    Date = 1,
    Press_val = 2,
    Depth = 3
  ) %>%
  dplyr::filter(!is.na(Date) & !is.na(Depth)) %>%
  dplyr::mutate(
    Date = as.POSIXct(Date, tz = "UTC"),
    Depth = abs(as.numeric(Depth)) # Depth in positive meters
  )

# Temperature
t_raw <- readxl::read_excel(excel_file, sheet = "Temp Data", skip = 1)
temp_ts <- t_raw %>%
  dplyr::select(
    Date = 1,
    Temp_val = 2,
    Temperature = 3
  ) %>%
  dplyr::filter(!is.na(Date) & !is.na(Temperature)) %>%
  dplyr::mutate(
    Date = as.POSIXct(Date, tz = "UTC"),
    Temperature = as.numeric(Temperature)
  )

# Combine into synchronized TS
ts <- dplyr::inner_join(depth_ts, temp_ts, by = "Date") %>%
  dplyr::select(Date, Depth, Temperature) %>%
  dplyr::arrange(Date)

message("Extracted ", nrow(ts), " matched depth-temperature observations.")

# 3. Bin Time Series to PDT-like profiles using HMMoce::bin_TempTS -------------
# Create daily target dates spanning tag deployment
out_dates <- seq.POSIXt(
  from = tag_meta$time_coverage_start,
  to = tag_meta$time_coverage_end,
  by = "day"
)

# Bin resolution (e.g. 8 m standard, or 10/25 m)
pdt <- bin_TempTS(ts, out_dates = out_dates, bin_res = 8)

message("Generated PDT profiles for ", length(unique(pdt$Date)), " unique days.")
print(head(pdt))

# 4. Extract Daily MinMax Temperature Data --------------------------------------
t_minmax <- readxl::read_excel(excel_file, sheet = "Temp Data (MinMax)", skip = 1) %>%
  dplyr::select(
    Date = 1,
    MinTemp = 4,
    MaxTemp = 5
  ) %>%
  dplyr::filter(!is.na(Date)) %>%
  dplyr::mutate(
    Date = as.POSIXct(Date, tz = "UTC"),
    MinTemp = as.numeric(MinTemp),
    MaxTemp = as.numeric(MaxTemp)
  )

message("Extracted ", nrow(t_minmax), " daily Min/Max SST records.")
