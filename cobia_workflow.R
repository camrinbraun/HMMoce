# ==============================================================================
# Cobia PSAT Tag 36419 Workflow Script
# Direct processing from Microwave Telemetry XLS export to HMMoce
# ==============================================================================

library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tags2etuff)
source("c:/Users/benjamin.galuardi/Documents/GitHub/HMMoce/microwave_to_etuff.R")
devtools::load_all("c:/Users/benjamin.galuardi/Documents/GitHub/HMMoce")
# invisible(lapply(list.files("c:/Users/benjamin.galuardi/Documents/GitHub/HMMoce/R", 
#                             pattern = "\\.[Rr]$", full.names = TRUE), source))

base_dir <- "G:/.shortcut-targets-by-id/1jN4NVXfNcxVqIPxVvl3k3e3e4W3J5HgU/Cobia Shared Files"
rerun_dir <- file.path(base_dir, "hmmoce-rerun")

get_cobia_meta <- function(tag_id, dir = base_dir) {
  tag_id <- as.character(tag_id)
  tag_dir <- file.path(dir, "Updated PSAT Reports", paste0(tag_id, "_updated(2013)"))
  excel_file <- file.path(tag_dir, paste0(tag_id, "_updated(2013).xls"))
  
  if (!file.exists(excel_file)) {
    stop(paste("Excel file not found for tag:", tag_id, "at", excel_file))
  }
  
  tracks_file <- file.path(dir, "Data files", "Cobia tracks.csv")
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

# 2. Generate and Load eTUFF File ----------------------------------------------
out_dir <- file.path(rerun_dir, instrument_id)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

etuff_out_file <- file.path(out_dir, paste0(instrument_id, "_eTUFF.txt"))
if (!file.exists(etuff_out_file)) {
  message("Generating eTUFF file for tag ", instrument_id, " from: ", excel_file)
  etuff <- microwave_to_etuff(excel_file, tag_meta, out_file = etuff_out_file)
} else {
  message("Loading existing eTUFF file: ", etuff_out_file)
  etuff <- tags2etuff::read_etuff(etuff_out_file)
}

# 3. Setup Start / End Bounds & Daily Timesteps --------------------------------
t_start <- as.POSIXct(tag_meta$time_coverage_start, tz = "UTC")
t_end   <- as.POSIXct(tag_meta$time_coverage_end, tz = "UTC")

iniloc <- data.frame(
  day   = c(lubridate::day(t_start), lubridate::day(t_end)),
  month = c(lubridate::month(t_start), lubridate::month(t_end)),
  year  = c(lubridate::year(t_start), lubridate::year(t_end)),
  lat   = c(as.numeric(tag_meta$geospatial_lat_start), as.numeric(tag_meta$geospatial_lat_end)),
  lon   = c(as.numeric(tag_meta$geospatial_lon_start), as.numeric(tag_meta$geospatial_lon_end))
)
rownames(iniloc) <- c("start", "end")

tag <- as.POSIXct(sprintf("%04d-%02d-%02d", iniloc["start", "year"], iniloc["start", "month"], iniloc["start", "day"]), tz = "UTC")
pop <- as.POSIXct(sprintf("%04d-%02d-%02d", iniloc["end", "year"], iniloc["end", "month"], iniloc["end", "day"]), tz = "UTC")
dateVec <- seq.POSIXt(tag, pop, by = "24 hours")

message("Deployment spans: ", format(tag, "%Y-%m-%d"), " to ", format(pop, "%Y-%m-%d"), " (", length(dateVec), " daily steps)")
print(iniloc)

# 4. Extract Observations from eTUFF -------------------------------------------

## 4.1 Daily SST (Sea Surface Temperature)
if ("tempMax" %in% names(etuff$etuff)) {
  tag.sst <- etuff$etuff %>%
    dplyr::filter(!is.na(tempMax)) %>%
    dplyr::mutate(Date = as.Date(DateTime)) %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(
      Depth = 0,
      Temperature = mean(tempMax, na.rm = TRUE)
    ) %>%
    dplyr::mutate(Date = as.POSIXct(Date, tz = "UTC")) %>%
    as.data.frame()
  message("Extracted ", nrow(tag.sst), " daily SST records from eTUFF.")
} else if ("temperature" %in% names(etuff$etuff)) {
  tag.sst <- etuff$etuff %>%
    dplyr::filter(!is.na(temperature) & (is.na(depth) | depth <= 5)) %>%
    dplyr::mutate(Date = as.Date(DateTime)) %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(
      Depth = 0,
      Temperature = mean(temperature, na.rm = TRUE)
    ) %>%
    dplyr::mutate(Date = as.POSIXct(Date, tz = "UTC")) %>%
    as.data.frame()
  message("Extracted ", nrow(tag.sst), " surface temperature records from eTUFF.")
} else {
  tag.sst <- NULL
}

## 4.2 Depth-Temperature Profile Data (PDT)
if (all(c("depth", "temperature") %in% names(etuff$etuff))) {
  ts_df <- etuff$etuff %>%
    dplyr::filter(!is.na(depth) & !is.na(temperature)) %>%
    dplyr::select(Date = DateTime, Depth = depth, Temperature = temperature)
  
  if (nrow(ts_df) > 0) {
    pdt <- bin_TempTS(ts_df, out_dates = dateVec, bin_res = 8)
    message("Generated PDT profiles for ", length(unique(pdt$Date)), " unique days.")
  } else {
    pdt <- NULL
  }
} else {
  pdt <- NULL
}

## 4.3 Light / Geolocation Estimates
if (all(c("latitude", "longitude") %in% names(etuff$etuff))) {
  lightloc <- etuff$etuff %>%
    dplyr::filter(!is.na(latitude) & !is.na(longitude)) %>%
    dplyr::select(Date = DateTime, Latitude = latitude, Longitude = longitude) %>%
    dplyr::mutate(
      Error.Semi.minor.axis = 50,
      Error.Semi.major.axis = 100,
      Offset = 0,
      Offset.orientation = 0
    ) %>%
    as.data.frame()
  message("Extracted ", nrow(lightloc), " light-based location estimates.")
} else {
  lightloc <- NULL
}

## 4.4 Maximum Daily Depth (for Bathymetry masking)
if ("depthMax" %in% names(etuff$etuff)) {
  mmd <- etuff$etuff %>%
    dplyr::filter(!is.na(depthMax)) %>%
    dplyr::select(DateTime, depthMax) %>%
    as.data.frame()
} else if ("depth" %in% names(etuff$etuff)) {
  mmd <- etuff$etuff %>%
    dplyr::filter(!is.na(depth)) %>%
    dplyr::mutate(Date = as.Date(DateTime)) %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(depthMax = max(depth, na.rm = TRUE)) %>%
    dplyr::mutate(DateTime = as.POSIXct(Date, tz = "UTC")) %>%
    dplyr::select(DateTime, depthMax) %>%
    as.data.frame()
} else {
  mmd <- NULL
}
if (!is.null(mmd)) message("Extracted ", nrow(mmd), " daily maximum depth records.")

# 5. Spatial Domain & Grid Setup -----------------------------------------------

# Define study area bounding box around Gulf of Mexico & US SE Atlantic
# Bounds dynamically envelope deployment coords with buffer
pad <- 3.0
lon_coords <- c(iniloc$lon, if (!is.null(lightloc)) lightloc$Longitude else NULL)
lat_coords <- c(iniloc$lat, if (!is.null(lightloc)) lightloc$Latitude else NULL)

sp.lim <- list(
  lonmin = max(-98, floor(min(lon_coords, na.rm = TRUE) - pad)),
  lonmax = min(-65, ceiling(max(lon_coords, na.rm = TRUE) + pad)),
  latmin = max(18,  floor(min(lat_coords, na.rm = TRUE) - pad)),
  latmax = min(40,  ceiling(max(lat_coords, na.rm = TRUE) + pad))
)

message("Spatial Limits: Lon [", sp.lim$lonmin, ", ", sp.lim$lonmax, "] | Lat [", sp.lim$latmin, ", ", sp.lim$latmax, "]")

# Setup spatial grid for likelihood computations
locs.grid <- setup.locs.grid(sp.lim, res = "hycom")

# Bathymetry grid
bathy.dir <- file.path(rerun_dir, "EnvData", "bathy")
if (!dir.exists(bathy.dir)) dir.create(bathy.dir, recursive = TRUE)

bathy_file <- file.path(bathy.dir, "bathy.nc")
if (file.exists(bathy_file)) {
  message("Loading cached bathymetry: ", bathy_file)
  bathy <- raster::raster(bathy_file)
} else {
  message("Fetching NOAA bathymetry raster for study area...")
  bathy <- tryCatch(
    get.bath.data(spatLim = sp.lim, save.dir = bathy.dir, res = 0.5),
    error = function(e) {
      message("Bathymetry remote download fallback: ", e$message)
      NULL
    }
  )
}
# 6. Observation Likelihood Calculations ---------------------------------------

## 6.1 Light / Geolocation Likelihood
if (!is.null(lightloc) && nrow(lightloc) > 0) {
  light_rds <- file.path(out_dir, "L.light.rds")
  if (!file.exists(light_rds)) {
    message("Calculating light likelihood (L.light)...")
    L.light <- calc.lightloc(
      lightloc,
      locs.grid = locs.grid,
      dateVec = dateVec,
      errEll = FALSE,
      lon_only = TRUE
    )
    saveRDS(L.light, light_rds)
  } else {
    message("Loading cached L.light: ", light_rds)
    L.light <- readRDS(light_rds)
  }
} else {
  L.light <- NULL
}

## 6.2 Bathymetry Likelihood
if (!is.null(bathy) && !is.null(mmd) && nrow(mmd) > 0) {
  bathy_rds <- file.path(out_dir, "L.bathy.rds")
  if (!file.exists(bathy_rds)) {
    message("Calculating bathymetry likelihood (L.bathy)...")
    bathy_crop <- raster::crop(bathy, raster::extent(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax))
    L.bathy <- calc.bathy(
      mmd,
      bathy_crop,
      dateVec = dateVec,
      focalDim = 5,
      sens.err = 5,
      lik.type = "max"
    )
    saveRDS(L.bathy, bathy_rds)
  } else {
    message("Loading cached L.bathy: ", bathy_rds)
    L.bathy <- readRDS(bathy_rds)
  }
} else {
  L.bathy <- NULL
}

## 6.3 SST Likelihood
# Note: Requires local or downloaded SST rasters (e.g. NOAA OISST / GHRSST / GLORYS / HYCOM)
# L.sst <- calc.sst(tag.sst, sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
L.sst <- NULL

## 6.4 Combine Available Likelihood Rasters
lik_list <- list(light = L.light, bathy = L.bathy, sst = L.sst)
lik_list <- lik_list[!sapply(lik_list, is.null)]
message("Available observation likelihoods: ", paste(names(lik_list), collapse = ", "))
