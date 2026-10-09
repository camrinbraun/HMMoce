# ==============================================================================
# Cobia PSAT Tag 36419 Workflow Script
# Direct processing from Microwave Telemetry XLS export to HMMoce
# ==============================================================================

library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tags2etuff)
# source("c:/Users/benjamin.galuardi/Documents/GitHub/HMMoce/microwave_to_etuff.R")
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
    dplyr::select(Date = DateTime, MaxDepth = depthMax) %>%
    as.data.frame()
} else if ("depth" %in% names(etuff$etuff)) {
  mmd <- etuff$etuff %>%
    dplyr::filter(!is.na(depth)) %>%
    dplyr::mutate(d = as.Date(DateTime)) %>%
    dplyr::group_by(d) %>%
    dplyr::summarise(MaxDepth = max(depth, na.rm = TRUE)) %>%
    dplyr::mutate(Date = as.POSIXct(d, tz = "UTC")) %>%
    dplyr::select(Date, MaxDepth) %>%
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

## 6.3 GLORYS Ocean Model Data (3D Profile OHC & SST) -----------------------
glorys.dir <- file.path(rerun_dir, "EnvData", "glorys")
if (!dir.exists(glorys.dir)) dir.create(glorys.dir, recursive = TRUE)

# Interactive download function for GLORYS daily NetCDF files via copernicusmarine
# To run interactively in R:
# 1) Make sure `copernicusmarine` CLI is logged in: `copernicusmarine login`
# 2) Set download_glorys <- TRUE to fetch missing days
download_glorys <- TRUE

# Detect copernicusmarine executable path
cm_cli <- Sys.which("copernicusmarine")
if (cm_cli == "") {
  user_cm <- file.path(Sys.getenv("APPDATA"), "Python", "Python313", "Scripts", "copernicusmarine.exe")
  if (file.exists(user_cm)) cm_cli <- paste0('"', user_cm, '"') else cm_cli <- "copernicusmarine"
}

if (download_glorys) {
  message("Checking/downloading GLORYS daily files for date range...")
  temp_product_id <- "cmems_mod_glo_phy_my_0.083deg_P1D-m"
  
  # Identify missing dates only
  all_nc <- paste0(temp_product_id, "_", format(as.Date(dateVec), "%Y%m%d"), ".nc")
  missing_idx <- which(!file.exists(file.path(glorys.dir, all_nc)))
  
  if (length(missing_idx) > 0) {
    n_workers <- min(4, parallel::detectCores() - 1, length(missing_idx))
    message("Downloading ", length(missing_idx), " files using ", n_workers, " parallel workers...")
    
    cl <- parallel::makeCluster(n_workers)
    parallel::clusterExport(cl, varlist = c("cm_cli", "temp_product_id", "sp.lim", "glorys.dir", "dateVec"), envir = environment())
    
    parallel::parLapply(cl, missing_idx, function(i) {
      d_str <- format(as.Date(dateVec[i]), "%Y%m%d")
      d_dash <- format(as.Date(dateVec[i]), "%Y-%m-%d")
      out_nc <- paste0(temp_product_id, "_", d_str, ".nc")
      cmd <- sprintf(
        '%s subset -i %s -x %f -X %f -y %f -Y %f -t %s -T %s -z 0 -Z 2000 --variable thetao -o "%s" -f "%s" --force-download',
        cm_cli, temp_product_id, sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax,
        d_dash, d_dash, glorys.dir, out_nc
      )
      system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    })
    parallel::stopCluster(cl)
    message("Parallel download complete.")
  } else {
    message("All GLORYS daily files already downloaded.")
  }
}

# 6.3.1 GLORYS OHC Profile Likelihood (3D)
glorys_files <- list.files(glorys.dir, pattern = "\\.nc$", full.names = TRUE)
if (length(glorys_files) > 0 && !is.null(pdt) && nrow(pdt) > 0) {
  ohc_rds <- file.path(out_dir, "L.ohc.glorys.rds")
  if (!file.exists(ohc_rds)) {
    message("Calculating GLORYS OHC profile likelihood...")
    L.ohc <- calc.ohc.glorys(
      pdt = pdt,
      filename = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
      ohc.dir = glorys.dir,
      dateVec = dateVec
    )
    saveRDS(L.ohc, ohc_rds)
  } else {
    message("Loading cached L.ohc: ", ohc_rds)
    L.ohc <- readRDS(ohc_rds)
  }
} else {
  L.ohc <- NULL
}

# 6.3.2 GLORYS SST Likelihood
if (length(glorys_files) > 0 && !is.null(tag.sst) && nrow(tag.sst) > 0) {
  sst_rds <- file.path(out_dir, "L.sst.glorys.rds")
  if (!file.exists(sst_rds)) {
    message("Calculating GLORYS SST likelihood...")
    L.sst <- calc.sst.par.glorys(
      tag.sst,
      filename = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
      sst.dir = glorys.dir,
      dateVec = dateVec,
      sens.err = 1
    )
    saveRDS(L.sst, sst_rds)
  } else {
    message("Loading cached L.sst: ", sst_rds)
    L.sst <- readRDS(sst_rds)
  }
} else {
  L.sst <- NULL
}

## 6.4 Combine Available Likelihood Rasters -----------------------------------
L.rasters <- list(bathy = L.bathy, light = L.light, ohc = L.ohc, sst = L.sst)
L.rasters <- L.rasters[!sapply(L.rasters, is.null)]
message("Combining observation likelihoods: ", paste(names(L.rasters), collapse = ", "))

# Resample to common finest grid resolution
resamp.idx <- which.min(lapply(L.rasters, function(x) raster::res(x)[1]))
L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])

# Prepare daily maxDepth vector aligned with dateVec (zero-padded)
mmd_match <- match(as.Date(dateVec), as.Date(mmd$Date))
maxDepth_vec <- ifelse(!is.na(mmd_match), mmd$MaxDepth[mmd_match], 0)

# Form master observation likelihood array
L <- make.L(
  ras.list = L.res$L.rasters,
  iniloc = iniloc,
  dateVec = dateVec,
  maxDepth = maxDepth_vec,
  bathy = bathy
)
saveRDS(L, file.path(out_dir, "L.master.rds"))
message("Master observation likelihood L constructed: dim = ", paste(dim(L), collapse = " x "))

# 7. Model Comparison with Parameter Optimization -----------------------------

# 7.1 Define Likelihood Combinations to Test
combos <- list(
  all             = c("bathy", "light", "ohc", "sst"),
  bathy_light_sst = c("bathy", "light", "sst"),
  bathy_light_ohc = c("bathy", "light", "ohc"),
  bathy_light     = c("bathy", "light"),
  bathy_ohc       = c("bathy", "ohc"),
  bathy_sst       = c("bathy", "sst")
)
# Keep only combos where all requested likelihoods exist in L.res$L.rasters
combos <- combos[sapply(combos, function(x) all(x %in% names(L.res$L.rasters)))]

# 7.2 Run Each Combination
model_results <- list()
comp_rows <- list()

for (m_name in names(combos)) {
  sel_liks <- combos[[m_name]]
  message("\n==================================================")
  message("Running model: ", m_name, " (", paste(sel_liks, collapse = " + "), ")")
  message("==================================================")
  
  # Form observation likelihood for this subset
  L <- make.L(
    ras.list = L.res$L.rasters[sel_liks],
    iniloc = iniloc,
    dateVec = dateVec,
    maxDepth = maxDepth_vec,
    bathy = bathy
  )
  
  # Parameter optimization
  message("Optimizing movement parameters with HMMoce:::opt.params...")
  ncores <- min(parallel::detectCores() - 1, 8)
  pars_fit <- tryCatch({
    HMMoce:::opt.params(
      pars.init = c(2, 0.2, 0.6, 0.8),
      lower.bounds = c(0.1, 0.001, 0.1, 0.1),
      upper.bounds = c(6, 0.6, 0.9, 0.9),
      g = L.res$g,
      L = L,
      alg.opt = "ga",
      max_iter = 15,
      run = 10,
      p_size = 50,
      write.results = FALSE,
      ncores = ncores
    )
  }, error = function(e) {
    message("GA optimization failed (", e$message, "), falling back to optim (L-BFGS-B)...")
    HMMoce:::opt.params(
      pars.init = c(2, 0.2, 0.6, 0.8),
      lower.bounds = c(0.1, 0.001, 0.1, 0.1),
      upper.bounds = c(6, 0.6, 0.9, 0.9),
      g = L.res$g,
      L = L,
      alg.opt = "optim",
      write.results = FALSE
    )
  })
  
  pars <- as.numeric(pars_fit$par)
  if (length(pars) == 4) {
    sigmas <- pars[1:2]
sizes <- rep(ceiling(sigmas[1] * 4), 2)
    pb <- pars[3:4]
    muadvs <- c(0, 0)
    P <- matrix(c(pb[1], 1 - pb[1], 1 - pb[2], pb[2]), nrow = 2, ncol = 2, byrow = TRUE)
  } else {
    sigmas <- pars[1]
    sizes <- rep(ceiling(sigmas[1] * 4), 2)
    pb <- NULL
    muadvs <- 0
    P <- NULL
  }
  
  # Movement kernels
  if (sizes[1] %% 2 == 0) sizes[1] <- sizes[1] + 1
  K1 <- HMMoce:::gausskern.pg(sizes[1], sigmas[1], muadv = muadvs[1])
  K1 <- HMMoce:::mask.K(K1)
  
  if (!is.null(pb)) {
    if (sizes[2] %% 2 == 0) sizes[2] <- sizes[2] + 1
    K2 <- HMMoce:::gausskern.pg(sizes[2], sigmas[2], muadv = muadvs[2])
    K2 <- HMMoce:::mask.K(K2)
K <- list(K1, K2)
    m_states <- 2
  } else {
    K <- list(K1)
    m_states <- 1
  }
  
  # Forward filter
message("Running HMM filter...")
  f <- hmm.filter(g = L.res$g, L = L, K = K, P = P, m = m_states)

  # Likelihood metrics & AIC
  nllf <- -sum(log(f$psi[f$psi > 0]))
  aic <- 2 * nllf + 2 * length(pars)
  
  # Backward smoother
message("Running HMM smoother...")
s <- hmm.smoother(f, K = K, L = L, P = P)

  # Most probable track
tr <- calc.track(s, g = L.res$g, dateVec = dateVec, iniloc = iniloc, method = "mean")
  tr$Date <- as.POSIXct(dateVec, tz = "UTC")
  
  # Save individual model outputs
  fit_file <- file.path(out_dir, paste0(instrument_id, "_", m_name, "_fit.rds"))
  track_file <- file.path(out_dir, paste0(instrument_id, "_", m_name, "_track.csv"))
  saveRDS(list(model = m_name, pars = pars, s = s, tr = tr, aic = aic, nll = nllf), fit_file)
  write.csv(tr, track_file, row.names = FALSE)

  # Diagnostic PNG
  plot_file <- file.path(out_dir, paste0(instrument_id, "_", m_name, "_plotHMM.png"))
png(plot_file, width = 8, height = 9, units = "in", res = 300)
  plotHMM(s = s, track = tr, dateVec = dateVec, ptt = paste(instrument_id, m_name), behav.pts = TRUE, save.plot = FALSE)
dev.off()
  
  # Endpoint error (distance to popoff in km)
  last_idx <- nrow(tr)
  end_err_km <- raster::pointDistance(
    c(tr$lon[last_idx], tr$lat[last_idx]),
    c(iniloc$lon[2], iniloc$lat[2]),
    lonlat = TRUE
  ) / 1000
  
  comp_rows[[m_name]] <- data.frame(
    model = m_name,
    likelihoods = paste(sel_liks, collapse = "+"),
    sigma1 = round(sigmas[1], 3),
    sigma2 = if (length(sigmas) > 1) round(sigmas[2], 3) else NA,
    p11 = if (!is.null(pb)) round(pb[1], 3) else NA,
    p22 = if (!is.null(pb)) round(pb[2], 3) else NA,
    nll = round(nllf, 2),
    aic = round(aic, 2),
    end_error_km = round(end_err_km, 1),
    stringsAsFactors = FALSE
  )
  
  tr$model <- m_name
  model_results[[m_name]] <- tr
}

# 7.3 Model Comparison Summary Table
comp_table <- do.call(rbind, comp_rows)
comp_table <- comp_table[order(comp_table$aic), ]
write.csv(comp_table, file.path(out_dir, paste0(instrument_id, "_model_comparison.csv")), row.names = FALSE)
message("\nModel Comparison Table:")
print(comp_table)

# 8. Multi-Model Track Comparison Plot -----------------------------------------
library(ggplot2)

world_map <- map_data("world")
all_tracks <- do.call(rbind, model_results)

p_comp <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey80", color = "grey60") +
  coord_fixed(xlim = c(sp.lim$lonmin, sp.lim$lonmax), ylim = c(sp.lim$latmin, sp.lim$latmax)) +
  geom_path(data = all_tracks, aes(x = lon, y = lat, color = model), linewidth = 0.9, alpha = 0.85) +
  geom_point(data = iniloc[1, ], aes(x = lon, y = lat), fill = "green", color = "black", shape = 21, size = 4) +
  geom_point(data = iniloc[2, ], aes(x = lon, y = lat), fill = "red", color = "black", shape = 21, size = 4) +
  theme_bw() +
  labs(
    title = paste("Track Comparison across Likelihood Combinations - Tag", instrument_id),
    subtitle = "Green = Release, Red = Popoff",
    x = "Longitude", y = "Latitude", color = "Model"
  )

if (!is.null(lightloc) && nrow(lightloc) > 0) {
  p_comp <- p_comp +
    geom_point(data = lightloc, aes(x = Longitude, y = Latitude), color = "orange", alpha = 0.3, size = 1)
}

comp_plot_file <- file.path(out_dir, paste0(instrument_id, "_track_comparison.png"))
ggsave(comp_plot_file, plot = p_comp, width = 10, height = 8, dpi = 300)
message("Saved comparison track map to: ", comp_plot_file)

# 8.2 Detailed ggplot Map with Raw Light Locations & Endpoints
library(ggplot2)

world_map <- map_data("world")
tr$Date <- as.POSIXct(dateVec, tz = "UTC")

p_track <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey75", color = "grey50") +
  coord_fixed(xlim = c(sp.lim$lonmin, sp.lim$lonmax), ylim = c(sp.lim$latmin, sp.lim$latmax)) +
  theme_bw() +
  labs(title = paste("HMMoce Estimated Track - Cobia", instrument_id), x = "Longitude", y = "Latitude")

# Add raw light geolocations if available
if (!is.null(lightloc) && nrow(lightloc) > 0) {
  p_track <- p_track +
    geom_point(data = lightloc, aes(x = Longitude, y = Latitude), color = "orange", alpha = 0.4, size = 1.2)
}

# Add estimated trajectory
p_track <- p_track +
  geom_path(data = tr, aes(x = lon, y = lat), color = "darkblue", linewidth = 0.8) +
  geom_point(data = tr, aes(x = lon, y = lat, color = as.numeric(Date)), size = 2) +
  scale_color_viridis_c(name = "Date", breaks = as.numeric(pretty(tr$Date, n = 5)), labels = function(x) format(as.Date(as.POSIXct(x, origin = "1970-01-01")), "%b %d")) +
  geom_point(data = iniloc[1, ], aes(x = lon, y = lat), fill = "green", color = "black", shape = 21, size = 4) +
  geom_point(data = iniloc[2, ], aes(x = lon, y = lat), fill = "red", color = "black", shape = 21, size = 4)

ggsave(file.path(out_dir, paste0(instrument_id, "_track_map.png")), plot = p_track, width = 9, height = 7, dpi = 300)
message("Saved track map to: ", file.path(out_dir, paste0(instrument_id, "_track_map.png")))
print(p_track)
