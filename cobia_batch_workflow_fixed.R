# ==============================================================================
# Cobia PSAT Batch Workflow Script (All Tags in Folder)
# Runs the full cobia_workflow_fixed_params.R pipeline across all Cobia tags
# ==============================================================================

library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(sf)
library(terra)
library(tags2etuff)
devtools::load_all(".")

save_rast_rds <- function(r, file) {
  saveRDS(if (inherits(r, "SpatRaster")) terra::wrap(r) else r, file)
}
load_rast_rds <- function(file) {
  x <- readRDS(file)
  if (inherits(x, "PackedSpatRaster")) terra::unwrap(x) else if (inherits(x, "SpatRaster")) x else terra::rast(x)
}

# 1. Google Drive / Local Directory Setup --------------------------------------
gdrive_root_id <- "1jN4NVXfNcxVqIPxVvl3k3e3e4W3J5HgU" # "Cobia Shared Files"
local_gdrive   <- "G:/.shortcut-targets-by-id/1jN4NVXfNcxVqIPxVvl3k3e3e4W3J5HgU/Cobia Shared Files"
use_gdrive_api <- !dir.exists(local_gdrive)

if (!use_gdrive_api) {
  base_dir <- local_gdrive
} else {
  library(googledrive)
  base_dir <- file.path(getwd(), "cobia_gdrive_cache")

  get_or_create_gdrive_folder <- function(parent_id, folder_name) {
    items <- drive_ls(as_id(parent_id))
    hit <- items[items$name == folder_name, ]
    if (nrow(hit) > 0) return(hit$id[1])
    drive_mkdir(folder_name, path = as_id(parent_id))$id[1]
  }

  sync_gdrive_tag <- function(tag_id, base_dir, root_id = gdrive_root_id) {
    root_ls <- drive_ls(as_id(root_id))

    tracks_local <- file.path(base_dir, "Data files", "Cobia tracks.csv")
    if (!file.exists(tracks_local)) {
      dir.create(dirname(tracks_local), recursive = TRUE, showWarnings = FALSE)
      df_id <- root_ls$id[root_ls$name == "Data files"][1]
      tr_f  <- drive_ls(as_id(df_id), pattern = "Cobia tracks\\.csv")
      drive_download(as_id(tr_f$id[1]), path = tracks_local, overwrite = TRUE)
    }

    folder_nm <- paste0(tag_id, "_updated(2013)")
    xls_nm    <- paste0(folder_nm, ".xls")
    xls_local <- file.path(base_dir, "Updated PSAT Reports", folder_nm, xls_nm)
    if (!file.exists(xls_local)) {
      dir.create(dirname(xls_local), recursive = TRUE, showWarnings = FALSE)
      psat_id <- root_ls$id[root_ls$name == "Updated PSAT Reports"][1]
      tag_f   <- drive_ls(as_id(psat_id), pattern = paste0("^", tag_id, "_updated"))
      xls_f   <- drive_ls(as_id(tag_f$id[1]), pattern = "\\.xls$")
      drive_download(as_id(xls_f$id[1]), path = xls_local, overwrite = TRUE)
    }

    rerun_hit <- root_ls[root_ls$name == "hmmoce-rerun", ]
    if (nrow(rerun_hit) > 0) {
      rerun_ls <- drive_ls(as_id(rerun_hit$id[1]))
      tag_out  <- rerun_ls[rerun_ls$name == as.character(tag_id), ]
      if (nrow(tag_out) > 0) {
        local_tag_out <- file.path(base_dir, "hmmoce-rerun", tag_id)
        dir.create(local_tag_out, recursive = TRUE, showWarnings = FALSE)
        remote_files <- drive_ls(as_id(tag_out$id[1]))
        for (k in seq_len(nrow(remote_files))) {
          dest <- file.path(local_tag_out, remote_files$name[k])
          if (!file.exists(dest)) {
            drive_download(as_id(remote_files$id[k]), path = dest, overwrite = TRUE)
          }
        }
      }
    }
  }

  push_gdrive_outputs <- function(tag_id, out_dir, files = NULL, root_id = gdrive_root_id) {
    if (!use_gdrive_api) return(invisible(NULL))
    rerun_id <- get_or_create_gdrive_folder(root_id, "hmmoce-rerun")
    tag_id_drive <- get_or_create_gdrive_folder(rerun_id, as.character(tag_id))
    if (is.null(files)) files <- list.files(out_dir, full.names = TRUE)
    files <- files[file.exists(files) & !dir.exists(files)]
    for (f in files) {
      message("Uploading to Google Drive: ", basename(f))
      drive_put(media = f, path = as_id(tag_id_drive), name = basename(f))
    }
  }
}

rerun_dir <- file.path(base_dir, "hmmoce-rerun")
if (!dir.exists(rerun_dir)) dir.create(rerun_dir, recursive = TRUE)

get_cobia_meta <- function(tag_id, dir = base_dir) {
  tag_id <- as.character(tag_id)
  if (use_gdrive_api) sync_gdrive_tag(tag_id, dir)
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

  data.frame(
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
}

# Discover all Cobia tag IDs in "Updated PSAT Reports"
psat_reports_dir <- file.path(base_dir, "Updated PSAT Reports")
if (dir.exists(psat_reports_dir)) {
  tag_dirs <- list.dirs(psat_reports_dir, full.names = FALSE, recursive = FALSE)
  tag_ids  <- sort(sub("_updated.*$", "", tag_dirs[grepl("_updated", tag_dirs)]))
} else {
  tag_ids <- c("36419", "36420", "36421", "36422", "40302", "49157", "55550", "55552")
}
message("Found ", length(tag_ids), " Cobia tags to process: ", paste(tag_ids, collapse = ", "))

# Shared Copernicus CLI and map objects
download_glorys <- TRUE
cm_cli <- Sys.which("copernicusmarine")
if (cm_cli == "") {
  cm_candidates <- c(
    path.expand("~/miniforge3/bin/copernicusmarine"),
    path.expand("~/mambaforge/bin/copernicusmarine"),
    path.expand("~/miniconda3/bin/copernicusmarine"),
    path.expand("~/.local/bin/copernicusmarine"),
    file.path(Sys.getenv("APPDATA"), "Python", "Python313", "Scripts", "copernicusmarine.exe")
  )
  cm_hit <- Find(file.exists, cm_candidates)
  cm_cli <- if (!is.null(cm_hit)) shQuote(cm_hit) else "copernicusmarine"
} else {
  cm_cli <- shQuote(unname(cm_cli))
}
cm_help <- tryCatch(system(paste(cm_cli, "subset --help"), intern = TRUE, ignore.stderr = TRUE), error = function(e) character())
cm_force_flag <- if (any(grepl("--force-download", cm_help, fixed = TRUE))) "--force-download" else "--overwrite"

world_sf <- sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
bathy.dir <- file.path(rerun_dir, "EnvData", "bathy")
if (!dir.exists(bathy.dir)) dir.create(bathy.dir, recursive = TRUE)
glorys.dir <- file.path(rerun_dir, "EnvData", "glorys")
if (!dir.exists(glorys.dir)) dir.create(glorys.dir, recursive = TRUE)

# ==============================================================================
# Single-Tag Processing Function
# ==============================================================================
process_cobia_tag <- function(instrument_id) {
  message("\n############################################################")
  message("PROCESSING COBIA TAG: ", instrument_id)
  message("############################################################")

  tag_meta <- get_cobia_meta(instrument_id, base_dir)
  excel_file <- tag_meta$excel_file

  # 2. Generate and Load eTUFF File --------------------------------------------
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

  # 3. Setup Start / End Bounds & Daily Timesteps ------------------------------
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

  # 4. Extract Observations from eTUFF -----------------------------------------

  ## 4.1 Daily SST
  if ("tempMax" %in% names(etuff$etuff)) {
    tag.sst <- etuff$etuff %>%
      dplyr::filter(!is.na(tempMax)) %>%
      dplyr::mutate(Date = as.Date(DateTime)) %>%
      dplyr::group_by(Date) %>%
      dplyr::summarise(Depth = 0, Temperature = mean(tempMax, na.rm = TRUE)) %>%
      dplyr::mutate(Date = as.POSIXct(Date, tz = "UTC")) %>%
      as.data.frame()
    message("Extracted ", nrow(tag.sst), " daily SST records from eTUFF.")
  } else if ("temperature" %in% names(etuff$etuff)) {
    tag.sst <- etuff$etuff %>%
      dplyr::filter(!is.na(temperature) & (is.na(depth) | depth <= 5)) %>%
      dplyr::mutate(Date = as.Date(DateTime)) %>%
      dplyr::group_by(Date) %>%
      dplyr::summarise(Depth = 0, Temperature = mean(temperature, na.rm = TRUE)) %>%
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
      dplyr::filter(!is.na(longitude) & longitude >= -100 & longitude <= -75) %>%
      dplyr::select(Date = DateTime, Latitude = latitude, Longitude = longitude) %>%
      dplyr::mutate(
        Latitude = ifelse(!is.na(Latitude) & Latitude >= 18 & Latitude <= 31, Latitude, NA_real_),
        Error.Semi.minor.axis = 70000,
        Error.Semi.major.axis = 100000 * 3,
        Offset = 0,
        Offset.orientation = 0
      ) %>%
      as.data.frame()
    message("Extracted ", nrow(lightloc), " light-based location estimates.")
  } else {
    lightloc <- NULL
  }

  ## 4.4 Maximum Daily Depth
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

  # 5. Spatial Domain & Grid Setup ---------------------------------------------
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
  locs.grid <- setup.locs.grid(sp.lim, res = "hycom")

  bathy_file <- file.path(bathy.dir, "bathy.nc")
  bathy <- NULL
  if (file.exists(bathy_file)) {
    bathy_cand <- terra::rast(bathy_file)
    ext_b <- terra::ext(bathy_cand)
    if (terra::xmin(ext_b) <= sp.lim$lonmin && terra::xmax(ext_b) >= sp.lim$lonmax &&
        terra::ymin(ext_b) <= sp.lim$latmin && terra::ymax(ext_b) >= sp.lim$latmax) {
      message("Loading cached bathymetry: ", bathy_file)
      bathy <- bathy_cand
    }
  }
  if (is.null(bathy)) {
    message("Fetching NOAA bathymetry raster covering study area...")
    bathy_dl_lim <- list(lonmin = -98, lonmax = -65, latmin = 18, latmax = 40)
    bathy <- tryCatch(
      get.bath.data(spatLim = bathy_dl_lim, save.dir = bathy.dir, res = 0.5),
      error = function(e) {
        message("Bathymetry remote download fallback: ", e$message)
        if (file.exists(bathy_file)) terra::rast(bathy_file) else NULL
      }
    )
  }

  # 6. Observation Likelihood Calculations -------------------------------------

  ## 6.1 Light / Geolocation Likelihood
  if (!is.null(lightloc) && nrow(lightloc) > 0) {
    light_rds <- file.path(out_dir, "L.light.rds")
    if (!file.exists(light_rds)) {
      message("Calculating light likelihood (L.light)...")
      L.light <- calc.lightloc(
        lightloc = lightloc,
        locs.grid = locs.grid,
        dateVec = dateVec,
        errEll = TRUE,
        lon_only = FALSE
      )
      save_rast_rds(L.light, light_rds)
    } else {
      message("Loading cached L.light: ", light_rds)
      L.light <- load_rast_rds(light_rds)
    }
  } else {
    L.light <- NULL
  }

  ## 6.2 Bathymetry Likelihood
  if (!is.null(bathy) && !is.null(mmd) && nrow(mmd) > 0) {
    bathy_rds <- file.path(out_dir, "L.bathy.rds")
    if (!file.exists(bathy_rds)) {
      message("Calculating bathymetry likelihood (L.bathy)...")
      bathy_crop <- terra::crop(bathy, terra::ext(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax))
      L.bathy <- calc.bathy(
        mmd,
        bathy_crop,
        dateVec = dateVec,
        focalDim = 5,
        sens.err = 5,
        lik.type = "max"
      )
      save_rast_rds(L.bathy, bathy_rds)
    } else {
      message("Loading cached L.bathy: ", bathy_rds)
      L.bathy <- load_rast_rds(bathy_rds)
    }
  } else {
    L.bathy <- NULL
  }

  ## 6.3 GLORYS Ocean Model Data (3D Profile OHC & SST) ------------------------
  temp_product_id <- "cmems_mod_glo_phy_my_0.083deg_P1D-m"
  if (download_glorys) {
    all_nc <- paste0(temp_product_id, "_", format(as.Date(dateVec), "%Y%m%d"), ".nc")
    nc_paths <- file.path(glorys.dir, all_nc)
    # Re-download any file that is missing OR whose spatial extent does not cover sp.lim
    needs_dl <- vapply(nc_paths, function(fp) {
      if (!file.exists(fp)) return(TRUE)
      ex <- tryCatch(terra::ext(terra::rast(fp, lyrs = 1)), error = function(e) NULL)
      if (is.null(ex)) return(TRUE)
      (terra::xmin(ex) > sp.lim$lonmin + 0.1) || (terra::xmax(ex) < sp.lim$lonmax - 0.1) ||
        (terra::ymin(ex) > sp.lim$latmin + 0.1) || (terra::ymax(ex) < sp.lim$latmax - 0.1)
    }, logical(1))
    missing_idx <- which(needs_dl)

    if (length(missing_idx) > 0) {
      # Download full Gulf/SE Atlantic Cobia domain [-98, -75] x [18, 38] so shared files cover all tags
      dl_lim <- list(
        lonmin = min(-98, sp.lim$lonmin),
        lonmax = max(-75, sp.lim$lonmax),
        latmin = min(18,  sp.lim$latmin),
        latmax = max(38,  sp.lim$latmax)
      )
      n_workers <- min(4, parallel::detectCores() - 1, length(missing_idx))
      message("Downloading ", length(missing_idx), " GLORYS daily files using ", n_workers, " parallel workers...")

      cl <- parallel::makeCluster(n_workers)
      parallel::clusterExport(cl, varlist = c("cm_cli", "cm_force_flag", "temp_product_id", "dl_lim", "glorys.dir", "dateVec"), envir = environment())

      dl_res <- parallel::parLapply(cl, missing_idx, function(i) {
        d_str <- format(as.Date(dateVec[i]), "%Y%m%d")
        d_dash <- format(as.Date(dateVec[i]), "%Y-%m-%d")
        out_nc <- paste0(temp_product_id, "_", d_str, ".nc")
        cmd <- sprintf(
          '%s subset -i %s -x %f -X %f -y %f -Y %f -t %s -T %s -z 0 -Z 2000 --variable thetao -o "%s" -f "%s" %s',
          cm_cli, temp_product_id, dl_lim$lonmin, dl_lim$lonmax, dl_lim$latmin, dl_lim$latmax,
          d_dash, d_dash, glorys.dir, out_nc, cm_force_flag
        )
        system(cmd, intern = TRUE)
      })
      parallel::stopCluster(cl)
      message("Parallel GLORYS download complete.")
    } else {
      message("All GLORYS daily files for tag ", instrument_id, " already cached.")
    }
  }

  glorys_files <- list.files(glorys.dir, pattern = "\\.nc$", full.names = TRUE)
  ext_tag <- terra::ext(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax)

  # 6.3.1 GLORYS OHC Profile Likelihood (3D)
  if (length(glorys_files) > 0 && !is.null(pdt) && nrow(pdt) > 0) {
    ohc_rds <- file.path(out_dir, "L.ohc.glorys.rds")
    if (!file.exists(ohc_rds)) {
      message("Calculating GLORYS OHC profile likelihood...")
      L.ohc <- calc.ohc.glorys(
        pdt = pdt,
        filename = temp_product_id,
        ohc.dir = glorys.dir,
        dateVec = dateVec
      )
      L.ohc <- terra::crop(L.ohc, ext_tag)
      save_rast_rds(L.ohc, ohc_rds)
    } else {
      message("Loading cached L.ohc: ", ohc_rds)
      L.ohc <- terra::crop(load_rast_rds(ohc_rds), ext_tag)
    }
  } else {
    L.ohc <- NULL
  }

  # 6.3.2 GLORYS Level-Specific 3D Profile Likelihood
  if (length(glorys_files) > 0 && !is.null(pdt) && nrow(pdt) > 0) {
    glorys_prof_rds <- file.path(out_dir, "L.glorys.rds")
    if (!file.exists(glorys_prof_rds)) {
      message("Calculating GLORYS level-specific 3D profile likelihood...")
      L.glorys <- calc.glorys.par(
        pdt = pdt,
        filename = temp_product_id,
        glorys.dir = glorys.dir,
        dateVec = dateVec,
        use.se = TRUE
      )
      L.glorys <- terra::crop(L.glorys, ext_tag)
      save_rast_rds(L.glorys, glorys_prof_rds)
    } else {
      message("Loading cached L.glorys: ", glorys_prof_rds)
      L.glorys <- terra::crop(load_rast_rds(glorys_prof_rds), ext_tag)
    }
  } else {
    L.glorys <- NULL
  }

  # 6.3.3 Combined PDT Likelihood (L.glorys + L.ohc via make.L.pdt)
  if (!is.null(L.glorys) && !is.null(L.ohc)) {
    pdt_rds <- file.path(out_dir, "L.pdt.rds")
    if (!file.exists(pdt_rds)) {
      message("Combining L.glorys and L.ohc into L.pdt via make.L.pdt...")
      pdt_ras <- list(glorys = L.glorys, ohc = L.ohc)
      pdt_res_idx <- which.min(lapply(pdt_ras, function(x) terra::res(x)[1]))
      pdt_res <- resample.grid(pdt_ras, pdt_ras[[pdt_res_idx]])
      L.pdt <- make.L.pdt(pdt_res$L.rasters, dateVec = dateVec)
      save_rast_rds(L.pdt, pdt_rds)
    } else {
      message("Loading cached L.pdt: ", pdt_rds)
      L.pdt <- terra::crop(load_rast_rds(pdt_rds), ext_tag)
    }
  } else {
    L.pdt <- NULL
  }

  # 6.3.4 GLORYS SST Likelihood
  if (length(glorys_files) > 0 && !is.null(tag.sst) && nrow(tag.sst) > 0) {
    sst_rds <- file.path(out_dir, "L.sst.glorys.rds")
    if (!file.exists(sst_rds)) {
      message("Calculating GLORYS SST likelihood...")
      L.sst <- calc.sst.par.glorys(
        tag.sst,
        filename = temp_product_id,
        sst.dir = glorys.dir,
        dateVec = dateVec,
        sens.err = 1
      )
      L.sst <- terra::crop(L.sst, ext_tag)
      save_rast_rds(L.sst, sst_rds)
    } else {
      message("Loading cached L.sst: ", sst_rds)
      L.sst <- terra::crop(load_rast_rds(sst_rds), ext_tag)
    }
  } else {
    L.sst <- NULL
  }

  ## 6.4 Combine Available Likelihood Rasters ----------------------------------
  L.rasters <- list(light = L.light, ohc = L.ohc, glorys = L.glorys, pdt = L.pdt, bathy = L.bathy, sst = L.sst)
  L.rasters <- L.rasters[!sapply(L.rasters, is.null)]
  message("Combining observation likelihoods: ", paste(names(L.rasters), collapse = ", "))

  resamp.idx <- which.min(lapply(L.rasters, function(x) terra::res(x)[1]))
  L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])

  if ("light" %in% names(L.res$L.rasters)) {
    ocean_src <- intersect(c("glorys", "sst", "ohc"), names(L.res$L.rasters))
    if (length(ocean_src) > 0) {
      ocean_mask <- max(L.res$L.rasters[[ocean_src[1]]], na.rm = TRUE) > 0
    } else if (!is.null(bathy)) {
      b_res <- if (terra::compareGeom(bathy, L.res$L.rasters$light, stopOnError = FALSE)) bathy else terra::resample(bathy, L.res$L.rasters$light)
      ocean_mask <- !is.na(b_res) & (b_res < 0)
    } else {
      ocean_mask <- 1
    }
    light_max <- terra::global(L.res$L.rasters$light, "max", na.rm = TRUE)[, 1]
    light_active <- which(is.finite(light_max) & light_max > 0)
    for (ii in light_active) {
      lyr <- L.res$L.rasters$light[[ii]]
      lyr[is.na(lyr) | lyr < 1e-3] <- 1e-3
      L.res$L.rasters$light[[ii]] <- lyr * ocean_mask
    }
  }

  maxDepth_vec <- rep(0, length(dateVec))
  weights <- c(light = 1.0, sst = 1.0, ohc = 1.0, glorys = 1.0, pdt = 1.0, bathy = 1.0)
  ras_weighted <- Map(function(r, nm) r ^ weights[[nm]], L.res$L.rasters, names(L.res$L.rasters))

  L <- make.L(
    ras.list = ras_weighted,
    iniloc = iniloc,
    dateVec = dateVec,
    maxDepth = maxDepth_vec,
    bathy = if (any(c("glorys", "sst", "ohc") %in% names(L.res$L.rasters))) NULL else bathy
  )
  saveRDS(L, file.path(out_dir, "L.master.rds"))
  message("Master observation likelihood L constructed: dim = ", paste(dim(L), collapse = " x "))

  # 7. Model Comparison with Movement Parameters -------------------------------
  fit <- tryCatch(
    opt.params(
      pars.init    = c(4, 1, 0.85, 0.85),
      lower.bounds = c(2, 0.2, 0.2, 0.2),
      upper.bounds = c(10, 4.0, 0.95, 0.95),
      g = L.res$g, L = L, alg.opt = "rtmb"
    ),
    error = function(e) {
      message("opt.params fallback to default parameters for tag ", instrument_id, ": ", e$message)
      list(par = c(4, 1, 0.85, 0.85))
    }
  )

  sigmas <- fit$par[1:2]
  sizes  <- rep(ceiling(sigmas[1] * 6), 2)
  if (sizes[1] %% 2 == 0) sizes[1] <- sizes[1] + 1
  if (sizes[2] %% 2 == 0) sizes[2] <- sizes[2] + 1
  muadvs <- c(0, 0)

  pb <- fit$par[3:4]
  P  <- matrix(c(pb[1], 1 - pb[1], 1 - pb[2], pb[2]), nrow = 2, ncol = 2, byrow = TRUE)

  K1 <- HMMoce:::gausskern.pg(sizes[1], sigmas[1], muadv = muadvs[1])
  K1 <- HMMoce:::mask.K(K1)
  K2 <- HMMoce:::gausskern.pg(sizes[2], sigmas[2], muadv = muadvs[2])
  K2 <- HMMoce:::mask.K(K2)
  K  <- list(K1, K2)
  m_states <- 2

  combos <- list(
    ohc_light        = c("ohc", "light"),
    glorys_light     = c("glorys", "light"),
    glorys_light_sst = c("glorys", "light", "sst"),
    ohc_light_sst    = c("ohc", "light", "sst"),
    pdt_light        = c("pdt", "light"),
    ohc              = c("ohc"),
    glorys           = c("glorys"),
    light            = c("light"),
    light_bathy      = c("light", "bathy"),
    light_bathy_sst  = c("light", "bathy", "sst"),
    light_ohc_bathy  = c("light", "bathy", "ohc"),
    light_pdt_bathy  = c("light", "bathy", "pdt")
  )
  combos <- combos[sapply(combos, function(x) all(x %in% names(L.res$L.rasters)))]

  model_results <- list()
  model_fits <- list()
  comp_rows <- list()

  for (m_name in names(combos)) {
    sel_liks <- combos[[m_name]]
    message("\n--------------------------------------------------")
    message("Tag ", instrument_id, " | Model: ", m_name, " (", paste(sel_liks, collapse = " + "), ")")
    message("--------------------------------------------------")

    use_bathy <- "bathy" %in% sel_liks
    L_sub <- make.L(
      ras.list = ras_weighted[sel_liks],
      iniloc = iniloc,
      dateVec = dateVec,
      maxDepth = if (use_bathy) maxDepth_vec else rep(0, length(dateVec)),
      bathy = if (use_bathy) bathy else NULL
    )

    f <- hmm.filter(g = L.res$g, L = L_sub, K = K, P = P, m = m_states)
    nllf <- -sum(log(f$psi[f$psi > 0]))
    aic  <- 2 * nllf + 2 * length(sigmas)

    s <- hmm.smoother(f, K = K, L = L_sub, P = P)

    tr <- calc.track(s, g = L.res$g, dateVec = dateVec, iniloc = iniloc, method = "mean")
    tr$Date <- as.POSIXct(dateVec, tz = "UTC")
    tr$model <- m_name
    tr$tag_id <- instrument_id

    fit_obj <- list(model = m_name, likelihoods = sel_liks, pars = c(sigmas, pb), s = s, tr = tr, aic = aic, nll = nllf)
    fit_file <- file.path(out_dir, paste0(instrument_id, "_", m_name, "_fixed_fit.rds"))
    track_file <- file.path(out_dir, paste0(instrument_id, "_", m_name, "_fixed_track.csv"))
    saveRDS(fit_obj, fit_file)
    write.csv(tr, track_file, row.names = FALSE)

    last_idx <- nrow(tr)
    end_err_km <- terra::distance(
      matrix(c(tr$lon[last_idx], tr$lat[last_idx]), ncol = 2),
      matrix(c(iniloc$lon[2], iniloc$lat[2]), ncol = 2),
      lonlat = TRUE
    )[1] / 1000

    comp_rows[[m_name]] <- data.frame(
      tag_id = instrument_id,
      model = m_name,
      likelihoods = paste(sel_liks, collapse = "+"),
      sigma1 = round(sigmas[1], 3),
      sigma2 = round(sigmas[2], 3),
      p11 = round(pb[1], 3),
      p22 = round(pb[2], 3),
      nll = round(nllf, 2),
      aic = round(aic, 2),
      end_error_km = round(end_err_km, 1),
      stringsAsFactors = FALSE
    )

    model_fits[[m_name]] <- fit_obj
    model_results[[m_name]] <- tr
  }

  comp_table <- do.call(rbind, comp_rows)
  comp_table <- comp_table[order(comp_table$aic), ]
  write.csv(comp_table, file.path(out_dir, paste0(instrument_id, "_fixed_model_comparison.csv")), row.names = FALSE)
  print(comp_table)

  # 8. Post-Run Plotting for All Models ----------------------------------------
  all_tracks <- do.call(rbind, model_results)

  p_comp <- ggplot() +
    geom_sf(data = world_sf, fill = "grey80", color = "grey50", linewidth = 0.3) +
    geom_path(data = all_tracks, aes(x = lon, y = lat, color = model), linewidth = 0.9, alpha = 0.85) +
    geom_point(data = all_tracks, aes(x = lon, y = lat), fill = "black", color = "black", shape = 19, size = 0.4) +
    geom_point(data = iniloc[1, ], aes(x = lon, y = lat), fill = "green", color = "black", shape = 24, size = 4) +
    geom_point(data = iniloc[2, ], aes(x = lon, y = lat), fill = "red", color = "black", shape = 25, size = 4) +
    coord_sf(xlim = c(sp.lim$lonmin, sp.lim$lonmax), ylim = c(sp.lim$latmin, sp.lim$latmax), expand = FALSE) +
    theme_bw(base_size = 13) +
    theme(axis.text = element_text(size = 11, color = "black")) +
    labs(
      title = paste("Track Comparison across Likelihood Combinations - Tag", instrument_id),
      subtitle = "Green = Release, Red = Popoff",
      x = "Longitude", y = "Latitude", color = "Model"
    )

  if (!is.null(lightloc) && nrow(lightloc) > 0) {
    p_comp <- p_comp +
      geom_point(data = lightloc, aes(x = Longitude, y = Latitude), color = "orange", alpha = 0.3, size = 1)
  }

  comp_plot_file <- file.path(out_dir, paste0(instrument_id, "_fixed_track_comparison.png"))
  ggsave(comp_plot_file, plot = p_comp, width = 10, height = 8, dpi = 300)

  for (m_name in names(model_fits)) {
    fit_m <- model_fits[[m_name]]
    tr_m <- fit_m$tr
    tr_m$Month <- factor(format(tr_m$Date, "%b %Y"), levels = unique(format(tr_m$Date, "%b %Y")))
    sel_liks_str <- paste(fit_m$likelihoods, collapse = " + ")

    plot_file <- file.path(out_dir, paste0(instrument_id, "_", m_name, "_fixed_plotHMM.png"))
    png(plot_file, width = 8, height = 9, units = "in", res = 300)
    tryCatch(
      plotHMM(distr = fit_m$s, track = tr_m, dateVec = dateVec, ptt = paste(instrument_id, m_name, "fixed"), behav.pts = TRUE, save.plot = FALSE),
      finally = dev.off()
    )

    p_track <- ggplot() +
      geom_sf(data = world_sf, fill = "grey80", color = "grey50", linewidth = 0.3)

    if (!is.null(lightloc) && nrow(lightloc) > 0) {
      p_track <- p_track +
        geom_point(data = lightloc, aes(x = Longitude, y = Latitude), color = "orange", alpha = 0.6, size = 1.4)
    }

    p_track <- p_track +
      geom_path(data = tr_m, aes(x = lon, y = lat), color = "grey35", linewidth = 0.8) +
      geom_point(
        data = tr_m,
        aes(x = lon, y = lat, fill = Month),
        shape = 21, color = "black", stroke = 0.35, size = 3
      ) +
      scale_fill_viridis_d(option = "turbo", name = "Month", begin = 0.1, end = 0.9) +
      geom_point(data = iniloc[1, ], aes(x = lon, y = lat), fill = "green", color = "black", shape = 24, size = 4) +
      geom_point(data = iniloc[2, ], aes(x = lon, y = lat), fill = "red", color = "black", shape = 25, size = 4) +
      coord_sf(xlim = c(sp.lim$lonmin, sp.lim$lonmax), ylim = c(sp.lim$latmin, sp.lim$latmax), expand = FALSE) +
      theme_bw(base_size = 13) +
      theme(
        axis.text = element_text(size = 11, color = "black"),
        plot.subtitle = element_text(size = 11, face = "italic")
      ) +
      labs(
        title = paste("HMMoce Estimated Track - Cobia", instrument_id),
        subtitle = paste0("Model: ", m_name, " (", sel_liks_str, ")"),
        x = "Longitude", y = "Latitude"
      )

    track_plot_file <- file.path(out_dir, paste0(instrument_id, "_", m_name, "_track_map.png"))
    ggsave(track_plot_file, plot = p_track, width = 9, height = 7, dpi = 300)
  }

  if (use_gdrive_api) push_gdrive_outputs(instrument_id, out_dir)
  return(comp_table)
}

# ==============================================================================
# Run Batch Loop Across All Tags
# ==============================================================================
batch_summary <- list()
batch_errors  <- list()

for (tid in tag_ids) {
  res <- tryCatch(
    process_cobia_tag(tid),
    error = function(e) {
      message("ERROR processing tag ", tid, ": ", e$message)
      batch_errors[[tid]] <<- e$message
      NULL
    }
  )
  if (!is.null(res)) batch_summary[[tid]] <- res
}

if (length(batch_summary) > 0) {
  master_comp <- do.call(rbind, batch_summary)
  master_file <- file.path(rerun_dir, "all_tags_fixed_model_comparison.csv")
  write.csv(master_comp, master_file, row.names = FALSE)
  message("\n============================================================")
  message("Batch complete! Saved master comparison table to: ", master_file)
  message("============================================================")
  print(master_comp)
}

if (length(batch_errors) > 0) {
  message("\nTags with errors:")
  print(batch_errors)
}
