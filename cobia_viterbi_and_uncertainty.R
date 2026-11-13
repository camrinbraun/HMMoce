# ==============================================================================
# Cobia PSAT Posterior Viterbi Track Re-Estimation & Posterior Uncertainty Script
# Model: glorys_light (3D GLORYS temperature-depth profile + Light geolocation)
# Produces:
#   1. Viterbi most probable track via hmm.posteriorviterbi()
#   2. Daily posterior uncertainty table (SD, 50% & 95% credible intervals, 95% HPD area)
#   3. Comparison map: Mean vs. Viterbi track + 50%/95% HPD contours + daily 95% ellipses
#   4. Posterior Utilization Distribution (UD) heat map with 50%/95% HPD volume contours
#   5. Lat/Lon time-series plot with 50% & 95% posterior credible ribbons + Mean vs. Viterbi
# ==============================================================================

library(dplyr)
library(lubridate)
library(ggplot2)
library(sf)
library(terra)
devtools::load_all(".")

load_rast_rds <- function(file) {
  x <- readRDS(file)
  if (inherits(x, "PackedSpatRaster")) terra::unwrap(x) else if (inherits(x, "SpatRaster")) x else terra::rast(x)
}

# 1. Directory & Metadata Setup ------------------------------------------------
gdrive_root_id <- "1jN4NVXfNcxVqIPxVvl3k3e3e4W3J5HgU"
local_gdrive   <- "G:/.shortcut-targets-by-id/1jN4NVXfNcxVqIPxVvl3k3e3e4W3J5HgU/Cobia Shared Files"
base_dir       <- if (dir.exists(local_gdrive)) local_gdrive else file.path(getwd(), "cobia_gdrive_cache")
rerun_dir      <- file.path(base_dir, "hmmoce-rerun")

get_cobia_meta <- function(tag_id, dir = base_dir) {
  tag_id <- as.character(tag_id)
  tracks_file <- file.path(dir, "Data files", "Cobia tracks.csv")
  tracks <- read.csv(tracks_file)
  tag_track <- tracks[tracks[, 1] == as.numeric(tag_id), ]
  tag_track$Date <- sprintf("%04d-%02d-%02d", tag_track$Year, tag_track$Month, tag_track$Day)
  tag_track <- tag_track[order(tag_track$Date), ]

  data.frame(
    instrument_name      = tag_id,
    time_coverage_start  = as.POSIXct(paste(tag_track$Date[1], "00:00:00"), tz = "UTC"),
    time_coverage_end    = as.POSIXct(paste(tag_track$Date[nrow(tag_track)], "00:00:00"), tz = "UTC"),
    geospatial_lat_start = tag_track$Lat_N[1],
    geospatial_lon_start = tag_track$Lon_E[1],
    geospatial_lat_end   = tag_track$Lat_N[nrow(tag_track)],
    geospatial_lon_end   = tag_track$Lon_E[nrow(tag_track)],
    stringsAsFactors     = FALSE
  )
}

# Helper: weighted empirical quantile from a 1D marginal probability vector
weighted_quantile <- function(x, w, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  w[is.na(w) | w < 0] <- 0
  s <- sum(w)
  if (s <= 0) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  cdf <- cumsum(w[ord]) / s
  vapply(probs, function(p) x[which(cdf >= p)[1]], numeric(1))
}

# Helper: convert a 2D probability matrix (lon x lat) to an HPD volume raster
# where cell values represent the cumulative probability volume (0 = peak, 0.95 = 95% HPD contour)
calc_hpd_volume_matrix <- function(prob_mat) {
  prob_mat[is.na(prob_mat) | prob_mat < 0] <- 0
  tot <- sum(prob_mat)
  if (tot <= 0) return(prob_mat * NA_real_)
  p_vec <- as.vector(prob_mat) / tot
  ord   <- order(p_vec, decreasing = TRUE)
  v_vec <- numeric(length(p_vec))
  v_vec[ord] <- cumsum(p_vec[ord])
  matrix(v_vec, nrow = nrow(prob_mat), ncol = ncol(prob_mat))
}

# Helper: build 2D covariance ellipse polygon points for a given day
make_ellipse_df <- function(mu_lon, mu_lat, cov_mat, date_val, level = 0.95, n_pts = 60) {
  if (any(!is.finite(c(mu_lon, mu_lat, cov_mat)))) return(NULL)
  eig <- eigen(cov_mat, symmetric = TRUE)
  eig$values[eig$values < 0] <- 0
  if (max(eig$values) < 1e-6) return(NULL)
  r_chi <- sqrt(qchisq(level, df = 2))
  theta <- seq(0, 2 * pi, length.out = n_pts)
  circle <- rbind(cos(theta), sin(theta))
  pts <- r_chi * (eig$vectors %*% (sqrt(eig$values) * circle))
  data.frame(
    Date = date_val,
    lon  = mu_lon + pts[1, ],
    lat  = mu_lat + pts[2, ]
  )
}

# Discover all Cobia tag IDs
psat_reports_dir <- file.path(base_dir, "Updated PSAT Reports")
if (dir.exists(psat_reports_dir)) {
  tag_dirs <- list.dirs(psat_reports_dir, full.names = FALSE, recursive = FALSE)
  tag_ids  <- sort(sub("_updated.*$", "", tag_dirs[grepl("_updated", tag_dirs)]))
} else {
  tag_ids <- c("36419", "36420", "36421", "36422", "40302", "49157", "55550", "55552")
}

world_sf <- sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))

# ==============================================================================
# Per-Tag Viterbi & Posterior Uncertainty Function
# ==============================================================================
run_viterbi_and_uncertainty <- function(instrument_id, model_name = "glorys_light", percentile_v = 99.9, approach_v = 'joint') {
  out_dir  <- file.path(rerun_dir, instrument_id)
  fit_file <- file.path(out_dir, paste0(instrument_id, "_", model_name, "_fixed_fit.rds"))
  if (!file.exists(fit_file)) {
    message("Skipping tag ", instrument_id, ": fit file not found at ", fit_file)
    return(NULL)
  }

  message("\n============================================================")
  message("Tag ", instrument_id, " | Running Viterbi & Uncertainty for ", model_name)
  message("============================================================")

  fit_obj <- readRDS(fit_file)
  s       <- fit_obj$s
  pars    <- fit_obj$pars
  tr_mean <- fit_obj$tr

  tag_meta <- get_cobia_meta(instrument_id, base_dir)
  t_start  <- tag_meta$time_coverage_start
  t_end    <- tag_meta$time_coverage_end

  iniloc <- data.frame(
    day   = c(lubridate::day(t_start), lubridate::day(t_end)),
    month = c(lubridate::month(t_start), lubridate::month(t_end)),
    year  = c(lubridate::year(t_start), lubridate::year(t_end)),
    lat   = c(as.numeric(tag_meta$geospatial_lat_start), as.numeric(tag_meta$geospatial_lat_end)),
    lon   = c(as.numeric(tag_meta$geospatial_lon_start), as.numeric(tag_meta$geospatial_lon_end))
  )
  rownames(iniloc) <- c("start", "end")

  tag     <- as.POSIXct(sprintf("%04d-%02d-%02d", iniloc["start", "year"], iniloc["start", "month"], iniloc["start", "day"]), tz = "UTC")
  pop     <- as.POSIXct(sprintf("%04d-%02d-%02d", iniloc["end", "year"], iniloc["end", "month"], iniloc["end", "day"]), tz = "UTC")
  dateVec <- seq.POSIXt(tag, pop, by = "24 hours")

  # Reconstruct spatial grid g matching s dimensions (m x T x ncol x nrow)
  glorys_rds <- file.path(out_dir, "L.glorys.rds")
  if (file.exists(glorys_rds)) {
    r_tmpl <- load_rast_rds(glorys_rds)[[1]]
    cand_files <- file.path(out_dir, c("L.light.rds", "L.bathy.rds", "L.ohc.glorys.rds", "L.sst.glorys.rds"))
    cand_files <- cand_files[file.exists(cand_files)]
    if ((ncol(r_tmpl) != dim(s)[3] || nrow(r_tmpl) != dim(s)[4]) && length(cand_files) > 0) {
      exts <- lapply(cand_files, function(f) terra::ext(load_rast_rds(f)))
      common_ext <- Reduce(terra::intersect, c(list(terra::ext(r_tmpl)), exts))
      r_tmpl <- terra::crop(r_tmpl, common_ext)
    }
    if (ncol(r_tmpl) != dim(s)[3] || nrow(r_tmpl) != dim(s)[4]) {
      r_tmpl <- terra::rast(terra::ext(r_tmpl), ncols = dim(s)[3], nrows = dim(s)[4], crs = terra::crs(r_tmpl))
    }
    g <- setup.grid.raster(r_tmpl)
  } else {
    stop("Template raster L.glorys.rds not found for tag ", instrument_id)
  }

  lon_vec <- g$lon[1, ]
  lat_vec <- g$lat[, 1]
  sp.lim  <- list(lonmin = min(lon_vec), lonmax = max(lon_vec), latmin = min(lat_vec), latmax = max(lat_vec))

  # 1. Run Posterior Viterbi ---------------------------------------------------
  message("Running hmm.posteriorviterbi(approach = 'joint')...")
  tr_vit <- hmm.posteriorviterbi(
    pars       = pars,
    s          = s,
    g          = g,
    iniloc     = iniloc,
    dateVec    = dateVec,
    lims.lat   = range(lat_vec),
    lims.lon   = range(lon_vec),
    percentile = percentile_v,
    approach   =  approach_v
  )
  tr_vit$Date   <- as.POSIXct(tr_vit$date, tz = "UTC")
  tr_vit$method <- "Posterior Viterbi"
  tr_mean$Date  <- as.POSIXct(tr_mean$Date, tz = "UTC")
  tr_mean$method <- "Posterior Mean"

  vit_file <- file.path(out_dir, paste0(instrument_id, "_", model_name, "_viterbi_track.csv"))
  write.csv(tr_vit, vit_file, row.names = FALSE)
  message("Saved Viterbi track to: ", vit_file)

  # 2. Quantify Daily & Cumulative Posterior Uncertainty -----------------------
  T_steps <- dim(s)[2]
  nx      <- length(lon_vec)
  ny      <- length(lat_vec)
  dx_deg  <- abs(mean(diff(lon_vec)))
  dy_deg  <- abs(mean(diff(lat_vec)))

  # Daily state-combined posterior P_t(lon, lat)
  post_daily <- apply(s, c(2, 3, 4), sum, na.rm = TRUE) # dim: (T, nx, ny)
  ud_sum     <- matrix(0, nrow = nx, ncol = ny)

  unc_rows    <- vector("list", T_steps)
  ellipse_lst <- vector("list", T_steps)
  lon_grid_mat <- matrix(lon_vec, nrow = nx, ncol = ny)
  lat_grid_mat <- matrix(lat_vec, nrow = nx, ncol = ny, byrow = TRUE)

  for (t in seq_len(T_steps)) {
    P_t <- post_daily[t, , ]
    P_t[is.na(P_t) | P_t < 0] <- 0
    tot_t <- sum(P_t)
    if (tot_t > 0) P_t <- P_t / tot_t
    ud_sum <- ud_sum + P_t

    # Marginal distributions over lon and lat
    p_lon <- rowSums(P_t)
    p_lat <- colSums(P_t)

    mu_lon <- sum(lon_vec * p_lon)
    mu_lat <- sum(lat_vec * p_lat)
    sd_lon <- sqrt(max(0, sum(((lon_vec - mu_lon) ^ 2) * p_lon)))
    sd_lat <- sqrt(max(0, sum(((lat_vec - mu_lat) ^ 2) * p_lat)))
    cov_ll <- sum((lon_grid_mat - mu_lon) * (lat_grid_mat - mu_lat) * P_t)

    q_lon <- weighted_quantile(lon_vec, p_lon, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    q_lat <- weighted_quantile(lat_vec, p_lat, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

    # Daily 95% HPD area in km^2
    hpd_t <- calc_hpd_volume_matrix(P_t)
    cell_area_km2 <- (dx_deg * 111.32 * cos(mu_lat * pi / 180)) * (dy_deg * 110.57)
    n_cells_95 <- sum(!is.na(hpd_t) & hpd_t <= 0.95)
    n_cells_50 <- sum(!is.na(hpd_t) & hpd_t <= 0.50)

    unc_rows[[t]] <- data.frame(
      tag_id        = instrument_id,
      Date          = dateVec[t],
      mean_lon      = mu_lon,
      mean_lat      = mu_lat,
      viterbi_lon   = tr_vit$lon[t],
      viterbi_lat   = tr_vit$lat[t],
      sd_lon_deg    = sd_lon,
      sd_lat_deg    = sd_lat,
      lon_q025      = q_lon[1],
      lon_q250      = q_lon[2],
      lon_q500      = q_lon[3],
      lon_q750      = q_lon[4],
      lon_q975      = q_lon[5],
      lat_q025      = q_lat[1],
      lat_q250      = q_lat[2],
      lat_q500      = q_lat[3],
      lat_q750      = q_lat[4],
      lat_q975      = q_lat[5],
      hpd50_area_km2 = round(n_cells_50 * cell_area_km2, 1),
      hpd95_area_km2 = round(n_cells_95 * cell_area_km2, 1),
      p_resid       = tr_vit$p.resid[t],
      stringsAsFactors = FALSE
    )

    # Sample weekly 95% credible ellipses (every 5 days) to keep map legible
    if (t %% 5 == 0 && t > 1 && t < T_steps) {
      cov_mat <- matrix(c(sd_lon ^ 2, cov_ll, cov_ll, sd_lat ^ 2), 2, 2)
      ellipse_lst[[t]] <- make_ellipse_df(mu_lon, mu_lat, cov_mat, dateVec[t], level = 0.95)
    }
  }

  unc_df <- do.call(rbind, unc_rows)
  unc_file <- file.path(out_dir, paste0(instrument_id, "_", model_name, "_posterior_uncertainty.csv"))
  write.csv(unc_df, unc_file, row.names = FALSE)
  message("Saved daily posterior uncertainty metrics to: ", unc_file)

  # Build cumulative Posterior Utilization Distribution (UD) & HPD volume SpatRaster
  ud_norm  <- ud_sum / sum(ud_sum, na.rm = TRUE)
  hpd_vol  <- calc_hpd_volume_matrix(ud_norm) # 0 at peak density -> 1 at tail

  ext_g <- terra::ext(min(lon_vec) - dx_deg / 2, max(lon_vec) + dx_deg / 2,
                      min(lat_vec) - dy_deg / 2, max(lat_vec) + dy_deg / 2)
  # Convert (nx, ny) South-up to (ny, nx) North-up for terra::rast
  ud_rast  <- terra::rast(t(ud_norm / max(ud_norm, na.rm = TRUE))[ny:1, , drop = FALSE], extent = ext_g, crs = "EPSG:4326")
  hpd_rast <- terra::rast(t(hpd_vol)[ny:1, , drop = FALSE], extent = ext_g, crs = "EPSG:4326")

  # Extract 50% and 95% HPD contour lines as sf
  hpd_contours_sf <- tryCatch({
    cl <- terra::as.contour(hpd_rast, levels = c(0.50, 0.95))
    sf_obj <- sf::st_as_sf(cl)
    sf_obj$Isopleth <- factor(paste0(round(sf_obj$level * 100), "% HPD"), levels = c("50% HPD", "95% HPD"))
    sf_obj
  }, error = function(e) NULL)

  # Convert UD raster to data.frame for ggplot
  ud_df <- as.data.frame(ud_rast, xy = TRUE)
  names(ud_df) <- c("lon", "lat", "rel_prob")
  hpd_df <- as.data.frame(hpd_rast, xy = TRUE)
  names(hpd_df) <- c("lon", "lat", "hpd_vol")
  ud_df$hpd_vol <- hpd_df$hpd_vol
  ud_df <- ud_df[!is.na(ud_df$hpd_vol) & ud_df$hpd_vol <= 0.99, ]

  # Combine Mean & Viterbi tracks for plotting
  both_tracks <- rbind(
    data.frame(Date = tr_mean$Date, lon = tr_mean$lon, lat = tr_mean$lat, method = "Posterior Mean", stringsAsFactors = FALSE),
    data.frame(Date = tr_vit$Date,  lon = tr_vit$lon,  lat = tr_vit$lat,  method = "Posterior Viterbi", stringsAsFactors = FALSE)
  )
  both_tracks$Month <- factor(format(both_tracks$Date, "%b %Y"), levels = unique(format(both_tracks$Date, "%b %Y")))
  ellipses_df <- do.call(rbind, ellipse_lst)

  # Endpoint errors (km)
  err_mean_km <- terra::distance(
    matrix(c(tr_mean$lon[nrow(tr_mean)], tr_mean$lat[nrow(tr_mean)]), ncol = 2),
    matrix(c(iniloc$lon[2], iniloc$lat[2]), ncol = 2), lonlat = TRUE
  )[1] / 1000
  err_vit_km <- terra::distance(
    matrix(c(tr_vit$lon[nrow(tr_vit)], tr_vit$lat[nrow(tr_vit)]), ncol = 2),
    matrix(c(iniloc$lon[2], iniloc$lat[2]), ncol = 2), lonlat = TRUE
  )[1] / 1000

  # 3. Plot 1: Mean vs. Viterbi Track + Daily 95% Posterior Ellipses & HPD -----
  p_comp <- ggplot() +
    geom_tile(
      data = ud_df[ud_df$hpd_vol <= 0.95, ],
      aes(x = lon, y = lat, fill = ifelse(hpd_vol <= 0.50, "50% Core HPD", "95% Credible HPD")),
      alpha = 0.35
    ) +
    scale_fill_manual(
      name   = "Posterior Region",
      values = c("50% Core HPD" = "#3182bd", "95% Credible HPD" = "#bdd7e7")
    )

  if (!is.null(ellipses_df) && nrow(ellipses_df) > 0) {
    p_comp <- p_comp +
      geom_path(
        data = ellipses_df,
        aes(x = lon, y = lat, group = Date),
        color = "steelblue4", linetype = "dotted", linewidth = 0.4, alpha = 0.65
      )
  }

  p_comp <- p_comp +
    geom_sf(data = world_sf, fill = "grey80", color = "grey50", linewidth = 0.3) +
    geom_path(
      data = both_tracks,
      aes(x = lon, y = lat, color = method, linetype = method),
      linewidth = 1.0
    ) +
    geom_point(
      data = both_tracks,
      aes(x = lon, y = lat, color = method, shape = method),
      size = 2.0, alpha = 0.85
    ) +
    scale_color_manual(
      name   = "Track Estimator",
      values = c("Posterior Mean" = "#0072B2", "Posterior Viterbi" = "#D55E00")
    ) +
    scale_linetype_manual(
      name   = "Track Estimator",
      values = c("Posterior Mean" = "solid", "Posterior Viterbi" = "twodash")
    ) +
    scale_shape_manual(
      name   = "Track Estimator",
      values = c("Posterior Mean" = 16, "Posterior Viterbi" = 17)
    ) +
    geom_point(data = iniloc[1, ], aes(x = lon, y = lat), fill = "green", color = "black", shape = 24, size = 4) +
    geom_point(data = iniloc[2, ], aes(x = lon, y = lat), fill = "red", color = "black", shape = 25, size = 4) +
    coord_sf(xlim = c(sp.lim$lonmin, sp.lim$lonmax), ylim = c(sp.lim$latmin, sp.lim$latmax), expand = FALSE) +
    theme_bw(base_size = 13) +
    theme(
      axis.text = element_text(size = 11, color = "black"),
      plot.subtitle = element_text(size = 11, face = "italic")
    ) +
    labs(
      title    = paste("Mean vs. Viterbi Track & Posterior Uncertainty - Cobia", instrument_id),
      subtitle = sprintf("Model: glorys_light | Endpoint Error: Mean = %.1f km, Viterbi = %.1f km", err_mean_km, err_vit_km),
      x = "Longitude", y = "Latitude"
    )

  map_comp_file <- file.path(out_dir, paste0(instrument_id, "_", model_name, "_mean_vs_viterbi_map.png"))
  ggsave(map_comp_file, plot = p_comp, width = 9.5, height = 7.5, dpi = 300)
  message("Saved Mean vs. Viterbi comparison map: ", map_comp_file)

  # 4. Plot 2: Continuous Posterior UD Heatmap + 50%/95% HPD Contours ---------
  p_ud <- ggplot() +
    geom_tile(data = ud_df, aes(x = lon, y = lat, fill = 100 * (1 - hpd_vol))) +
    scale_fill_viridis_c(
      option = "magma", direction = -1,
      name = "Posterior\nDensity (%)"
    )

  if (!is.null(hpd_contours_sf) && nrow(hpd_contours_sf) > 0) {
    p_ud <- p_ud +
      geom_sf(
        data = hpd_contours_sf,
        aes(linetype = Isopleth),
        color = "cyan3", linewidth = 0.65
      ) +
      scale_linetype_manual(name = "Credible Contour", values = c("50% HPD" = "solid", "95% HPD" = "dashed"))
  }

  p_ud <- p_ud +
    geom_sf(data = world_sf, fill = "grey80", color = "grey50", linewidth = 0.3) +
    geom_path(
      data = both_tracks,
      aes(x = lon, y = lat, color = method),
      linewidth = 0.9
    ) +
    scale_color_manual(
      name   = "Track Estimator",
      values = c("Posterior Mean" = "white", "Posterior Viterbi" = "#00FFCC")
    ) +
    geom_point(data = iniloc[1, ], aes(x = lon, y = lat), fill = "green", color = "black", shape = 24, size = 4) +
    geom_point(data = iniloc[2, ], aes(x = lon, y = lat), fill = "red", color = "black", shape = 25, size = 4) +
    coord_sf(xlim = c(sp.lim$lonmin, sp.lim$lonmax), ylim = c(sp.lim$latmin, sp.lim$latmax), expand = FALSE) +
    theme_bw(base_size = 13) +
    theme(
      axis.text = element_text(size = 11, color = "black"),
      plot.subtitle = element_text(size = 11, face = "italic")
    ) +
    labs(
      title    = paste("Cumulative Posterior Utilization Distribution (UD) - Cobia", instrument_id),
      subtitle = paste0("Model: ", model_name, " (50% Core & 95% Credible Contours)"),
      x = "Longitude", y = "Latitude"
    )

  ud_plot_file <- file.path(out_dir, paste0(instrument_id, "_", model_name, "_posterior_UD_map.png"))
  ggsave(ud_plot_file, plot = p_ud, width = 9.5, height = 7.5, dpi = 300)
  message("Saved Posterior UD map: ", ud_plot_file)

  # 5. Plot 3: Daily Lat/Lon Posterior 50% & 95% Credible Ribbons --------------
  ts_lon <- data.frame(
    Date       = unc_df$Date,
    Coordinate = "Longitude (°W)",
    mean_est   = unc_df$mean_lon,
    vit_est    = unc_df$viterbi_lon,
    q025       = unc_df$lon_q025,
    q250       = unc_df$lon_q250,
    q750       = unc_df$lon_q750,
    q975       = unc_df$lon_q975
  )
  ts_lat <- data.frame(
    Date       = unc_df$Date,
    Coordinate = "Latitude (°N)",
    mean_est   = unc_df$mean_lat,
    vit_est    = unc_df$viterbi_lat,
    q025       = unc_df$lat_q025,
    q250       = unc_df$lat_q250,
    q750       = unc_df$lat_q750,
    q975       = unc_df$lat_q975
  )
  ts_both <- rbind(ts_lon, ts_lat)
  ts_lines <- rbind(
    data.frame(Date = ts_both$Date, Coordinate = ts_both$Coordinate, val = ts_both$mean_est, method = "Posterior Mean"),
    data.frame(Date = ts_both$Date, Coordinate = ts_both$Coordinate, val = ts_both$vit_est,  method = "Posterior Viterbi")
  )

  p_ts <- ggplot(ts_both, aes(x = Date)) +
    geom_ribbon(aes(ymin = q025, ymax = q975, fill = "95% Credible Interval"), alpha = 0.30) +
    geom_ribbon(aes(ymin = q250, ymax = q750, fill = "50% Credible Interval"), alpha = 0.45) +
    geom_line(data = ts_lines, aes(x = Date, y = val, color = method, linetype = method), linewidth = 0.9) +
    facet_wrap(~Coordinate, ncol = 1, scales = "free_y") +
    scale_fill_manual(
      name   = "Posterior Uncertainty",
      values = c("50% Credible Interval" = "#2b8cbe", "95% Credible Interval" = "#a6bddb")
    ) +
    scale_color_manual(
      name   = "Track Estimator",
      values = c("Posterior Mean" = "#0072B2", "Posterior Viterbi" = "#D55E00")
    ) +
    scale_linetype_manual(
      name   = "Track Estimator",
      values = c("Posterior Mean" = "solid", "Posterior Viterbi" = "twodash")
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text = element_text(size = 11, color = "black"),
      strip.text = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 11, face = "italic")
    ) +
    labs(
      title    = paste("Daily Posterior Credible Intervals (50% & 95%) - Cobia", instrument_id),
      subtitle = paste0("Model: ", model_name, " (glorys + light)"),
      x = "Date", y = "Coordinate (Degrees)"
    )

  ts_plot_file <- file.path(out_dir, paste0(instrument_id, "_", model_name, "_uncertainty_timeseries.png"))
  ggsave(ts_plot_file, plot = p_ts, width = 10, height = 7, dpi = 300)
  message("Saved Posterior Uncertainty time-series plot: ", ts_plot_file)

  data.frame(
    tag_id                  = instrument_id,
    model                   = model_name,
    n_days                  = T_steps,
    mean_end_error_km       = round(err_mean_km, 1),
    viterbi_end_error_km    = round(err_vit_km, 1),
    median_sd_lon_deg       = round(median(unc_df$sd_lon_deg, na.rm = TRUE), 3),
    median_sd_lat_deg       = round(median(unc_df$sd_lat_deg, na.rm = TRUE), 3),
    median_hpd50_area_km2   = round(median(unc_df$hpd50_area_km2, na.rm = TRUE), 1),
    median_hpd95_area_km2   = round(median(unc_df$hpd95_area_km2, na.rm = TRUE), 1),
    stringsAsFactors        = FALSE
  )
}

# ==============================================================================
# Run Across All Cobia Tags
# ==============================================================================
summary_list <- list()
for (tid in tag_ids) {
  res <- tryCatch(
    run_viterbi_and_uncertainty(tid, model_name = "glorys_light", percentile_v = 95, approach_v = 'joint'),
    error = function(e) {
      message("ERROR on tag ", tid, ": ", e$message)
      NULL
    }
  )
  if (!is.null(res)) summary_list[[tid]] <- res
}

if (length(summary_list) > 0) {
  viterbi_summary <- do.call(rbind, summary_list)
  summary_file <- file.path(rerun_dir, "all_tags_glorys_light_mean_vs_viterbi_summary.csv")
  write.csv(viterbi_summary, summary_file, row.names = FALSE)
  message("\n============================================================")
  message("Saved Mean vs. Viterbi & Uncertainty summary to: ", summary_file)
  message("============================================================")
  print(viterbi_summary)
}
