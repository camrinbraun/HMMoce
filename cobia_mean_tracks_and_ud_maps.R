# ==============================================================================
# Cobia PSAT Final Posterior Mean Tracks & Posterior UD Maps (glorys_light only)
# Produces:
#   1. Per-tag Posterior UD + 50%/95% HPD map (Posterior Mean track only)
#   2. Per-tag Posterior Mean track map with monthly color scheme
#   3. Combined All-Tags Posterior Mean track map (monthly color scheme, no UD)
#   4. Combined daily mean tracks + uncertainty CSV & per-tag summary CSV
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
model_name     <- "glorys_light"

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
    geospatial_lat_start = as.numeric(tag_track$Lat_N[1]),
    geospatial_lon_start = as.numeric(tag_track$Lon_E[1]),
    geospatial_lat_end   = as.numeric(tag_track$Lat_N[nrow(tag_track)]),
    geospatial_lon_end   = as.numeric(tag_track$Lon_E[nrow(tag_track)]),
    stringsAsFactors     = FALSE
  )
}

weighted_quantile <- function(x, w, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  w[is.na(w) | w < 0] <- 0
  s <- sum(w)
  if (s <= 0) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  cdf <- cumsum(w[ord]) / s
  vapply(probs, function(p) x[which(cdf >= p)[1]], numeric(1))
}

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

psat_reports_dir <- file.path(base_dir, "Updated PSAT Reports")
if (dir.exists(psat_reports_dir)) {
  tag_dirs <- list.dirs(psat_reports_dir, full.names = FALSE, recursive = FALSE)
  tag_ids  <- sort(sub("_updated.*$", "", tag_dirs[grepl("_updated", tag_dirs)]))
} else {
  tag_ids <- c("36419", "36420", "36421", "36422", "40302", "49157", "55550", "55552")
}

world_sf <- sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))

# ==============================================================================
# Per-Tag Mean Track & Posterior UD Map Function
# ==============================================================================
process_tag_mean_maps <- function(instrument_id, model_name = "glorys_light") {
  out_dir  <- file.path(rerun_dir, instrument_id)
  fit_file <- file.path(out_dir, paste0(instrument_id, "_", model_name, "_fixed_fit.rds"))
  if (!file.exists(fit_file)) {
    message("Skipping tag ", instrument_id, ": fit file not found at ", fit_file)
    return(NULL)
  }

  message("\n============================================================")
  message("Tag ", instrument_id, " | Generating Mean Track & UD Maps (", model_name, ")")
  message("============================================================")

  fit_obj <- readRDS(fit_file)
  s       <- fit_obj$s
  tr_mean <- fit_obj$tr
  tr_mean$Date <- as.POSIXct(tr_mean$Date, tz = "UTC")
  tr_mean$Month <- factor(format(tr_mean$Date, "%b %Y"), levels = unique(format(tr_mean$Date, "%b %Y")))

  tag_meta <- get_cobia_meta(instrument_id, base_dir)
  t_start  <- tag_meta$time_coverage_start
  t_end    <- tag_meta$time_coverage_end

  iniloc <- data.frame(
    tag_id = instrument_id,
    day    = c(lubridate::day(t_start), lubridate::day(t_end)),
    month  = c(lubridate::month(t_start), lubridate::month(t_end)),
    year   = c(lubridate::year(t_start), lubridate::year(t_end)),
    lat    = c(tag_meta$geospatial_lat_start, tag_meta$geospatial_lat_end),
    lon    = c(tag_meta$geospatial_lon_start, tag_meta$geospatial_lon_end),
    type   = c("Release", "Popoff"),
    stringsAsFactors = FALSE
  )
  rownames(iniloc) <- c("start", "end")

  tag_dt  <- as.POSIXct(sprintf("%04d-%02d-%02d", iniloc["start", "year"], iniloc["start", "month"], iniloc["start", "day"]), tz = "UTC")
  pop_dt  <- as.POSIXct(sprintf("%04d-%02d-%02d", iniloc["end", "year"], iniloc["end", "month"], iniloc["end", "day"]), tz = "UTC")
  dateVec <- seq.POSIXt(tag_dt, pop_dt, by = "24 hours")

  # Reconstruct spatial grid g matching s dimensions (m x T x ncol x nrow)
  glorys_rds <- file.path(out_dir, "L.glorys.rds")
  if (!file.exists(glorys_rds)) stop("Template raster L.glorys.rds not found for tag ", instrument_id)
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

  lon_vec <- g$lon[1, ]
  lat_vec <- g$lat[, 1]
  sp.lim  <- list(lonmin = min(lon_vec), lonmax = max(lon_vec), latmin = min(lat_vec), latmax = max(lat_vec))

  # Daily & Cumulative Posterior Uncertainty
  T_steps <- dim(s)[2]
  nx      <- length(lon_vec)
  ny      <- length(lat_vec)
  dx_deg  <- abs(mean(diff(lon_vec)))
  dy_deg  <- abs(mean(diff(lat_vec)))

  p_resid_vec <- if (dim(s)[1] > 1) rowSums(s, dims = 2)[2, ] else rep(NA_real_, T_steps)
  post_daily  <- apply(s, c(2, 3, 4), sum, na.rm = TRUE)
  ud_sum      <- matrix(0, nrow = nx, ncol = ny)
  unc_rows    <- vector("list", T_steps)

  for (t in seq_len(T_steps)) {
    P_t <- post_daily[t, , ]
    P_t[is.na(P_t) | P_t < 0] <- 0
    tot_t <- sum(P_t)
    if (tot_t > 0) P_t <- P_t / tot_t
    ud_sum <- ud_sum + P_t

    p_lon  <- rowSums(P_t)
    p_lat  <- colSums(P_t)
    mu_lon <- sum(lon_vec * p_lon)
    mu_lat <- sum(lat_vec * p_lat)
    sd_lon <- sqrt(max(0, sum(((lon_vec - mu_lon)^2) * p_lon)))
    sd_lat <- sqrt(max(0, sum(((lat_vec - mu_lat)^2) * p_lat)))
    q_lon  <- weighted_quantile(lon_vec, p_lon)
    q_lat  <- weighted_quantile(lat_vec, p_lat)

    hpd_t <- calc_hpd_volume_matrix(P_t)
    cell_area_km2 <- (dx_deg * 111.32 * cos(mu_lat * pi / 180)) * (dy_deg * 110.57)
    n_cells_50 <- sum(!is.na(hpd_t) & hpd_t <= 0.50)
    n_cells_95 <- sum(!is.na(hpd_t) & hpd_t <= 0.95)

    unc_rows[[t]] <- data.frame(
      tag_id         = instrument_id,
      Date           = tr_mean$Date[t],
      lon            = tr_mean$lon[t],
      lat            = tr_mean$lat[t],
      sd_lon_deg     = sd_lon,
      sd_lat_deg     = sd_lat,
      lon_q025       = q_lon[1],
      lon_q250       = q_lon[2],
      lon_q500       = q_lon[3],
      lon_q750       = q_lon[4],
      lon_q975       = q_lon[5],
      lat_q025       = q_lat[1],
      lat_q250       = q_lat[2],
      lat_q500       = q_lat[3],
      lat_q750       = q_lat[4],
      lat_q975       = q_lat[5],
      hpd50_area_km2 = round(n_cells_50 * cell_area_km2, 1),
      hpd95_area_km2 = round(n_cells_95 * cell_area_km2, 1),
      p_resid        = p_resid_vec[t],
      stringsAsFactors = FALSE
    )
  }
  track_unc_df <- do.call(rbind, unc_rows)

  # Cumulative Posterior UD & 50%/95% HPD Contours
  ud_norm <- ud_sum / sum(ud_sum, na.rm = TRUE)
  hpd_vol <- calc_hpd_volume_matrix(ud_norm)
  ext_g   <- terra::ext(min(lon_vec) - dx_deg / 2, max(lon_vec) + dx_deg / 2,
                        min(lat_vec) - dy_deg / 2, max(lat_vec) + dy_deg / 2)
  ud_rast  <- terra::rast(t(ud_norm / max(ud_norm, na.rm = TRUE))[ny:1, , drop = FALSE], extent = ext_g, crs = "EPSG:4326")
  hpd_rast <- terra::rast(t(hpd_vol)[ny:1, , drop = FALSE], extent = ext_g, crs = "EPSG:4326")

  hpd_contours_sf <- tryCatch({
    cl <- terra::as.contour(hpd_rast, levels = c(0.50, 0.95))
    sf_obj <- sf::st_as_sf(cl)
    sf_obj$Isopleth <- factor(paste0(round(sf_obj$level * 100), "% HPD"), levels = c("50% HPD", "95% HPD"))
    sf_obj
  }, error = function(e) NULL)

  ud_df <- as.data.frame(ud_rast, xy = TRUE)
  names(ud_df) <- c("lon", "lat", "rel_prob")
  hpd_df <- as.data.frame(hpd_rast, xy = TRUE)
  names(hpd_df) <- c("lon", "lat", "hpd_vol")
  ud_df$hpd_vol <- hpd_df$hpd_vol
  ud_df <- ud_df[!is.na(ud_df$hpd_vol) & ud_df$hpd_vol <= 0.99, ]

  err_mean_km <- terra::distance(
    matrix(c(tr_mean$lon[nrow(tr_mean)], tr_mean$lat[nrow(tr_mean)]), ncol = 2),
    matrix(c(iniloc$lon[2], iniloc$lat[2]), ncol = 2), lonlat = TRUE
  )[1] / 1000

  # 1. Per-Tag Posterior UD Heatmap (Posterior Mean Track Only) ----------------
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
    geom_path(data = tr_mean, aes(x = lon, y = lat), color = "black", linewidth = 1.6) +
    geom_path(
      data = tr_mean,
      aes(x = lon, y = lat, color = "Posterior Mean"),
      linewidth = 0.95
    ) +
    scale_color_manual(name = "Track", values = c("Posterior Mean" = "white")) +
    geom_point(data = iniloc[1, ], aes(x = lon, y = lat), fill = "green", color = "black", shape = 24, size = 4) +
    geom_point(data = iniloc[2, ], aes(x = lon, y = lat), fill = "red", color = "black", shape = 25, size = 4) +
    coord_sf(xlim = c(sp.lim$lonmin, sp.lim$lonmax), ylim = c(sp.lim$latmin, sp.lim$latmax), expand = FALSE) +
    theme_bw(base_size = 13) +
    theme(
      axis.text = element_text(size = 11, color = "black"),
      plot.subtitle = element_text(size = 11, face = "italic"),
      legend.key = element_rect(fill = "grey25", color = NA)
    ) +
    labs(
      title    = paste("Cumulative Posterior Utilization Distribution (UD) - Cobia", instrument_id),
      subtitle = sprintf("Model: %s (50%% Core & 95%% Credible Contours | Endpoint Error: %.1f km)", model_name, err_mean_km),
      x = "Longitude", y = "Latitude"
    )

  ud_plot_file <- file.path(out_dir, paste0(instrument_id, "_", model_name, "_mean_posterior_UD_map.png"))
  ggsave(ud_plot_file, plot = p_ud, width = 9.5, height = 7.5, dpi = 300)
  message("Saved Mean Posterior UD map: ", ud_plot_file)

  # 2. Per-Tag Monthly Color Scheme Mean Track Map -----------------------------
  p_track <- ggplot() +
    geom_sf(data = world_sf, fill = "grey80", color = "grey50", linewidth = 0.3) +
    geom_path(data = tr_mean, aes(x = lon, y = lat), color = "grey35", linewidth = 0.8) +
    geom_point(
      data = tr_mean,
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
      title    = paste("HMMoce Posterior Mean Track - Cobia", instrument_id),
      subtitle = paste0("Model: ", model_name, " (glorys + light)"),
      x = "Longitude", y = "Latitude"
    )

  track_plot_file <- file.path(out_dir, paste0(instrument_id, "_", model_name, "_mean_track_map.png"))
  ggsave(track_plot_file, plot = p_track, width = 9, height = 7, dpi = 300)
  message("Saved Monthly Mean Track map: ", track_plot_file)

  summary_row <- data.frame(
    tag_id                = instrument_id,
    model                 = model_name,
    n_days                = T_steps,
    mean_end_error_km     = round(err_mean_km, 1),
    median_sd_lon_deg     = round(median(track_unc_df$sd_lon_deg, na.rm = TRUE), 3),
    median_sd_lat_deg     = round(median(track_unc_df$sd_lat_deg, na.rm = TRUE), 3),
    median_hpd50_area_km2 = round(median(track_unc_df$hpd50_area_km2, na.rm = TRUE), 1),
    median_hpd95_area_km2 = round(median(track_unc_df$hpd95_area_km2, na.rm = TRUE), 1),
    stringsAsFactors      = FALSE
  )

  list(track_df = track_unc_df, iniloc = iniloc, summary = summary_row)
}

# ==============================================================================
# Run Across All Cobia Tags
# ==============================================================================
all_tracks_list  <- list()
all_iniloc_list  <- list()
all_summary_list <- list()

for (tid in tag_ids) {
  res <- tryCatch(
    process_tag_mean_maps(tid, model_name = model_name),
    error = function(e) {
      message("ERROR on tag ", tid, ": ", e$message)
      NULL
    }
  )
  if (!is.null(res)) {
    all_tracks_list[[tid]]  <- res$track_df
    all_iniloc_list[[tid]]  <- res$iniloc
    all_summary_list[[tid]] <- res$summary
  }
}

# ==============================================================================
# 3. Combined All-Tags Mean Track Map (Month Colors, No Posterior) & CSVs
# ==============================================================================
if (length(all_tracks_list) > 0) {
  all_tracks_df  <- do.call(rbind, all_tracks_list)
  all_iniloc_df  <- do.call(rbind, all_iniloc_list)
  all_summary_df <- do.call(rbind, all_summary_list)

  # Calendar Month factor (Jan-Dec) for unified seasonal color scale across years
  obs_months <- month.abb[month.abb %in% unique(format(all_tracks_df$Date, "%b"))]
  all_tracks_df$Month <- factor(format(all_tracks_df$Date, "%b"), levels = obs_months)

  lon_rng <- range(c(all_tracks_df$lon, all_iniloc_df$lon), na.rm = TRUE)
  lat_rng <- range(c(all_tracks_df$lat, all_iniloc_df$lat), na.rm = TRUE)
  pad_deg <- 1.5

  start_pts <- all_iniloc_df[all_iniloc_df$type == "Release", ]
  end_pts   <- all_iniloc_df[all_iniloc_df$type == "Popoff", ]

  p_all <- ggplot() +
    geom_sf(data = world_sf, fill = "grey80", color = "grey50", linewidth = 0.3) +
    geom_path(
      data = all_tracks_df,
      aes(x = lon, y = lat, group = tag_id),
      color = "grey35", linewidth = 0.7, alpha = 0.85
    ) +
    geom_point(
      data = all_tracks_df,
      aes(x = lon, y = lat, fill = Month),
      shape = 21, color = "black", stroke = 0.25, size = 2.5
    ) +
    scale_fill_viridis_d(option = "turbo", name = "Month", begin = 0.1, end = 0.9) +
    geom_point(data = start_pts, aes(x = lon, y = lat), fill = "green", color = "black", shape = 24, size = 3.6) +
    geom_point(data = end_pts,   aes(x = lon, y = lat), fill = "red",   color = "black", shape = 25, size = 3.6) +
    coord_sf(
      xlim = c(lon_rng[1] - pad_deg, lon_rng[2] + pad_deg),
      ylim = c(lat_rng[1] - pad_deg, lat_rng[2] + pad_deg),
      expand = FALSE
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text = element_text(size = 11, color = "black"),
      plot.subtitle = element_text(size = 11, face = "italic")
    ) +
    labs(
      title    = "HMMoce Final Posterior Mean Tracks - All Cobia Tags",
      subtitle = sprintf("Model: %s (n = %d tags | Green = Release, Red = Popoff)", model_name, length(all_tracks_list)),
      x = "Longitude", y = "Latitude"
    )

  all_map_file <- file.path(rerun_dir, paste0("all_tags_", model_name, "_mean_tracks_map.png"))
  ggsave(all_map_file, plot = p_all, width = 10, height = 8, dpi = 300)
  message("\nSaved combined all-tags Mean Track map: ", all_map_file)

  tracks_csv  <- file.path(rerun_dir, paste0("all_tags_", model_name, "_mean_tracks.csv"))
  summary_csv <- file.path(rerun_dir, paste0("all_tags_", model_name, "_mean_summary.csv"))
  write.csv(all_tracks_df, tracks_csv, row.names = FALSE)
  write.csv(all_summary_df, summary_csv, row.names = FALSE)
  message("Saved combined daily tracks CSV: ", tracks_csv)
  message("Saved summary CSV: ", summary_csv)
  print(all_summary_df)
}
