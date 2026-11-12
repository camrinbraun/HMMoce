library(terra)

plot_L_series <- function(L_raster, dateVec, days_to_plot = c(1, 10, 20, 30, 40, 50, 60, 70, 93), main_title = "Likelihood") {
  # Keep only indices present in dateVec and within raster bounds

dateVec <- dateVec[days_to_plot]

  # max_len <- if (inherits(L_raster, "list")) length(L_raster) else raster::nlayers(L_raster)

  
  if (inherits(L_raster, "list")) {
    r_terra <- terra::rast(lapply(days_to_plot, function(i) terra::rast(L_raster[[i]])))
  } else {
    r_terra <- terra::rast(L_raster)[[days_to_plot]]
  }
  
  # Names strictly match the subsetted days_to_plot indices
  names(r_terra) <- format(as.Date(dateVec[days_to_plot]), "%Y-%m-%d")
  
  terra::plot(r_terra, 
              col = viridis::viridis(100), 
              main = paste(main_title, names(r_terra)),
              axes = TRUE)
}

array_to_rast <- function(L, g, dateVec, days_to_plot = c(1, 10, 20, 30, 40, 50, 60, 70, 93), main_title = "Combined L") {
  # valid_idx <- days_to_plot[days_to_plot >= 1 & days_to_plot <= min(length(dateVec), dim(L)[1])]
  
  sub_list <- lapply(days_to_plot, function(i) {
    mat <- t(L[i, , ])
    r <- terra::rast(mat, 
                     extent = terra::ext(min(g$lon), max(g$lon), min(g$lat), max(g$lat)),
                     crs = "EPSG:4326")
    terra::flip(r, direction = "vertical")
  })
  
  r_stack <- terra::rast(sub_list)
  names(r_stack) <- format(as.Date(dateVec[days_to_plot]), "%Y-%m-%d")
  
  terra::plot(r_stack, 
              col = viridis::viridis(100), 
              main = paste(main_title, names(r_stack)),
              axes = TRUE)
}

plot_L_series(L.res$L.rasters$ohc, dateVec, days_to_plot = c(1,3,25,38), main_title = "OHC Likelihood")

array_to_rast(L.res$L.rasters$ohc, dateVec, days_to_plot = c(1,3,25,38), main_title = "OHC Likelihood")


# plot the raster time series for inspection ----

library(terra)

# Convert 3D [time, lon, lat] (e.g. L) or 4D [state, time, lon, lat] (f$phi or s) to SpatRaster
hmm_array_to_rast <- function(arr, g, dateVec) {
  if (length(dim(arr)) == 4) {
    # Sum posterior across behavioral states -> [time, lon, lat]
    arr <- apply(arr, c(2, 3, 4), sum, na.rm = TRUE)
  }
  lon_v <- g$lon[1, ]; lat_v <- g$lat[, 1]
  dx <- abs( lon_v[2] - lon_v[1] ) / 2
  dy <- abs( lat_v[2] - lat_v[1] ) / 2
  ex <- ext(min(lon_v) - dx, max(lon_v) + dx, min(lat_v) - dy, max(lat_v) + dy)

  r <- flip(rast(aperm(arr, c(3, 2, 1)), extent = ex, crs = "EPSG:4326"), direction = "vertical")
  names(r) <- format(as.Date(dateVec), "%Y-%m-%d")
  r
}

# Plot selected days with Release (green) and Popoff (red) overlaid
plot_hmm_days <- function(x, dateVec, iniloc, days = c(1, 2, 5, 15, 30, 60, 90, 92, 93), nr = 3, nc = 3) {
  r <- if (inherits(x, "SpatRaster")) x else rast(x)
  names(r) <- format(as.Date(dateVec), "%Y-%m-%d")
  par(mfrow = c(nr, nc), oma = c(0, 0, 0, 0))
  invisible(lapply(days, function(i) {
    plot(r[[i]], main = paste0("Day ", i, ": ", names(r)[i]), mar = c(1.2, 1.2, 1.4, 2.2))
    points(iniloc$lon[1], iniloc$lat[1], pch = 21, bg = "green", cex = 1.3)
    points(iniloc$lon[2], iniloc$lat[2], pch = 21, bg = "red",   cex = 1.3)
  }))
}

# 1) Inspect individual likelihoods (e.g. glorys vs ohc vs light)
plot_hmm_days(L.res$L.rasters$glorys, dateVec, iniloc, days = c(1,5,10,20,30))

plot_hmm_days(L.res$L.rasters$light,  dateVec, iniloc, days = c(1,5,10,20,30))

plot_hmm_days(L.light,  dateVec, iniloc, days = c(1,5,10,20,30))


# 2) Inspect combined L, forward filter (f$phi), and smoothed posterior (s)
L_rast    <- hmm_array_to_rast(L,     L.res$g, dateVec)
phi_rast  <- hmm_array_to_rast(f$phi, L.res$g, dateVec)
post_rast <- hmm_array_to_rast(s,     L.res$g, dateVec)

plot_hmm_days(L_rast,    dateVec, iniloc)
plot_hmm_days(phi_rast,  dateVec, iniloc)
plot_hmm_days(post_rast, dateVec, iniloc)
