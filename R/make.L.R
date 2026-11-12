#' Combine individual source likelihoods
#' 
#' \code{make.L} combines individual likelihoods from various data sources (e.g.
#' SST, OHC) to make overall combined likelihoods for each time point
#' 
#' @param ras.list a list of likelihood rasters
#' @param iniloc is data.frame of tag and pop locations. Required columns are 'date' (POSIX), 'lon', 'lat'.
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param maxDepth a vector of the daily max depths (must be positive and be zero-padded for NA days)
#' @param bathy is the original bathy raster pre-resampling
#' @param known.locs is data frame of known locations containing named columns 
#'   of 'date' (POSIX), 'lon', 'lat'. Default is NULL.
#' @return an overall likelihood array, L
#' 
#' @examples
#' \dontrun{
#' L <- make.L(ras.list, iniloc, dateVec)
#' }
#' @export
#'   

make.L <- function(ras.list, iniloc, dateVec, maxDepth = NULL, bathy = NULL, known.locs = NULL){
  
  r_tmpl <- if (inherits(ras.list[[1]], "SpatRaster")) ras.list[[1]] else terra::rast(ras.list[[1]])
  nr <- terra::nrow(r_tmpl)
  nc <- terra::ncol(r_tmpl)
  T_len <- length(dateVec)
  
  # Convert all input rasters to in-memory 3D arrays (nrow x ncol x T) once
  arr_list <- lapply(ras.list, function(r) {
    sr <- if (inherits(r, "SpatRaster")) r else terra::rast(r)
    terra::as.array(sr)
  })
  n_ras <- length(arr_list)
  
  # Prepare 2D bathymetry matrix once (skip resample if geometries already match)
  bathy_mat <- NULL
  if (!is.null(bathy)) {
    if (!is.null(maxDepth) && any(maxDepth > 0) && length(maxDepth) != T_len){
      stop('Error: Vector of daily maximum depths is not same length at date vector. Ensure mmd was zero-padded for days with NA data.')
    }
    bathy_sr <- if (inherits(bathy, "SpatRaster")) bathy else terra::rast(bathy)
    b_min <- suppressWarnings(min(terra::values(bathy_sr), na.rm = TRUE))
    if (is.finite(b_min) && b_min < 0) {
      bathy_sr <- bathy_sr * -1
    }
    bathy_sr[bathy_sr <= 0] <- NA
    if (!terra::compareGeom(bathy_sr, r_tmpl, stopOnError = FALSE)) {
      bathy_sr <- terra::resample(bathy_sr, r_tmpl)
    }
    bathy_mat <- terra::as.matrix(bathy_sr, wide = TRUE)
  }
  
  L_arr <- array(0, dim = c(nr, nc, T_len))
  
  ## COMBINE THE LIKELIHOOD LAYERS IN MEMORY
  for (i in seq_len(T_len)){
    active_mats <- list()
    for (bb in seq_len(n_ras)){
      mat <- arr_list[[bb]][,,i]
      mv <- suppressWarnings(max(mat, na.rm = TRUE))
      if (is.finite(mv) && mv > 0) {
        active_mats[[length(active_mats) + 1L]] <- mat / mv
      }
    }
    
    if (length(active_mats) == 0L) next
    
    if (length(active_mats) == 1L){
      day_mat <- active_mats[[1L]]
    } else {
      s_prod <- Reduce(`*`, active_mats)
      m_prod <- suppressWarnings(max(s_prod, na.rm = TRUE))
      if (is.finite(m_prod) && m_prod > 0) {
        day_mat <- s_prod / m_prod
      } else {
        s_sum <- Reduce(function(a, b) {
          a[is.na(a)] <- 0
          b[is.na(b)] <- 0
          a + b
        }, active_mats)
        m_sum <- suppressWarnings(max(s_sum, na.rm = TRUE))
        day_mat <- if (is.finite(m_sum) && m_sum > 0) s_sum / m_sum else matrix(0, nr, nc)
      }
    }
    
    ## daily bathy mask (always zero out land; filter by depth when maxDepth[i] > 0)
    if (!is.null(bathy_mat)) {
      day_mat[is.na(bathy_mat)] <- 0
      if (!is.null(maxDepth) && maxDepth[i] > 0) {
        bathy.i <- maxDepth[i] * (1 - 20 / 100)
        day_mat[bathy_mat < bathy.i] <- 0
      }
    }
    
    L_arr[,,i] <- day_mat
  }
  
  #----------------------------------------------------------------------------------#
  # MAKE ALL NA'S VERY TINY FOR THE CONVOLUTION ON INTERIOR DAYS
  # Run BEFORE iniloc / known.locs so fixed start/end/known days stay exact point masses (1 vs 0)
  #----------------------------------------------------------------------------------#
  L_arr[is.na(L_arr) | L_arr <= 1e-15] <- 1e-15

  ## ADD START/END LOCATIONS
  if(!is.null(iniloc)){
    print('Adding start and end locations from iniloc...')
    iniloc$date <- as.POSIXct(
      sprintf('%04d-%02d-%02d', as.numeric(iniloc$year), as.numeric(iniloc$month), as.numeric(iniloc$day)),
      tz = 'UTC'
    )
    if (class(iniloc$date)[1] != class(dateVec)[1]) stop('dateVec and known.locs$date both need to be of class POSIXct.')
    iniloc$dateVec <- findInterval(iniloc$date, dateVec)
    
    for(i in unique(iniloc$dateVec)){
      known.locs.i <- iniloc[which(iniloc$dateVec == i),]
      L_arr[,,i] <- 0
      cells <- terra::cellFromXY(r_tmpl, as.matrix(known.locs.i[, c('lon', 'lat')]))
      rc <- terra::rowColFromCell(r_tmpl, cells)
      valid_rc <- !is.na(rc[,1]) & !is.na(rc[,2])
      if (any(valid_rc)) {
        L_arr[cbind(rc[valid_rc, 1], rc[valid_rc, 2], i)] <- 1
      }
    }
  }
  
  ## ADD KNOWN LOCATIONS
  if(!is.null(known.locs)){
    print('Known locations are being added to the likelihoods...')
    if (class(known.locs$date)[1] != class(dateVec)[1]) stop('dateVec and known.locs$date both need to be of class POSIXct.')
    known.locs <- known.locs[which(known.locs$date <= max(dateVec)),]
    known.locs$dateVec <- findInterval(known.locs$date, dateVec)
    
    for(i in unique(known.locs$dateVec)){
      known.locs.i <- known.locs[which(known.locs$dateVec == i),]
      if(nrow(known.locs.i) > 1){
        warning(paste0('Multiple locations supplied at time step ', dateVec[i], '. Only the first one is being used.'))
        known.locs.i <- known.locs.i[1,]
      }
      L_arr[,,i] <- 0
      cells <- terra::cellFromXY(r_tmpl, as.matrix(known.locs.i[, c('lon', 'lat')]))
      rc <- terra::rowColFromCell(r_tmpl, cells)
      valid_rc <- !is.na(rc[,1]) & !is.na(rc[,2])
      if (any(valid_rc)) {
        L_arr[cbind(rc[valid_rc, 1], rc[valid_rc, 2], i)] <- 1
      }
    }
  }
  
  # CONVERT NORTH-UP (nrow, ncol, T) ARRAY TO SOUTH-UP (T, ncol, nrow) ARRAY
  L <- aperm(L_arr[nr:1, , , drop = FALSE], c(3, 2, 1))
  
  print('Finishing make.L...', sep='')
  
  return(L)
}
