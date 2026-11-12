#' Calculate Bathymetry-based likelihood
#' 
#' \code{calc.bathy} compares tag max depth to bathymetry grid and calculates 
#' likelihoods
#' 
#' @param mmd is dataframe containing at least cols: "Date" (POSIXct) and "MaxDepth" in meters (positive values are wet)
#' @param bathy.grid is raster of bathymetry. If the minimum bathymetry values are < 0, the function automatically converts in-water values to positive and masks land.
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param focalDim is integer for dimensions of raster::focal used to calculate 
#'   sd() of env grid cell. If left to NULL (default), this dimension will
#'   approximately incorporate 0.25 degrees.
#' @param sens.err is numeric indicating the percent sensor error in the tag-measured max depth. This allows some uncertainty when calculating the integral for the likelihood and doesnt have to necessarily reflect the actual sensor error. Default is 5 percent.
#' @param lik.type is character indicating which likelihood type to use in the bathymetry calculation. Options are dnorm (a traditional likelihood bounded by tag measurement +/- sens.err; experimental) or max ("one-sided" likelihood >= tag-measured max depth; DEFAULT). The latter choice acts more like a mask in that it doesnt allow likelihoods in water shallower than the max depth for each time point.
#' 
#' @return likelihood is raster brick of likelihood surfaces representing
#' matches between tag-based sst and remotely sensed sst maps
#'   
#' @export
#' 
#' @author Originally written by Paul Gatti

calc.bathy <- function(mmd, bathy.grid, dateVec, focalDim = NULL, sens.err = 5, lik.type = 'max'){
  
  if (lik.type == 'dnorm') warning('Bathymetry likelihood calculation with lik.type = normal is experimental. If you use it, please send feedback on your experience as we work to improve it.')
  if (!inherits(bathy.grid, "SpatRaster")) bathy.grid <- terra::rast(bathy.grid)
  
  ## convert a negative bathy grid to positive to match expectations and mask land
  b_min <- suppressWarnings(min(terra::values(bathy.grid), na.rm = TRUE))
  if (is.finite(b_min) && b_min < 0){
    bathy.grid <- bathy.grid * -1
    bathy.grid[bathy.grid < 0] <- NA
  }
  
  if (class(mmd$Date)[1] != class(dateVec)[1]) stop('dateVec and mmd$Date both need to be of class POSIXct.')
  mmd <- mmd[which(mmd$Date <= max(dateVec)),]
  mmd$dateVec <- findInterval(mmd$Date, dateVec)
  mmd <- data.frame(mmd %>% group_by(dateVec) %>% 
                          dplyr::summarise(maxDepth = max(MaxDepth, na.rm=TRUE), .groups = 'drop_last'))
  
  print(paste("Starting bathymetry likelihood calculation..."))
  t0 <- Sys.time()
  
  if (is.null(focalDim)){
    focalDim <- round(0.25 / terra::res(bathy.grid)[1], 0)
    if (focalDim %% 2 == 0) focalDim <- focalDim - 1
    if (focalDim < 3) focalDim <- 3
  }
  
  dat <- terra::as.matrix(bathy.grid, wide = TRUE)
  if (lik.type == 'dnorm') {
    sd_r <- terra::focal(bathy.grid, w = focalDim, fun = "sd", na.rm = TRUE)
    sdx  <- terra::as.matrix(sd_r, wide = TRUE)
  }
  
  T <- length(dateVec)
  print(paste("Starting iterations through deployment period ", "..."))
  
  L_arr <- array(0, dim = c(terra::nrow(bathy.grid), terra::ncol(bathy.grid), T))
  
  for (i in 1:T) {
    print(dateVec[i])
    mmd.i <- mmd[which(mmd$dateVec == i),]
    if (nrow(mmd.i) == 0) next
    
    bathy.i <- c(mmd.i$maxDepth * (1 - sens.err / 100), mmd.i$maxDepth * (1 + sens.err / 100))
    
    if (lik.type == 'max'){
      lik.bathy <- matrix(0, nrow = nrow(dat), ncol = ncol(dat))
      lik.bathy[!is.na(dat) & dat >= bathy.i[1]] <- 1
    } else if (lik.type == 'dnorm'){
      lik.bathy <- likint3(dat, sdx, bathy.i[1], bathy.i[2])
      m <- suppressWarnings(max(lik.bathy, na.rm = TRUE))
      if (is.finite(m) && m > 0) lik.bathy <- lik.bathy / m
    } else{
      stop('Argument lik.type can only be one of max or dnorm.')
    }
    
    lik.bathy[is.na(lik.bathy) | is.infinite(lik.bathy) | lik.bathy < 0] <- 0
    L_arr[,,i] <- lik.bathy
  }
  
  print(paste("Making final likelihood raster..."))
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  L.bathy <- terra::rast(L_arr, extent = terra::ext(bathy.grid), crs = crs)
  names(L.bathy) <- as.character(dateVec)
  
  t1 <- Sys.time()
  print(paste("Bathymetric calculations took ", round(as.numeric(difftime(t1, t0, units = "mins")), 2), "minutes..."))
  return(L.bathy)
}