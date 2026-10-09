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
  
  ## convert a negative bathy grid to positive to match expectations and mask land
  if (raster::cellStats(bathy.grid, 'min', na.rm=T) < 0){
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
  #if (class(mmd$Date) != 'Date') stop('mmd$Date must be of class Date')
  
  # compute bathy.grid sd
  sdx = raster::focal(bathy.grid, w = matrix(1, nrow = focalDim, ncol = focalDim), fun = function(x) stats::sd(x,na.rm = TRUE), pad = TRUE)
  ## sdx to matrix for likint3
  sdx = t(raster::as.matrix(raster::flip(sdx, 2)))
  ## bathy grid to matrix for likint3
  dat <- base::t(raster::as.matrix(raster::flip(bathy.grid, 2)))
  
  T <- length(dateVec)
  print(paste("Starting iterations through deployment period ", "..."))
  
  L.bathy <- array(0, dim = c(dim(bathy.grid)[2:1], length(dateVec)))
  lon.agg <- seq(raster::extent(bathy.grid)[1], raster::extent(bathy.grid)[2], 
                 length.out = dim(bathy.grid)[2])
  lat.agg <- seq(raster::extent(bathy.grid)[3], raster::extent(bathy.grid)[4], 
                 length.out = dim(bathy.grid)[1])
  
  for (i in 1:T) {

    print(dateVec[i])
    
    mmd.i <- mmd[which(mmd$dateVec == i),]
    if (nrow(mmd.i) == 0) next
    
    bathy.i <- c(mmd.i$maxDepth * (1 - sens.err / 100), mmd.i$maxDepth * (1 + sens.err / 100)) # sensor error
    
    if (lik.type == 'max'){
      ## create bathy max from max depth of tag to max of bathy grid allowed (=1), otherwise 0
      lik.bathy <- bathy.grid
      lik.bathy[lik.bathy < bathy.i[1]] <- 0
      lik.bathy[lik.bathy >= bathy.i[1]] <- 1
      lik.bathy <- base::t(raster::as.matrix(raster::flip(lik.bathy, 2)))
      
    } else if (lik.type == 'dnorm'){
      
      ## the actual likelihood calculation
      lik.bathy <- likint3(dat, sdx, bathy.i[1], bathy.i[2])
      lik.bathy = as.matrix(lik.bathy) / max(as.matrix(lik.bathy), na.rm = TRUE)
    } else{
      stop('Argument lik.type can only be one of max or dnorm.')
    }
    
    lik.bathy[is.na(lik.bathy) | is.infinite(lik.bathy)] <- 0
    L.bathy[,,i] <- lik.bathy
  }
  
  L.bathy <- aperm(L.bathy,c(2,1,3))
  L.bathy <- L.bathy[,dim(L.bathy)[2]:1,]
  
  print(paste("Making final likelihood raster..."))
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.bathy <- list(x = lon.agg, y = lat.agg, z = L.bathy)

  ex <- raster::extent(list.bathy)
  L.bathy <- raster::brick(list.bathy$z, xmn = ex[1], xmx = ex[2], 
                           ymn = ex[3], ymx = ex[4], transpose = F, crs)
  L.bathy <- raster::flip(raster::flip(L.bathy, direction = "y"), direction = 'x')
  L.bathy[L.bathy < 0] <- 0
  
  t1 <- Sys.time()
  print(paste("Bathymetric calculations took ", round(as.numeric(difftime(t1, t0, units = "mins")), 2), "minutes..."))
  return(L.bathy)
}