#' Calculates error ellipses around input locations
#' 
#' Error ellipses are generated around WC tag based positions using
#' tag-calculated error structure.
#' 
#' @param locs is data frame of light-based location estimates containing columns for Date (POSIXct), Longitude, Latitude, Error.Semi.minor.axis, Error.Semi.major.axis, Offset, and Offset.orientation.
#' @param locs.grid is list output from \code{setup.locs.grid}
#' 
#' @return an array of error ellipses for input light-based locations  

calc.errEll <- function(locs, locs.grid){

  # arithmetic converts from meters to degrees
  # add transformation due to projection?
  
  # set up the offset
  locs$Offset[which(is.na(locs$Offset))] <- 0
  locs$Offset.orientation[which(is.na(locs$Offset.orientation))] <- 0
  
  #then do shift based on:
  if(locs$Offset.orientation == 0){
    
    shiftDist <- locs$Offset / 1000 / 111
    
  } else if(locs$Offset.orientation == 180){
    
    shiftDist <- (-1 * (locs$Offset / 1000 / 111))
    
  }
  
  # if offset is really large, we revert to longitude only. if not, we use elliptical error.
  if(shiftDist >= -10 & shiftDist <= 10){
    # set up a larger grid to base ellipse on and to shift that error, if necessary (GPE only)
    ngrid <- rev(dim(locs.grid$lon))
    lon1 <- seq(min(locs.grid$lon[1,]) - 10, max(locs.grid$lon[1,]) + 10, by = locs.grid$dlo)
    lat1 <- seq(min(locs.grid$lat[,1]) - 10, max(locs.grid$lat[,1]) + 10, by = locs.grid$dla)
    g1 <- meshgrid(lon1, lat1)
    
    ## longitude error
    if ('Error.Semi.minor.axis' %in% names(locs)){
      # create longitude likelihood based on GPE data
      slon.sd <- locs$Error.Semi.minor.axis[1] / 1000 / 111 #semi minor axis
    } else if ('longitudeError' %in% names(locs)){
      # create longitude likelihood based on GPE data
      slon.sd <- locs$longitudeError #semi minor axis
    }
    
    ## latitude error
    if ('Error.Semi.major.axis' %in% names(locs)){
      # create lat likelihood based on GPE data
      slat.sd <- locs$Error.Semi.major.axis[1] / 1000 / 111 #semi minor axis
    } else if ('latitudeError' %in% names(locs)){
      # create lat likelihood based on GPE data
      slat.sd <- locs$latitudeError #semi minor axis
    }
    
    # calc semi minor axis based on longitude error
    #slon.sd <- locs$Error.Semi.minor.axis / 1000 / 111 #semi minor axis
    L.light.lon <- stats::dnorm(t(g1$X), locs$Longitude, slon.sd) # Longitude data
    #slat.sd <- locs$Error.Semi.major.axis / 1000 / 111 #semi major axis
    L.light.lat <- stats::dnorm(t(g1$Y), locs$Latitude, slat.sd)
    
    #image.plot(g$lon[1,],g$lat[,1],L.light.lat*L.light.lon)
    
    # create the ellipse by multiplying lat * lon error
    L <- raster::flip(raster::raster(t(L.light.lat * L.light.lon), xmn = min(lon1), 
                                     xmx = max(lon1), ymn = min(lat1), ymx = max(lat1)), direction = 'y')
    
    # then shift the ellipse based on shiftDist calculated above
    Ls <- raster::shift(L, y = shiftDist)
    
    # create an extent raster to resample our current raster back to the uniform extent/resolution
    L.ext <- raster::flip(raster::raster(locs.grid$lon, xmn = min(locs.grid$lon[1,]), 
                                         xmx = max(locs.grid$lon[1,]), ymn = min(locs.grid$lat[,1]),
                                         ymx = max(locs.grid$lat[,1])), direction = 'y')
    # fill with 1s
    L.ext[L.ext <= 0] = 1
    
    # then crop our shifted raster
    Lsx <- raster::crop(Ls, L.ext)
    rr <- raster::resample(Lsx, L.ext)
    L.light <- t(raster::as.matrix(raster::flip(rr, direction = 'y')))
    
  } else if(shiftDist < -10 | shiftDist > 10){
    
    # if supposed shift in error ellipse is >10 degrees, we revert to longitude only
    if ('Error.Semi.minor.axis' %in% names(locs)){
      # create longitude likelihood based on GPE data
      slon.sd <- locs$Error.Semi.minor.axis[1] / 1000 / 111 #semi minor axis
    } else if ('longitudeError' %in% names(locs)){
      # create longitude likelihood based on GPE data
      slon.sd <- locs$longitudeError #semi minor axis
    }
    
    L.light <- stats::dnorm(t(locs.grid$lon), locs$Longitude, slon.sd)
    
  }
  
  return(L.light)
  
}
