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
  
  lat_valid <- !is.na(locs$Latitude[1]) &&
               locs$Latitude[1] >= min(locs.grid$lat[,1]) &&
               locs$Latitude[1] <= max(locs.grid$lat[,1])

  # if offset is really large or latitude is out of bounds/NA, revert to longitude only. if not, use elliptical error.
  if(lat_valid && shiftDist >= -10 & shiftDist <= 10){
    ## longitude error
    if ('Error.Semi.minor.axis' %in% names(locs)){
      slon.sd <- locs$Error.Semi.minor.axis[1] / 1000 / 111
    } else if ('longitudeError' %in% names(locs)){
      slon.sd <- locs$longitudeError[1]
    }
    
    ## latitude error
    if ('Error.Semi.major.axis' %in% names(locs)){
      slat.sd <- locs$Error.Semi.major.axis[1] / 1000 / 111
    } else if ('latitudeError' %in% names(locs)){
      slat.sd <- locs$latitudeError[1]
    }
    
    L.light.lon <- stats::dnorm(t(locs.grid$lon), locs$Longitude[1], slon.sd)
    L.light.lat <- stats::dnorm(t(locs.grid$lat), locs$Latitude[1] + shiftDist, slat.sd)
    L.light <- L.light.lon * L.light.lat
    
  } else {
    
    # if supposed shift in error ellipse is >10 degrees or latitude is invalid, we revert to longitude only
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
