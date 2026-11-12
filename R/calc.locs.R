#' Calculate Position-based Likelihood
#' 
#' \code{calc.locs} calculates likelihood estimates for each day of animal tag 
#' data using tag-based locations.
#' 
#' GPS and Argos positions are given a "likelihood" using this function but are
#' currently both considered to be "known" positions without error.
#' 
#' @param locs is data frame from -Locations file output from DAP/Tag Portal for WC tags and 
#'   contains GPS, Argos, and GPE locations as applicable.
#' @param gps is data frame from -FastGPS file output from WC Tag Portal
#' @param iniloc is 2 x 5 dataframe containing day, month, year, lat, lon for 
#'   both tag and pop locations
#' @param dateVec is vector of dates from tag to pop-up in 1 day increments.
#' @param locs.grid is list output from \code{setup.locs.grid}
#'   
#' @return L is an array of lon x lat likelihood surfaces (matrices) for each 
#'   time point (3rd dimension)
#'   

calc.locs <- function(locs, gps = NULL, iniloc, locs.grid, dateVec){
  
  print(paste('Starting likelihood calculation...'))
  
  # get rid of GPE locations and check length of resulting locs file
  
  locs <- locs[which(locs$Type != 'GPE'),]
  if(length(locs[,1]) < 1){
    stop('No non-GPE positions available in input locs file. Make sure you have non-GPE positions in your -Locations.csv file.')
  }
  
  # set some date vectors for our input data
  locDates <- as.Date(locs$Date, format = findDateFormat(locs$Date))
  if(!is.null(gps)){
    gpsDates <- as.Date(gps$InitTime, format = findDateFormat(gps$InitTime))  
  } else{
    gpsDates = NULL
  }
  
  if(any(duplicated(locDates))){
    # run a simplify function
    locList <- simplifyLocs(locs, locDates)
    locs <- locList$locs; locDates <- locList$locDates
  }
  
  if(any(duplicated(gpsDates))){
    # run a simplify function
    locList <- simplifyLocs(gps, gpsDates)
    gps <- locList$locs; gpsDates <- locList$locDates
    
  }
  
  # set up results array
  row <- dim(locs.grid$lon)[1]
  col <- dim(locs.grid$lon)[2]
  L.locs <- array(0, dim = c(col, row, length(dateVec)))
  
  # add tag/pop locations as known
  ilo <- which.min(abs(locs.grid$lon[1,] - iniloc$lon[1]))
  ila <- which.min(abs(locs.grid$lat[,1] - iniloc$lat[1]))
  L.locs[ilo, ila, 1] <- 1   # Initial location is known
  
  elo <- which.min(abs(locs.grid$lon[1,] - iniloc$lon[2]))
  ela <- which.min(abs(locs.grid$lat[,1] - iniloc$lat[2]))
  L.locs[elo, ela, length(dateVec)] <- 1  # End location is known

  print(paste('Starting iterations through deployment period...'))
  
  for(t in 2:(length(dateVec)) - 1){
    
    if(!is.null(gps) & dateVec[t] %in% gpsDates){
      
      # if GPS exists then other forms of data for that time point are obsolete
      idx <- which(gpsDates == dateVec[t]) # set index to identify position in gps file
      glo <- which.min(abs(locs.grid$lon[1,] - gps$InitLon[idx]))
      gla <- which.min(abs(locs.grid$lat[,1] - gps$InitLat[idx]))
      L.locs[glo, gla, t] <- 1
      
    } else if(!is.null(locs) & dateVec[t] %in% locDates){
      # set index to identify position in locs file
      idx <- which(locDates == dateVec[t])
      
      if(locs$Type[idx] == 'GPS'){ #locs includes GPS
        # if GPS exists then other forms of data for that time point are obsolete
        glo <- which.min(abs(locs.grid$lon[1,] - locs$Longitude[idx]))
        gla <- which.min(abs(locs.grid$lat[,1] - locs$Latitude[idx]))
        L.locs[glo, gla, t] <- 1
        
      } else if(locs$Type[idx] == 'Argos'){ #locs includes Argos
        # if Argos exists, GPE positions are obsolete
        alo <- which.min(abs(locs.grid$lon[1,] - locs$Longitude[idx]))
        ala <- which.min(abs(locs.grid$lat[,1] - locs$Latitude[idx]))
        L.locs[alo, ala, t] <- 1
        
      } else{
        
        stop('Error: Error ellipse unspecified.')
        
      }
      
    } else{
      
      # no data so we skip this day
      
    }
    
  }
  
  print(paste('Making final likelihood raster...'))
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  L_arr <- aperm(L.locs, c(2, 1, 3))
  L_arr <- L_arr[dim(L_arr)[1]:1, , , drop = FALSE]
  L.locs <- terra::rast(
    L_arr,
    extent = terra::ext(min(locs.grid$lon[1,]), max(locs.grid$lon[1,]), min(locs.grid$lat[,1]), max(locs.grid$lat[,1])),
    crs = crs
  )
  names(L.locs) <- as.character(dateVec)
  
  return(list(L.locs = L.locs, gpsIdx = which(dateVec %in% gpsDates)))
  
}
