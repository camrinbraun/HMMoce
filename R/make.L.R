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
  
  ## generate blank results likelihood
  L <- ras.list[[1]] * 0
  L[is.na(L)] <- 0
  
  # check bathy mask requirements
  if (!is.null(bathy) && !is.null(maxDepth) && any(maxDepth > 0)) {
    if (length(maxDepth) != length(dateVec)){
      stop('Error: Vector of daily maximum depths is not same length at date vector. Ensure mmd was zero-padded for days with NA data.')
    }
    
    bathy <- raster::readAll(bathy)
    ## convert a negative bathy grid to positive to match expectations and mask land
    if (raster::cellStats(bathy, 'min', na.rm=T) < 0){
      bathy <- bathy * -1
      bathy[bathy < 0] <- NA
    }
    bathy <- raster::resample(bathy, L)
  } else {
    bathy <- NULL
  }
  
  ## COMBINE THE LIKELIHOOD LAYERS
  ## for each day:
  for (i in 1:length(dateVec)){
    
    #print(i)
    
    ## get relevant likelihoods across the list
    for (bb in 1:length(ras.list)){
      if (bb == 1){
        s <- raster::stack(ras.list[[bb]][[i]])
      } else{
        s <- raster::stack(s, ras.list[[bb]][[i]])
      }
    }
    
    ## check for layers with all NA or all <= 0
    sum_NA <- raster::cellStats(!is.na(s), sum, na.rm = TRUE) == 0
    max_vals <- suppressWarnings(raster::cellStats(s, 'max', na.rm = TRUE))
    empty_layer <- sum_NA | is.na(max_vals) | is.infinite(max_vals) | (max_vals <= 0)
    
    ## if all layers are empty, skip to next day
    if (all(empty_layer)){
      next
    } else if (any(empty_layer)){
      s <- s[[-which(empty_layer)]]
    }
    
    ## sum & normalize whatever layers remain
    if (raster::nlayers(s) == 1){
      L[[i]] <- s / raster::cellStats(s, 'max', na.rm = TRUE)
    } else{
      s_sum <- sum(s, na.rm = TRUE)
      L[[i]] <- s_sum / raster::cellStats(s_sum, 'max', na.rm = TRUE)
    }
    
    ## daily bathy mask
    if (!is.null(bathy) && !is.null(maxDepth) && maxDepth[i] > 0) {
      bathy.i <- maxDepth[i] * (1 - 20 / 100) # give 20% of max depth buffer for likelihood flexibility
      L[[i]][bathy < bathy.i] <- 0
    }
    
  }
  
  ## ADD START/END LOCATIONS
  if(!is.null(iniloc)){
    # convert input date, lat, lon to likelihood surfaces
    print('Adding start and end locations from iniloc...')
    
    # get lat/lon vectors
    #lon <- seq(raster::extent(L)[1], raster::extent(L)[2], length.out=dim(L)[2])
    #lat <- seq(raster::extent(L)[3], raster::extent(L)[4], length.out=dim(L)[1])
    iniloc$date <- as.POSIXct(
      sprintf('%04d-%02d-%02d', as.numeric(iniloc$year), as.numeric(iniloc$month), as.numeric(iniloc$day)),
      tz = 'UTC'
    )
    if (class(iniloc$date)[1] != class(dateVec)[1]) stop('dateVec and known.locs$date both need to be of class POSIXct.')
    iniloc$dateVec <- findInterval(iniloc$date, dateVec)
    
    for(i in unique(iniloc$dateVec)){
      known.locs.i <- iniloc[which(iniloc$dateVec == i),]
      
      #x = which.min((known.locs.i$lon - lon) ^ 2)
      #y = which.min((known.locs.i$lat - lat) ^ 2)
      
      ## erase other likelihood results for this day
      L[[i]] <- L[[i]] * 0
      
      # assign the known location for this day, i, as 1 (known) in likelihood raster
      L[[i]][raster::cellFromXY(L[[i]], known.locs.i[,c('lon','lat')])] <- 1
      
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
      
      if(length(known.locs.i[,1]) > 1){
        # if multiple known locations are provided for a given day, only the first is used
        warning(paste0('Multiple locations supplied at time step ', dateVec[i], '. Only the first one is being used.'))
        known.locs.i <- known.locs.i[1,]
      }
      
      ## erase other likelihood results for this day
      L[[i]] <- L[[i]] * 0
      
      ## assign the known location for this day, i, as 1 (known) in likelihood raster
      L[[i]][raster::cellFromXY(L[[i]], known.locs.i[,c('lon','lat')])] <- 1
      
    }
    
  }
  
  #----------------------------------------------------------------------------------#
  # MAKE ALL NA'S VERY TINY FOR THE CONVOLUTION
  # the previous steps may have taken care of this...
  #----------------------------------------------------------------------------------#
  L[L <= 1e-15] <- 1e-15
  L[is.na(L)] <- 1e-15
  
  # CONVERT OUTPUT RASTER INTO AN ARRAY
  L <- aperm(raster::as.array(raster::flip(L, direction = 'y')), c(3, 2, 1))
  
  print('Finishing make.L...', sep='')
  
  return(L)
}
