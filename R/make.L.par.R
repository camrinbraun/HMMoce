#' Combine individual source likelihoods
#' 
#' \code{make.L} combines individual likelihoods from various data sources (e.g.
#' SST, OHC) to make overall combined likelihoods for each time point
#' 
#' @param ras.list a list of likelihood rasters
#' @param iniloc is data.frame of tag and pop locations. Required columns are 'date' (POSIX), 'lon', 'lat'.
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param known.locs is data frame of known locations containing named columns 
#'   of 'date' (POSIX), 'lon', 'lat'. Default is NULL.
#' @param ncores specify number of cores, or leave blank and use whatever you 
#'   have!
#' @return an overall likelihood array, L
#' 
#' @examples
#' \dontrun{
#' L <- make.L(ras.list, iniloc, dateVec)
#' }
#' @export
#' @importFrom foreach %dopar%

make.L.par <- function(ras.list, iniloc, dateVec, known.locs = NULL, ncores = NULL){
  
  ## generate blank results likelihood
  L <- ras.list[[1]] * 0
  L[is.na(L)] <- 0
  
  if (is.null(ncores)) ncores <- ceiling(parallel::detectCores() * .9)
  if (is.na(ncores) | ncores < 0) ncores <- ceiling(as.numeric(system('nproc', intern=T)) * .9)
  
  # BEGIN PARALLEL STUFF  
  
  # ncores = detectCores()  # should be an input argument
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  
  print(paste0('Processing in parallel using ', ncores, ' cores... '))
  
  ans = foreach::foreach(i = 1:length(dateVec)) %dopar%{
    
  ## COMBINE THE LIKELIHOOD LAYERS
  ## for each day:
  #for (i in 1:length(dateVec)){
    
    #print(i)
    
    ## get relevant likelihoods across the list
    for (bb in 1:length(ras.list)){
      if (bb == 1){
        s <- raster::stack(ras.list[[bb]][[i]])
      } else{
        s <- raster::stack(s, ras.list[[bb]][[i]])
      }
    }
    
    ## check for layers with all NA
    sum_NA <- raster::cellStats(!is.na(s), sum, na.rm=T) == 0
    
    ## check for layers that sum to 0
    sum_zero <- rep(NA, raster::nlayers(s))
    for (bb in 1:raster::nlayers(s)){
      s_bb <- s[[bb]]
      s_bb[s_bb == 0] <- 1
      
      sum_zero[bb] <- ifelse(raster::cellStats(s_bb, 'sum') == raster::ncell(s), TRUE, FALSE)
      
    }
    
    ## if all layers are NA, just fill with a 0 raster
    if (length(which(sum_NA)) == raster::nlayers(s)){
      ## nothing left for this iteration
      s <- s * 0
      
      ## if a layer is all NA, drop it
    } else if (length(which(sum_NA)) < raster::nlayers(s) & length(which(sum_NA)) > 0){
      s <- s[[-which(sum_NA)]]
    } else{
      ## do nothing
    }
    
    ## sum & normalize whatever layers remain
    if (raster::nlayers(s) == 1){
      l.ans <- s / raster::cellStats(s, 'max') ## do not remove NA yet
    } else{
      sum_s <- sum(s)
      l.ans <- sum_s / raster::cellStats(sum_s, 'max') ## do not remove NA yet  
    }
    
    l.ans
  }
  
  parallel::stopCluster(cl)
  
  
  for(i in 1:length(dateVec)){
    L[[i]] = ans[[i]]
  }
  
  ## ADD START/END LOCATIONS
  if(!is.null(iniloc)){
    # convert input date, lat, lon to likelihood surfaces
    print('Adding start and end locations from iniloc...')
    
    # get lat/lon vectors
    #lon <- seq(raster::extent(L)[1], raster::extent(L)[2], length.out=dim(L)[2])
    #lat <- seq(raster::extent(L)[3], raster::extent(L)[4], length.out=dim(L)[1])
    iniloc$date <- as.POSIXct(paste(iniloc$year, iniloc$month, iniloc$day, sep='-'), format='%Y-%m-%d', tz='UTC')
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

