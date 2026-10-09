#' Calculate Light Likelihood from lightloc Output
#' 
#' \code{calc.lightloc} calculates likelihood estimates for each day of animal tag 
#' data.
#' 
#' Light errors are parameterized using elliptical error values output in 
#' '-Locations.csv' (WC tags).
#' 
#' @param lightloc is data frame of light-based location estimates. If errEll is FALSE, only Date (POSIXct), Longitude, and Error.Semi.minor.axis (in meters, default output from WC tags) are required. longitudeError (in decimal degrees) can be supplied in place of Error.Semi.minor.axis. If errEll is TRUE, additional required columns are Latitude and Error.Semi.major.axis (meters). These are default outputs from WC tags. Offset (meters) and Offset.orientation (degrees of rotation) are also default outputs from WC tags and are optional columns to include in input. latitudeError (decimal degrees) can be supplied in place of Error.Semi.major.axis. In that case, offset variables are ignored. 
#' @param locs.grid is list output from \code{setup.locs.grid}
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param errEll is logical indicating whether error ellipses should be 
#'   generated for light-based likelihoods. If FALSE (default), only longitude is used to generate likelihoods. If TRUE, both latitude and longitude are used resulting in an ellipse-shaped likelihood. 
#' @param lon_only is logical indicating whether the likelihood should be generated for longitude only (default is TRUE). If FALSE, latitude only will be used to generate the likelihood. This argument is only valid is errEll = FALSE.
#' @references Musyl MK, Domeier ML, Nasby-Lucas N, Brill RW, McNaughton LM, 
#'   Swimmer JY, Lutcavage MS, Wilson SG, Galuardi B, Liddle JB (2011) 
#'   Performance of pop-up satellite archival tags. Mar Ecol Prog Ser
#' @export
#' @return L is an array of lon x lat likelihood surfaces (matrices) for each 
#'   time point (3rd dimension)
#' @seealso \code{\link{calc.srss}}
#' @examples
#' \dontrun{
#' ## Setup for calculating light likelihood
#' # Read the data
#' locsFile <- system.file("extdata", "141259-Locations-GPE2.csv", package = "HMMoce")
#' locs <- read.table(locsFile, sep = ',', header = TRUE, blank.lines.skip = FALSE)
#' 
#' # Set spatial and temporal limits
#' sp.lim <- list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50)
#' locs.grid <- setup.locs.grid(sp.lim)
#' iniloc <- data.frame(matrix(c(13, 10, 2015, 41.3, -69.27, 10, 4, 2016, 40.251, -36.061),
#'  nrow = 2, ncol = 5, byrow = TRUE))
#' names(iniloc) <- list('day','month','year','lat','lon')
#' tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), 
#' format = '%d/%m/%Y', tz='UTC')
#' pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), 
#' format = '%d/%m/%Y', tz='UTC')
#' dateVec <- as.Date(seq(tag, pop, by = 'day')) 
#' 
#' # Full example light calculation
#' L.light <- calc.lightloc(locs, iniloc = iniloc, locs.grid = locs.grid,
#'                      dateVec = dateVec, errEll = TRUE, gpeOnly = TRUE)
#' }

calc.lightloc <- function(lightloc, locs.grid, dateVec, errEll = TRUE, lon_only = TRUE){
  
  print(paste('Starting light likelihood calculation...'))
  t0 <- Sys.time()
  
  # check date formats match
  lightloc <- lightloc[which(lightloc$Date <= max(dateVec)),]
  lightloc$dateVec <- findInterval(lightloc$Date, dateVec)
  T <- length(dateVec)
  
  # set up results array
  row <- dim(locs.grid$lon)[1]
  col <- dim(locs.grid$lon)[2]
  lat <- locs.grid$lat[,1]
  lon <- locs.grid$lon[1,]
  L.lightloc <- array(0, dim = c(col, row, length(dateVec)))
  
  print(paste('Starting iterations through deployment period...'))
  
  for(t in 1:T){
    
    # data for this time step T
    light.t <- lightloc[which(lightloc$dateVec == t),]
    L.lightloc.try <- array(0, dim = c(col, row, nrow(light.t)))
    
    if (nrow(light.t) == 0) next
    
    for (ii in 1:nrow(light.t)){
      
      locs.ii <- light.t[ii,]
      
      if(errEll == FALSE){ ## then we only care about longitude
        
        if (lon_only){
          if ('Error.Semi.minor.axis' %in% names(locs.ii)){
            if (locs.ii$Error.Semi.minor.axis[1] < 70000) warning('Some values of Error.Semi.minor.axis are < 70km which usually does not represent the actual error in these measurements.')
            # create longitude likelihood based on GPE data
            slon.sd <- locs.ii$Error.Semi.minor.axis[1] / 1000 / 111 #semi minor axis
            
          } else if ('longitudeError' %in% names(locs.ii)){
            if (locs.ii$longitudeError[1] < 0.7) warning('Some values of longitudeError are < 0.7deg which usually does not represent the actual error in these measurements.')
            # create longitude likelihood based on GPE data
            slon.sd <- locs.ii$longitudeError #semi minor axis
            
          } else{
            warning('No longitude error specified. Fixed at 0.7deg based on Musyl et al 2011. Although note that this error is likely higher for many species.')
            slon.sd <- 0.7 ## musyl et al 2011
          }
          
          # use normally distributed error from position using fixed std dev
          L.lightloc.try[,,ii] <- stats::dnorm(t(locs.grid$lon), locs.ii$Longitude[1], slon.sd)
          
        } else{
          ## lat only
          if ('Error.Semi.major.axis' %in% names(locs.ii)){
            if (locs.ii$Error.Semi.major.axis[1] < 70000) warning('Some values of Error.Semi.major.axis are < 70km which usually does not represent the actual error in these measurements.')
            # create longitude likelihood based on GPE data
            slat.sd <- locs.ii$Error.Semi.major.axis[1] / 1000 / 111 #semi minor axis
            
          } else if ('latitudeError' %in% names(locs.ii)){
            if (locs.ii$latitudeError[1] < 0.7) warning('Some values of latitudeError are < 0.7deg which usually does not represent the actual error in these measurements.')
            # create longitude likelihood based on GPE data
            slat.sd <- locs.ii$latitudeError #semi minor axis
            
          } else{
            warning('No longitude error specified. Fixed at 3.5deg based on Doherty et al 2017 and Biais et al 2017.')
            slat.sd <- 3.5
          }
          
          # use normally distributed error from position using fixed std dev
          L.lightloc.try[,,ii] <- stats::dnorm(t(locs.grid$lat), locs.ii$Latitude[1], slat.sd)
          
        }
       
      } else if(errEll == TRUE){
        #if (locs.ii$Error.Semi.minor.axis[1] < 100000) warning('Some values of Error.Semi.minor.axis are < 100km which usually does not represent the actual error in these measurements.')
       
        if ('latitudeError' %in% names(locs.ii)){
         print('Offset variables are being set to 0 because latitudeError was supplied')
         locs.ii$Offset <- 0
         locs.ii$Offset.orientation <- 0
        }
        
        L.lightloc.try[,,ii] <- calc.errEll(locs.ii[1,], locs.grid)
        
      } 
      
    }
    
    ## sum over ii's in case there is more than 1 row of data per day
    L.lightloc[,,t] <- apply(L.lightloc.try, 1:2, sum, na.rm = T)
    
    ## normalize
    L.lightloc[,,t] = L.lightloc[,,t] / max(L.lightloc[,,t], na.rm=T)
   
  }
  
  print(paste('Making final likelihood raster...'))
  
  # this performs some transformations to the likelihood array to convert to useable raster
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.locs <- list(x = locs.grid$lon[1,], y = locs.grid$lat[,1], z = L.lightloc)
  ex <- raster::extent(list.locs)
  L.lightloc <- raster::brick(list.locs$z, xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], transpose = T, crs)
  L.lightloc <- raster::flip(L.lightloc, direction = 'y')
  
  L.lightloc[L.lightloc < 0] <- 0
  
  t1 <- Sys.time()
  print(paste('Light calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  return(L.lightloc)
  
}
