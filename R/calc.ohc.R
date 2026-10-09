#' Calculate Ocean Heat Content (OHC) likelihood surface
#' 
#' Compare tag data to OHC grid and calculate likelihoods
#' 
#' @param pdt input depth-temperature profile data. Need at least cols: Date (POSIXct), Depth, MinTemp, MaxTemp
#' @param filename is the first part of the filename specified to the download 
#'   function \code{\link{get.env}}. For example, if downloaded files were 
#'   specific to a particular dataset, you may want to identify that with a name
#'   like 'tuna' or 'shark1'. This results in a downloaded filename of, for 
#'   example, 'tuna_date.nc'. This filename is required here so the calc 
#'   function knows where to get the env data.
#' @param isotherm default '' in which isotherm is calculated on the fly based 
#'   on daily tag data. Otherwise, numeric isotherm constraint can be specified 
#'   (e.g. 20 deg C).
#' @param ohc.dir directory of downloaded hycom (or other) data
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param bathy is logical indicating whether or not a bathymetric mask should 
#'   be applied. This does NOT use an input bathymetry dataset but rather masks the depth-temperature likelihood based on the possible locations given the maximum recorded tag-depth (i.e. areas with bottom depth shallower than tag-recorded max depth are masked out).
#' @param use.se is logical indicating whether or not to use SE when using
#'   regression to predict temperature at specific depth levels.
#' @param flip_y is logical indicating whether or not to flip the resulting likelihood in the y. Set this to true if output likelihoods are upside down.
#'   
#' @return likelihood is raster brick of likelihood surfaces representing 
#'   estimated position based on tag-based OHC compared to calculated OHC using 
#'   HYCOM
#'   
#' @export
#' @seealso \code{\link{calc.woa}}
#' @references Luo J, Ault JS, Shay LK, Hoolihan JP, Prince ED, Brown C a.,
#'   Rooker JR (2015) Ocean Heat Content Reveals Secrets of Fish Migrations.
#'   PLoS One 10:e0141101
#'   

calc.ohc <- function(pdt, filename, isotherm = '', ohc.dir, dateVec, bathy = TRUE, use.se = TRUE, flip_y = FALSE){

  names(pdt) <- tolower(names(pdt))
  #options(warn=1)
  
  if (class(pdt$date)[1] != class(dateVec)[1]) stop('dateVec and pdt$date both need to be of class POSIXct.')
  
  t0 <- Sys.time()
  print(paste('Starting OHC likelihood calculation...'))
  
  # constants for OHC calc
  cp <- 3.993 # kJ/kg*C <- heat capacity of seawater
  rho <- 1025 # kg/m3 <- assumed density of seawater
  
  # calculate midpoint of tag-based min/max temps
  if(length(grep('mean', names(pdt))) > 1){
    pdt$useTemp <- pdt[,grep('mean', names(pdt))[1]]
  } else if(length(grep('mean', names(pdt))) == 1){
    pdt$useTemp <- pdt[,grep('mean', names(pdt))]
  } else{
    pdt$useTemp <- (pdt$maxtemp + pdt$mintemp) / 2
  }
  
  pdt <- pdt[which(pdt$date <= max(dateVec)),]
  pdt$dateVec <- findInterval(pdt$date, dateVec)
  T <- length(dateVec)
  
  if(isotherm != '') iso.def <- TRUE else iso.def <- FALSE
  
  # open nc and get the indices for the vars
  nc1 =  RNetCDF::open.nc(dir(ohc.dir, full.names = T)[1])
  ncnames = NULL
  nmax <- RNetCDF::file.inq.nc(nc1)$nvars - 1
  for(ii in 0:nmax) ncnames[ii + 1] <- RNetCDF::var.inq.nc(nc1, ii)$name
  temp.idx <- grep('temp', ncnames, ignore.case=TRUE) - 1
  lat.idx <- grep('lat', ncnames, ignore.case=TRUE) - 1
  lon.idx <- grep('lon', ncnames, ignore.case=TRUE) - 1
  if (any(ncnames == 'z')) ncnames[which(ncnames %in% 'z')] <- 'depth'
  dep.idx <- grep('dep', ncnames, ignore.case=TRUE) - 1
  
  ## better handling of depth
  if (length(dep.idx) == 0){
    depth <- c(0, 2, 4, 6, 8, 10, 12, 15, 20, 25,
               30, 35, 40, 45, 50, 60, 70, 80, 90,
               100, 125, 150, 200, 250, 300, 350, 
               400, 500, 600, 700, 800, 900, 1000,
               1250, 1500, 2000, 2500, 3000, 4000, 5000)
    
  } else{
    depth <- RNetCDF::var.get.nc(nc1, dep.idx)
    if (max(depth) < 100) depth <- c(0, 2, 4, 6, 8, 10, 12, 15, 20, 25,
                                     30, 35, 40, 45, 50, 60, 70, 80, 90,
                                     100, 125, 150, 200, 250, 300, 350, 
                                     400, 500, 600, 700, 800, 900, 1000,
                                     1250, 1500, 2000, 2500, 3000, 4000, 5000)
  }
    
  # get attributes, if they exist
  ncatts <- NULL
  nmax <- RNetCDF::var.inq.nc(nc1, temp.idx)$natts - 1
  for(ii in 0:nmax) ncatts[ii + 1] <- RNetCDF::att.inq.nc(nc1, temp.idx, ii)$name
  scale.idx <- grep('scale', ncatts, ignore.case=TRUE) - 1
  if(length(scale.idx) != 0){
    scale <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=scale.idx)
  } else{
    scale <- 1
  }
  off.idx <- grep('off', ncatts, ignore.case=TRUE) - 1
  if(length(off.idx) != 0){
    offset <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=off.idx)
  } else{
    offset <- 0
  }
  
  # get and check the vars
  lon <- RNetCDF::var.get.nc(nc1, lon.idx)
  if(length(dim(lon)) == 2) lon <- lon[,1]
  if(!any(lon < 180)) lon <- lon - 360
  lat <- RNetCDF::var.get.nc(nc1, lat.idx)
  if(length(dim(lat)) == 2) lat <- lat[1,]
  
  # result will be array of likelihood surfaces
  L.ohc <- array(0, dim = c(length(lon), length(lat), length(dateVec)))
  
  print(paste('Starting iterations through deployment period ', '...'))
  
  for(i in 1:T){
    
    print(dateVec[i])
    
    pdt.i <- pdt[which(pdt$dateVec == i),]
    if (nrow(pdt.i) == 0) next
    
    # open day's hycom data
    nc <- RNetCDF::open.nc(paste(ohc.dir, filename,'_', as.Date(dateVec[i]), '.nc', sep=''))
    dat <- RNetCDF::var.get.nc(nc, temp.idx) * scale + offset
    
    #extracts depth from tag data for day i
    y <- pdt.i$depth[!is.na(pdt.i$depth)] 
    y[y < 0] <- 0
    
    #extract temperature from tag data for day i
    x <- pdt.i$useTemp[!is.na(pdt.i$depth)]  
    
    # use the which.min
    depIdx = unique(apply(as.data.frame(pdt.i$depth), 1, FUN=function(x) which.min((x - depth) ^ 2)))
    hycomDep <- depth[depIdx]
    
    if(bathy){
      mask <- dat[,,max(depIdx)]
      mask[is.na(mask)] <- NA
      mask[!is.na(mask)] <- 1
      for(bb in 1:length(depth)){
        dat[,,bb] <- dat[,,bb] * mask
      }
    }
    
    # make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
    suppressWarnings(
    fit.low <- locfit::locfit(pdt.i$mintemp ~ pdt.i$depth, maxk=500)
    )
    suppressWarnings(
    fit.high <- locfit::locfit(pdt.i$maxtemp ~ pdt.i$depth, maxk=500)
    )
    n = length(hycomDep)
      
    #suppressWarnings(
    pred.low = stats::predict(fit.low, newdata = hycomDep, se = TRUE, get.data = TRUE)
    #suppressWarnings(
    pred.high = stats::predict(fit.high, newdata = hycomDep, se = TRUE, get.data = TRUE)
    
    if (use.se){
      # data frame for next step
      df = data.frame(low = pred.low$fit - pred.low$se.fit * sqrt(n),
                      high = pred.high$fit + pred.high$se.fit * sqrt(n),
                      depth = hycomDep)
    } else{
      # data frame for next step
      df = data.frame(low = pred.low$fit,# - pred.low$se.fit * sqrt(n),
                      high = pred.high$fit,# + pred.high$se.fit * sqrt(n),
                      depth = hycomDep)
    }

    # isotherm is minimum temperature recorded for that time point
    if(iso.def == FALSE) isotherm <- min(df$low, na.rm = TRUE)
    
    # perform tag data integration at limits of model fits
    minT.ohc <- cp * rho * sum(df$low - isotherm, na.rm = TRUE) / 10000
    maxT.ohc <- cp * rho * sum(df$high - isotherm, na.rm = TRUE) / 10000
    
    # Perform hycom integration
    #dat[dat < isotherm] <- NA
    dat <- dat - isotherm
    ohc <- cp * rho * apply(dat[,,depIdx], 1:2, sum, na.rm = TRUE) / 10000 
    ohc[ohc == 0] <- NA
    
    # calc sd of OHC
    # focal calc on mean temp and write to sd var
    r = raster::flip(raster::raster(t(ohc)), 2)
    sdx = raster::focal(r, w = matrix(1, nrow = 9, ncol = 9),
                        fun = function(x) stats::sd(x, na.rm = TRUE), pad=TRUE)
    sdx = t(raster::as.matrix(raster::flip(sdx, 2)))

    # compare hycom to that day's tag-based ohc
    #lik.ohc <- likint3(ohc, sdx, minT.ohc, maxT.ohc)
    #lik.ohc <- lik.ohc / max(lik.ohc, na.rm = T)
    
    lik.try <- try(likint3(ohc, sdx, minT.ohc, maxT.ohc), TRUE)
    
    if(class(lik.try)[1] == 'try-error' & use.se == FALSE){
      
      # try ohc again with use.se = T
      df = data.frame(low = pred.low$fit - pred.low$se.fit * sqrt(n),
                      high = pred.high$fit + pred.high$se.fit * sqrt(n),
                      depth = hycomDep)
      
      minT.ohc <- cp * rho * sum(df$low - isotherm, na.rm = TRUE) / 10000
      maxT.ohc <- cp * rho * sum(df$high - isotherm, na.rm = TRUE) / 10000
      
      lik.try <- try(likint3(ohc, sdx, minT.ohc, maxT.ohc), TRUE)
      
      if (class(lik.try)[1] == 'try-error'){
        lik.try <- ohc * 0
        warning(paste('Warning: likint3 failed after trying with and without SE prediction of depth-temp profiles. This is most likely a divergent integral for ', dateVec[i], '...', sep=''))
      }
      
    } else if (class(lik.try)[1] == 'try-error' & use.se == TRUE){
      lik.try <- ohc * 0
      warning(paste('Warning: likint3 failed after trying with and without SE prediction of depth-temp profiles. This is most likely a divergent integral for ', dateVec[i], '...', sep=''))
    }
    
    lik.ohc <- lik.try / max(lik.try, na.rm = TRUE)
    L.ohc[,,i] = lik.ohc

  }

  print(paste('Making final likelihood raster...'))
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  if(!any(lon < 180)) lon <- lon - 360
  list.ohc <- list(x = lon, y = lat, z = L.ohc)
  ex <- raster::extent(list.ohc)
  L.ohc <- raster::brick(list.ohc$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=TRUE, crs)
  if (flip_y){
    L.ohc <- raster::flip(L.ohc, direction = 'y')
    warning('Output raster is being flipped in the y. If this is not desired, use need_flip=FALSE.')
  }

  L.ohc[L.ohc < 0] <- 0
  
  names(L.ohc) = as.character(dateVec)
  
  t1 <- Sys.time()
  print(paste('OHC calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  # return ohc likelihood surfaces
  return(L.ohc)
  
}

