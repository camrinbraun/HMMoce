#' Hycom Profile LIkelihood in Parallel
#' 
#' Calculate Hycom profile likelihood surface in parallel
#' 
#' @param pdt input depth-temperature profile data. Need at least cols: Date (POSIXct), Depth, MinTemp, MaxTemp 
#' @param filename is the first part of the filename specified to the download 
#'   function \code{\link{get.env}}. For example, if downloaded files were 
#'   specific to a particular dataset, you may want to identify that with a name
#'   like 'tuna' or 'shark1'. This results in a downloaded filename of, for 
#'   example, 'tuna_date.nc'. This filename is required here so the calc 
#'   function knows where to get the env data.
#' @param hycom.dir directory of downloaded hycom (or other) data
#' @param focalDim is integer for dimensions of raster::focal used to calculate 
#'   sd() of temperature grid cell. Recommend focalDim = 9 for Hycom data at 
#'   0.08deg resolution.
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param use.se is logical indicating whether or not to use SE when using 
#'   regression to predict temperature at specific depth levels.
#' @param ncores specify number of cores, or leave blank and use whatever you 
#'   have!
#' @param flip_y is logical indicating whether or not to flip the resulting likelihood in the y. Set this to true if output likelihoods are upside down.
#'   
#' @return a raster brick of Hycom profile likelihood
#' @export
#' @importFrom foreach %dopar%
#'

calc.hycom.par <- function(pdt, filename, hycom.dir, focalDim = 9, dateVec, use.se = TRUE, ncores = NULL, flip_y = FALSE){
  
  names(pdt) <- tolower(names(pdt))
  #options(warn=-1)
  if (class(pdt$date)[1] != class(dateVec)[1]) stop('dateVec and pdt$date both need to be of class POSIXct.')
  
  t0 <- Sys.time()
  print(paste('Starting Hycom profile likelihood calculation...'))
  
  if (is.null(ncores)) ncores <- ceiling(parallel::detectCores() * .9)
  if (is.na(ncores) | ncores < 0) ncores <- ceiling(as.numeric(system('nproc', intern=T)) * .9)
  
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
  
  print(paste0('Generating profile likelihood for ', dateVec[1], ' through ', dateVec[length(dateVec)]))
  
  # open nc and get the indices for the vars
  nc1 =  RNetCDF::open.nc(dir(hycom.dir, full.names = TRUE)[1])
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
  #depth <- RNetCDF::var.get.nc(nc1, dep.idx)
  lon <- RNetCDF::var.get.nc(nc1, lon.idx)
  if(length(dim(lon)) == 2) lon <- lon[,1]
  if(!any(lon < 180)) lon <- lon - 360
  lat <- RNetCDF::var.get.nc(nc1, lat.idx)
  if(length(dim(lat)) == 2) lat <- lat[1,]
  
  # result will be array of likelihood surfaces
  L.hycom <- array(0, dim = c(length(lon), length(lat), length(dateVec)))
  
  # BEGIN PARALLEL STUFF  
  
  
  # ncores = detectCores()  # should be an input argument
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  
  print(paste0('Processing in parallel using ', ncores, ' cores... '))
  
  ans = foreach::foreach(i = 1:T, .export = 'likint3', .packages = c('raster')) %dopar%{
    
    pdt.i <- pdt[which(pdt$dateVec == i),]
    if (nrow(pdt.i) == 0) return(NA)
    
    # open day's hycom data
    nc <- RNetCDF::open.nc(paste(hycom.dir, filename, '_', as.Date(dateVec[i]), '.nc', sep=''))
    dat <- RNetCDF::var.get.nc(nc, temp.idx) * scale + offset
    
    #extracts depth from tag data for day i
    y <- pdt.i$depth[!is.na(pdt.i$depth)] 
    y[y < 0] <- 0
    
    #extract temperature from tag data for day i
    x <- pdt.i$useTemp[!is.na(pdt.i$depth)]  
    
    # use the which.min
    depIdx = unique(apply(as.data.frame(pdt.i$depth), 1, FUN = function(x) which.min((x - depth) ^ 2)))
    hycomDep <- depth[depIdx]
    
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
    
    # calculate sd using Le Bris neighbor method and focal()
    sd.i = array(NA, dim = c(dim(dat)[1:2], length(depIdx)))
    
    for(ii in 1:length(depIdx)){
      r = raster::flip(raster::raster(t(dat[,,depIdx[ii]])), 2)
      f1 = raster::focal(r, w = matrix(1, nrow = focalDim, ncol = focalDim), fun = function(x) stats::sd(x, na.rm = TRUE), pad = TRUE)
      f1 = t(raster::as.matrix(raster::flip(f1, 2)))
      sd.i[,,ii] = f1
    }
    
    # setup the likelihood array for each day. Will have length (dim[3]) = n depths
    lik.pdt = array(NA, dim = c(dim(dat)[1], dim(dat)[2], length(depIdx)))
    
    for (b in 1:length(depIdx)) {
      #calculate the likelihood for each depth level, b
      lik.try <- try(likint3(dat[,,depIdx[b]], sd.i[,,b], df[b, 1], df[b, 2]), TRUE)
      class.try <- class(lik.try)
      
      if(!any(which(lik.try > 0))) class.try <- 'try-error'
      
      if(class.try[1] == 'try-error' & use.se == FALSE){
        df[b,1] <- pred.low$fit[b] - pred.low$se.fit[b] * sqrt(n)
        df[b,2] <- pred.high$fit[b] - pred.high$se.fit[b] * sqrt(n)
        
        lik.try <- try(likint3(dat[,,depIdx[b]], sd.i[,,b], df[b, 1], df[b, 2]), TRUE)
        class.try <- class(lik.try)
        
        if(!any(which(lik.try > 0))) class.try <- 'try-error'
        
        if (class.try[1] == 'try-error'){
          lik.try <- dat[,,depIdx[b]] * 0
          warning(paste('Warning: likint3 failed after trying with and without SE prediction of depth-temp profiles. This is most likely a divergent integral for ', dateVec[i], '...', sep=''))
        }
        
      } else if (class.try[1] == 'try-error' & use.se == TRUE){
        lik.try <- dat[,,depIdx[b]] * 0
        warning(paste('Warning: likint3 failed after trying with and without SE prediction of depth-temp profiles. This is most likely a divergent integral for ', dateVec[i], '...', sep=''))
      }
      
      lik.pdt[,,b] <- lik.try
      
    }
    
    lik.pdt0 <- lik.pdt
    lik.pdt0[is.na(lik.pdt0)] <- 0
    use.idx <- unique(which(lik.pdt0 != 0, arr.ind=TRUE)[,3])
    
    # multiply likelihood across depth levels for each day
    lik.pdt <- apply(lik.pdt[,,use.idx], 1:2, FUN=function(x) prod(x, na.rm=FALSE))
    
    lik.pdt
    
  }
  
  parallel::stopCluster(cl)
  
  # lapply to put parallel answers back together
  lik.pdt = lapply(ans, function(x) x / max(x, na.rm = TRUE))
  
  for(i in 1:T){
    L.hycom[,,i] = lik.pdt[[i]]
  }
  
  print(paste('Making final likelihood raster...'))
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  L.hycom <- raster::brick(L.hycom, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), transpose=TRUE, crs)
  if (flip_y){
    L.hycom <- raster::flip(L.hycom, direction = 'y')
    warning('Output raster is being flipped in the y. If this is not desired, use need_flip=FALSE.')
  }
  
  names(L.hycom) = as.character(dateVec)
  
  t1 <- Sys.time()
  print(paste('Hycom profile calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  #options(warn=2)
  
  # return hycom likelihood surfaces
  return(L.hycom)
  
}
