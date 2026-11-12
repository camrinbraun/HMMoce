#' Calculate bt-based likelihood in parallel
#' 
#' \code{calc.bottomTemp.par} compares tag-measured bottom temperature to a bottom temperature grid and calculates 
#' likelihoods in parallel
#' 
#' @param tag.bt is data frame containing tag-collected bottom temperature data. Requires at least cols: Date (POSIXct), Temperature
#' @param filename is the first part of the filename specified to the download 
#'   function \code{\link{get.env}}. For example, if downloaded files were 
#'   specific to a particular dataset, you may want to identify that with a name
#'   like 'tuna' or 'shark1'. This results in a downloaded filename of, for 
#'   example, 'tuna_date.nc'. This filename is required here so the calc 
#'   function knows where to get the env data.
#' @param bt.dir local directory where bottom temp nc grids are stored
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param focalDim is integer for dimensions of raster::focal used to calculate 
#'   sd() of env grid cell. If left to NULL (default), this dimension will
#'   approximately incorporate 0.25 degrees.
#' @param sens.err is numeric indicating the percent sensor error in the tag bottom temp
#'   sensor. Default is 1.
#' @param varName is name of temperature variable in the bottom temperature nc grids. Default is "Temperature".
#' @param ncores is integer indicating number of cores used in this parallel
#'  computation. Defaults to using a detection function that chooses cores for
#'  you.
#'   
#' @return likelihood is raster brick of likelihood surfaces representing
#'   matches between tag-based bottom temperature and environmental bottom temperature grids
#'   
#' @export
#' 
#' @seealso \code{\link{calc.sst.par}}
#'   

calc.bottomTemp.par <- function(tag.bt, filename, bt.dir, dateVec, focalDim = NULL, sens.err = 1, varName = 'Temperature', ncores = NULL){
  
  print(paste('Starting bottom temperature likelihood calculation...'))
  
  if (is.null(ncores)) ncores <- ceiling(parallel::detectCores() * .9)
  if (is.na(ncores) | ncores < 0) ncores <- ceiling(as.numeric(system('nproc', intern=T)) * .9)
  
  t0 <- Sys.time()
  
  if (class(tag.bt$Date)[1] != class(dateVec)[1]) stop('dateVec and tag.sst$Date both need to be of class POSIXct.')

  tag.bt <- tag.bt[which(tag.bt$Date <= max(dateVec)),]
  tag.bt$dateVec <- findInterval(tag.bt$Date, dateVec)
  tag.bt <- data.frame(tag.bt %>% group_by(dateVec) %>% 
                         dplyr::summarise(minT = min(Temperature, na.rm=T), maxT = max(Temperature, na.rm=T), .groups = 'drop_last'))
  
  T <- length(dateVec)
  
  print(paste('Starting iterations through deployment period ', '...'))
  
  # GET EVERYTHING SETUP BEFORE PARALLEL
  time1 <- tag.bt$time[1]
  bt1 <- c(tag.bt$minT[1] * (1 - sens.err), tag.bt$maxT[1] * (1 + sens.err)) # sensor error
  
  # open day's bt data
  nc1 <- RNetCDF::open.nc(paste(bt.dir, filename, '_', as.Date(dateVec[tag.bt$dateVec[1]]), '.nc', sep='')) #add lat lon in filename '.nc', sep=''))
  
  # get correct name in bt data
  ncnames = NULL
  nmax <- RNetCDF::file.inq.nc(nc1)$nvars - 1
  for(ii in 0:nmax) ncnames[ii + 1] <- RNetCDF::var.inq.nc(nc1, ii)$name
  nameidx <- grep(varName, ncnames, ignore.case = TRUE) - 1
  dat <- RNetCDF::var.get.nc(nc1, nameidx)
  lon <- RNetCDF::var.get.nc(nc1, 'longitude')
  lat <- RNetCDF::var.get.nc(nc1, 'latitude')
  
  # result will be array of likelihood surfaces
  L.bt <- array(0, dim = c(length(lon), length(lat), length(dateVec)))
  
  # calc sd of bt via terra
  r = terra::flip(terra::rast(t(dat), extent = terra::ext(min(lon), max(lon), min(lat), max(lat)), crs = "EPSG:4326"), direction = "vertical")
  
  if(round(terra::res(r)[1], 2) < 0.1){
    aggFact <- round(0.1 / round(terra::res(r)[1], 2), 0)
    if(aggFact > 1) r <- terra::aggregate(r, fact = aggFact)
  } 
  
  if(is.null(focalDim)){
    focalDim <- round(0.25 / terra::res(r)[1], 0)
    if(focalDim %% 2 == 0) focalDim <- focalDim - 1
    if (focalDim == 1) focalDim <- 3
  }
  
  tmpl_ext <- as.vector(terra::ext(r))
  L.bt <- array(0, dim = c(terra::ncol(r), terra::nrow(r), length(dateVec)))
  
  ## END of SETUP RUN
  
  print('Processing in parallel... ')
  
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  
  ans = foreach::foreach(i = 1:T, .export = 'likint3', .packages = c('terra')) %dopar%{
    
    tag.bt.i <- tag.bt[which(tag.bt$dateVec == i),]
    if (nrow(tag.bt.i) == 0) return(NA)
    
    bt.i <- c(tag.bt.i$minT * (1 - sens.err / 100), tag.bt.i$maxT * (1 + sens.err / 100))
    
    nc <- RNetCDF::open.nc(paste(bt.dir, filename, '_', as.Date(dateVec[i]), '.nc', sep=''))
    dat <- RNetCDF::var.get.nc(nc, nameidx)
    
    r <- terra::flip(terra::rast(t(dat), extent = terra::ext(min(lon), max(lon), min(lat), max(lat)), crs = "EPSG:4326"), direction = "vertical")
    
    if(round(terra::res(r)[1], 2) < 0.1){
      aggFact <- round(0.1 / round(terra::res(r)[1], 2), 0)
      if(aggFact != 1) r <- terra::aggregate(r, fact = aggFact)
    } 
    
    sd_r <- terra::focal(r, w = focalDim, fun = "sd", na.rm = TRUE)
    sdx  <- t(terra::as.matrix(terra::flip(sd_r, direction = "vertical"), wide = TRUE))
    dat  <- t(terra::as.matrix(terra::flip(r, direction = "vertical"), wide = TRUE))
    
    lik.bt <- likint3(dat, sdx, bt.i[1], bt.i[2])
    lik.bt
  }
  
  parallel::stopCluster(cl)
  
  lik.bt <- lapply(ans, function(x) x / max(x, na.rm = T))
  
  for (i in 1:T){
    L.bt[,,i] <- lik.bt[[i]]
  }
  
  print(paste('Making final likelihood raster...'))
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  L.bt[is.na(L.bt) | L.bt < 0] <- 0
  L.bt <- terra::flip(terra::rast(aperm(L.bt, c(2, 1, 3)), extent = terra::ext(tmpl_ext), crs = crs), direction = "vertical")
  names(L.bt) <- as.character(dateVec)
  
  t1 <- Sys.time()
  print(paste('Bottom temperature calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  # return bt likelihood surfaces
  return(L.bt)
  
}
