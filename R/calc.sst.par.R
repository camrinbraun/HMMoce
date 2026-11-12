#'Calculate SST-based likelihood in parallel
#'
#'\code{calc.sst.par} compares tag SST to remotely sensed SST and calculates 
#'likelihoods in parallel
#'
#' @param tag.sst is data frame containing tag-collected SST data. Requires at least cols: Date (POSIXct), Temperature
#'@param filename is the first part of the filename specified to the download 
#'  function \code{\link{get.env}}. For example, if downloaded files were 
#'  specific to a particular dataset, you may want to identify that with a name 
#'  like 'tuna' or 'shark1'. This results in a downloaded filename of, for 
#'  example, 'tuna_date.nc'. This filename is required here so the calc function
#'  knows where to get the env data.
#'@param sst.dir local directory where remote sensing SST downloads are stored
#'@param dateVec is vector of POSIXct dates for each time step of the likelihood
#'@param focalDim is integer for dimensions of raster::focal used to calculate 
#'  sd() of env grid cell. If left to NULL (default), this dimension will 
#'  approximately incorporate 0.25 degrees.
#'@param sens.err is numeric indicating the percent sensor error in the tag sst
#'  sensor. Default is 1.
#'@param auto.aggr is logical indicating whether to automatically aggregate the calculated likelihood raster if the resolution is higher than 0.1.
#'@param ncores is integer indicating number of cores used in this parallel
#'  computation. Defaults to using a detection function that chooses cores for
#'  you.
#'  
#'@return likelihood is raster brick of likelihood surfaces representing matches
#'  between tag-based sst and remotely sensed sst maps
#'  
#'@export
#'@importFrom foreach "%dopar%"
#'@seealso \code{\link{calc.sst}}
#'  

calc.sst.par <- function(tag.sst, filename, sst.dir, dateVec, focalDim = NULL, sens.err = 1, auto.aggr = TRUE, ncores = NULL){
  
  print(paste('Starting SST likelihood calculation...'))
  
  if (is.null(ncores)) ncores <- ceiling(parallel::detectCores() * .9)
  if (is.na(ncores) | ncores < 0) ncores <- ceiling(as.numeric(system('nproc', intern=T)) * .9)
  
  t0 <- Sys.time()
  
  tag.sst <- tag.sst[which(tag.sst$Date <= max(dateVec)),]
  tag.sst$dateVec <- findInterval(tag.sst$Date, dateVec)
  tag.sst <- data.frame(tag.sst %>% group_by(dateVec) %>% 
                          dplyr::summarise(minT = min(Temperature, na.rm=T), maxT = max(Temperature, na.rm=T), .groups = 'drop_last'))
  T <- length(dateVec)

  # GET EVERYTHING SETUP BEFORE PARALLEL
  nc <- RNetCDF::open.nc(paste(sst.dir, filename, '_', as.Date(dateVec[min(tag.sst$dateVec)]), '.nc', sep=''))
  ncnames = NULL
  nmax <- RNetCDF::file.inq.nc(nc)$nvars - 1
  for(ii in 0:nmax) ncnames[ii + 1] <- RNetCDF::var.inq.nc(nc, ii)$name
  nameidx <- grep('sst', ncnames, ignore.case=TRUE) - 1
  
  lon <- c(RNetCDF::var.get.nc(nc, 'longitude'))
  lat <- c(RNetCDF::var.get.nc(nc, 'latitude'))
  dat <- RNetCDF::var.get.nc(nc, nameidx)
  
  r <- terra::flip(terra::rast(t(dat), extent = terra::ext(min(lon), max(lon), min(lat), max(lat)), crs = "EPSG:4326"), direction = "vertical")
  
  if (is.null(focalDim)){
    focalDim <- round(0.25 / terra::res(r)[1], 0)
    if (focalDim %% 2 == 0) focalDim <- focalDim - 1
    if (focalDim == 1) focalDim <- 3
  }
  
  if(auto.aggr & round(terra::res(r)[1], 2) < 0.1){
    print('Raster is very high resolution. Automatically coarsening using terra::aggregate. If you do not want this behavior, set auto.aggr = FALSE')
    aggFact <- round(0.1 / round(terra::res(r)[1], 2), 0)
    if(aggFact > 1) r <- terra::aggregate(r, fact = aggFact)
  }
  
  tmpl_nr  <- terra::nrow(r)
  tmpl_nc  <- terra::ncol(r)
  tmpl_ext <- as.vector(terra::ext(r))
  L_arr <- array(0, dim = c(tmpl_nr, tmpl_nc, T))
  
  print('Processing in parallel... ')
  
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
  ans = foreach::foreach(i = 1:T, .export = 'likint3', .packages = c('terra')) %dopar%{
    tag.sst.i <- tag.sst[which(tag.sst$dateVec == i),]
    if (nrow(tag.sst.i) == 0) return(NULL)
    sst.i <- c(tag.sst.i$minT * (1 - sens.err / 100), tag.sst.i$maxT * (1 + sens.err / 100))
    
    nc <- RNetCDF::open.nc(paste(sst.dir, filename, '_', as.Date(dateVec[i]), '.nc', sep=''))
    dat <- RNetCDF::var.get.nc(nc, nameidx)
    
    r <- terra::flip(terra::rast(t(dat), extent = terra::ext(min(lon), max(lon), min(lat), max(lat)), crs = "EPSG:4326"), direction = "vertical")
    if(auto.aggr & round(terra::res(r)[1], 2) < 0.1){
      aggFact <- round(0.1 / round(terra::res(r)[1], 2), 0)
      if(aggFact > 1) r <- terra::aggregate(r, fact = aggFact)
    }
    
    sd_r <- terra::focal(r, w = focalDim, fun = "sd", na.rm = TRUE)
    dat  <- terra::as.matrix(r, wide = TRUE)
    sdx  <- terra::as.matrix(sd_r, wide = TRUE)
    
    lik.sst <- likint3(dat, sdx, sst.i[1], sst.i[2])
    lik.sst
  }
  
  parallel::stopCluster(cl)
  on.exit(NULL)
  
  for (i in 1:T){
    if (is.matrix(ans[[i]])) {
      m <- suppressWarnings(max(ans[[i]], na.rm = TRUE))
      if (is.finite(m) && m > 0) {
        slice <- ans[[i]] / m
        slice[is.na(slice) | slice < 0] <- 0
        L_arr[,,i] <- slice
      }
    }
  }
  
  print(paste('Making final likelihood raster...'))
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  L.sst <- terra::rast(L_arr, extent = terra::ext(tmpl_ext), crs = crs)
  names(L.sst) <- as.character(dateVec)
  
  t1 <- Sys.time()
  print(paste('SST calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  # return sst likelihood surfaces
  return(L.sst)
  
}
