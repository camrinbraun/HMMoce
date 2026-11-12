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

calc.sst.par.glorys <- function(tag.sst, filename, sst.dir, dateVec, focalDim = NULL, sens.err = 1, auto.aggr = FALSE, ncores = NULL){
  
  print(paste('Starting SST likelihood calculation...'))
  
  if (is.null(ncores)) ncores <- ceiling(parallel::detectCores() * .9)
  if (is.na(ncores) | ncores < 0) ncores <- ceiling(as.numeric(system('nproc', intern=T)) * .9)
  if (class(tag.sst$Date)[1] != class(dateVec)[1]) stop('dateVec and tag.sst$Date both need to be of class POSIXct.')
  
  t0 <- Sys.time()
  
  tag.sst <- tag.sst[which(tag.sst$Date <= max(dateVec)),]
  tag.sst$dateVec <- findInterval(tag.sst$Date, dateVec)
  tag.sst <- data.frame(tag.sst %>% group_by(dateVec) %>% 
                          dplyr::summarise(minT = min(Temperature, na.rm=T), maxT = max(Temperature, na.rm=T), .groups = 'drop_last'))
  T <- length(dateVec)

  # open template nc using terra to get grid dimensions and extent
  fl_nc <- list.files(sst.dir, pattern = "\\.nc$", full.names = TRUE)
  if (length(fl_nc) == 0) stop("Source raster file does not exist in sst.dir.")
  first_nc <- NULL
  for (d_idx in seq_along(dateVec)) {
    cand <- c(
      file.path(sst.dir, paste0(filename, '_', format(dateVec[d_idx], '%Y-%m-%d'), '.nc')),
      file.path(sst.dir, paste0(filename, '_', format(dateVec[d_idx], '%Y%m%d'), '.nc'))
    )
    cand <- cand[file.exists(cand)]
    if (length(cand) > 0) {
      first_nc <- cand[1]
      break
    }
  }
  if (is.null(first_nc)) first_nc <- fl_nc[1]
  r_tmpl   <- terra::rast(first_nc, lyrs = 1)
  tmpl_nr  <- terra::nrow(r_tmpl)
  tmpl_nc  <- terra::ncol(r_tmpl)
  tmpl_ext <- as.vector(terra::ext(r_tmpl))
  if (!any(tmpl_ext[1:2] < 180)) tmpl_ext[1:2] <- tmpl_ext[1:2] - 360
  
  ## deal with focalDim if NULL
  if (is.null(focalDim)){
    focalDim <- round(0.25 / terra::res(r_tmpl)[1], 0)
    if (focalDim %% 2 == 0) focalDim <- focalDim - 1
    if (focalDim == 1) focalDim <- 3
  }
  
  # result will be 3D array (nrow x ncol x T) in terra North-up order
  L_arr <- array(0, dim = c(tmpl_nr, tmpl_nc, T))
  
  print('Processing in parallel... ')
  
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
  ans = foreach::foreach(i = 1:T, .export = 'likint3', .packages = c('terra')) %dopar%{
      
    tag.sst.i <- tag.sst[which(tag.sst$dateVec == i),]
    if (nrow(tag.sst.i) == 0) return(NULL)
    sst.i <- c(tag.sst.i$minT * (1 - sens.err / 100), tag.sst.i$maxT * (1 + sens.err / 100))
    
    glorys_file <- file.path(sst.dir, paste0(filename, '_', format(dateVec[i], '%Y-%m-%d'), '.nc'))
    if (!file.exists(glorys_file)) glorys_file <- file.path(sst.dir, paste0(filename, '_', format(dateVec[i], '%Y%m%d'), '.nc'))
    if (!file.exists(glorys_file)){
      glorys_flist <- list.files(sst.dir)
      m <- grep(format(dateVec[i], '%Y%m%d'), glorys_flist)
      if (length(m) > 0) glorys_file <- file.path(sst.dir, glorys_flist[m[1]])
    }
    if (!file.exists(glorys_file)) return(NULL)
    
    r <- terra::rast(glorys_file, lyrs = 1)
    if (terra::nrow(r) != tmpl_nr || terra::ncol(r) != tmpl_nc) {
      target_g <- terra::rast(nrows = tmpl_nr, ncols = tmpl_nc, extent = terra::ext(tmpl_ext), crs = "EPSG:4326")
      r <- terra::resample(r, target_g)
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
    val <- ans[[i]]
    if (is.matrix(val)) {
      m <- suppressWarnings(max(val, na.rm = TRUE))
      if (is.finite(m) && m > 0) {
        slice <- val / m
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
