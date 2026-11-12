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
#' @param ohc.dir directory of downloaded glorys (or other) data
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param bathy is logical indicating whether or not a bathymetric mask should 
#'   be applied. This does NOT use an input bathymetry dataset but rather masks the depth-temperature likelihood based on the possible locations given the maximum recorded tag-depth (i.e. areas with bottom depth shallower than tag-recorded max depth are masked out).
#' @param use.se is logical indicating whether or not to use SE when using
#'   regression to predict temperature at specific depth levels.
#' @param ncores specify number of cores, or leave blank and use whatever you 
#'   have
#' @param flip_y is logical indicating whether or not to flip the resulting likelihood in the y. Set this to true if output likelihoods are upside down.
#'   
#' @return likelihood is raster brick of likelihood surfaces representing 
#'   estimated position based on tag-based OHC compared to calculated OHC using 
#'   glorys
#'   
#' @export
#' @seealso \code{\link{calc.woa}}
#' @references Luo J, Ault JS, Shay LK, Hoolihan JP, Prince ED, Brown C a.,
#'   Rooker JR (2015) Ocean Heat Content Reveals Secrets of Fish Migrations.
#'   PLoS One 10:e0141101
#'   

calc.ohc.glorys.par <- function(pdt, filename, isotherm = '', ohc.dir, dateVec, depthVec = NULL, bathy = TRUE, use.se = TRUE, flip_y = FALSE, ncores=NULL){

  names(pdt) <- tolower(names(pdt))
  
  if (is.null(ncores)) ncores <- ceiling(parallel::detectCores() * .9)
  if (is.na(ncores) | ncores < 0) ncores <- ceiling(as.numeric(system('nproc', intern=T)) * .9)
  
  if (class(pdt$date)[1] != class(dateVec)[1]) stop('dateVec and pdt$date both need to be of class POSIXct.')
  
  t0 <- Sys.time()
  print(paste('Starting OHC likelihood calculation...'))
  
  # constants for OHC calc
  cp <- 3.993 # kJ/kg*C <- heat capacity of seawater
  rho <- 1025 # kg/m3 <- assumed density of seawater
  
  # calculate midpoint of tag-based min/max temps
  #if(length(grep('mean', names(pdt))) > 1){
  #  pdt$useTemp <- pdt[,grep('mean', names(pdt))[1]]
  #} else if(length(grep('mean', names(pdt))) == 1){
  #  pdt$useTemp <- pdt[,grep('mean', names(pdt))]
  #} else{
  #  pdt$useTemp <- (pdt$maxtemp + pdt$mintemp) / 2
  #}
  
  pdt <- pdt[which(pdt$date <= max(dateVec)),]
  pdt$dateVec <- findInterval(pdt$date, dateVec)
  T <- length(dateVec)
  
  if(isotherm != '') iso.def <- TRUE else iso.def <- FALSE
  
  # open template nc using terra to get grid dimensions and extent
  fl_nc <- list.files(ohc.dir, pattern = "\\.nc$", full.names = TRUE)
  if (length(fl_nc) == 0) stop("Source raster file does not exist in ohc.dir.")
  first_nc <- NULL
  for (d_idx in seq_along(dateVec)) {
    cand <- c(
      file.path(ohc.dir, paste0(filename, '_', format(dateVec[d_idx], '%Y-%m-%d'), '.nc')),
      file.path(ohc.dir, paste0(filename, '_', format(dateVec[d_idx], '%Y%m%d'), '.nc'))
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
  
  ## better handling of depth
  if (is.null(depthVec)){
    depth <- c(0.494025, 1.541375,    2.645669,    3.819495,    5.078224,    6.440614,    7.929560,    9.572997,   11.405000,
               13.467140,   15.810070,   18.495560,   21.598820,   25.211411,   29.444731,   34.434151,   40.344051,   47.373692,
               55.764290,   65.807266,   77.853851,   92.326073,  109.729301,  130.666000,  155.850693,  186.125595,  222.475204,
               266.040314,  318.127411,  380.213013,  453.937714,  541.088928,  643.566772,  763.333130,  902.339294, 1062.439941,
               1245.291016, 1452.250977, 1684.284058, 1941.892944, 2225.077881, 2533.335938, 2865.702881, 3220.820068, 3597.031982,
               3992.483887, 4405.224121, 4833.291016, 5274.784180, 5727.916992)
  } else {
    depth <- depthVec
  }
  
  # result will be 3D array (nrow x ncol x T) in terra North-up order
  L_arr <- array(0, dim = c(tmpl_nr, tmpl_nc, T))
  
  # BEGIN PARALLEL STUFF  
  print('Processing in parallel... ')
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
  ans = foreach::foreach(i = 1:T, .export = 'likint3', .packages = c('terra', 'locfit')) %dopar%{
    
    pdt.i <- pdt[which(pdt$dateVec == i),]
    if (nrow(pdt.i) < 3) return(NULL)
    
    depIdx <- sort(unique(vapply(pdt.i$depth, function(x) which.min((x - depth)^2), integer(1))))
    if (length(depIdx) < 2) return(NULL)
    glorysDep <- depth[depIdx]
    
    glorys_file <- file.path(ohc.dir, paste0(filename, '_', format(dateVec[i], '%Y-%m-%d'), '.nc'))
    if (!file.exists(glorys_file)) glorys_file <- file.path(ohc.dir, paste0(filename, '_', format(dateVec[i], '%Y%m%d'), '.nc'))
    if (!file.exists(glorys_file)){
      glorys_flist <- list.files(ohc.dir)
      m <- grep(format(dateVec[i], '%Y%m%d'), glorys_flist)
      if (length(m) > 0) glorys_file <- file.path(ohc.dir, glorys_flist[m[1]])
    }
    if (!file.exists(glorys_file)) return(NULL)
    
    r_sub <- terra::rast(glorys_file, lyrs = depIdx)
    if (terra::nrow(r_sub) != tmpl_nr || terra::ncol(r_sub) != tmpl_nc) {
      target_g <- terra::rast(nrows = tmpl_nr, ncols = tmpl_nc, extent = terra::ext(tmpl_ext), crs = "EPSG:4326")
      r_sub <- terra::resample(r_sub, target_g)
    }
    dat <- terra::as.array(r_sub)
    deep_ocean <- is.finite(dat[,,length(depIdx)])
    
    suppressWarnings(fit.low  <- locfit::locfit(pdt.i$mintemp ~ pdt.i$depth, maxk = 500))
    suppressWarnings(fit.high <- locfit::locfit(pdt.i$maxtemp ~ pdt.i$depth, maxk = 500))
    
    pred.low  <- stats::predict(fit.low,  newdata = glorysDep, se = TRUE, get.data = TRUE)
    pred.high <- stats::predict(fit.high, newdata = glorysDep, se = TRUE, get.data = TRUE)
    
    if (use.se && all(is.finite(pred.low$se.fit)) && all(is.finite(pred.high$se.fit))){
      df <- data.frame(low = pred.low$fit - pred.low$se.fit,
                       high = pred.high$fit + pred.high$se.fit,
                       depth = glorysDep)
      if (any(!is.finite(df$low)) || any(!is.finite(df$high)) || min(df$low, na.rm = TRUE) < -2) {
        df <- data.frame(low = pred.low$fit, high = pred.high$fit, depth = glorysDep)
      }
    } else{
      df <- data.frame(low = pred.low$fit, high = pred.high$fit, depth = glorysDep)
    }
    df$high <- pmax(df$high, df$low + 0.25)

    if(iso.def == FALSE) isotherm <- max(0, min(df$low, na.rm = TRUE))
    
    minT.ohc <- cp * rho * sum(df$low - isotherm, na.rm = TRUE) / 10000
    maxT.ohc <- cp * rho * sum(df$high - isotherm, na.rm = TRUE) / 10000
    
    ohc <- cp * rho * rowSums(dat - isotherm, na.rm = TRUE, dims = 2) / 10000
    ohc[ohc == 0] <- NA_real_
    if (bathy) ohc[!deep_ocean] <- NA_real_
    
    r_ohc  <- terra::rast(ohc, extent = terra::ext(tmpl_ext), crs = "EPSG:4326")
    sd_ohc <- terra::focal(r_ohc, w = 9, fun = "sd", na.rm = TRUE)
    sdx    <- terra::as.matrix(sd_ohc, wide = TRUE)

    lik.try <- try(likint3(ohc, sdx, minT.ohc, maxT.ohc), TRUE)
    
    if (inherits(lik.try, 'try-error') || max(lik.try, na.rm = TRUE) == 0){
      df <- data.frame(low = pred.low$fit,
                       high = pmax(pred.high$fit, pred.low$fit + 0.25),
                       depth = glorysDep)
      if(iso.def == FALSE) isotherm <- max(0, min(df$low, na.rm = TRUE))
      minT.ohc <- cp * rho * sum(df$low - isotherm, na.rm = TRUE) / 10000
      maxT.ohc <- cp * rho * sum(df$high - isotherm, na.rm = TRUE) / 10000
      ohc <- cp * rho * rowSums(dat - isotherm, na.rm = TRUE, dims = 2) / 10000
      ohc[ohc == 0] <- NA_real_
      if (bathy) ohc[!deep_ocean] <- NA_real_
      r_ohc  <- terra::rast(ohc, extent = terra::ext(tmpl_ext), crs = "EPSG:4326")
      sd_ohc <- terra::focal(r_ohc, w = 9, fun = "sd", na.rm = TRUE)
      sdx    <- terra::as.matrix(sd_ohc, wide = TRUE)
      lik.try <- try(likint3(ohc, sdx, minT.ohc, maxT.ohc), TRUE)
    }
    if (bathy && is.matrix(lik.try)) lik.try[!deep_ocean] <- NA_real_
    lik.try
  }

  parallel::stopCluster(cl)
  on.exit(NULL)
  
  for(i in 1:T){
    ans_i <- ans[[i]]
    if (is.matrix(ans_i)) {
      m <- suppressWarnings(max(ans_i, na.rm = TRUE))
      if (is.finite(m) && m > 0) {
        slice <- ans_i / m
        slice[is.na(slice) | slice < 0] <- 0
        L_arr[,,i] <- slice
      }
    }
  }
  
  print(paste('Making final likelihood raster...'))
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  L.ohc <- terra::rast(L_arr, extent = terra::ext(tmpl_ext), crs = crs)
  if (flip_y){
    L.ohc <- terra::flip(L.ohc, direction = "vertical")
    warning('Output raster is being flipped in the y. If this is not desired, use flip_y=FALSE.')
  }
  names(L.ohc) <- as.character(dateVec)
  
  t1 <- Sys.time()
  print(paste('OHC calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  # return ohc likelihood surfaces
  return(L.ohc)
  
}

