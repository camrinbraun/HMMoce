#' glorys Profile LIkelihood in Parallel
#' 
#' Calculate glorys profile likelihood surface in parallel
#' 
#' @param pdt input depth-temperature profile data. Need at least cols: Date (POSIXct), Depth, MinTemp, MaxTemp 
#' @param filename is the first part of the filename specified to the download 
#'   function \code{\link{get.env}}. For example, if downloaded files were 
#'   specific to a particular dataset, you may want to identify that with a name
#'   like 'tuna' or 'shark1'. This results in a downloaded filename of, for 
#'   example, 'tuna_date.nc'. This filename is required here so the calc 
#'   function knows where to get the env data.
#' @param glorys.dir directory of downloaded glorys (or other) data
#' @param focalDim is integer for dimensions of raster::focal used to calculate 
#'   sd() of temperature grid cell. Recommend focalDim = 9 for glorys data at 
#'   0.08deg resolution.
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param use.se is logical indicating whether or not to use SE when using 
#'   regression to predict temperature at specific depth levels.
#' @param flip_y is logical indicating whether or not to flip the resulting likelihood in the y. Set this to true if output likelihoods are upside down.
#'   
#' @return a raster brick of glorys profile likelihood
#' @export
#' @importFrom foreach %dopar%
#'

calc.glorys <- function(pdt, filename, glorys.dir, focalDim = 9, dateVec, depthVec = NULL, use.se = TRUE, flip_y = FALSE){
  
  names(pdt) <- tolower(names(pdt))
  #options(warn=-1)
  if (class(pdt$date)[1] != class(dateVec)[1]) stop('dateVec and pdt$date both need to be of class POSIXct.')
  
  t0 <- Sys.time()
  print(paste('Starting glorys profile likelihood calculation...'))

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
  
  print(paste0('Generating profile likelihood for ', dateVec[1], ' through ', dateVec[length(dateVec)]))
  
  # open template nc using terra to get grid dimensions and extent
  first_nc <- NULL
  for (d_idx in seq_along(dateVec)) {
    cand <- c(
      file.path(glorys.dir, paste0(filename, '_', format(dateVec[d_idx], '%Y-%m-%d'), '.nc')),
      file.path(glorys.dir, paste0(filename, '_', format(dateVec[d_idx], '%Y%m%d'), '.nc'))
    )
    cand <- cand[file.exists(cand)]
    if (length(cand) > 0) {
      first_nc <- cand[1]
      break
    }
  }
  if (is.null(first_nc)) first_nc <- dir(glorys.dir, pattern = "\\.nc$", full.names = TRUE)[1]
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
  
  for(i in 1:T){
    print(dateVec[i])
    
    pdt.i <- pdt[which(pdt$dateVec == i),]
    if (nrow(pdt.i) < 3) next
    
    depIdx <- sort(unique(vapply(pdt.i$depth, function(x) which.min((x - depth)^2), integer(1))))
    if (length(depIdx) < 2) next
    glorysDep <- depth[depIdx]
    
    glorys_file <- file.path(glorys.dir, paste0(filename, '_', format(dateVec[i], '%Y-%m-%d'), '.nc'))
    if (!file.exists(glorys_file)) glorys_file <- file.path(glorys.dir, paste0(filename, '_', format(dateVec[i], '%Y%m%d'), '.nc'))
    if (!file.exists(glorys_file)){
      glorys_flist <- list.files(glorys.dir)
      m <- grep(format(dateVec[i], '%Y%m%d'), glorys_flist)
      if (length(m) > 0) glorys_file <- file.path(glorys.dir, glorys_flist[m[1]])
    }
    if (!file.exists(glorys_file)) next
    
    r_sub <- terra::rast(glorys_file, lyrs = depIdx)
    if (terra::nrow(r_sub) != tmpl_nr || terra::ncol(r_sub) != tmpl_nc) {
      target_g <- terra::rast(nrows = tmpl_nr, ncols = tmpl_nc, extent = terra::ext(tmpl_ext), crs = "EPSG:4326")
      r_sub <- terra::resample(r_sub, target_g)
    }
    
    sd_sub <- terra::focal(r_sub, w = focalDim, fun = "sd", na.rm = TRUE)
    dat    <- terra::as.array(r_sub)
    sd.i   <- terra::as.array(sd_sub)
    
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
    df$high <- pmax(df$high, df$low + 1.0)
    
    lik.pdt <- array(NA_real_, dim = c(dim(dat)[1], dim(dat)[2], length(depIdx)))
    active_b <- logical(length(depIdx))
    deep_ocean <- is.finite(dat[,,length(depIdx)]) & is.finite(sd.i[,,length(depIdx)])
    
    for (b in seq_along(depIdx)) {
      lik.try <- try(likint3(dat[,,b], sd.i[,,b], df$low[b], df$high[b]), TRUE)
      if (inherits(lik.try, 'try-error') || !any(lik.try > 0, na.rm = TRUE)){
        low_b  <- pred.low$fit[b]
        high_b <- max(pred.high$fit[b], low_b + 1.0)
        if (is.finite(pred.low$se.fit[b]) && is.finite(pred.high$se.fit[b])) {
          low_b  <- max(-2, pred.low$fit[b] - pred.low$se.fit[b])
          high_b <- max(low_b + 1.0, pred.high$fit[b] + pred.high$se.fit[b])
        }
        lik.try <- try(likint3(dat[,,b], sd.i[,,b], low_b, high_b), TRUE)
        if (inherits(lik.try, 'try-error') || !any(lik.try > 0, na.rm = TRUE)){
          lik.try <- dat[,,b] * 0
        }
      }
      if (any(lik.try > 0, na.rm = TRUE)) {
        active_b[b] <- TRUE
        lik.try[deep_ocean & (is.na(lik.try) | lik.try < 1e-3)] <- 1e-3
        lik.try[!deep_ocean] <- NA_real_
      }
      lik.pdt[,,b] <- lik.try
    }
    
    use.idx <- which(active_b)
    if (length(use.idx) == 0) next
    
    log_sum <- matrix(0, nrow = dim(dat)[1], ncol = dim(dat)[2])
    n_valid <- matrix(0L, nrow = dim(dat)[1], ncol = dim(dat)[2])
    for (b in use.idx) {
      lb <- lik.pdt[,,b]
      ok <- is.finite(lb) & lb > 0
      log_sum[ok] <- log_sum[ok] + log(lb[ok])
      n_valid[ok] <- n_valid[ok] + 1L
    }
    out <- matrix(NA_real_, nrow = dim(dat)[1], ncol = dim(dat)[2])
    has_obs <- n_valid > 0L
    out[has_obs] <- exp(log_sum[has_obs] / n_valid[has_obs])
    m <- suppressWarnings(max(out, na.rm = TRUE))
    if (is.finite(m) && m > 0) {
      slice <- out / m
      slice[is.na(slice) | slice < 0] <- 0
      L_arr[,,i] <- slice
    }
  }
  
  print(paste('Making final likelihood raster...'))
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  L.glorys <- terra::rast(L_arr, extent = terra::ext(tmpl_ext), crs = crs)
  if (flip_y){
    L.glorys <- terra::flip(L.glorys, direction = "vertical")
    warning('Output raster is being flipped in the y. If this is not desired, use flip_y=FALSE.')
  }
  names(L.glorys) <- as.character(dateVec)
  
  t1 <- Sys.time()
  print(paste('glorys profile calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  #options(warn=2)
  
  # return glorys likelihood surfaces
  return(L.glorys)
  
}
