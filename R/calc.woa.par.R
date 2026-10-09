#' Calculate WOA profile likelihood in parallel
#' 
#' Calculate Depth-temperature profile based likelihood
#' 
#' \code{calc.woa.par} calculates likelihood of animal position based on 
#' summarized depth-temperature profiles
#' 
#' Tag-based depth-temperature profile summaries are compared to climatological 
#' profiles from the World Ocean Atlas (WOA) "matched" to generate position 
#' likelihoods. This essentially attempts to estimate animal position based on 
#' the water mass it is in, particularly if extensive diving performs thorough 
#' sampling of the environment. However, remember the in situ data is being 
#' compared to climatological means (WOA) or the results of an oceanographic 
#' model (HYCOM).
#' 
#' @param pdt input PDT data output from \code{\link{read.wc}} and 
#'   \code{\link{extract.pdt}}
#' @param ptt is unique tag identifier
#' @param woa.data is (typically) a list of monthly global 1/4deg climatology
#'   data from WOA13. See \code{\link{get.env}}.
#' @param focalDim is integer for dimensions of raster::focal used to calculate 
#'   sd() of temperature grid cell. Recommend focalDim = 3 if woa.data = woa.one
#'   and 9 if using woa.quarter.
#' @param dateVec is vector of dates from tag to pop-up in 1 day increments.
#' @param sp.lim is list of limits as \code{list(xmin, xmax, ymin, ymax)}
#' @param use.se is logical indicating whether or not to use SE when using 
#'   regression to predict temperature at specific depth levels.
#' @param ncores is integer indicating number of cores used in this parallel 
#'   computation. Defaults to using a detection function that chooses cores for 
#'   you.
#'   
#' @export
#' @return raster brick of likelihood
#' @importFrom foreach "%dopar%"
#' @seealso \code{\link{calc.ohc}}

calc.woa.par <- function(pdt, ptt, woa.data = NULL, dateVec, sp.lim = NULL, focalDim = NULL, use.se = TRUE, ncores = NULL){
  
  .Defunct('calc.hycom', package = 'HMMoce')
  
  options(warn=-1)

  t0 <- Sys.time()
  print(paste('Starting WOA likelihood calculation...'))
  
  if (is.null(ncores)) ncores <- ceiling(parallel::detectCores() * .9)
  if (is.na(ncores) | ncores < 0) ncores <- ceiling(as.numeric(system('nproc', intern=T)) * .9)
  
  if(is.null(woa.data)){
    stop('Error: data must be specified')
  }
  
  if (is.null(focalDim)){
    stop('Error: focalDim must be specified.')
  }
  
  if (!is.null(sp.lim)){
    lon.idx <- c(which.min(abs(woa.data$lon - sp.lim[[1]])):
                   which.min(abs(woa.data$lon - sp.lim[[2]])))
    lat.idx <- c(which.min(abs(woa.data$lat - sp.lim[[3]])):
                   which.min(abs(woa.data$lat - sp.lim[[4]])))
    
    woa.data[[1]] <- woa.data[[1]][lon.idx, lat.idx,,]
    woa.data$lon <- woa.data$lon[lon.idx]
    woa.data$lat <- woa.data$lat[lat.idx]
    
  }
  
  depth <- c(0, seq(2.5, 97.5, by=5), seq(112.5, 487.5, by=25), seq(525, 1475, by=50))
  
  # get unique time points
  dateVec = lubridate::parse_date_time(dateVec, '%Y-%m-%d')
  
  udates <- unique(lubridate::parse_date_time(pdt$Date, orders = '%Y-%m-%d %H%:%M:%S'))
  T <- length(udates)
  
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  print(paste0('Generating WOA profile likelihood for ', udates[1], ' through ', udates[length(udates)]))
  
  L.prof <- array(0, dim = c(length(woa.data$lon), length(woa.data$lat), length(dateVec)))
  
  # BEGIN PARALLEL STUFF
  
  print('Processing in parallel... ')
  
  # ncores <- detectCores() # as input argument
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  
  
ans = foreach::foreach(i = 1:T) %dopar%{
  
  # define time based on tag data
  time <- as.Date(udates[i])
  pdt.i <- pdt[which(pdt$Date == time), ]
  
  #extracts depth from tag data for day i
  y <- pdt.i$Depth[!is.na(pdt.i$Depth)]
  y[y < 0] <- 0
  
  #extract temperature from tag data for day i
  x <- pdt.i$MidTemp[!is.na(pdt.i$Depth)]
  
  # use the which.min
  depIdx = apply(as.data.frame(pdt.i$Depth), 1, FUN = function(x) which.min((x - depth) ^ 2))
  
  # make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
  suppressWarnings(
  fit.low <- locfit::locfit(pdt.i$MinTemp ~ pdt.i$Depth)
  )
  suppressWarnings(
  fit.high <- locfit::locfit(pdt.i$MaxTemp ~ pdt.i$Depth)
  )
  
  n = length(depth[depIdx])
  
  pred.low <- stats::predict(fit.low, newdata = depth[depIdx], se = T, get.data = T)
  pred.high <- stats::predict(fit.high, newdata = depth[depIdx], se = T, get.data = T)
  
  if (use.se){
    # data frame for next step
    df = data.frame(low = pred.low$fit - pred.low$se.fit * sqrt(n),
                    high = pred.high$fit + pred.high$se.fit * sqrt(n),
                    depth = depth[depIdx])
  } else{
    # data frame for next step
    df = data.frame(low = pred.low$fit,# - pred.low$se.fit * sqrt(n),
                    high = pred.high$fit,# + pred.high$se.fit * sqrt(n),
                    depth = depth[depIdx])
  }
  
  pdtMonth <-
    as.numeric(format(as.Date(pdt.i$Date), format = '%m'))[1]
  
  wdat = woa.data[[1]]
  
  dat.i = wdat[, , , pdtMonth] #extract months climatology
  
  # calculate sd using Le Bris neighbor method and focal()
  # sd.i = array(NA, dim = c(dim(dat.i)[1:2], length(depth)))
  
  sd.i = array(NA, dim = c(dim(dat.i)[1:2], length(depIdx)))
  
  for(ii in 1:length(depIdx)){
    r = raster::flip(raster::raster(t(dat.i[,,depIdx[ii]])), 2)
    f1 = raster::focal(r, w = matrix(1, nrow = focalDim, ncol = focalDim), fun = function(x) stats::sd(x, na.rm = T))
    f1 = t(raster::as.matrix(raster::flip(f1, 2)))
    sd.i[,,ii] = f1
  } 
  
  # make index of dates for filling in lik.prof
  didx <- base::match(udates, dateVec)
  
  # setup the likelihood array for each day. Will have length (dim[3]) = n depths
  lik.pdt = array(NA, dim = c(dim(dat.i)[1], dim(dat.i)[2], length(depIdx)))
  
  for (b in 1:length(depIdx)) {
    #calculate the likelihood for each depth level, b
    lik.try <- try(likint3(dat.i[,,depIdx[b]], sd.i[,,b], df[b, 1], df[b, 2]), TRUE)
    
    if(class(lik.try) == 'try-error' & use.se == FALSE){
      df[b,1] <- pred.low$fit[b] - pred.low$se.fit[b] * sqrt(n)
      df[b,2] <- pred.high$fit[b] - pred.high$se.fit[b] * sqrt(n)
      
      lik.try <- try(likint3(dat.i[,,depIdx[b]], sd.i[,,b], df[b, 1], df[b, 2]), TRUE)
      
      if (class(lik.try) == 'try-error'){
        lik.try <- dat.i[,,depIdx[b]] * 0
        warning(paste('Warning: likint3 failed after trying with and without SE prediction of depth-temp profiles. This is most likely a divergent integral for ', time, '...', sep=''))
      }
      
    } else if (class(lik.try) == 'try-error' & use.se == TRUE){
      lik.try <- dat.i[,,depIdx[b]] * 0
      warning(paste('Warning: likint3 failed after trying with and without SE prediction of depth-temp profiles. This is most likely a divergent integral for ', time, '...', sep=''))
    }
    
    lik.pdt[,,b] <- lik.try
    
  }
  
  # control for full surfaces of NA vals
  lik.pdt0 <- lik.pdt
  lik.pdt0[is.na(lik.pdt0)] <- 0
  use.idx <- unique(which(lik.pdt0 != 0, arr.ind=T)[,3])
  lik.pdt[lik.pdt == 1] <- NA
  
  # multiply likelihood across depth levels for each day
  lik.pdt <- apply(lik.pdt[,,use.idx], 1:2, FUN=function(x) prod(x, na.rm=F))
  
}

parallel::stopCluster(cl)

# make index of dates for filling in lik.prof
didx <- match(udates, dateVec)

# lapply to normalize
lik.pdt <- lapply(ans, function(x) (x / max(x, na.rm = T)))

# Fill in the L.prof from the list output
ii = 1
for (i in didx) {
  L.prof[,,i] = lik.pdt[[ii]]
  ii = ii + 1
}

print(paste('Making final likelihood raster...'))

crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
L.ras <- raster::brick(L.prof, xmn = min(woa.data$lon), xmx = max(woa.data$lon), ymn = min(woa.data$lat), ymx = max(woa.data$lat), transpose = T, crs)
L.ras <- raster::flip(L.ras, direction = 'y')

t1 <- Sys.time()
print(paste('WOA calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))

options(warn = 2)

return(L.ras)
  
}
