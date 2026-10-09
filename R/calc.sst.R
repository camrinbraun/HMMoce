#' Calculate SST-based likelihood
#' 
#' \code{calc.sst} compares tag SST to remotely sensed SST and calculates 
#' likelihoods
#' 
#' @param tag.sst is data frame containing tag-collected SST data. Requires at least cols: Date (POSIXct), Temperature
#' @param filename is the first part of the filename specified to the download 
#'   function \code{\link{get.env}}. For example, if downloaded files were 
#'   specific to a particular dataset, you may want to identify that with a name
#'   like 'tuna' or 'shark1'. This results in a downloaded filename of, for 
#'   example, 'tuna_date.nc'. This filename is required here so the calc 
#'   function knows where to get the env data.
#' @param sst.dir local directory where remote sensing SST downloads are stored
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param focalDim is integer for dimensions of raster::focal used to calculate 
#'   sd() of env grid cell. If left to NULL (default), this dimension will
#'   approximately incorporate 0.25 degrees.
#' @param sens.err is numeric indicating the percent sensor error in the tag sst
#'   sensor. Default is 1.
#' @param auto.aggr is logical indicating whether to automatically aggregate the calculated likelihood raster if the resolution is higher than 0.1.
#'   
#' @return likelihood is raster brick of likelihood surfaces representing
#'   matches between tag-based sst and remotely sensed sst maps
#'   
#' @export
#' 
#' @seealso \code{\link{calc.ohc}}
#'   

calc.sst <- function(tag.sst, filename, sst.dir, dateVec, focalDim = NULL, sens.err = 1, auto.aggr = TRUE){
  
  print(paste('Starting SST likelihood calculation...'))
  t0 <- Sys.time()
  
  if (class(tag.sst$Date)[1] != class(dateVec)[1]) stop('dateVec and tag.sst$Date both need to be of class POSIXct.')
  
  tag.sst <- tag.sst[which(tag.sst$Date <= max(dateVec)),]
  tag.sst$dateVec <- findInterval(tag.sst$Date, dateVec)
  tag.sst <- data.frame(tag.sst %>% group_by(dateVec) %>% 
                          dplyr::summarise(minT = min(Temperature, na.rm=T), maxT = max(Temperature, na.rm=T), .groups = 'drop_last'))
  T <- length(dateVec)
  
  ## setup
  nc1 <- RNetCDF::open.nc(dir(sst.dir, full.names = T)[1])
  ncnames = NULL
  nmax <- RNetCDF::file.inq.nc(nc1)$nvars - 1
  for(ii in 0:nmax) ncnames[ii + 1] <- RNetCDF::var.inq.nc(nc1, ii)$name
  nameidx <- grep('sst', ncnames, ignore.case=TRUE) - 1
  lon <- RNetCDF::var.get.nc(nc1, 'longitude')
  lat <- RNetCDF::var.get.nc(nc1, 'latitude')
  # result will be array of likelihood surfaces
  L.sst <- array(0, dim = c(length(lon), length(lat), T))
  print(paste('Starting iterations through deployment period ', '...'))

  for(i in 1:T){
    
    print(dateVec[i])
  
    # tag data for this time step T
    tag.sst.i <- tag.sst[which(tag.sst$dateVec == i),]
    if (nrow(tag.sst.i) == 0) next
    
    sst.i <- c(tag.sst.i$minT * (1 - sens.err / 100), tag.sst.i$maxT * (1 + sens.err / 100)) # sensor error
    
    # open day's sst data
    nc <- RNetCDF::open.nc(paste(sst.dir, filename, '_', as.Date(dateVec[i]), '.nc', sep=''))
   
    ## get grid
    dat <- RNetCDF::var.get.nc(nc, nameidx) # for OI SST
    # calc sd of SST
    # focal calc on mean temp and write to sd var
    r <- raster::flip(raster::raster(t(dat), xmn=min(lon), xmx=max(lon),
                                    ymn=min(lat), ymx=max(lat)), 2)
    
    if(auto.aggr & round(raster::res(r)[1], 2) < 0.1){
      aggFact <- round(0.1 / round(raster::res(r)[1], 2), 0)
      if(aggFact != 1) r <- raster::aggregate(r, fact = aggFact)
      if (i == 1 | !exists('lik.sst')){
        print('Raster is very high resolution. Automatically coarsening using raster::aggregate. If you do not want this behavior, set auto.aggr = FALSE')
        L.sst <- array(0, dim = c(dim(r)[2], dim(r)[1], T))
      }
    } 
    
    if (is.null(focalDim)){
      focalDim <- round(0.25 / raster::res(r)[1], 0)
      if (focalDim %% 2 == 0) focalDim <- focalDim - 1
      if (focalDim == 1) focalDim <- 3 ## if this equals 1, weird things happen
    }

    sdx = raster::focal(r, w = matrix(1, nrow = focalDim, ncol = focalDim), fun = function(x) stats::sd(x, na.rm = T), pad = TRUE)
    sdx = t(raster::as.matrix(raster::flip(sdx, 2)))
    dat <- base::t(raster::as.matrix(raster::flip(r, 2)))

    # compare sst to that day's tag-based ohc
    lik.sst <- likint3(dat, sdx, sst.i[1], sst.i[2])
    #idx <- which(dateVec == time)
    L.sst[,,i] = lik.sst / max(lik.sst, na.rm=T)
    
  }
    
  print(paste('Making final likelihood raster...'))
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  lon.agg <- seq(raster::extent(r)[1], raster::extent(r)[2], length.out=dim(r)[2])
  lat.agg <- seq(raster::extent(r)[3], raster::extent(r)[4], length.out=dim(r)[1])
  list.sst <- list(x = lon.agg, y = lat.agg, z = L.sst)
  ex <- raster::extent(list.sst)
  L.sst <- raster::brick(list.sst$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
  L.sst <- raster::flip(L.sst, direction = 'y')
    
  L.sst[L.sst < 0] <- 0
  
  t1 <- Sys.time()
  print(paste('SST calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  # return sst likelihood surfaces
  return(L.sst)
    
}
  