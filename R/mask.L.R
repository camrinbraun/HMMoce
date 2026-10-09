#' Mask L likelihood
#' 
#' \code{mask.L} masks L likelihood based on a certain extent from previous days
#' filter prediction
#' 
#' L likelihood is masked based on a percent size (in addition to) the extent of
#' the previous days prediction kernel. If that kernel is smaller than the 
#' migratory kernel being used, the migratory kernel size is defaulted to. User 
#' can also specify a minimum bound size.
#' 
#' @param pred.t is prediction at t-1 that is used for 1) the size of the mask 
#'   and 2) multiplying with L for resulting post likelihood surface
#' @param L.t is data based likelihood layer, L[t], for the time step of
#'   interest (usually a day)
#' @param lon vector of longitude values corresponding to dims of the previous 2
#'   layers
#' @param lat vector of latitude values corresponding to dims of the previous 2
#'   layers
#' @param par0 is vector of movement parameter values, likely output from
#'   \code{opt.params}.
#' @param bound.thr is numeric indicating the percent threshold that is added 
#'   and subtracted from the bounding box of the filter output from the 
#'   previous day before masking. Default is .05 (5 percent).
#' @param minBounds is size (in degrees) of the minimum bounding box around the 
#'   previous days filter prediction that L data within that box will be 
#'   included. Outside this box (centered on t-1 filter prediction), L will be 
#'   masked out.
#' @importFrom methods as
#' @return post matrix as a product of the input prediction and L values
#' 

mask.L <- function(pred.t, L.t, lon, lat, par0, bound.thr = .05, minBounds=NULL){
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.pred <- list(x = lon, y = lat, z = pred.t)
  ex <- raster::extent(list.pred)
  r <- raster::flip(raster::raster(t(list.pred$z), xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], crs), 'y')
  r.m <- r
  lwr.thr <- bound.thr * raster::cellStats(r.m, 'max')
  r.m[r.m <= lwr.thr] <- NA
  r.m <- raster::trim(r.m)
  #r.t <- trim(r.m)
  #plot(r.t)
  
  ex <- raster::extent(r.m)
  #p1 <- as(ex, 'SpatialPolygons')  
  x.diff <- .1 * (ex[2] - ex[1])
  y.diff <- .1 * (ex[4] - ex[3])
  ex[1] <- ex[1] - x.diff
  ex[2] <- ex[2] + x.diff
  ex[3] <- ex[3] - y.diff
  ex[4] <- ex[4] + y.diff
  
  if(is.null(minBounds)) minBounds <- par0$migr * raster::res(r)[1]
  x.diff <- ex[2] - ex[1]
  y.diff <- ex[4] - ex[3]
  
  if(x.diff < minBounds | y.diff < minBounds){
    ex <- raster::extent(r.m)
    mid.x <- (ex[2] + ex[1]) / 2
    ex[1] <- mid.x - minBounds / 2
    ex[2] <- mid.x + minBounds / 2
    
    mid.y <- (ex[4] + ex[3]) / 2
    ex[3] <- mid.y - minBounds / 2
    ex[4] <- mid.y + minBounds / 2
    
  }
  
  new.ex <- as(ex, 'SpatialPolygons')  
  
  p.mask <- raster::mask(r, new.ex, updatevalue=1, inverse=T)
  p.mask[p.mask < 1] <- 1e-15
  #r.mask <- p.mask * r
  #plot(r.mask)
  #plot(p1, add=T, col='green')
  #plot(p2, add=T, col='red')
  
  p.mask <- t(raster::as.matrix(raster::flip(p.mask, 'y')))
  
  new.L <- p.mask * L.t
  if(max(new.L, na.rm=T) > 1e-15){
    post <- pred.t * new.L
  } else{
    post <- pred.t
  }
  
  return(post)
}