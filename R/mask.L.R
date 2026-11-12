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
  lwr.thr <- bound.thr * max(pred.t, na.rm = TRUE)
  idx <- which(pred.t > lwr.thr, arr.ind = TRUE)
  if (nrow(idx) == 0) return(pred.t)
  
  ex <- c(min(lon[idx[,1]]), max(lon[idx[,1]]), min(lat[idx[,2]]), max(lat[idx[,2]]))
  x.diff <- .1 * (ex[2] - ex[1])
  y.diff <- .1 * (ex[4] - ex[3])
  ex <- c(ex[1] - x.diff, ex[2] + x.diff, ex[3] - y.diff, ex[4] + y.diff)
  
  res_x <- if (length(lon) > 1) abs(lon[2] - lon[1]) else 0.08
  if(is.null(minBounds)) minBounds <- (if (is.list(par0)) par0$migr else par0[1]) * res_x
  if((ex[2] - ex[1]) < minBounds || (ex[4] - ex[3]) < minBounds){
    mid.x <- mean(range(lon[idx[,1]]))
    mid.y <- mean(range(lat[idx[,2]]))
    ex <- c(mid.x - minBounds / 2, mid.x + minBounds / 2, mid.y - minBounds / 2, mid.y + minBounds / 2)
  }
  
  in_x <- lon >= ex[1] & lon <= ex[2]
  in_y <- lat >= ex[3] & lat <= ex[4]
  p.mask <- outer(in_x, in_y, function(a, b) ifelse(a & b, 1, 1e-15))
  
  new.L <- p.mask * L.t
  if(max(new.L, na.rm=T) > 1e-15){
    post <- pred.t * new.L
  } else{
    post <- pred.t
  }
  
  return(post)
}