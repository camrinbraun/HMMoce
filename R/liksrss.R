#' Calculate density function and integrate between limits
#' 
#' #' \code{liksrss} calculates density function for a normal distribution and
#' integrates between limits of SRSS times (dawn/dusk)
#' 
#' @param obs input light observation. See \code{\link{calc.srss}}
#' @param srss is raster of possible SR or SS times within study area
#' @param srsd is raster of SD of SRSS times in study area calculated with \code{raster::focal}
#'   
#' @return an array of SRSS-based likelihoods

liksrss <- function(obs, srss, srsd){
  if(is.na(obs)){
    srssout <- srss * 0
  } else{
    sr_v <- as.vector(terra::values(srss))
    sd_v <- as.vector(terra::values(srsd))
    sd_v[is.na(sd_v)] <- 0
    res <- stats::dnorm(obs, sr_v, sd_v)
    srssout <- srss
    terra::values(srssout) <- res
  }
  srssout
}  