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
      
    srssout = srss
    raster::values(srssout) = 0
    
  } else{
    # midT = (maxT + minT) / 2
    # Tsd = (maxT - minT) / 4
    #d = obs-srss
    # widx = w >= minT & w <= maxT & !is.na(w)
    sdf = data.frame(sr = as.vector(srss), srsd = as.vector(srsd))
    sdf$srsd[is.na(sdf$srsd)] = 0
    # wint = apply(wdf, 1, function(x) pracma::integral(dnorm, minT, maxT, mean = x[1], sd = x[2]))
    # wint = apply(wdf, 1, function(x) integrate(dnorm, x[1]-x[2], x[1]+x[2], mean = midT, sd = Tsd * 2)$value) 
    res = stats::dnorm(obs, sdf$sr, sdf$srsd)
    srssout = srss
    raster::values(srssout) = res
    # w = w * 0
    # w[widx] = wint
    # w
  }
  srssout
}  