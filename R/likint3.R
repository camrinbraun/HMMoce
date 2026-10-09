#' Calculate density function and integrate between limits
#' 
#' #' \code{likint3} calculates density function for a normal distribution and integrates between limits
#'  
#' @param w is an array of grid values like sea surface temperature
#' @param wsd is an array containing the sd of the grid values, usually from \code{raster::focal}
#' @param minT is an integer representing the lower limit of the tag-measured variable (e.g. SST)
#' @param maxT is an integer representing the upper limit of the tag-measured variable
#' 
#' @return an array of dim(w) that represents the likelihood of the tag-measured variable as compared to the input grid
#' @export 
#' @examples 
#' # Dummy example environmental grid
#' env.grid <- matrix(seq(1,100,by=1), ncol=10)
#'
#' # Generates likelihood of measurement
#' # as compared to environmental grid
#' likint3(env.grid, 2, 55, 70)

likint3 <- function(w, wsd, minT, maxT){
  #lwr <- minT - .75 * (maxT - minT)
  #upr <- maxT + .75 * (maxT - minT)
  #widx = w >= minT-1 & w <= maxT+1 & !is.na(w)
  #widx = w >= minT-1 & w <= maxT+1 & 
  widx <- !is.na(w)
  wdf = data.frame(w = as.vector(w[widx]), wsd = as.vector(wsd[widx]))
  wdf$wsd[is.na(wdf$wsd)] = 1e-3
  # wint = apply(wdf, 1, function(x) pracma::integral(dnorm, minT, maxT, mean = x[1], sd = x[2]))
  wint = apply(wdf, 1, function(x) stats::integrate(stats::dnorm, 
                                                    # !!!!!! *2 ????
                                                    # minT, maxT, mean = x[1], sd = x[2]*2 )$value)
                                                    minT, maxT, mean = x[1], sd = x[2])$value)
  w = w * 0
  w[widx] = wint
  w
} 