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
  widx <- !is.na(w)
  if (!any(widx)) return(w * 0)
  wsd_v <- if (length(wsd) == 1L) rep(wsd, sum(widx)) else as.vector(wsd[widx])
  wsd_v[is.na(wsd_v) | wsd_v <= 0] <- 1e-3
  w_out <- w * 0
  w_out[widx] <- stats::pnorm(maxT, mean = as.vector(w[widx]), sd = wsd_v) -
                 stats::pnorm(minT, mean = as.vector(w[widx]), sd = wsd_v)
  w_out
} 