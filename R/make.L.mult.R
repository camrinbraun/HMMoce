#' Combine individual source likelihoods
#' 
#' \code{make.L} combines individual likelihoods from various data sources (e.g. SST, OHC) to make 
#' overall combined likelihoods for each time point. This multiplies (i.e., does not average) likelihoods
#' 
#' @param ras.list a list of likelihood rasters
#' @param iniloc is data.frame of tag and pop locations. Required columns are 'date' (POSIX), 'lon', 'lat'.
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param maxDepth a vector of the daily max depths (must be positive and be zero-padded for NA days)
#' @param bathy is the original bathy raster pre-resampling
#' @param known.locs is data frame of known locations containing named columns 
#'   of 'date' (POSIX), 'lon', 'lat'. Default is NULL.
#' @return an overall likelihood array, L
#' 
#' @examples
#' \dontrun{
#' L <- make.L(ras.list, iniloc, dateVec)
#' }
#' @export
#'   

make.L.mult <- function(ras.list, iniloc, dateVec, maxDepth = NULL, bathy = NULL, known.locs = NULL){
  make.L(ras.list = ras.list, iniloc = iniloc, dateVec = dateVec, maxDepth = maxDepth, bathy = bathy, known.locs = known.locs)
}
