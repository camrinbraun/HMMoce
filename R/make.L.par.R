#' Combine individual source likelihoods
#' 
#' \code{make.L} combines individual likelihoods from various data sources (e.g.
#' SST, OHC) to make overall combined likelihoods for each time point
#' 
#' @param ras.list a list of likelihood rasters
#' @param iniloc is data.frame of tag and pop locations. Required columns are 'date' (POSIX), 'lon', 'lat'.
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param known.locs is data frame of known locations containing named columns 
#'   of 'date' (POSIX), 'lon', 'lat'. Default is NULL.
#' @param ncores specify number of cores, or leave blank and use whatever you 
#'   have!
#' @return an overall likelihood array, L
#' 
#' @examples
#' \dontrun{
#' L <- make.L(ras.list, iniloc, dateVec)
#' }
#' @export
#' @importFrom foreach %dopar%

make.L.par <- function(ras.list, iniloc, dateVec, known.locs = NULL, ncores = NULL){
  make.L(ras.list = ras.list, iniloc = iniloc, dateVec = dateVec, known.locs = known.locs)
}

