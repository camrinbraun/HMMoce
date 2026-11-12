#' Calculate Bathymetry-based likelihood in parallel
#'
#' \code{calc.bathy.par} compares tag max depth to bathymetry grid and calculates
#' likelihoods in parallel
#'
#' @param mmd is dataframe containing at least cols: "Date" (POSIXct) and "MaxDepth" in meters (positive values are wet)
#' @param bathy.grid is raster of bathymetry. If the minimum bathymetry values are < 0, the function automatically converts in-water values to positive and masks land.
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param focalDim is integer for dimensions of raster::focal used to calculate
#'   sd() of env grid cell. If left to NULL (default), this dimension will
#'   approximately incorporate 0.25 degrees.
#' @param sens.err is numeric indicating the percent sensor error in the tag-measured max depth. This allows some uncertainty when calculating the integral for the likelihood and doesnt have to necessarily reflect the actual sensor error. Default is 5 percent.
#' @param lik.type is character indicating which likelihood type to use in the bathymetry calculation. Options are dnorm (a traditional likelihood bounded by tag measurement +/- sens.err; experimental) or max ("one-sided" likelihood >= tag-measured max depth; DEFAULT). The latter choice acts more like a mask in that it doesnt allow likelihoods in water shallower than the max depth for each time point.
#' @param ncores is integer indicating number of cores used in this parallel
#'  computation. Defaults to using a detection function that chooses cores for
#'  you.
#'
#' @return likelihood is raster brick of likelihood surfaces representing
#' matches between tag-based depth and a bathymetry field
#'
#' @export
#' @importFrom foreach "%dopar%"
#' @author Originally written by Paul Gatti

calc.bathy.par <- function(mmd, bathy.grid, dateVec, focalDim = NULL, sens.err = 5, lik.type = 'max', ncores = NULL){
  calc.bathy(mmd = mmd, bathy.grid = bathy.grid, dateVec = dateVec, focalDim = focalDim, sens.err = sens.err, lik.type = lik.type)
}
