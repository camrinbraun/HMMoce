#' Create coarse overall likelihood for efficient parameter estimation
#' 
#' \code{coarse.L} makes a more coarse (aggregated) version of the overall likelihood array (output from \code{make.L}). This results in a coarse version of the overall likelihood grid used for more efficient parameter estimation.
#' @param L is overall likelihood array output from \code{make.L}
#' @param ras.list a list of likelihood rasters as input to \code{make.L}. Can also be just one of these raster layers rather than of class "list".
#' @param aggregate_fact is \code{fact} input to \code{raster::aggregate}. 
#' @export

coarse.L <- function(L, ras.list, aggregate_fact = NULL){
  r <- if (is.list(ras.list)) ras.list[[1]] else ras.list
  if (!inherits(r, "SpatRaster")) r <- terra::rast(r)
  
  ## convert L (T, ncol, nrow South-up) back to North-up SpatRaster (nrow, ncol, T)
  L_arr <- aperm(L, c(3, 2, 1))
  L_arr <- L_arr[dim(L_arr)[1]:1, , , drop = FALSE]
  L.mle <- terra::rast(L_arr, extent = terra::ext(r), crs = terra::crs(r))
  
  if (is.null(aggregate_fact)) aggregate_fact <- 5
  L.mle <- terra::aggregate(L.mle, fact = aggregate_fact, fun = "max", na.rm = TRUE)
  g.mle <- setup.grid.raster(L.mle)
  
  arr_mle <- terra::as.array(L.mle)
  L.mle <- aperm(arr_mle[dim(arr_mle)[1]:1, , , drop = FALSE], c(3, 2, 1))
  
  return(list(L.mle = L.mle, g.mle = g.mle))
}


