#' Setup the discrete spatial grid for the HMM
#' 
#' @param grid.ras is a raster for which a grid is desired
#'   
#' @return a list
#' @export
#'   

setup.grid.raster <- function(grid.ras){
  if (!inherits(grid.ras, "SpatRaster")) grid.ras <- terra::rast(grid.ras)
  
  lo <- terra::xFromCol(grid.ras, seq_len(terra::ncol(grid.ras)))
  la <- rev(terra::yFromRow(grid.ras, seq_len(terra::nrow(grid.ras))))
  g <- meshgrid(lo, la)
  dlo <- terra::xres(grid.ras)
  dla <- terra::yres(grid.ras)
  
  list(lon = g$X, lat = g$Y, dlo = dlo, dla = dla)
}