#' Coerce irregular grid as nc to raster
#' 
#' Main purpose of this function is for converting bathymetry data that is typically acquired as an irregular grid
#' 
#' @param fname is input filename of the .nc file
#' @param varid is character indicating the name of the variable of interest in the input nc variable
#' @examples
#' \dontrun{
#' bathy <- irregular_ncToRaster(paste0(bathy.dir, 'bathy.nc'), varid = 'topo')
#' }
#' @export
#' @return a raster version of the input .nc file coerced to regular grid

irregular_ncToRaster <- function(fname, varid){
  
  nc <- RNetCDF::open.nc(fname)
  lon <- as.numeric(RNetCDF::var.get.nc(nc, variable = "longitude"))
  lat <- as.numeric(RNetCDF::var.get.nc(nc, variable = "latitude"))
  bdata <- RNetCDF::var.get.nc(nc, variable = varid)
  
  if (lat[2] - lat[1] > 0){
    flip = TRUE
  } else{
    lat = lat[order(lat)]
    flip = FALSE
  }
  
  bathy = list(x = lon, y = lat, data = bdata)
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  bathy <- terra::rast(t(bathy$data), extent = terra::ext(min(lon), max(lon), min(lat), max(lat)), crs = crs)
  if (flip) bathy <- terra::flip(bathy, direction = "vertical")
  
  bathy
  
}