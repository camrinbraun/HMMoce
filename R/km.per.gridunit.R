#' Input Speed Parameter Guidance by Martin C. Arostegui aka Dr. Mola
#' km.per.gridunit() calculates the min and max longitudinal width (in km) of a single grid unit
#' from the spatial limits of the grid used to model the track ('limits') and grid resolution ('res')

#' The calculation uses Great Circle distance over the surface of the ellipsoid Earth and
#' relies upon the WGS84 coordinate system with longitude [-180,180] and latitude [-90,90]

#' If the spatial limits of the grid include coordinates in both North & South hemispheres,
#' the max longitudinal width occurs at the equator and the min at the farthest lat from the equator

#' If the spatial limits of the grid include coordinates in only one hemisphere (North or South),
#' the max longitudinal width occurs at the nearest lat to the equator and 
#' the min at the farthest lat from the equator

#' @param limits is designed to accept the 'sp.lim' list generated when defining spatial limits,
#' which includes columns; lonmin, lonmax, latmin, latmax

#' @param res must be given in decimal degrees (e.g., 0.08 for hycom)

#' Output is modified to always yield units of km

km.per.gridunit <- function(limits, res){
  
  if (limits$latmax > 0 & limits$latmin < 0){
    # Min and max latitude of spatial grid
    lats <- c(limits$latmax,limits$latmin)
    # Identifies latitude farthest from equator
    lat.farthest <- lats[which.max(abs(lats))]
    # Matrix of grid resolution longitude differences at farthest lat from equator and at the equator
    # Used when the spatial grid occurs in both hemispheres
    londist <- matrix(data = c(limits$lonmin, lat.farthest,
                               limits$lonmin + res, lat.farthest,
                               limits$lonmin, 0,
                               limits$lonmin + res, 0),
                      nrow=4,byrow=T)
  } else {
    # Matrix of grid resolution longitude differences at max and min lat 
    # Used when the spatial grid occurs in only one hemisphere
    londist <- matrix(data = c(limits$lonmin, limits$latmax,
                               limits$lonmin + res, limits$latmax,
                               limits$lonmin, limits$latmin,
                               limits$lonmin + res, limits$latmin),
                      nrow=4,byrow=T)
  }
  
  # Empty vector
  rangedist <- vector(length=2)
  
  # Vector of grid resolution longitude differences
   rangedist[1] <- st_length(st_sfc(st_linestring(londist[1:2,]), crs = 4326))/1000
   rangedist[2] <- st_length(st_sfc(st_linestring(londist[3:4,]), crs = 4326))/1000
   
  return(sort(round(rangedist,digits=2)))
}
