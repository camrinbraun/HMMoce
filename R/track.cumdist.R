#' Track Metrics Series by Martin C. Arostegui aka Dr. Mola
#' track.cumdist() calculates the cumulative distance travelled over a tag deployment track, 
#' using the coordinates at each time step from the model output of HMMoce

#' The calculation uses Great Circle distance over the surface of the ellipsoid Earth and
#' relies upon the WGS84 coordinate system with longitude [-180,180] and latitude [-90,90]

#' @param track must include paired columns containing "lon" and "lat" coordinates in WGS84 format

#' Output is modified to always yield units of km

track.cumdist <- function(track){
  # Builds trackline in WGS84 coordinate system
  line = st_sfc(st_linestring(as.matrix(track[,c("lon","lat")])), crs = 4326)
  
  # Cumulative distance over track (km)
  cumdist <- round(as.vector(st_length(line)/1000), digits = 2)
  return(cumdist)
}