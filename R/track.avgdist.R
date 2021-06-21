#' Track Metrics Series by Martin C. Arostegui aka Dr. Mola
#' track.avgdist() calculates the average distance travelled per model time step of a track, 
#' using the coordinates at each time step from the model output of HMMoce

#' The calculation uses Great Circle distance over the surface of the ellipsoid Earth and
#' relies upon the WGS84 coordinate system with longitude [-180,180] and latitude [-90,90]

#' @param track must include paired columns containing "lon" and "lat" coordinates in WGS84 format

#' Output is modified to always yield units of km/*time* where time is the time step duration
#' If daily, the unit is km/d; if at 12 or 6 hr intervals, the unit is km/X hrs

track.avgdist <- function(track){
  # Builds trackline in WGS84 coordinate system
  line = st_sfc(st_linestring(as.matrix(track[,c("lon","lat")])), crs = 4326)
  
  # Avg distance travelled per model time step (km/d or km/X hrs, depending on model formulation)
  avgdist <- round(as.vector(st_length(line) / (length(dateVec)-1) / 1000), digits = 2)
  return(avgdist)
}