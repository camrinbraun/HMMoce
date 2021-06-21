#' Track Metrics Series by Martin C. Arostegui aka Dr. Mola
#' track.stepdist() calculates a vector of the distance travelled per time step of a tag deployment track, 
#' using the coordinates at each time step from the model output of HMMoce

#' The calculation uses Great Circle distance over the surface of the ellipsoid Earth and
#' relies upon the WGS84 coordinate system with longitude [-180,180] and latitude [-90,90]

#' @param track must include paired columns containing "lon" and "lat" coordinates in WGS84 format

#' Output is modified to always yield units of km

track.stepdist <- function(track){
  # Empty vector
  stepdist <- vector(length=(nrow(tr)-1))
  
  # Vector of distances travelled between each model time step (km)
  for (i in 1:(nrow(tr)-1)){
    line = st_sfc(st_linestring(as.matrix(tr[i:(i+1),c("lon","lat")])), crs = 4326)
    stepdist[i] <- st_length(line)/1000
  }
  return(stepdist)
}