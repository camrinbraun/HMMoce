calc.dist2shore <- function(bathy_temp, auto.aggr = TRUE){
  
  ## silence the PROJ warnings
  default_warn <- getOption("warn") 
  options(warn = -1) 
  
  newproj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  ex <- raster::extent(bathy_temp)
  proj4string(bathy_temp) <- newproj
  #bathy_temp <- raster::rotate(bathy_temp)
  shore <- rasterToContour(bathy_temp, levels=-1)
  proj4string(shore) <- proj4string(bathy_temp)
  shore <- raster::crop(shore, ex)
  
  r1 <- raster(ex, resolution=c(.05))
  dd <- gDistance(shore, as(r1,"SpatialPoints"), byid=TRUE)
  r1[] = apply(dd, 1, min)
  
  bathy_temp[bathy_temp >= 0] <- NA
  bathy_temp <- raster::resample(bathy_temp, r1)
  r2 <- raster::mask(r1, bathy_temp)
  
  ## un-silence warnings
  options(warn = default_warn)
  
  
  return(r2)
}
