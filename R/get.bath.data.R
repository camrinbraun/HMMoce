#' Download bathymetry data
#' 
#' Download ETOPO bathymetry. Resolution is either 30 second or one minute 
#' resolution
#' 
#' @param spatLim is a list of spatial limits as \code{list(lonmin, lonmax, latmin, latmax)}
#' @param save.dir is destination folder. Default is a temporary directory.
#' @param seaonly is logical indicating whether you want to download only
#'   bathymetry below sea level.
#' @param res is numeric indicating resolution in minutes. Choices currently are 0.5 or 1minute.
#' @return Downloads a NetCDF file containing ETopo bathymetry and also returns the output as a \code{raster} object.
#' @author Function originally written by Ben Galuardi.
#' @examples
#' \dontrun{
#' # Not run to prevent actual data download
#' sp.lim <- list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50)
#' bathy <- get.bath.data(sp.lim, save.dir = tempdir())
#' }
#' @export
#' @importFrom curl curl_download
#' @note Be patient! The download can take a few minutes!

get.bath.data <- function(spatLim, save.dir = tempdir(), seaonly = TRUE, res = c(.5)){
  
  if (res == 1){
    cat('ERDDAP downloading: Topography, Smith & Sandwell v11.1, 1/60-degree \n UCSD   (Dataset ID: usgsCeSS111)')
    opt = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSS91.nc?topo[(LATHIGH):(LATLOW)][(LONLOW):(LONHIGH)]"
    varid = 'topo'
  } else if(res == 0.5){
    cat('ERDDAP downloading: Topography, SRTM30+ Version 1.0, 30 arc second, Global \n 	Scripps   (Dataset ID: usgsCeSrtm30v1)')
    #opt ="https://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSrtm30v6.nc?topo[(LATHIGH):(LATLOW)][(LONLOW):(LONHIGH)]"
    opt = 'http://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.nc?altitude[(LATLOW):1:(LATHIGH)][(LONLOW):1:(LONHIGH)]'
    varid = 'altitude'
  }
  
  if (spatLim$lonmax > 180){ ## 0 to 360
    print('Detected input coordinates > 180, downloading multiple files.')
    tdir <- tempdir()
    
    fname1 = paste(tdir, "bathy1.nc", sep = "/")
    fname2 = paste(tdir, "bathy2.nc", sep = "/")
    
    ex180 <- raster::extent(-180, 180, -90, 90)
    ex360 <- raster::extent(180, 360, -90, 90)
    
    ## get part 1 - whatever is above 180
    ex1 <- raster::intersect(raster::extent(unlist(spatLim)), ex360)
    ex1 <- raster::extent(raster::rotate(raster::raster(ex1)))
    #if (ex1@xmin == -180) ex1@xmin <- -179.99
    original_dir <- getwd()
    tdir <- tempdir()
    
    lonlow <- ex1@xmin
    lonhigh <- ex1@xmax
    latlow <- ex1@ymin
    lathigh <- ex1@ymax
    
    opt_save <- opt
    #opt <- opt_save
    opt <- sub("LATLOW", latlow, opt)
    opt <- sub("LATHIGH", lathigh, opt)
    opt <- sub("LONLOW", lonlow, opt)
    opt <- sub("LONHIGH", lonhigh, opt)
    
    cat(opt)
    curl::curl_download(opt, fname1, quiet=FALSE)
    #utils::download.file(opt, fname)
    
    
    ## get part 2 - whatever is below 180, if any
    ex2 <- raster::intersect(raster::extent(unlist(spatLim)), ex180)
    if (is.null(ex2)){
      ## the first download got the full extent 
      
      bathy <- irregular_ncToRaster(fname1, varid)
      
      setwd(original_dir)
      raster::writeRaster(bathy, paste0(save.dir, '/bathy.nc'), format='CDF', overwrite=TRUE, varname = 'topo')
      if (seaonly == T) bathy[bathy >= 0] <- NA
      return(bathy)
      
    }
    #if (ex2@xmax == 180) ex2@xmax <- 179.995
    
    lonlow <- ex2@xmin
    lonhigh <- ex2@xmax
    latlow <- ex2@ymin
    lathigh <- ex2@ymax
    
    opt <- opt_save
    opt <- sub("LATLOW", latlow, opt)
    opt <- sub("LATHIGH", lathigh, opt)
    opt <- sub("LONLOW", lonlow, opt)
    opt <- sub("LONHIGH", lonhigh, opt)
    
    cat(opt)
    curl::curl_download(opt, fname2, quiet=FALSE)
    #utils::download.file(opt, fname)
    
    
    print('merging those files to a single output file')
    r1 <- irregular_ncToRaster(fname1, varid)
    r2 <- irregular_ncToRaster(fname2, varid)
    r1r <- raster::shift(raster::rotate(raster::shift(r1, 180)), 180)
    
    ## reset original directory after messing w temp files
    setwd(original_dir)
    
    e <- raster::extent(spatLim$lonmin, spatLim$lonmax, spatLim$latmin, spatLim$latmax)
    template <- raster::raster(e, ncols = ceiling((e@xmax - e@xmin) / raster::res(r1)[1]),
                               nrows = ceiling((e@ymax - e@ymin) / raster::res(r1)[1]))
    raster::projection(template) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    crs(r1) <- crs(template)
    crs(r2) <- crs(template)
    r1p <- try(raster::resample(r1, template), TRUE)
    if (class(r1p) == 'try-error') r1p <- raster::resample(r1r, template)
    r2p <- raster::resample(r2, template)
    r3 <- raster::merge(r1p, r2p)
    
    #print(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))
    #r1r <- raster::shift(raster::rotate(raster::shift(r1, 180)), 180)
    #r3 <- raster::merge(r1r, r2)
    raster::writeRaster(r3, paste0(save.dir, '/bathy.nc'), format='CDF', overwrite=TRUE, varname = 'topo')
    if (file.exists(paste0(save.dir, '/bathy.nc'))) print(paste0('File output to ', paste0(save.dir, '/bathy.nc')))
    bathy <- r3
    
  } else {
    
    original_dir <- getwd()
    tdir <- tempdir()
    fname <- paste(tdir, "bathy.nc", sep = "/")
    
    opt <- sub("LATLOW", spatLim$latmin, opt)
    opt <- sub("LATHIGH", spatLim$latmax, opt)
    opt <- sub("LONLOW", spatLim$lonmin, opt)
    opt <- sub("LONHIGH", spatLim$lonmax, opt)
    
    cat(opt)
    curl::curl_download(opt, fname, quiet=FALSE)
    #utils::download.file(opt, fname)
    
    ## bathy is on irregular grid so needs to be coerced
    bathy <- irregular_ncToRaster(fname, varid)
    
    setwd(original_dir)
    raster::writeRaster(bathy, paste0(save.dir, '/bathy.nc'), format='CDF', overwrite=TRUE, varname = 'topo')
    if (file.exists(paste0(save.dir, '/bathy.nc'))) print(paste0('File output to ', paste0(save.dir, '/bathy.nc')))
    
  }
  
 
  if (seaonly == T) bathy[bathy >= 0] <- NA
  
  return(bathy)
  
}

