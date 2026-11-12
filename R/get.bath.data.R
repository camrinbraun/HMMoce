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
    
    ex180 <- terra::ext(-180, 180, -90, 90)
    ex360 <- terra::ext(180, 360, -90, 90)
    sp_ext <- terra::ext(spatLim$lonmin, spatLim$lonmax, spatLim$latmin, spatLim$latmax)
    
    ## get part 1 - whatever is above 180
    ex1 <- terra::intersect(sp_ext, ex360)
    ex1 <- terra::ext(ex1[1] - 360, ex1[2] - 360, ex1[3], ex1[4])
    original_dir <- getwd()
    tdir <- tempdir()
    
    lonlow <- ex1[1]
    lonhigh <- ex1[2]
    latlow <- ex1[3]
    lathigh <- ex1[4]
    
    opt_save <- opt
    opt <- sub("LATLOW", latlow, opt)
    opt <- sub("LATHIGH", lathigh, opt)
    opt <- sub("LONLOW", lonlow, opt)
    opt <- sub("LONHIGH", lonhigh, opt)
    
    cat(opt)
    curl::curl_download(opt, fname1, quiet=FALSE)
    
    ## get part 2 - whatever is below 180, if any
    ex2 <- terra::intersect(sp_ext, ex180)
    if (is.null(ex2)){
      bathy <- irregular_ncToRaster(fname1, varid)
      setwd(original_dir)
      terra::writeCDF(bathy, paste0(save.dir, '/bathy.nc'), overwrite=TRUE, varname = 'topo')
      if (seaonly == T) bathy[bathy >= 0] <- NA
      return(bathy)
    }
    
    lonlow <- ex2[1]
    lonhigh <- ex2[2]
    latlow <- ex2[3]
    lathigh <- ex2[4]
    
    opt <- opt_save
    opt <- sub("LATLOW", latlow, opt)
    opt <- sub("LATHIGH", lathigh, opt)
    opt <- sub("LONLOW", lonlow, opt)
    opt <- sub("LONHIGH", lonhigh, opt)
    
    cat(opt)
    curl::curl_download(opt, fname2, quiet=FALSE)
    
    print('merging those files to a single output file')
    r1 <- irregular_ncToRaster(fname1, varid)
    r2 <- irregular_ncToRaster(fname2, varid)
    r1r <- terra::shift(terra::rotate(terra::shift(r1, dx = 180)), dx = 180)
    
    setwd(original_dir)
    
    e <- terra::ext(spatLim$lonmin, spatLim$lonmax, spatLim$latmin, spatLim$latmax)
    template <- terra::rast(e, ncols = ceiling((e[2] - e[1]) / terra::res(r1)[1]),
                            nrows = ceiling((e[4] - e[3]) / terra::res(r1)[1]),
                            crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    terra::crs(r1) <- terra::crs(template)
    terra::crs(r2) <- terra::crs(template)
    r1p <- try(terra::resample(r1, template), TRUE)
    if (inherits(r1p, 'try-error')) r1p <- terra::resample(r1r, template)
    r2p <- terra::resample(r2, template)
    r3 <- terra::merge(r1p, r2p)
    
    terra::writeCDF(r3, paste0(save.dir, '/bathy.nc'), overwrite=TRUE, varname = 'topo')
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
    
    bathy <- irregular_ncToRaster(fname, varid)
    
    setwd(original_dir)
    terra::writeCDF(bathy, paste0(save.dir, '/bathy.nc'), overwrite=TRUE, varname = 'topo')
    if (file.exists(paste0(save.dir, '/bathy.nc'))) print(paste0('File output to ', paste0(save.dir, '/bathy.nc')))
    
  }
  
 
  if (seaonly == T) bathy[bathy >= 0] <- NA
  
  return(bathy)
  
}

