#' Download and Read Oceanographic Data
#' 
#' \code{get.env} accesses oceanographic data like sea surface temperature from 
#' a remote server and downloads the temporal and spatial extent of interest for
#' further use
#' 
#' @param uniqueDates is a POSIXct vector of desired dates
#' @param filename is first part of the filename specified to the download 
#'   function. For example, if downloaded files were specific to a particular 
#'   dataset, you may want to identify that with a name like 'tuna' or 'shark1'.
#'   This results in a downloaded filename of, for example, 'tuna_date.nc'.
#' @param type is a character string indicating whether you're after sea surface
#'   temperature 'sst', hybrid coordinate ocean model 'hycom', or world ocean 
#'   atlas 'woa' data
#' @param spatLim is a list of spatial limits as \code{list(xmin, xmax, ymin, 
#'   ymax)}
#' @param resol is character describing the desired resolution in degrees if 
#'   type = 'woa', otherwise NULL. Choices are 'one' or 'quarter'.
#' @param save.dir is the directory to save the downloaded data to
#' @param sst.type is character indicating type of desired SST product. Choices 
#'   are currently Optimum Interpolation ('oi')
#'   \url{https://www.ncdc.noaa.gov/oisst}, MUR ('mur')\url{https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.html} or a high-resolution composite
#'   ('ghr') \url{https://www.ghrsst.org/}.
#' @param depLevels is an integer describing which depth levels to download from
#'   Hycom (e.g. 1=surface). Default is NULL and all levels are downloaded.
#'   
#' @return nothing, just downloads the data to your local machine
#' @examples
#' \dontrun{
#' # Not run to prevent actual data download
#' sp.lim <- list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50)
#' # FOR OI SST DATA
#' get.env(as.Date('2015-10-01'), filename='oisst', type = 'sst', 
#' sst.type='oi', spatLim = sp.lim, save.dir = tempdir())
#' 
#' # FOR HYCOM DATA
#' get.env(as.Date('2015-10-01'), filename='hycom', type = 'hycom', 
#' spatLim = sp.lim, save.dir = tempdir())
#' 
#' # FOR WORLD OCEAN ATLAS DATA
#' get.env(type = 'woa', resol = 'quarter', save.dir = woa.dir)
#' }
#' @export

get.env <- function(uniqueDates = NULL, filename = NULL, type = NULL, spatLim = NULL, resol = NULL, save.dir = getwd(), sst.type=NULL, depLevels=NULL, template=NULL, ...){
  
  original_dir <- getwd()
  
  if(is.null(type)){
    
    stop('Type of environmental data desired not specified.')
    
  } else if(type == 'sst'){
    if(is.null(sst.type)){
      warning('Warning: if type=sst then sst.type should be specified. Using default GHRsst.')
      sst.type = 'ghr'
    }
    
    if(sst.type == 'oi'){
      for(i in 1:length(uniqueDates)){
        time <- as.Date(uniqueDates[i])
        repeat{
          get.oi.sst(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''), download.file = TRUE, dir = save.dir, ...) # filenames based on dates from above
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(save.dir, filename, '_', time, '.nc', sep = '')), silent = T)
          }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
          if(class(err) != 'try-error') break
        }
      }
      
      ## GHR
    } else if(sst.type == 'ghr'){
      for(i in 1:length(uniqueDates)){
        time <- as.Date(uniqueDates[i])
        
        if (spatLim$lonmax > 180){ ## 0 to 360
          cat('Detected input coordinates > 180, downloading multiple files.','\n')
          ex180 <- raster::extent(-180,180,-90,90)
          ex360 <- raster::extent(180,360,-90,90)
          
          ## get part 1 - whatever is above 180
          ex1 <- raster::intersect(raster::extent(unlist(spatLim)), ex360)
          ex1 <- raster::extent(raster::rotate(raster::raster(ex1)))
          if (ex1@xmin == -180) ex1@xmin <- -179.995
          tdir <- tempdir()
          repeat{
            get.ghr.sst(c(ex1@xmin, ex1@xmax, ex1@ymin, ex1@ymax), time, filename = paste(filename, '_', time, '_1.nc', sep = ''), download.file = TRUE, dir = tdir)#, ...) # filenames based on dates from above
            tryCatch({
              err <- try(RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_1.nc', sep = '')), silent = T)
            }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
            if(class(err) != 'try-error') break
          }
          
          ## get part 2 - whatever is below 180, if any
          ex2 <- raster::intersect(raster::extent(unlist(spatLim)), ex180)
          if (ex2@xmax == 180) ex2@xmax <- 179.995
          
          repeat{
            get.ghr.sst(c(ex2@xmin, ex2@xmax, ex2@ymin, ex2@ymax), time, filename = paste(filename, '_', time, '_2.nc', sep = ''), download.file = TRUE, dir = tdir)#, ...) # filenames based on dates from above
            tryCatch({
              err <- try(RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_2.nc', sep = '')), silent = T)
            }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
            if(class(err) != 'try-error') break
          }
          
          cat('merging those files to a single output file','\n')
          nc1 <- RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_1.nc', sep = ''))
          lon1 <- RNetCDF::var.get.nc(nc1, 'longitude')
          lat1 <- RNetCDF::var.get.nc(nc1, 'latitude')
          dat1 <- RNetCDF::var.get.nc(nc1, 'SST') # for OI SST
          r1 <- raster::flip(raster::raster(t(dat1), xmn=min(lon1), xmx=max(lon1),
                                            ymn=min(lat1), ymx=max(lat1)), 2)
          
          nc2 <- RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_2.nc', sep = ''))
          lon2 <- RNetCDF::var.get.nc(nc2, 'longitude')
          lat2 <- RNetCDF::var.get.nc(nc2, 'latitude')
          dat2 <- RNetCDF::var.get.nc(nc2, 'SST') # for OI SST
          r2 <- raster::flip(raster::raster(t(dat2), xmn=min(lon2), xmx=max(lon2),
                                            ymn=min(lat2), ymx=max(lat2)), 2)
          
          ## reset original directory after messing w temp files
          setwd(original_dir)
          
          #print(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))
          r1r <- raster::shift(raster::rotate(raster::shift(r1, 180)), 180)
          r3 <- raster::merge(r1r, r2)
          raster::writeRaster(r3, paste(save.dir, '/', filename, '_', time, '.nc', sep = ''), format='CDF', varname = 'SST', overwrite=TRUE)
          if (file.exists(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))) cat(paste0('File output to ', paste(save.dir, '/', filename, '_', time, '.nc', sep = '')),'\n')
          
        } else{ ## -180 to 180
          repeat{
            get.ghr.sst(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''), download.file = TRUE, dir = save.dir, ...) # filenames based on dates from above
            tryCatch({
              err <- try(RNetCDF::open.nc(paste(save.dir, filename, '_', time, '.nc', sep = '')), silent = T)
            }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
            if(class(err) != 'try-error') break
          }
        }
        
      }
      
      ## MUR
    } else if(sst.type == 'mur'){
      for(i in 1:length(uniqueDates)){
        time <- as.Date(uniqueDates[i])
        
        if (spatLim$lonmax > 180){ ## 0 to 360
          cat('Detected input coordinates > 180, downloading multiple files.','\n')
          ex180 <- raster::extent(-180,180,-90,90)
          ex360 <- raster::extent(180,360,-90,90)
          
          ## get part 1 - whatever is above 180
          ex1 <- raster::intersect(raster::extent(unlist(spatLim)), ex360)
          ex1 <- raster::extent(raster::rotate(raster::raster(ex1)))
          if (ex1@xmin == -180) ex1@xmin <- -179.99
          original_dir <- getwd()
          tdir <- tempdir()
          repeat{
            get.mur.sst(c(ex1@xmin, ex1@xmax, ex1@ymin, ex1@ymax), time, filename = paste(filename, '_', time, '_1.nc', sep = ''), download.file = TRUE, dir = tdir)#, ...) # filenames based on dates from above
            tryCatch({
              err <- try(RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_1.nc', sep = '')), silent = T)
            }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
            if(class(err) != 'try-error') break
          }
          
          ## get part 2 - whatever is below 180, if any
          ex2 <- raster::intersect(raster::extent(unlist(spatLim)), ex180)
          #if (ex1@xmax == 180) ex1@xmax <- 179.99
          
          repeat{
            get.mur.sst(c(ex2@xmin, ex2@xmax, ex2@ymin, ex2@ymax), time, filename = paste(filename, '_', time, '_2.nc', sep = ''), download.file = TRUE, dir = tdir)#, ...) # filenames based on dates from above
            tryCatch({
              err <- try(RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_2.nc', sep = '')), silent = T)
            }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
            if(class(err) != 'try-error') break
          }
          
          cat('merging those files to a single output file','\n')
          nc1 <- RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_1.nc', sep = ''))
          lon1 <- RNetCDF::var.get.nc(nc1, 'longitude')
          lat1 <- RNetCDF::var.get.nc(nc1, 'latitude')
          dat1 <- RNetCDF::var.get.nc(nc1, 'analysed_sst') # for OI SST
          r1 <- raster::flip(raster::raster(t(dat1), xmn=min(lon1), xmx=max(lon1),
                                           ymn=min(lat1), ymx=max(lat1)), 2)
          
          nc2 <- RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_2.nc', sep = ''))
          lon2 <- RNetCDF::var.get.nc(nc2, 'longitude')
          lat2 <- RNetCDF::var.get.nc(nc2, 'latitude')
          dat2 <- RNetCDF::var.get.nc(nc2, 'analysed_sst') # for OI SST
          r2 <- raster::flip(raster::raster(t(dat2), xmn=min(lon2), xmx=max(lon2),
                                            ymn=min(lat2), ymx=max(lat2)), 2)
          
          ## reset original directory after messing w temp files
          setwd(original_dir)
          
          #print(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))
          r1r <- raster::shift(raster::rotate(raster::shift(r1, 180)), 180)
          r3 <- raster::merge(r1r, r2)
          raster::writeRaster(r3, paste(save.dir, '/', filename, '_', time, '.nc', sep = ''), format='CDF', varname = 'SST', overwrite=TRUE)
          if (file.exists(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))) cat(paste0('File output to ', paste(save.dir, '/', filename, '_', time, '.nc', sep = '')),'\n')
          
        } else{ ## -180 to 180
          repeat{
            get.mur.sst(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''), download.file = TRUE, dir = save.dir, ...) # filenames based on dates from above
            tryCatch({
              err <- try(RNetCDF::open.nc(paste(save.dir, filename, '_', time, '.nc', sep = '')), silent = T)
            }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
            if(class(err) != 'try-error') break
          }
        }
        
      }
    }
    
    
  } else if(type == 'hycom'){
    
    
    for(i in 1:length(uniqueDates)){
      time <- as.Date(uniqueDates[i])
      
      if (spatLim$lonmin >= -180 & spatLim$lonmax <= 180 & time < as.Date('2013-01-01')){
        # 180 coords, before 2013 -> normal
        print(paste('case 1'))
        repeat{
          get.hycom(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''),
                    download.file = TRUE, dir = save.dir, depLevels=depLevels, ...) 
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(save.dir,'/', filename, '_', time, '.nc', sep = '')), silent = T)
          }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
          if(class(err) != 'try-error') break
        }
        
      } else if (spatLim$lonmin >= 0 & spatLim$lonmax <= 180 & time >= as.Date('2013-01-01')){
        # 180 coords after 2013 -> hycom native is 360 so:
        #  - if you have 0 to 180 coords, you're good
        print(paste('case 2'))
        
        repeat{
          get.hycom(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''),
                    download.file = TRUE, dir = save.dir, depLevels=depLevels, ...) 
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(save.dir,'/', filename, '_', time, '.nc', sep = '')), silent = T)
          }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
          if(class(err) != 'try-error') break
        }
        
      } else if(spatLim$lonmin < 0 & spatLim$lonmax <= 0 & time >= as.Date('2013-01-01')){
        # 180 coords after 2013 -> hycom native is 360 so:
        #  - if you have -180 to 0, convert to get hycom then convert back
        print(paste('case 3'))
        
        spatLim.temp <- spatLim
        spatLim.temp$lonmin <- make360(spatLim$lonmin)
        spatLim.temp$lonmax <- make360(spatLim$lonmax)
        repeat{
          get.hycom(spatLim.temp, time, filename = paste(filename, '_', time, '.nc', sep = ''),
                    download.file = TRUE, dir = save.dir, depLevels=depLevels, ...) 
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(save.dir,'/', filename, '_', time, '.nc', sep = '')), silent = T)
          }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
          if(class(err) != 'try-error') break
        }
        
        ## get var names
        nc1 <- RNetCDF::open.nc(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))
        ncnames = NULL
        nmax <- RNetCDF::file.inq.nc(nc1)$nvars - 1
        for(ii in 0:nmax) ncnames[ii + 1] <- RNetCDF::var.inq.nc(nc1, ii)$name
        temp.idx <- grep('temp', ncnames, ignore.case=TRUE) - 1
        lat.idx <- grep('lat', ncnames, ignore.case=TRUE) - 1
        lon.idx <- grep('lon', ncnames, ignore.case=TRUE) - 1
        #dep.idx <- grep('dep', ncnames, ignore.case=TRUE) - 1
        
        # get attributes, if they exist
        ncatts <- NULL
        nmax <- RNetCDF::var.inq.nc(nc1, temp.idx)$natts - 1
        for(ii in 0:nmax) ncatts[ii + 1] <- RNetCDF::att.inq.nc(nc1, temp.idx, ii)$name
        scale.idx <- grep('scale', ncatts, ignore.case=TRUE) - 1
        if(length(scale.idx) != 0){
          scale <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=scale.idx)
        } else{
          scale <- 1
        }
        off.idx <- grep('off', ncatts, ignore.case=TRUE) - 1
        if(length(off.idx) != 0){
          offset <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=off.idx)
        } else{
          offset <- 0
        }
        
        lon1 <- RNetCDF::var.get.nc(nc1, lon.idx)
        lat1 <- RNetCDF::var.get.nc(nc1, lat.idx)
        dat1 <- RNetCDF::var.get.nc(nc1, temp.idx) * scale + offset
        r1 <- raster::flip(raster::brick(aperm(dat1, c(2,1,3)), xmn=min(lon1), xmx=max(lon1),
                                         ymn=min(lat1), ymx=max(lat1)), 2)
        
        #r1 <- raster::brick(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))
        r1r <- raster::rotate(r1)
        
        if (!is.null(template)){
          print('Input raster template is being used to resample the downloaded hycom grid. This is an experimental feature.')
          r1r <- raster::resample(r1r, template)
        }
        
        raster::writeRaster(r1r, paste(save.dir, '/', filename, '_', time, '.nc', sep = ''), format='CDF', overwrite=TRUE, varname = ncnames[grep('temp', ncnames, ignore.case=TRUE)])
        if (file.exists(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))) cat(paste0('File output to ', paste(save.dir, '/', filename, '_', time, '.nc', sep = '')),'\n')
        
      } else if(spatLim$lonmin < 0 & spatLim$lonmax > 0 & spatLim$lonmax <= 180 & time >= as.Date('2013-01-01')){
        # 180 coords after 2013 -> hycom native is 360 so:
        #  - if you span 0, get both "normal" 0 to 180 and converted < 0 coords and combine
        print(paste('case 4'))
        
        cat('Detected input coordinates in 180 range but HYCOM is 360 starting in 2013, downloading multiple files to fix this.','\n')
        ex180 <- raster::extent(-180,-0.0001,-90,90)
        ex360 <- raster::extent(0,360,-90,90)
        
        ## get part 1 - whatever is below 0
        ex1 <- raster::intersect(raster::extent(unlist(spatLim)), ex180)
        ex1@xmin <- make360(ex1@xmin)
        ex1@xmax <- make360(ex1@xmax)
        
        tdir <- tempdir()
        repeat{
          get.hycom(c(ex1@xmin, ex1@xmax, ex1@ymin, ex1@ymax), time, filename = paste(filename, '_', time, '_1.nc', sep = ''),
                    download.file = TRUE, dir = tdir, depLevels=depLevels)#, ...) 
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_1.nc', sep = '')), silent = T)
          }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
          if(class(err) != 'try-error') break
        }
        
        ## get part 2 - whatever is >= 0, if any
        ex2 <- raster::intersect(raster::extent(unlist(spatLim)), ex360)
        if (is.null(ex2)){
          ## the first download got the full extent so just copy it over to final file
          file.copy(paste0(tdir, '/', filename, '_', time, '_1.nc'),
                    paste0(save.dir,'/', filename, '_', time, '.nc'))
          next
        } else{
          repeat{
            get.hycom(c(ex2@xmin, ex2@xmax, ex2@ymin, ex2@ymax), time, filename = paste(filename, '_', time, '_2.nc', sep = ''),
                      download.file = TRUE, dir = tdir, depLevels=depLevels)#, ...) 
            tryCatch({
              err <- try(RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_2.nc', sep = '')), silent = T)
            }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
            if(class(err) != 'try-error') break
          }
        }
        
        cat('merging those files to a single output file','\n')
        nc1 <- RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_1.nc', sep = ''))
        
        ## get var names
        ncnames = NULL
        nmax <- RNetCDF::file.inq.nc(nc1)$nvars - 1
        for(ii in 0:nmax) ncnames[ii + 1] <- RNetCDF::var.inq.nc(nc1, ii)$name
        temp.idx <- grep('temp', ncnames, ignore.case=TRUE) - 1
        lat.idx <- grep('lat', ncnames, ignore.case=TRUE) - 1
        lon.idx <- grep('lon', ncnames, ignore.case=TRUE) - 1
        #dep.idx <- grep('dep', ncnames, ignore.case=TRUE) - 1
        
        # get attributes, if they exist
        ncatts <- NULL
        nmax <- RNetCDF::var.inq.nc(nc1, temp.idx)$natts - 1
        for(ii in 0:nmax) ncatts[ii + 1] <- RNetCDF::att.inq.nc(nc1, temp.idx, ii)$name
        scale.idx <- grep('scale', ncatts, ignore.case=TRUE) - 1
        if(length(scale.idx) != 0){
          scale <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=scale.idx)
        } else{
          scale <- 1
        }
        off.idx <- grep('off', ncatts, ignore.case=TRUE) - 1
        if(length(off.idx) != 0){
          offset <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=off.idx)
        } else{
          offset <- 0
        }
        
        lon1 <- RNetCDF::var.get.nc(nc1, lon.idx)
        lat1 <- RNetCDF::var.get.nc(nc1, lat.idx)
        dat1 <- RNetCDF::var.get.nc(nc1, temp.idx) * scale + offset
        r1 <- raster::flip(raster::brick(aperm(dat1, c(2,1,3)), xmn=min(lon1), xmx=max(lon1),
                                         ymn=min(lat1), ymx=max(lat1)), 2)
        
        nc2 <- RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_2.nc', sep = ''))
        lon2 <- RNetCDF::var.get.nc(nc2, lon.idx)
        lat2 <- RNetCDF::var.get.nc(nc2, lat.idx)
        dat2 <- RNetCDF::var.get.nc(nc2, temp.idx) * scale + offset
        r2 <- raster::flip(raster::brick(aperm(dat2, c(2,1,3)), xmn=min(lon2), xmx=max(lon2),
                                         ymn=min(lat2), ymx=max(lat2)), 2)
        
        ## reset original directory after messing w temp files
        setwd(original_dir)
        
        #print(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))
        r1r <- raster::rotate(r1)
        cr1 <- try(compareRaster(r1,r2), TRUE)
        cr2 <- try(compareRaster(r1r,r2), TRUE)
        if (class(cr1) != 'try-error'){
          r3 <- raster::merge(r1, r2)
        } else if (class(cr2) != 'try-error'){
          r3 <- raster::merge(r1r, r2)
        } else{
          e <- raster::extent(spatLim$lonmin, spatLim$lonmax, spatLim$latmin, spatLim$latmax)
          template <- raster::raster(e, ncols = ceiling((e@xmax - e@xmin) / raster::res(r1)[1]),
                                     nrows = ceiling((e@ymax - e@ymin) / raster::res(r1)[1]))
          raster::projection(template) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
          raster::crs(r1) <- raster::crs(template)
          raster::crs(r2) <- raster::crs(template)
          r1p <- try(raster::resample(r1, template), TRUE)
          if (class(r1p) == 'try-error') r1p <- raster::resample(r1r, template)
          r2p <- raster::resample(r2, template)
          r3 <- raster::merge(r1p, r2p)
          
        }
        #r3 <- raster::flip(r3, 'y')
        raster::writeRaster(r3, paste(save.dir, '/', filename, '_', time, '.nc', sep = ''), format='CDF', overwrite=TRUE, varname = ncnames[grep('temp', ncnames, ignore.case=TRUE)])
        if (file.exists(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))) cat(paste0('File output to ', paste(save.dir, '/', filename, '_', time, '.nc', sep = '')),'\n')
        
      } else if(spatLim$lonmin > 180 & spatLim$lonmax > 180 & time < as.Date('2013-01-01')){
        # you have 360 coords, before 2013 -> hycom native is 180 so:
        #  - if > 180, convert your coords 
        print(paste('case 5'))
        
        spatLim.temp <- spatLim
        spatLim.temp$lonmin <- make180(spatLim$lonmin)
        spatLim.temp$lonmax <- make180(spatLim$lonmax)
        repeat{
          get.hycom(spatLim.temp, time, filename = paste(filename, '_', time, '.nc', sep = ''),
                    download.file = TRUE, dir = save.dir, depLevels=depLevels, ...) 
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(save.dir,'/', filename, '_', time, '.nc', sep = '')), silent = T)
          }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
          if(class(err) != 'try-error') break
        }
        
      } else if(spatLim$lonmin < 180 & spatLim$lonmax > 180 & time < as.Date('2013-01-01')){
        # you have 360 coords, before 2013 -> hycom native is 180 so:
        #  - if span 180, do both and combine
        print(paste('case 6'))
        
        cat('Detected input coordinates > 180 that also span the 180 line, downloading multiple files.','\n')
        ex180 <- raster::extent(-180,180,-90,90)
        ex360 <- raster::extent(180,360,-90,90)
        
        ## get part 1 - whatever is above 180
        ex1 <- raster::intersect(raster::extent(unlist(spatLim)), ex360)
        ex1 <- raster::extent(raster::rotate(raster::raster(ex1)))
        #if (ex1@xmin == -180) ex1@xmin <- -179.99
        #original_dir <- getwd()
        tdir <- tempdir()
        repeat{
          get.hycom(c(ex1@xmin, ex1@xmax, ex1@ymin, ex1@ymax), time, filename = paste(filename, '_', time, '_1.nc', sep = ''),
                    download.file = TRUE, dir = tdir, depLevels=depLevels)#, ...) 
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_1.nc', sep = '')), silent = T)
          }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
          if(class(err) != 'try-error') break
        }
        
        ## get part 2 - whatever is below 180, if any
        ex2 <- raster::intersect(raster::extent(unlist(spatLim)), ex180)
        if (is.null(ex2)){
          ## the first download got the full extent so just copy it over to final file
          file.copy(paste0(tdir, '/', filename, '_', time, '_1.nc'),
                    paste0(save.dir,'/', filename, '_', time, '.nc'))
          next
        } else{
          repeat{
            get.hycom(c(ex2@xmin, ex2@xmax, ex2@ymin, ex2@ymax), time, filename = paste(filename, '_', time, '_2.nc', sep = ''),
                      download.file = TRUE, dir = tdir, depLevels=depLevels)#, ...) 
            tryCatch({
              err <- try(RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_2.nc', sep = '')), silent = T)
            }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
            if(class(err) != 'try-error') break
          }
        }
        
        cat('merging those files to a single output file','\n')
        nc1 <- RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_1.nc', sep = ''))
        
        ## get var names
        ncnames = NULL
        nmax <- RNetCDF::file.inq.nc(nc1)$nvars - 1
        for(ii in 0:nmax) ncnames[ii + 1] <- RNetCDF::var.inq.nc(nc1, ii)$name
        temp.idx <- grep('temp', ncnames, ignore.case=TRUE) - 1
        lat.idx <- grep('lat', ncnames, ignore.case=TRUE) - 1
        lon.idx <- grep('lon', ncnames, ignore.case=TRUE) - 1
        #dep.idx <- grep('dep', ncnames, ignore.case=TRUE) - 1
        
        # get attributes, if they exist
        ncatts <- NULL
        nmax <- RNetCDF::var.inq.nc(nc1, temp.idx)$natts - 1
        for(ii in 0:nmax) ncatts[ii + 1] <- RNetCDF::att.inq.nc(nc1, temp.idx, ii)$name
        scale.idx <- grep('scale', ncatts, ignore.case=TRUE) - 1
        if(length(scale.idx) != 0){
          scale <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=scale.idx)
        } else{
          scale <- 1
        }
        off.idx <- grep('off', ncatts, ignore.case=TRUE) - 1
        if(length(off.idx) != 0){
          offset <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=off.idx)
        } else{
          offset <- 0
        }
        
        lon1 <- RNetCDF::var.get.nc(nc1, lon.idx)
        lat1 <- RNetCDF::var.get.nc(nc1, lat.idx)
        dat1 <- RNetCDF::var.get.nc(nc1, temp.idx) * scale + offset
        r1 <- raster::flip(raster::brick(aperm(dat1, c(2,1,3)), xmn=min(lon1), xmx=max(lon1),
                                         ymn=min(lat1), ymx=max(lat1)), 2)
        
        nc2 <- RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_2.nc', sep = ''))
        lon2 <- RNetCDF::var.get.nc(nc2, lon.idx)
        lat2 <- RNetCDF::var.get.nc(nc2, lat.idx)
        dat2 <- RNetCDF::var.get.nc(nc2, temp.idx) * scale + offset
        r2 <- raster::flip(raster::brick(aperm(dat2, c(2,1,3)), xmn=min(lon2), xmx=max(lon2),
                                         ymn=min(lat2), ymx=max(lat2)), 2)
        
        ## reset original directory after messing w temp files
        setwd(original_dir)
        
        #print(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))
        r1r <- raster::shift(raster::rotate(raster::shift(r1, 180)), 180)
        cr1 <- try(compareRaster(r1,r2), TRUE)
        cr2 <- try(compareRaster(r1r,r2), TRUE)
        if (class(cr1) != 'try-error'){
          r3 <- raster::merge(r1, r2)
        } else if (class(cr2) != 'try-error'){
          r3 <- raster::merge(r1r, r2)
        } else{
          e <- raster::extent(spatLim$lonmin, spatLim$lonmax, spatLim$latmin, spatLim$latmax)
          template <- raster::raster(e, ncols = ceiling((e@xmax - e@xmin) / raster::res(r1)[1]),
                                     nrows = ceiling((e@ymax - e@ymin) / raster::res(r1)[1]))
          raster::projection(template) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
          raster::crs(r1) <- raster::crs(template)
          raster::crs(r2) <- raster::crs(template)
          r1p <- try(raster::resample(r1, template), TRUE)
          if (class(r1p) == 'try-error') r1p <- raster::resample(r1r, template)
          r2p <- raster::resample(r2, template)
          r3 <- raster::merge(r1p, r2p)
          
        }
        #r3 <- raster::flip(r3, 'y')
        raster::writeRaster(r3, paste(save.dir, '/', filename, '_', time, '.nc', sep = ''), format='CDF', overwrite=TRUE, varname = ncnames[grep('temp', ncnames, ignore.case=TRUE)])
        if (file.exists(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))) cat(paste0('File output to ', paste(save.dir, '/', filename, '_', time, '.nc', sep = '')),'\n')
        
      } else if(spatLim$lonmin > 180 & spatLim$lonmax > 180 & time >= as.Date('2013-01-01')){
        # you have 360 coords, after 2013 -> hycom native is 360 so:
        #  - if you have > 180, your request is good but need to convert output
        cat('HYCOM outputs are 0-360 starting in 2013. Converting your output HYCOM data to the 180 coordinate system.','\n')
        
        print(paste('case 7'))
        
        tdir <- tempdir()
        repeat{
          get.hycom(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''),
                    download.file = TRUE, dir = tdir, depLevels=depLevels)#, ...) 
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '.nc', sep = '')), silent = T)
          }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
          if(class(err) != 'try-error') break
        }
        
        ## get var names
        nc1 <- RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '.nc', sep = ''))
        ncnames = NULL
        nmax <- RNetCDF::file.inq.nc(nc1)$nvars - 1
        for(ii in 0:nmax) ncnames[ii + 1] <- RNetCDF::var.inq.nc(nc1, ii)$name
        temp.idx <- grep('temp', ncnames, ignore.case=TRUE) - 1
        lat.idx <- grep('lat', ncnames, ignore.case=TRUE) - 1
        lon.idx <- grep('lon', ncnames, ignore.case=TRUE) - 1
        #dep.idx <- grep('dep', ncnames, ignore.case=TRUE) - 1
        
        # get attributes, if they exist
        ncatts <- NULL
        nmax <- RNetCDF::var.inq.nc(nc1, temp.idx)$natts - 1
        for(ii in 0:nmax) ncatts[ii + 1] <- RNetCDF::att.inq.nc(nc1, temp.idx, ii)$name
        scale.idx <- grep('scale', ncatts, ignore.case=TRUE) - 1
        if(length(scale.idx) != 0){
          scale <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=scale.idx)
        } else{
          scale <- 1
        }
        off.idx <- grep('off', ncatts, ignore.case=TRUE) - 1
        if(length(off.idx) != 0){
          offset <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=off.idx)
        } else{
          offset <- 0
        }
        
        lon1 <- RNetCDF::var.get.nc(nc1, lon.idx)
        lat1 <- RNetCDF::var.get.nc(nc1, lat.idx)
        dat1 <- RNetCDF::var.get.nc(nc1, temp.idx) * scale + offset
        r1 <- raster::flip(raster::brick(aperm(dat1, c(2,1,3)), xmn=min(lon1), xmx=max(lon1),
                                         ymn=min(lat1), ymx=max(lat1)), 2)
        
        #r1 <- raster::brick(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))
        r1r <- raster::rotate(r1)
        
        raster::writeRaster(r1r, paste(save.dir, '/', filename, '_', time, '.nc', sep = ''), format='CDF', overwrite=TRUE, varname = ncnames[grep('temp', ncnames, ignore.case=TRUE)])
        if (file.exists(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))) cat(paste0('File output to ', paste(save.dir, '/', filename, '_', time, '.nc', sep = '')),'\n')
        
      } else if(spatLim$lonmin < 180 & spatLim$lonmax > 180 & time >= as.Date('2013-01-01')){
        # you have 360 coords, after 2013 -> hycom native is 360 so:
        #  - if you span 180, get hycom as normal
        print(paste('case 8'))
        warning('Getting HYCOM data for longitudes > 180 that spans the 180 line, therefore this data cannot be reasonably converted to 180 coordinate system (i.e. you need to run HMMoce using 360 coords). This requires an advanced HMMoce skillset and extreme caution!')
        
        repeat{
          get.hycom(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''),
                    download.file = TRUE, dir = save.dir, depLevels=depLevels, ...) 
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(save.dir,'/', filename, '_', time, '.nc', sep = '')), silent = T)
          }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
          if(class(err) != 'try-error') break
        }
        
      } else {
        stop('No path identified to collect and re-format (if applicable) the appropriate HYCOM data for your combinaton of spatial limits and time.')
      }
   
    } # end udates
    
  } else if (type == 'woa'){
    
    if(is.null(resol)){
      stop('Error: If type = woa then resol must be specified. See ?get.env for help.')
    }
    
    filename <- get.woa(save.dir = save.dir, resol = resol, ...)
    print(paste('WOA data downloaded to ', filename,'...', sep=''))
    
  }
  
  setwd(original_dir)
  
  
}
