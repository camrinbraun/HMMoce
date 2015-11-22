extract.woa <- function(nc.dir,bbox,resolution, pdt){
  # Simple function to extract the desired temperature data from a global
  # dataset derived from monthly, 1/4 deg gridded climatology data 
  # contained in the 2013 World Ocean Atlas
  
  # Inputs (both required):
  #   NC.DIR is the directory to load the global nc file from; make sure it's
  #       the only .nc file in the given directory
  #   BBOX is a bounding box of form (long min, long max, lat min, lat max)
  #   RESOLUTION indicates whether oceanographic data is gridded at 'quarter'
  #       or 'one' degree resolution
  #   PDT is output dataframe from convert.pdt
  
  # Output:
  #   RETURNWOA is a list containing:
  #     DAT is an array of temperature data with dimensions (long, lat, depth, time)
  #       depth contains 57 standard depth levels by default and levels are defined
  #       in variable 'depth' contained here. time dimension covers the months of
  #       tag deployment as gathered from querying month data in variable 'pdt'
  #     LON/LAT are vectors of lon/lat bounds
  
  #--------------------------------------------------------------------------
  # - Written by Camrin Braun in R, December 20, 2014. cbraun@whoi.edu
  #--------------------------------------------------------------------------
  
  # load global nc
  ncfiles=dir(nc.dir, pattern='.nc')
  nc = open.ncdf(paste(nc.dir,ncfiles,sep=''))
  
  # retrieve var bounds from global nc
  lon = get.var.ncdf(nc, 'Longitude')
  lat = get.var.ncdf(nc, 'Latitude')
  depth = get.var.ncdf(nc, 'Depth')
  
  # set bounds for extracting data
  xmin = which.min((bbox[1]-lon)^2); xmax = which.min((bbox[2]-lon)^2)
  ymin = which.min((bbox[3]-lat)^2); ymax = which.min((bbox[4]-lat)^2)
  
  if(resolution == 'quarter'){
    xlen = 4*(bbox[2]-bbox[1]) # for quarter degree
    ylen = 4*(bbox[4]-bbox[3]) 
  } else if(resolution == 'one'){
    xlen = bbox[2]-bbox[1] # for one degree
    ylen = bbox[4]-bbox[3]
  } else{
    stop('Resolution of input oceanographic data not defined.')
  }  
  
  # define time bounds using tag data
  #month <- as.numeric(format(pdt$Date, format='%m'))
  #tmin = min(month); tmax = max(month); tlen = tmax - tmin + 1
  
  #if (tlen <= 1){
  #  tlen = 2
  #}
  
  dat = get.var.ncdf(nc, 'temp', start = c(xmin,ymin,1,1), count = c(xlen+1,ylen+1,57,12))
  
  returnWOA = list(dat = dat, lon = lon[xmin:xmax], lat = lat[ymin:ymax], depth = depth)
  
}