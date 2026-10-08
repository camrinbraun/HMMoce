#' Calculate Position-based Likelihood from SRSS
#' 
#' \code{calc.srss} calculates likelihood estimates for each day's estimated
#' dawn dusk times from the tag
#' 
#' Tag-measured sunrise/sunset (SRSS) times are first filtered according to the
#' min/max times possible. possible values are computed based on spatial limits
#' included in locs.grid input. after filtering, each yearday has a grid of
#' georeferenced SRSS times that the tag-measured times are compared to in order
#' to generate a likelihood.
#' 
#' @param light is data frame of Date (POSIXct) and Type where type is either "Dawn" or "Dusk". This is a standard product from -LightLoc file output from DAP/Tag Portal for WC tags.
#' @param locs.grid is list output from \code{setup.locs.grid}
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param res is resolution of light grid in degrees. default is 1 deg. higher
#'   resolution (e.g. res = .25 for 1/4 deg) takes considerably longer to
#'   compute and may not be appropriate for most scenarios
#' @param focalDim is integer for dimensions of raster::focal used to calculate 
#'   sd() of SRSS grid cell. Recommend focalDim = 3 if res=1 and focalDim = 9 if
#'   res=0.25.
#' @export
#' @return L is a raster of dim(lon x lat x dateVec) containing likelihood
#'   surfaces for each time point
#' @seealso \code{\link{calc.lightloc}}
#' 

calc.srss <- function(light = NULL, locs.grid, dateVec, res = 1, focalDim = 3){
  
  #options(warn=2)
  t0 <- Sys.time()
  print(paste('Starting light likelihood calculation using SRSS times...'))
  
  #=====
  # BUILD A SPOT FOR THE RESULTS TO GO
  #=====
  
  # need lat/lon vectors from locs.grid
  lon <- locs.grid$lon[1,]
  if (any(lon > 180)) stop('calc.srss cannot handle longitude values > 180 at this time')
  lat <- locs.grid$lat[,1]
  # then rebuild these vectors based on input resolution to this particular function, default is 1 deg
  lon <- seq(min(lon), max(lon), res)
  lat <- seq(min(lat), max(lat), res)
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  L.grid = numeric(length = c(length(lon) * length(lat) * length(dateVec)))
  dim(L.grid) = c(length(lon), length(lat), length(dateVec))
  list.ras <- list(x = lon, y = lat, z = L.grid)
  ex <- raster::extent(list.ras)
  L.light <- raster::brick(list.ras$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
  L.light <- raster::flip(L.light, direction = 'y')
  
  #==================
  # build the SRSS grids
  #==================
  print(paste('Building SRSS grids', '...'))
  
  # expand.grid and SpatialPoints establishes a grid
  xy = as.matrix(expand.grid(lon,lat))
  xy = sp::SpatialPoints(xy, proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
  
  # now do the building and rasterize
  sr.grid = numeric(length = c(length(lon) * length(lat) * 365))
  dim(sr.grid) = c(length(lon), length(lat), 365)
  ss.grid = sr.grid
  t <- Sys.time()
  fyear = seq(ISOdate(lubridate::year(dateVec[1]), 1, 1, tz = 'UTC'), ISOdate(lubridate::year(dateVec[1]), 12, 31, tz = 'UTC'), 'day')
  
  # NOAA solar calculation for sunrise/sunset fraction of day (replaces retired maptools::sunriset)
  calc_sun_day <- function(lon_vec, lat_vec, date_val, direction = "sunrise") {
    n <- as.numeric(difftime(date_val, as.POSIXct("2000-01-01 12:00:00", tz = "UTC"), units = "days"))
    Jstar <- n - (lon_vec / 360)
    M <- (357.5291 + 0.98560028 * Jstar) %% 360
    C <- 1.9148 * sin(M * pi / 180) + 0.02 * sin(2 * M * pi / 180) + 0.0003 * sin(3 * M * pi / 180)
    lambda <- (M + C + 180 + 102.9372) %% 360
    Jtransit <- 2451545.0 + Jstar + 0.0053 * sin(M * pi / 180) - 0.0069 * sin(2 * lambda * pi / 180)
    delta <- asin(sin(lambda * pi / 180) * sin(23.44 * pi / 180))
    cos_omega <- (sin(-0.83 * pi / 180) - sin(lat_vec * pi / 180) * sin(delta)) / (cos(lat_vec * pi / 180) * cos(delta))
    cos_omega <- pmax(pmin(cos_omega, 1), -1)
    omega <- acos(cos_omega) * 180 / pi
    if (direction == "sunrise") {
      J_event <- Jtransit - (omega / 360)
    } else {
      J_event <- Jtransit + (omega / 360)
    }
    return((J_event + 0.5) %% 1)
  }

  pts_grid <- expand.grid(lon = lon, lat = lat)
  sr.grid[,,1:365] = sapply(1:365, function(i) matrix(calc_sun_day(pts_grid$lon, pts_grid$lat, fyear[i], direction = "sunrise"), length(lon), length(lat)))
  ss.grid[,,1:365] = sapply(1:365, function(i) matrix(calc_sun_day(pts_grid$lon, pts_grid$lat, fyear[i], direction = "sunset"), length(lon), length(lat)))
  
  
  list.ras <- list(x = lon, y = lat, z = sr.grid*24*60)
  ex <- raster::extent(list.ras)
  sr.ras <- raster::brick(list.ras$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
  sr.ras <- raster::flip(sr.ras, direction = 'y')
  
  list.ras <- list(x = lon, y = lat, z = ss.grid*24*60)
  ss.ras <- raster::brick(list.ras$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
  ss.ras <- raster::flip(ss.ras, direction = 'y')
  
  # need to be able to cut SRSS times from tag that aren't within limits of the grid
  min.sr <- sapply(1:365, function(i) raster::cellStats(sr.ras[[i]],stat='min',na.rm=T))
  max.sr <- sapply(1:365, function(i) raster::cellStats(sr.ras[[i]],stat='max',na.rm=T))
  min.ss <- sapply(1:365, function(i) raster::cellStats(ss.ras[[i]],stat='min',na.rm=T))
  max.ss <- sapply(1:365, function(i) raster::cellStats(ss.ras[[i]],stat='max',na.rm=T))
  
  # make some calculations on the tag data: yday, dtime, etc
  
  light <- light[which(light$Date <= max(dateVec)),]
  light$dateVec <- findInterval(light$Date, dateVec)
  #by_dte <- dplyr::group_by(light, as.factor(light$dateVec))  # group by unique time points
  #light <- data.frame(dplyr::summarise_(by_dte, "min(Temperature)", "max(Temperature)"))
  #colnames(light) <- list('time', 'minT', 'maxT')
  #light$time <- dateVec[as.numeric(as.character(light$time))]
  light$yday <- lubridate::yday(light$Date)
  light$daymins <- lubridate::minute(light$Date) + (lubridate::hour(light$Date) * 60)
  light <- light[which(light$Type != ''),]
  #light$Date <- as.Date(format(light$Date, '%Y-%m-%d'))
  if (class(light$Date)[1] != class(dateVec)[1]) stop('dateVec and light$Date both need to be of class POSIXct.')
  
  print(paste('Starting daily calculations', '...'))
  
  T <- length(dateVec)
  
  for(t in 1:T){
    
    # data for this time step T
    light.t <- light[which(light$dateVec == t),]
    
    if(length(light.t[,1]) == 0){
      
    } else{
      if(length(light.t[,1]) == 1 & any(light.t$Type == 'Dawn')){
        # if we just have a dawn measurement
        didx <- light.t$yday[1]
        sr <- light.t$daymins[which(light.t$Type == 'Dawn')]
        
        if(sr < min.sr[didx] | sr > max.sr[didx]){
          sr <- NA
        }
        
        #light[which(light$Date %in% dateVec[t] & light$Type == 'Dawn'), 6] <- sr
        
        # now for likelihood
        # get the SD for this day, T
        srf <- raster::focal(sr.ras[[didx]], w = matrix(1, nrow = focalDim, ncol = focalDim), fun = function(x) stats::sd(x, na.rm = T))
        # the SR likelihood
        srlik <- liksrss(sr, srss = sr.ras[[didx]], srsd = srf)
        
        L.light[[t]] <- srlik
        
      } else if(length(light.t[,1]) == 1 & any(light.t$Type == 'Dusk')){
        # if we just have a dusk measurement
        didx <- light.t$yday[1]
        ss <- light.t$daymins[which(light.t$Type == 'Dusk')]
        
        if(ss < min.ss[didx] | ss > max.ss[didx]){
          ss <- NA
        }
        
        #light[which(light$Date %in% dateVec[t] & light$Type == 'Dusk'), 6] <- ss
        
        # and sunset
        ssf <- raster::focal(ss.ras[[didx]], w = matrix(1, nrow = focalDim, ncol = focalDim), fun = function(x) stats::sd(x, na.rm = T))
        sslik <- liksrss(ss, srss = ss.ras[[didx]], srsd = ssf)
        
        L.light[[t]] <- sslik
        
      } else{
        # if we have both dawn and dusk measurements
        didx <- light.t$yday[1]
        sr <- light.t$daymins[which(light.t$Type == 'Dawn')]
        ss <- light.t$daymins[which(light.t$Type == 'Dusk')]
        
        if(length(sr) == 0){sr <- NA}
        if(length(ss) == 0){ss <- NA}
        
        if(length(sr) > 1){
          # we want the first SR time if there are multiple
          sr <- sr[1]
        }
        
        if(length(ss) > 1){
          # we want the last SS time if there are multiple
          ss <- ss[length(ss)]
        }
        
        # filter based on possible grid values
        if(is.na(sr)){
          
        } else if(sr < min.sr[didx] | sr > max.sr[didx]){
          sr <- NA
        }
        
        if(is.na(ss)){
          
        } else if(ss < min.ss[didx] | ss > max.ss[didx]){
          ss <- NA
        }
        
        #light[which(light$Date %in% dateVec[t] & light$Type == 'Dawn'), 6] <- sr
        #light[which(light$Date %in% dateVec[t] & light$Type == 'Dusk'), 6] <- ss
        
        # now for likelihood
        # get the SD for this day, T
        srf <- raster::focal(sr.ras[[didx]], w = matrix(1, nrow = focalDim, ncol = focalDim), fun = function(x) stats::sd(x, na.rm = T))
        # the SR likelihood
        srlik <- liksrss(sr, srss = sr.ras[[didx]], srsd = srf)
        
        # and sunset
        ssf <- raster::focal(ss.ras[[didx]], w = matrix(1, nrow = focalDim, ncol = focalDim), fun = function(x) stats::sd(x, na.rm = T))
        sslik <- liksrss(ss, srss = ss.ras[[didx]], srsd = ssf)
        
        if(any(srlik[] != 0) & any(sslik[] != 0)){
          r <- srlik * sslik
          if(any(r[] != 0)){
            max.lat <- raster::xyFromCell(r, raster::which.max(r))[2]
            cds <- rbind(c(min(lon), max.lat), c(max(lon), max.lat))#, c(40,5), c(15,-45), c(-10,-25))
            lines <- sp::SpatialLines(list(sp::Lines(list(sp::Line(cds)), "1")))
            r[] <- c(unlist(raster::extract(r, lines)))
          } else if(any(srlik[] != 0)){
            #just srlik
            r <- srlik
          } else if(any(sslik[] != 0)){
            #just sslik
            r <- sslik
          }
          
        } else{
          if(any(srlik[] != 0)){
            #just srlik
            r <- srlik
          } else{
            #just sslik
            r <- sslik
          }
        }
        
        L.light[[t]] <- r
        
      }
      
    }
    
    if(raster::cellStats(L.light[[t]], 'max', na.rm=T) != 0){
      L.light[[t]] <- L.light[[t]] / raster::cellStats(L.light[[t]], 'max')
    }
    
  } # end for loop
  
  L.light[L.light < 0] <- 0
  
  t1 <- Sys.time()
  print(paste('Light calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  L.light
  
} # end function
