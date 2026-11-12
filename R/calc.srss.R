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
#' @param focalDim is integer for dimensions of terra::focal used to calculate 
#'   sd() of SRSS grid cell. Recommend focalDim = 3 if res=1 and focalDim = 9 if
#'   res=0.25.
#' @export
#' @return L is a SpatRaster of dim(lon x lat x dateVec) containing likelihood
#'   surfaces for each time point
#' @seealso \code{\link{calc.lightloc}}
#' 

calc.srss <- function(light = NULL, locs.grid, dateVec, res = 1, focalDim = 3){
  
  t0 <- Sys.time()
  print(paste('Starting light likelihood calculation using SRSS times...'))
  
  lon <- locs.grid$lon[1,]
  if (any(lon > 180)) stop('calc.srss cannot handle longitude values > 180 at this time')
  lat <- locs.grid$lat[,1]
  lon <- seq(min(lon), max(lon), res)
  lat <- seq(min(lat), max(lat), res)
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  dx <- if (length(lon) > 1) (max(lon) - min(lon)) / (length(lon) - 1) else res
  dy <- if (length(lat) > 1) (max(lat) - min(lat)) / (length(lat) - 1) else res
  ext_ll <- terra::ext(min(lon) - dx / 2, max(lon) + dx / 2, min(lat) - dy / 2, max(lat) + dy / 2)
  
  L_arr <- array(0, dim = c(length(lat), length(lon), length(dateVec)))
  L.light <- terra::rast(L_arr, extent = ext_ll, crs = crs)
  
  print(paste('Building SRSS grids', '...'))
  
  sr.grid <- array(0, dim = c(length(lon), length(lat), 365))
  ss.grid <- sr.grid
  fyear <- seq(ISOdate(lubridate::year(dateVec[1]), 1, 1, tz = 'UTC'), ISOdate(lubridate::year(dateVec[1]), 12, 31, tz = 'UTC'), 'day')
  
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
  sr.grid[,,1:365] <- sapply(1:365, function(i) matrix(calc_sun_day(pts_grid$lon, pts_grid$lat, fyear[i], direction = "sunrise"), length(lon), length(lat)))
  ss.grid[,,1:365] <- sapply(1:365, function(i) matrix(calc_sun_day(pts_grid$lon, pts_grid$lat, fyear[i], direction = "sunset"), length(lon), length(lat)))
  
  sr_arr <- aperm(sr.grid * 24 * 60, c(2, 1, 3))[length(lat):1, , , drop = FALSE]
  ss_arr <- aperm(ss.grid * 24 * 60, c(2, 1, 3))[length(lat):1, , , drop = FALSE]
  sr.ras <- terra::rast(sr_arr, extent = ext_ll, crs = crs)
  ss.ras <- terra::rast(ss_arr, extent = ext_ll, crs = crs)
  
  min.sr <- apply(sr_arr, 3, min, na.rm = TRUE)
  max.sr <- apply(sr_arr, 3, max, na.rm = TRUE)
  min.ss <- apply(ss_arr, 3, min, na.rm = TRUE)
  max.ss <- apply(ss_arr, 3, max, na.rm = TRUE)
  
  light <- light[which(light$Date <= max(dateVec)),]
  light$dateVec <- findInterval(light$Date, dateVec)
  light$yday <- lubridate::yday(light$Date)
  light$daymins <- lubridate::minute(light$Date) + (lubridate::hour(light$Date) * 60)
  light <- light[which(light$Type != ''),]
  if (class(light$Date)[1] != class(dateVec)[1]) stop('dateVec and light$Date both need to be of class POSIXct.')
  
  print(paste('Starting daily calculations', '...'))
  
  T <- length(dateVec)
  
  for(t in 1:T){
    light.t <- light[which(light$dateVec == t),]
    
    if(nrow(light.t) > 0){
      if(nrow(light.t) == 1 & any(light.t$Type == 'Dawn')){
        didx <- light.t$yday[1]
        sr <- light.t$daymins[which(light.t$Type == 'Dawn')]
        if(sr < min.sr[didx] | sr > max.sr[didx]) sr <- NA
        
        srf <- terra::focal(sr.ras[[didx]], w = focalDim, fun = "sd", na.rm = TRUE)
        srlik <- liksrss(sr, srss = sr.ras[[didx]], srsd = srf)
        L.light[[t]] <- srlik
        
      } else if(nrow(light.t) == 1 & any(light.t$Type == 'Dusk')){
        didx <- light.t$yday[1]
        ss <- light.t$daymins[which(light.t$Type == 'Dusk')]
        if(ss < min.ss[didx] | ss > max.ss[didx]) ss <- NA
        
        ssf <- terra::focal(ss.ras[[didx]], w = focalDim, fun = "sd", na.rm = TRUE)
        sslik <- liksrss(ss, srss = ss.ras[[didx]], srsd = ssf)
        L.light[[t]] <- sslik
        
      } else{
        didx <- light.t$yday[1]
        sr <- light.t$daymins[which(light.t$Type == 'Dawn')]
        ss <- light.t$daymins[which(light.t$Type == 'Dusk')]
        
        if(length(sr) == 0) sr <- NA
        if(length(ss) == 0) ss <- NA
        if(length(sr) > 1) sr <- sr[1]
        if(length(ss) > 1) ss <- ss[length(ss)]
        
        if(!is.na(sr) && (sr < min.sr[didx] | sr > max.sr[didx])) sr <- NA
        if(!is.na(ss) && (ss < min.ss[didx] | ss > max.ss[didx])) ss <- NA
        
        srf <- terra::focal(sr.ras[[didx]], w = focalDim, fun = "sd", na.rm = TRUE)
        srlik <- liksrss(sr, srss = sr.ras[[didx]], srsd = srf)
        
        ssf <- terra::focal(ss.ras[[didx]], w = focalDim, fun = "sd", na.rm = TRUE)
        sslik <- liksrss(ss, srss = ss.ras[[didx]], srsd = ssf)
        
        sr_v <- terra::values(srlik)
        ss_v <- terra::values(sslik)
        if(any(sr_v != 0, na.rm = TRUE) & any(ss_v != 0, na.rm = TRUE)){
          r <- srlik * sslik
          r_v <- terra::values(r)
          if(any(r_v != 0, na.rm = TRUE)){
            max_cell <- which.max(r_v)
            max_row <- terra::rowFromCell(r, max_cell)
            r_mat <- matrix(r_v, nrow = nrow(r), ncol = ncol(r), byrow = TRUE)
            row_vals <- r_mat[max_row, ]
            terra::values(r) <- rep(row_vals, times = nrow(r))
          } else if(any(sr_v != 0, na.rm = TRUE)){
            r <- srlik
          } else if(any(ss_v != 0, na.rm = TRUE)){
            r <- sslik
          }
        } else{
          if(any(sr_v != 0, na.rm = TRUE)){
            r <- srlik
          } else{
            r <- sslik
          }
        }
        L.light[[t]] <- r
      }
    }
    
    v_t <- terra::values(L.light[[t]])
    mx_t <- max(v_t, na.rm = TRUE)
    if(is.finite(mx_t) && mx_t != 0){
      L.light[[t]] <- L.light[[t]] / mx_t
    }
  }
  
  L.light[L.light < 0] <- 0
  
  t1 <- Sys.time()
  print(paste('Light calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  L.light
}
