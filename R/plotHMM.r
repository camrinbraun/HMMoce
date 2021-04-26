#' Plot track results of HMMoce
#' 
#' \code{plotHMM} uses HMM output via \code{calc.track} to make simple plots of 
#' calculated track and behavior state
#' 
#' @param distr is output array from \code{hmm.smoother}
#' @param track is output dataframe from \code{calc.track}
#' @param dateVec is vector of dates from tag to pop-up location by day.
#' @param ptt is unique ID for individual tag dataset
#' @param known is 3 column data frame containing date, lat, lon of known 
#'   movement track. This is only useful for comparing HMMoce results to known 
#'   track collected by SPOT or GPS, for example. Default is NULL.
#' @param resid is logical indicating whether you want to include a residual
#'   plot. This is not yet functional.
#' @param behav.pts is logical indicating whether to plot points in the
#'   resulting map colored by behav state.
#' @param save.plot is logical indicating whether you want the plot written to
#'   disk using \code{pdf}.
#' @importFrom fields "world"
#'   
#' @return NULL. A plot is rendered on screen or written to disk.
#' @references Pedersen MW, Patterson TA, Thygesen UH, Madsen H (2011)
#'   Estimating animal behavior and residency from movement data. Oikos
#'   120:1281-1290. doi: 10.1111/j.1600-0706.2011.19044.x
#' @export

plotHMM <- function(distr, track, dateVec, ptt, known = NULL, resid = FALSE, behav.pts = F, save.plot = FALSE){

  ## Show movement as animation
  #if(show.movie) show.movie(s)
  
  #### One-state
  if(dim(distr)[1]==1){
    # calculated track
    xl <- c(floor(min(track$lon)), ceiling(max(track$lon)))
    yl <- c(floor(min(track$lat)), ceiling(max(track$lat)))
    graphics::plot(track$lon, track$lat, type = 'n', main = 'Estimated movements', ylab = 'Latitude', xlab = 'Longitude', xlim = xl, ylim = yl)
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = "steelblue1")
    fields::world(add = TRUE, fill = TRUE, col = 'grey60')
    
    if(!is.null(known)){
      graphics::lines(known$lon, known$lat, col = 'white')
    }
    
    graphics::lines(track$lon, track$lat, col = 'black')
    
    graphics::points(track$lon[1], track$lat[1], bg = 'green', pch = 21)
    TT <- length(track$lon)
    graphics::points(track$lon[TT], track$lat[TT], bg = 'red', pch = 21)
    
    if(save.plot) grDevices::dev.off()
    
    ## NOT YET FUNCTIONAL
    ## Simple diagnostics plot ###
    ## Resample SST
    #if(save.plot) grDevices::pdf('../plot/sphmmDiagn.pdf',width=7,height=7)
    #if(!save.plot)
    if (resid){
      # do nothing for now
      #grDevices::dev.new()
      #graphics::par(mfrow = c(2,2))
      #ssterr <- sst - lsst$sst
      #sdsst <- sqrt(stats::var(ssterr))
      #ssterr <- ssterr / sdsst
      #lonerr <- sphmm$meanlon - lsst$lon
      #sdlon <- sqrt(stats::var(lonerr))
      #lonerr <- lonerr / sdlon
      #graphics::plot(track$date[ind], ssterr, xlab = 'Date', ylab = 'SST residual', main = 'Residuals through time', ylim = c(-3,3), pch = ".", cex = 3)
      #graphics::abline(h = c(-2,0,2), col = 'grey', lwd = 2, lty = 2)
      #stats::qqnorm(ssterr, main = 'QQ-plot of residuals', pch = ".", cex = 2)
      #graphics::abline(a = 0, b = 1, col = 'grey', lwd = 2)
      #graphics::plot(track$date[ind], lonerr, xlab = 'Date', ylab = 'Longitude residual', ylim = c(-3,3), pch = ".", cex = 3)
      #graphics::abline(h = c(-2,0,2), col = 'grey', lwd = 2, lty = 2)
      #stats::qqnorm(lonerr, main = '', pch = ".", cex = 2)
      #graphics::abline(a = 0, b = 1, col = 'grey', lwd = 2)
    }
  }
  
  #### Two-state
  else{
  ## Calc behaviour
  # get states from s
  sv <- -(apply(distr[1,,,], 1, sum) > apply(distr[2,,,], 1, sum)) + 1
  sv[sv == 0] <- NA
  grDevices::graphics.off();

  #if (class(dateVec) != 'Date'){
  #  stop('Error: dateVec must be of class "Date". See ?as.Date.')
  #}

  if(save.plot) grDevices::pdf(paste(ptt, '_track_results.pdf', sep = ''), width = 7, height = 8)

  graphics::par(mfrow = c(2,1))
  # behavior state plot
  graphics::plot(I(sv) ~ dateVec, col = 'grey', type = 'h', lwd=7, ylab = 'Probability of resident', main = 'Estimated behaviour', xlab = '', ylim=c(0,1))
  graphics::lines(track$p ~ dateVec, lwd = 2, col = 'red')

  # calculated track
  xl <- c(floor(min(track$lon)), ceiling(max(track$lon)))
  yl <- c(floor(min(track$lat)), ceiling(max(track$lat)))
  graphics::plot(track$lon, track$lat, type = 'n', main = 'Estimated movements', ylab = 'Latitude', xlab = 'Longitude', xlim = xl, ylim = yl)
  graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = "steelblue1")
  fields::world(add = TRUE, fill = TRUE, col = 'grey60')

  if(!is.null(known)){
    graphics::lines(known$lon, known$lat, col = 'white')
  }

  graphics::lines(track$lon, track$lat, col = 'black')

  if(behav.pts){
    na.idx <- which(is.na(sv))
    idx <- which(!is.na(sv))
    #graphics::points(track$lon[na.idx], track$lat[na.idx], bg='white', pch=21)
    graphics::points(track$lon[idx], track$lat[idx], bg='grey60', pch=21)
  }

  graphics::points(track$lon[1], track$lat[1], bg = 'green', pch = 21)
  TT <- length(track$lon)
  graphics::points(track$lon[TT], track$lat[TT], bg = 'red', pch = 21)

  if(save.plot) grDevices::dev.off()

  ## NOT YET FUNCTIONAL
  ## Simple diagnostics plot ###
  ## Resample SST
  #if(save.plot) grDevices::pdf('../plot/sphmmDiagn.pdf',width=7,height=7)
  #if(!save.plot)
  if (resid){
    # do nothing for now
    #grDevices::dev.new()
    #graphics::par(mfrow = c(2,2))
    #ssterr <- sst - lsst$sst
    #sdsst <- sqrt(stats::var(ssterr))
    #ssterr <- ssterr / sdsst
    #lonerr <- sphmm$meanlon - lsst$lon
    #sdlon <- sqrt(stats::var(lonerr))
    #lonerr <- lonerr / sdlon
    #graphics::plot(track$date[ind], ssterr, xlab = 'Date', ylab = 'SST residual', main = 'Residuals through time', ylim = c(-3,3), pch = ".", cex = 3)
    #graphics::abline(h = c(-2,0,2), col = 'grey', lwd = 2, lty = 2)
    #stats::qqnorm(ssterr, main = 'QQ-plot of residuals', pch = ".", cex = 2)
    #graphics::abline(a = 0, b = 1, col = 'grey', lwd = 2)
    #graphics::plot(track$date[ind], lonerr, xlab = 'Date', ylab = 'Longitude residual', ylim = c(-3,3), pch = ".", cex = 3)
    #graphics::abline(h = c(-2,0,2), col = 'grey', lwd = 2, lty = 2)
    #stats::qqnorm(lonerr, main = '', pch = ".", cex = 2)
    #graphics::abline(a = 0, b = 1, col = 'grey', lwd = 2)
  }
  #if(save.plot) grDevices::dev.off()
  }
  
  }
