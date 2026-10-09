#' Plot RD results
#'
#' \code{plotRD} uses HMM output via \code{calc.track} to calculate and plot residency distributions for each behavior state and combined data
#'
#' @param distr is output array from \code{hmm.smoother}
#' @param track is output dataframe from \code{calc.track}
#' @param known is 3 column data frame containing date, lat, lon of known
#'   movement track. This is only useful for comparing HMMoce results to known
#'   track collected by SPOT or GPS, for example. Default is NULL.
#' @param g grid from \code{\link{setup.grid}}
#' @param xlims a vector of length 2 indicating longitude limits (-180 to 180).
#' @param ylims a vector of length 2 indicating latitude limits
#' @param makePlot is logical indicating whether to make a plot of the RD
#' @param save.plot is logical indicating whether you want the plot written to
#'   disk using \code{pdf}.
#' @param filename is output filename if save.plot is TRUE. Default is NULL.
#' @importFrom fields "world"
#'
#' @return a list of one raster layer for the combined RD and one raster brick (2 layers) for the individual behavior RDs.
#' @export

plotRD <- function(distr, track, known=NULL, g, xlims, ylims, makePlot = TRUE, save.plot=FALSE, filename = NULL){
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  rd.cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdYlGn')))
  
  if (is.null(xlims)) xlims <- c(min(g$lon) - 0.1 * min(g$lon), max(g$lon) + 0.1 * max(g$lon))
  if (is.null(ylims)) ylims <- c(min(g$lat) - 0.1 * min(g$lat), max(g$lat) + 0.1 * max(g$lat))

  #### One-state
  if (dim(distr)[1] == 1){
    # generate summed across time pdf
    norm <- array(NA, dim=dim(distr)[c(3,4,2,1)])
    for (i in 1:dim(norm)[3]){
      norm[,,i,1] <- (distr[1,i,,] / max(distr[1,i,,], na.rm=T))
    }
    all <- apply(norm, 1:2, FUN=function(x) sum(x, na.rm=T))
    all <- all / max(all, na.rm=T)
    # get pdf as raster
    all <- raster::flip(raster::raster(t(all), xmn=min(g$lon), xmx=max(g$lon),
                                       ymn=min(g$lat), ymx=max(g$lat), crs=crs), 2)
    
    if(makePlot){
      # set up plot brks and cols
      #norm.rnge <- cellStats(norm, 'range')
      norm.breaks = seq(.001, 1, length.out=201)
      norm.mid = norm.breaks[1:(length(norm.breaks)-1)]
      plot.rd.col <- rd.cols(length(norm.breaks)-1)
      
      # set axes
      x.at <- pretty(xlims, 5)
      x.labels <- parse(text=paste(abs(x.at), "*degree~W", sep=""))
      y.at <- pretty(ylims, 5)
      y.labels <- parse(text=paste(abs(y.at), "*degree~N", sep=""))
      
      # set margins
      old.mar <- par()$mar
      mar.default <- par('mar')
      mar.bottom <- mar.default
      mar.bottom[1] <- 0
      
      # get some plot indices
      TT <- length(track$lon)
      
      if(save.plot) grDevices::pdf(paste(filename, '_RD.pdf', sep = ''), width = 7, height = 12)
      
      # BUILD IT
      nf <- layout(matrix(c(1,2), 1, 2, byrow=T), widths=c(7,2), heights=4)
      #layout.show(nf)
      
      # COMBINED PANEL
      #par (mar=mar.bottom)
      par(mar=c(2,2,1,1))
      raster::image(all, maxpixels=raster::ncell(all), xlim=xlims, ylim=ylims, axes=F,
                    col=plot.rd.col, breaks=norm.breaks, xlab='', ylab='')#, main='All')#, zlim=c(.05,1))#, axes=F)#, main=dt.idx) #, breaks=zbreaks
      axis(1, at=x.at, labels=x.labels)
      axis(2, at=y.at, labels=y.labels);
      world(add=T, fill=T, col='grey60')
      box()
      if(!is.null(known)) lines(known$lon, known$lat)
      
      # COLORBAR
      image(1, norm.mid, t(as.matrix(norm.mid)), breaks=norm.breaks, col=plot.rd.col, axes=FALSE, xlab="",
            ylab='Expected Proportion of Time Spent')
      axis(2, at=c(.05, .25, .5, .75, .95), labels=paste(rev(c(5, 25, 50, 75, 95)), '%'));box();
      
      #specific to one-state
      behavRD_input <- NULL
    }
  }
  
  #### Two-state
  if (dim(distr)[1] == 2){
    p.1 <- apply(distr[1,,,], 1, sum)
    p.2 <- apply(distr[2,,,], 1, sum)
    sv <- -(apply(distr[1,,,], 1, sum) > apply(distr[2,,,], 1, sum)) + 1
    #sv[sv == 0] <- NA
    
    # normalize for each behav at each time point
    norm <- array(NA, dim=dim(distr)[c(3,4,2,1)])
    for (i in 1:dim(norm)[3]){
      norm[,,i,1] <- (distr[1,i,,] / max(distr[1,i,,], na.rm=T)) * p.1[i]
      norm[,,i,2] <- (distr[2,i,,] / max(distr[2,i,,], na.rm=T)) * p.2[i]
    }
    
    all <- apply(norm, 1:2, FUN=function(x) sum(x, na.rm=T))
    all <- all / max(all, na.rm=T)
    
    # sum across time and normalize final surface
    norm <- apply(norm, c(1:2,4), FUN=function(x) sum(x, na.rm=T))
    norm[,,1] <- norm[,,1] / max(norm[,,1], na.rm=T)
    norm[,,2] <- norm[,,2] / max(norm[,,2], na.rm=T)
    
    # get pdf as raster
    norm <- raster::flip(raster::brick(norm, xmn=min(g$lon), xmx=max(g$lon),
                                       ymn=min(g$lat), ymx=max(g$lat), crs=crs, transpose=T), 2)
    all <- raster::flip(raster::raster(t(all), xmn=min(g$lon), xmx=max(g$lon),
                                       ymn=min(g$lat), ymx=max(g$lat), crs=crs), 2)
    
    if(makePlot){
      # set up plot brks and cols
      #norm.rnge <- cellStats(norm, 'range')
      norm.breaks = seq(.001, 1, length.out=201)
      norm.mid = norm.breaks[1:(length(norm.breaks)-1)]
      plot.rd.col <- rd.cols(length(norm.breaks)-1)
      
      # set axes
      x.at <- pretty(xlims, 5)
      x.labels <- parse(text=paste(abs(x.at), "*degree~W", sep=""))
      y.at <- pretty(ylims, 5)
      y.labels <- parse(text=paste(abs(y.at), "*degree~N", sep=""))
      
      # set margins
      old.mar <- par()$mar
      mar.default <- par('mar')
      mar.bottom <- mar.default
      mar.bottom[1] <- 0
      
      # get some plot indices
      idx1 <- which((sv+1) == 1)
      idx2 <- which((sv+1) == 2)
      TT <- length(track$lon)
      
      if(save.plot) grDevices::pdf(paste(filename, '_RD.pdf', sep = ''), width = 7, height = 12)
      
      # BUILD IT
      nf <- layout(matrix(c(1,2,
                            3,4,
                            5,6), 3, 2, byrow=T), widths=c(7,2), heights=c(4,4,4))
      #layout.show(nf)
      
      # COMBINED PANEL
      par(mar = mar.bottom)
      
      raster::image(all, maxpixels=raster::ncell(all), xlim=xlims, ylim=ylims, axes=F,
                    col=plot.rd.col, breaks=norm.breaks, xlab='', ylab='')#, main='All')#, zlim=c(.05,1))#, axes=F)#, main=dt.idx) #, breaks=zbreaks
      axis(1, at=x.at, labels=FALSE)
      axis(2, at=y.at, labels=y.labels);
      world(add=T, fill=T, col='grey60')
      box()
      if(!is.null(known)) lines(known$lon, known$lat)
      
      # COLORBAR
      image(1, norm.mid, t(as.matrix(norm.mid)), breaks=norm.breaks, col=plot.rd.col, axes=FALSE, xlab="",
            ylab='Expected Proportion of Time Spent')
      axis(2, at=c(.05, .25, .5, .75, .95), labels=paste(rev(c(5, 25, 50, 75, 95)), '%'));box();
      
      # MIGRATORY PANEL
      par(mar = mar.bottom)
      raster::image(norm[[1]], maxpixels=raster::ncell(norm[[1]]), xlim=xlims, ylim=ylims, axes=F,
                    col=plot.rd.col, breaks=norm.breaks, xlab='', ylab='')#, main='Migratory')#, zlim=c(.05,1))#, axes=F)#, main=dt.idx) #, breaks=zbreaks
      axis(1, at=x.at, labels=FALSE)
      axis(2, at=y.at, labels=y.labels)
      world(add=T, fill=T, col='grey60')
      box()
      mtext('Migratory behavior', 4, line=2)
      graphics::points(track$lon[idx1], track$lat[idx1], bg='white', pch=21, cex=1.25)
      graphics::points(track$lon[1], track$lat[1], bg = 'green', pch = 21, cex=1.75)
      graphics::points(track$lon[TT], track$lat[TT], bg = 'red', pch = 21, cex=1.75)
      
      # FILL COLORBAR SPACE
      plot.new()
      
      # RESIDENT PANEL
      par (mar=mar.default)
      raster::image(norm[[2]], maxpixels=raster::ncell(norm[[2]]), xlim=xlims, ylim=ylims, axes=F,
                    col=plot.rd.col, breaks=norm.breaks, xlab='', ylab='')#, main='Resident')#, zlim=c(.05,1))#, axes=F)#, main=dt.idx) #, breaks=zbreaks
      axis(1, at=x.at, labels=x.labels)
      axis(2, at=y.at, labels=y.labels);
      world(add=T, fill=T, col='grey60')
      box()
      mtext('Resident behavior', 4, line=2)
      graphics::points(track$lon[idx2], track$lat[idx2], bg='grey60', pch=21, cex=1.25)
      graphics::points(track$lon[1], track$lat[1], bg = 'green', pch = 21, cex=1.75)
      graphics::points(track$lon[TT], track$lat[TT], bg = 'red', pch = 21, cex=1.75)
      
      # FILL COLORBAR SPACE
      plot.new()
      
      if(save.plot) grDevices::dev.off()
      
      #specific to two-state
      behavRD_input <- norm
    }
  }
  
  return(list(allRD = all, behavRD = behavRD_input))
}
