#' Calculate most probable track from state estimates
#' 
#' \code{calc.track} uses HMM output via \code{hmm.smoother} to calculate most 
#' probable track and behavior state
#' 
#' @param distr is output array from \code{hmm.smoother}
#' @param g is one of the outputs from \code{resample.grid} which denotes what 
#'   spatial scale and grid you're working on
#' @param dateVec is vector of dates from tag to pop-up in 1 day increments.
#' @param method is character indicating what method to use for track 
#'   calculation. Currently only 'mean' and 'max' are supported.
#' @param iniloc is matrix of tag and pop locations. Default is NULL because
#'   this should be taken care of elsewhere.
#'   
#' @return calculated track
#' @export
#' @examples 
#' \dontrun{
#' 
#' # GET THE MOST PROBABLE TRACK
#' tr <- calc.track(s, g, dateVec, iniloc)
#' 
#' }
#' 

calc.track <- function(distr, g, dateVec, iniloc, method = 'mean'){
  ## Calculate track from probability distribution of location
  
  if (method == 'mean'){
    T <- dim(distr)[2]
    # Track calculated from mean
    lat <- apply(apply(distr, c(2, 4), sum) * repmat(t(as.matrix(g$lat[,1])), T, 1), 1, sum)
    lon <- apply(apply(distr, c(2, 3), sum) * repmat(t(as.matrix(g$lon[1,])), T, 1), 1, sum)
    
  } else if (method == 'max'){
    warning('Max method is not recommended and can result in unrealistic movements. Use with caution.')
    
    lat.mx <- apply(apply(distr, c(2, 4), sum), 1, which.max)
    #lat.mx <- apply(apply(distr, c(2, 4), sum) * repmat(t(as.matrix(g$lat[,1])), T, 1), 1, which.max)
    lat <- g$lat[,1][lat.mx]
    lon.mx <- apply(apply(distr, c(2, 3), sum), 1, which.max)
    #lon.mx <- apply(apply(distr, c(2, 3), sum) * repmat(t(as.matrix(g$lon[1,])), T, 1), 1, which.max)
    lon <- g$lon[1,][lon.mx]
    
    
  } else if (method == 'mode'){
    
    stop('Mode is currently unsupported.')
    
    Mode <- function(x) {
      # function from SO at http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    
    lat.md <- apply(apply(distr, c(2, 4), sum) * repmat(t(as.matrix(g$lat[,1])), T, 1), 1, FUN=function(x) Mode(x))
    lon.md <- apply(apply(distr, c(2, 3), sum) * repmat(t(as.matrix(g$lon[1,])), T, 1), 1,FUN=function(x) Mode(x))
    
    
    # Track calculated from mode
    row <- dim(g$lon)[1]
    col <- dim(g$lon)[2]
    modelat <- rep(0, T)
    modelon <- rep(0, T)
    
    for(t in 1:T){
      asd <- apply(distr[,t,,], c(2,3), sum)
      ind <- which.max(asd)
      x <- ceiling(ind / col)
      y <- ind %% row
      modelat[t] <- g$lat[y,x]
      modelon[t] <- g$lon[y,x]
    }
  } else if (method == 'Viterbi'){
    
    stop('Viterbi is currently unsupported.')
    
  }
  
  # calculate the estimated behavior state
  if(dim(distr)[1] > 1){
    p.resid <- apply(distr, c(1,2), sum)[2,]
  } else{
    p.resid <- rep(NA, length.out = length(dateVec))
  }
  
  track <- data.frame(cbind(date = dateVec, lon = lon, lat = lat, p = p.resid))
  track$date <- dateVec
  
  if (any(!c(iniloc$lon[1] - 1 < round(track$lon[1],0) & iniloc$lon[1] + 1 > round(track$lon[1],0),
  iniloc$lat[1] - 1 < round(track$lat[1],0) & iniloc$lat[1] + 1 > round(track$lat[1],0),
  iniloc$lon[2] - 1 < round(track$lon[nrow(track)],0) & iniloc$lon[2] + 1 > round(track$lon[nrow(track)],0),
  iniloc$lat[2] - 1 < round(track$lat[nrow(track)],0) & iniloc$lat[2] + 1 > round(track$lat[nrow(track)],0)))) warning('Known tag and pop-up locations are not within +/- 1deg of final calculated start/stop positions of the trackack. This usually means the input likelihoods were too far from the known locations for the filter/smoother process to realistically get the tag from the likely location(s) to the known location given the movement kernels specified. tracky a different likelihood input.')
  
  return(track)

}
