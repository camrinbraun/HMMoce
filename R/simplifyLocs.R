#' Simplifies locations to one per day
#' 
#' @param locations is data frame from -Locations file output from DAP/Tag Portal for WC tags and 
#'   contains GPS, Argos, and GPE locations as applicable.
#' @param loc.dts is vector of dates from locations input
#'   
#' @return simplified (one per day) version of locations input
#'

simplifyLocs <- function(locations, loc.dts){
  
  udates <- unique(loc.dts)
  
  for(tt in 1:length(udates)){
    # identify duplicate dates
    idx <- which(loc.dts == udates[tt])
    locs.i <- locations[idx,]
    
    # reconcile duplicates
    if('GPE.MSD' %in% names(locs.i)){
      # reconcile duplicate light data by selecting the position estimate with minimum MSD
      locs.i <- locs.i[which(locs.i$GPE.MSD == min(locs.i$GPE.MSD)),]
      
    } else if('Fastloc.Power' %in% names(locs.i)){
      # reconcile duplicate GPS data by selecting the position with smallest residual
      locs.i <- locs.i[which(locs.i$Residual == min(locs.i$Residual)),]
      
    } else{
      stop('Error: unable to determine type of input location data containing duplicate dates and thus needs thinning. See simplifyLocs function.')
    }
    
    if(tt == 1){
      locs.final <- locs.i
    } else{
      locs.final <- rbind(locs.final, locs.i)
    }
    
  }
  
  # output simplified locations file
  return(list(locs = locs.final, locDates = udates))
  
}