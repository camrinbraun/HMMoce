#' Combine individual source PDT likelihoods
#' 
#' \code{make.L.pdt} combines individual PDT-based likelihoods from various data sources (e.g.
#' HYCOM OHC, HYCOM depth-levels, GLORYS OHC, GLORYS depth-levels) to make overall 
#' combined PDT-based likelihoods for each time point
#' 
#' @param ras.list a list of likelihood rasters
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @return an overall PDT-based likelihood RasterStack
#' 
#' @examples
#' \dontrun{
#' L <- make.L(ras.list,dateVec)
#' }
#' @export
#'   
#' @author Martin C. Arostegui

make.L.pdt <- function(ras.list,dateVec){
  
  ## generate blank results likelihood
  L <- ras.list[[1]] * 0
  L[is.na(L)] <- 0
  
  ## COMBINE THE LIKELIHOOD LAYERS
  ## for each day:
  for (i in 1:length(dateVec)){
    
    #print(i)
    
    ## get relevant likelihoods across the list
    for (bb in 1:length(ras.list)){
      if (bb == 1){
        s <- raster::stack(ras.list[[bb]][[i]])
      } else{
        s <- raster::stack(s, ras.list[[bb]][[i]])
      }
    }
    
    ## check for layers with all NA or all <= 0
    sum_NA <- raster::cellStats(!is.na(s), sum, na.rm = TRUE) == 0
    max_vals <- suppressWarnings(raster::cellStats(s, 'max', na.rm = TRUE))
    empty_layer <- sum_NA | is.na(max_vals) | is.infinite(max_vals) | (max_vals <= 0)
    
    ## if all layers are empty, skip to next day
    if (all(empty_layer)){
      next
    } else if (any(empty_layer)){
      s <- s[[-which(empty_layer)]]
    }
    
    ## sum & normalize whatever layers remain
    if (raster::nlayers(s) == 1){
      L[[i]] <- s / raster::cellStats(s, 'max', na.rm = TRUE)
    } else{
      s_sum <- sum(s, na.rm = TRUE)
      L[[i]] <- s_sum / raster::cellStats(s_sum, 'max', na.rm = TRUE)
    }
    
  }

  print('Finishing make.L.pdt...', sep='')
  
  return(L)
}

