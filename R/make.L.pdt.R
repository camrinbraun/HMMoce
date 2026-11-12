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

make.L.pdt <- function(ras.list, dateVec){
  r_tmpl <- if (inherits(ras.list[[1]], "SpatRaster")) ras.list[[1]] else terra::rast(ras.list[[1]])
  nr <- terra::nrow(r_tmpl)
  nc <- terra::ncol(r_tmpl)
  T_len <- length(dateVec)
  
  arr_list <- lapply(ras.list, function(r) {
    sr <- if (inherits(r, "SpatRaster")) r else terra::rast(r)
    terra::as.array(sr)
  })
  
  L_arr <- array(0, dim = c(nr, nc, T_len))
  for (i in seq_len(T_len)){
    active_mats <- list()
    for (bb in seq_along(arr_list)){
      mat <- arr_list[[bb]][,,i]
      mv <- suppressWarnings(max(mat, na.rm = TRUE))
      if (is.finite(mv) && mv > 0) {
        mat[is.na(mat)] <- 0
        active_mats[[length(active_mats) + 1L]] <- mat / mv
      }
    }
    if (length(active_mats) == 0L) next
    s_sum <- Reduce(`+`, active_mats)
    m_sum <- suppressWarnings(max(s_sum, na.rm = TRUE))
    if (is.finite(m_sum) && m_sum > 0) {
      L_arr[,,i] <- s_sum / m_sum
    }
  }
  
  L <- terra::rast(L_arr, extent = terra::ext(r_tmpl), crs = terra::crs(r_tmpl))
  names(L) <- as.character(dateVec)
  print('Finishing make.L.pdt...', sep='')
  return(L)
}

