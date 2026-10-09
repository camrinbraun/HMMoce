#' Create coarse overall likelihood for efficient parameter estimation
#' 
#' \code{coarse.L} makes a more coarse (aggregated) version of the overall likelihood array (output from \code{make.L}). This results in a coarse version of the overall likelihood grid used for more efficient parameter estimation.
#' @param L is overall likelihood array output from \code{make.L}
#' @param ras.list a list of likelihood rasters as input to \code{make.L}. Can also be just one of these raster layers rather than of class "list".
#' @param aggregate_fact is \code{fact} input to \code{raster::aggregate}. 
#' @export

coarse.L <- function(L, ras.list, aggregate_fact = NULL){
  
  ## get example grid
  if (class(ras.list) == 'list') r <- ras.list[[1]]
  if (class(ras.list) == "RasterBrick" | class(ras.list) == "RasterStack" | class(ras.list) == "RasterLayer") r <- ras.list
  
  ## convert L back to raster
  L.mle <- aperm(L, c(3,2,1))
  lon.agg <- seq(raster::extent(r)[1], raster::extent(r)[2], length.out=dim(r)[2])
  lat.agg <- seq(raster::extent(r)[3], raster::extent(r)[4], length.out=dim(r)[1])
  list.mle <- list(x = lon.agg, y = lat.agg, z = L.mle)
  ex <- raster::extent(list.mle)
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  L.mle <- raster::brick(list.mle$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=F, crs)
  L.mle <- raster::flip(L.mle, direction = 'y')
  
  ## coarsen / aggregate
  if (is.null(aggregate_fact)) aggregate_fact <- 5
  
  L.mle <- raster::aggregate(L.mle, fact = aggregate_fact, fun = 'max')
  g.mle <- setup.grid.raster(L.mle)
  
  ## output array for parameter optim
  L.mle <- aperm(raster::as.array(raster::flip(L.mle, direction = 'y')), c(3, 2, 1))
  
  return(list(L.mle = L.mle, g.mle = g.mle))

}


