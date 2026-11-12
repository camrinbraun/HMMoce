#' HMM filter functions (modified for mini-expansion method)
#' 
#' Original expanding kernel code by Dr. Julie Nielsen, integrated by Dr. Martin Arostegui (aka Dr. Mola) into pre-existing hmm.filter function in HMMoce
#' 
#' Cite: Nielsen JK, Bryan DR, Rand KM, Arostegui MC, Braun CD, Galuardi B, McDermott SF. YYYY. Geolocation of a demersal fish (Pacific cod) in a high-latitude island chain (Aleutian Islands, Alaska). Animal Biotelemetry VV: PP-PP.
#' 
#' @param g grid from \code{\link{setup.grid}}
#' @param L is likelihood array output from \code{make.L}
#' @param K is list of length m where each element represents the movement (diffusion) kernel for that behavior state \code{\link{gausskern}}
#' @param sigmas is vector of length m where each element represents the one-dimensional variance for diffusion of the behavior state(s)
#' @param n.div is the number of mini-expansion steps the original diffusion kernel is divided into (default = 10, but test carefully for each use case)
#' @param bathy is the bathymetry raster used in landmasking - must be the same resolution and size as your grid (recommend bathy_resamp)
#' @param mini.idx is a numeric vector that corresponds to the positions of dateVec for which the mini-expansion should be applied.
#' @param P a matrix of dimensions m x m representing the probability for transitions between states.
#' @param m is integer indicating the number of desired behavior states. Currently limited to max 2.
#' @param maskL is logical indicating whether to mask the input L layer. Default is FALSE. See
#'   \code{mask.L} for details.
#' @param bound.thr is numeric indicating the percent threshold that is added 
#'   and subtracted from the bounding box of the filter output from the 
#'   previous day before masking. Default is .05 (5 percent).
#' @param minBounds is size (in degrees) of the minimum bounding box around the 
#'   previous days filter prediction that L data within that box will be 
#'   included. Outside this box (centered on t-1 filter prediction), L will be 
#'   masked out.
#'   
#' @return a list: list(phi = phi, pred = pred, psi = psi) where \itemize{ \item
#'   phi. is the probability for each state at each time step \item pred. is
#'   .... \item psi. is.... }
#' @references Pedersen MW, Patterson TA, Thygesen UH, Madsen H (2011)
#'   Estimating animal behavior and residency from movement data. Oikos
#'   120:1281-1290. doi: 10.1111/j.1600-0706.2011.19044.x
#' @export
#' @examples 
#' \dontrun{
#' # Not run as function relies on large arrays of likelihoods
#' # RUN THE FILTER STEP
#' f <- hmm.filter(g, L, K = list(K1, K2), P = P.final, m = 2)
#' nllf <- -sum(log(f$psi[f$psi>0])) # negative log-likelihood
#' 
#' }
#' 

hmm.filterExp <- function(g, L, K, sigmas, n.div = 10, bathy, P, m, maskL = FALSE, bound.thr = 0.1, minBounds = 10, mini.idx = NULL){
  
  ## Filter data to estimate locations and behaviour
  
  T <- dim(L)[1] # dimension of time 
  row <- dim(g$lon)[1] # nrows
  col <- dim(g$lon)[2] # ncols
  #m <- 2 # Number of behavioural states
  
  pred <- array(0, dim = c(m, T, col, row)) # empty array for prediction step. ordering in col before row emulates lon before lat
  phi  <- array(0, dim = c(m, T, col, row)) # posterior (final) step array
  
  # Start at the known initial location
  phi[1,1,,]  <- L[1,,] # first position is known
  pred[1,1,,] <- L[1,,] # first position is known
  
  for (ii in 1:m){
    # convert movement kernels from matrix to cimg for convolution
    K[[ii]] <- imager::as.cimg(K[[ii]])
  }
  
  ## set up the kernel(s) for the mini expansion
  # Miniaturize the kernel(s)
  if (is.null(n.div)) n.div <- 10
  print(paste0('Using ', n.div, ' divisions for the mini expansion.'))
  if (m == 1) {
    D1exp <- (sigmas ^ 2) / n.div
    sigmas.exp <- sqrt(D1exp)
    sizes.exp = rep(ceiling(sigmas.exp[1] * 4),2)
    muadvs.exp = c(0)
    K1exp <- gausskern.pg(sizes.exp[1], sigmas.exp[1], muadv = muadvs.exp[1])
    Kexp <- list(K1exp)
  } else if (m == 2) {
    D1exp <- (sigmas[1] ^ 2) / n.div
    D2exp <- (sigmas[2] ^ 2) / n.div
    sigmas.exp <- sqrt(c(D1exp, D2exp))
    sizes.exp = rep(ceiling(sigmas.exp[1] * 4), 2)
    muadvs.exp = c(0, 0)
    K1exp <- gausskern.pg(sizes.exp[1], sigmas.exp[1], muadv = muadvs.exp[1])
    K2exp <- gausskern.pg(sizes.exp[2], sigmas.exp[2], muadv = muadvs.exp[2])
    Kexp <- list(K1exp, K2exp)
  }
  
  for (ii in 1:m){
    # convert movement mini-kernels from matrix to cimg for convolution
    Kexp[[ii]] <- imager::as.cimg(Kexp[[ii]])
  }
  
  psi <- rep(0, T - 1) # sum of the probability of states at each step
  
  ## create empty lists for upcoming vars
  p <- list(); q <- list(); post <- list()
  
  ## reassign bathymetry raster as matrix for landmasking
  if (!inherits(bathy, "SpatRaster")) bathy <- terra::rast(bathy)
  bathy <- reverse(terra::as.matrix(bathy, wide = TRUE)) # reverse to match q[[ii]] format below
  
  # Start filter iterations
  for(t in 2:T){
    
    for (ii in 1:m){
      
      # convolve previous day's likelihood with movement kernels
      p[[ii]] = imager::as.cimg(t(phi[ii, t-1,,]))
      
      if (t %in% mini.idx){
        # Mini-expansions within a timestep
        expansion <- 1
        while (expansion <= n.div) {
          q[[ii]] = t(as.matrix(imager::convolve(p[[ii]], Kexp[[ii]])))
          #q1 = imager::convolve(p1, K1)
          #q2 = imager::convolve(p2, K2)
          #q1 = t(as.matrix(q1))
          #q2 = t(as.matrix(q2))
          q[[ii]][bathy >= 0] <- 0 # landmask each mini-expansion
          p[[ii]] <- imager::as.cimg(t(q[[ii]]))
          expansion = expansion+1
        }
        
      } else{
        q[[ii]] = t(as.matrix(imager::convolve(p[[ii]], K[[ii]])))
      }
      
    }
    
    # multiply by transition probability 
    if (m == 1){
      pred[1,t,,] <- q[[1]]
      
    } else if (m == 2){
      pred[1,t,,] <- P[1,1] * q[[1]] + P[2,1] * q[[2]]
      pred[2,t,,] <- P[1,2] * q[[1]] + P[2,2] * q[[2]]
      
    } else if (m > 2){
      stop('More than 2 behavior states is not currently supported.')
    }
    
    for (ii in 1:m){
      # is there a data-based likelihood observation for this day, t?
      sumL = sum(L[t,,])  
      if(sumL > 1e-6){
        if(maskL){
          post[[ii]] <- mask.L(pred.t = pred[ii,t,,], L.t = L[t,,], lon = g$lon[1,], lat = g$lat[,1], par0 = dim(K[[ii]])[1], bound.thr = bound.thr, minBounds=minBounds)
          #post1 <- mask.L(pred.t = pred[1,t,,], L.t = L[t,,], lon = g$lon[1,], lat = g$lat[,1], par0 = dim(K1)[1], bound.thr = bound.thr, minBounds=minBounds)
          #post2 <- mask.L(pred.t = pred[2,t,,], L.t = L[t,,], lon = g$lon[1,], lat = g$lat[,1], par0 = dim(K1)[1], bound.thr = bound.thr, minBounds=minBounds)
        } else{
          post[[ii]] <- pred[ii,t,,] * L[t,,]
          #post1 <- pred[1,t,,] * L[t,,]
          #post2 <- pred[2,t,,] * L[t,,]
        }
      }else{
        post[[ii]] <- pred[ii,t,,]
        #post1 <- pred[1,t,,]
        #post2 <- pred[2,t,,]
      }
    }
    
    psi[t-1] <- sum(unlist(lapply(post, FUN = function(x) sum(x, na.rm = T))))
    #psi[t-1] <- sum(as.vector(post[[1]]), na.rm=T) + sum(as.vector(post[[2]]), na.rm=T)
    
    for (ii in 1:m){
      phi[ii,t,,] <- post[[ii]] / (psi[t-1] + 1e-15)
    }
    #phi[1,t,,] <- post1 / (psi[t-1] + 1e-15)
    #phi[2,t,,] <- post2 / (psi[t-1] + 1e-15)
    
  } # t loop
  
  list(phi = phi, pred = pred, psi = psi)
  
}
