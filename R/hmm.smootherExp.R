#' Smoother recursion over filtered state estimates
#' 
#' \code{hmm.smoother} provides backward (starting at end) recursion over 
#' filtered state estimates as output from \code{hmm.filter}. The product of 
#' this function an array containing final state estimates.
#' 
#' @param f is array output from \code{hmm.filter}
#' @param L is likelihood array output from \code{make.L}
#' @param K is list of length m where each element represents the movement (diffusion) kernel for that behavior state \code{\link{gausskern}}
#' @param P a matrix of dimensions m x m representing the probability for transitions between states.
#'   
#' @return an array of the final state estimates of dim(state, time, lon, lat)
#' @references Pedersen MW, Patterson TA, Thygesen UH, Madsen H (2011)
#'   Estimating animal behavior and residency from movement data. Oikos
#'   120:1281-1290. doi: 10.1111/j.1600-0706.2011.19044.x
#' @examples 
#' \dontrun{
#' # Not run as function relies on large arrays of likelihoods
#' # RUN THE SMOOTHING STEP
#' s <- hmm.smoother(f, K1, K2, L, P.final)
#' }
#' 
#' @export
#' 

hmm.smootherExp <- function(f, L, K, P, bathy, n.div, mini.idx, sigmas){
  
  m <- dim(f$phi)[1]
  
  if (m == 1){
    
    s <- hmm.smoother1.exp(f=f, K=K, L=L, P=NULL, bathy, n.div, mini.idx, m=m, sigmas=sigmas)
    
  } else if (m == 2){
    
    s <- hmm.smoother2.exp(f=f, K1=K[[1]], K2=K[[2]], L=L, P=P, bathy, n.div, mini.idx, m=m, sigmas=sigmas)
    
  } 
  
  return(s)
  
}

hmm.smoother1.exp <- function(f, L, K, P=NULL, bathy, n.div, mini.idx, m=m, sigmas=sigmas){
  ## Smoothing the filtered estimates
  ## The equations for smoothing are presented in Pedersen et al. 2011, Oikos, Appendix
  m <- dim(f$phi)[1]
  T <- dim(f$phi)[2]
  row <- dim(f$phi)[3]
  col <- dim(f$phi)[4]
  
  ## set up the kernel(s) for the mini expansion
  # Miniaturize the kernel(s)
  if (is.null(n.div)) n.div <- 10
  print(paste0('Using ', n.div, ' divisions for the mini expansion.'))
    
  D1exp <- (sigmas ^ 2) / n.div
  sigmas.exp <- sqrt(D1exp)
  sizes.exp = rep(ceiling(sigmas.exp[1] * 4),2)
  muadvs.exp = c(0)
  K1exp <- gausskern.pg(sizes.exp[1], sigmas.exp[1], muadv = muadvs.exp[1])
  Kexp <- list(K1exp)
  
  for (ii in 1:m){
    # convert movement kernels from matrix to cimg for convolution
    Kexp[[ii]] <- imager::as.cimg(Kexp[[ii]])
    K[[ii]] <- imager::as.cimg(K[[ii]])
  }

  ## reassign bathymetry raster as matrix for landmasking
  if (!inherits(bathy, "SpatRaster")) bathy <- terra::rast(bathy)
  bathy <- reverse(terra::as.matrix(bathy, wide = TRUE)) # reverse to match q[[ii]] format below
  
  smooth <- array(0, dim = dim(f$phi))
  smooth[,T,,] <- f$phi[,T,,]
  
  ## create empty lists for upcoming vars
  p <- list(); Rp <- list(); post <- list()
  
  for(t in T:2){
    
    RAT <- smooth[,t,,] / (f$pred[,t,,] + 1e-15)
    ## preserve the same array dims even if m == 1
    if (m == 1) RAT <- array(RAT, dim = c(1L, dim(RAT)))

    for (ii in 1:m){
      
      # convolve today's smoother prediction with movement kernel
      p[[ii]] = imager::as.cimg(t(RAT[ii,,]))
      
      if (t %in% mini.idx){
        # Mini-expansions within a timestep
        expansion <- 1
        while (expansion <= n.div) {
          Rp[[ii]] = t(as.matrix(imager::convolve(p[[ii]], Kexp[[ii]])))
          Rp[[ii]][bathy >= 0] <- 0 # landmask each mini-expansion
          p[[ii]] <- imager::as.cimg(t(Rp[[ii]]))
          expansion = expansion+1
        }
        
      } else{
        Rp[[ii]] <- t(as.matrix(imager::convolve(p[[ii]], K[[ii]])))
      }
          
    }

    # multiply by transition probability 
    if (m == 1){
      post[[ii]] <- matrix(Rp[[1]], row, col)
      
    } else if (m == 2){
      post[[1]] <- matrix(P[1,1] * Rp[[1]] + P[1,2] * Rp[[2]], row, col)
      post[[2]] <- matrix(P[2,1] * Rp[[1]] + P[2,2] * Rp[[2]], row, col)
       
    } else if (m > 2){
      stop('More than 2 behavior states is not currently supported.')
    } else{
      print('Skipping mult trans prob')
    }

    if (T == t){
      for (ii in 1:m){
        post[[ii]] <- f$phi[1,t,,] * 0
        if (ii == 1) post[[1]] <- L[t,,]
        fac <- sum(unlist(lapply(post, FUN = function(x) sum(x, na.rm = T))))
        if (fac < 1e-200) fac <- 1e-200 ## added to prevent fac == 0
        smooth[ii,t,,] <- post[[ii]] / fac
        post[[ii]] <- post[[ii]] * f$phi[ii,t-1,,]
        fac <- sum(unlist(lapply(post, FUN = function(x) sum(x, na.rm = T))))
        if (fac < 1e-200) fac <- 1e-200 ## added to prevent fac == 0
        smooth[ii,t-1,,] <- post[[ii]] / fac
        
      }
        
    } else{
      for (ii in 1:m){
        
        post[[ii]] <- post[[ii]] * f$phi[ii,t-1,,]
        fac <- sum(unlist(lapply(post, FUN = function(x) sum(x, na.rm = T))))
        smooth[ii,t-1,,] <- post[[ii]] / fac

      }
      
    }
  }
  
  return(smooth)
  
}


hmm.smoother2.exp <- function(f, K1, K2, L, P, bathy, n.div, mini.idx, m=m, sigmas=sigmas){
  ## Smoothing the filtered estimates
  ## The equations for smoothing are presented in Pedersen et al. 2011, Oikos, Appendix
  T <- dim(f$phi)[2]
  row <- dim(f$phi)[3]
  col <- dim(f$phi)[4]
  
  # convert movement kernel from matrix to cimg for convolution
  K1 <- imager::as.cimg(K1)
  K2 <- imager::as.cimg(K2)
  
  ## reassign bathymetry raster as matrix for landmasking
  if (!inherits(bathy, "SpatRaster")) bathy <- terra::rast(bathy)
  bathy <- reverse(terra::as.matrix(bathy, wide = TRUE)) # reverse to match q[[ii]] format below
  
    D1exp <- (sigmas[1] ^ 2) / n.div
    D2exp <- (sigmas[2] ^ 2) / n.div
    sigmas.exp <- sqrt(c(D1exp, D2exp))
    sizes.exp = rep(ceiling(sigmas.exp[1] * 4), 2)
    muadvs.exp = c(0, 0)
    Kexp1 <- gausskern.pg(sizes.exp[1], sigmas.exp[1], muadv = muadvs.exp[1])
    Kexp2 <- gausskern.pg(sizes.exp[2], sigmas.exp[2], muadv = muadvs.exp[2])
    Kexp1 <- imager::as.cimg(Kexp1)
    Kexp2 <- imager::as.cimg(Kexp2)
    #Kexp <- list(K1exp, K2exp)

    #for (ii in 1:m){
    #  # convert movement kernels from matrix to cimg for convolution
    #  Kexp[[ii]] <- imager::as.cimg(Kexp[[ii]])
    #  #K[[ii]] <- imager::as.cimg(K[[ii]])
    #}
    
  smooth <- array(0, dim = dim(f$phi))
  smooth[,T,,] <- f$phi[,T,,]
  
  #smooth <- f$phi  #default; fill in as the prediction step.
  
  for(t in T:2){
    RAT <- smooth[,t,,] / (f$pred[,t,,] + 1e-15)
    
    p1 = imager::as.cimg(t(RAT[1,,]))
    p2 = imager::as.cimg(t(RAT[2,,]))
    
    if (t %in% mini.idx){
      # Mini-expansions within a timestep
      expansion <- 1
      while (expansion <= n.div) {
        Rp1 = t(as.matrix(imager::convolve(p1, Kexp1)))
        Rp1[bathy >= 0] <- 0 # landmask each mini-expansion
        p1 <- imager::as.cimg(t(Rp1))
        
        Rp2 = t(as.matrix(imager::convolve(p2, Kexp2)))
        Rp2[bathy >= 0] <- 0 # landmask each mini-expansion
        p2 <- imager::as.cimg(t(Rp2))
        
        expansion = expansion+1
      }
      
    } else{
      # convolve today's smoother prediction with movement kernel
      Rp1 <- imager::convolve(p1, K1)
      Rp2 <- imager::convolve(p2, K2)
      Rp1 = t(as.matrix(Rp1))
      Rp2 = t(as.matrix(Rp2))
      
    }
    
    post1 <- matrix(P[1,1] * Rp1 + P[1,2] * Rp2, row, col)
    post2 <- matrix(P[2,1] * Rp1 + P[2,2] * Rp2, row, col)
    
    if(T == t){
      post1 <- f$phi[1,t,,] * 0
      post2 <- L[t,,]
      fac <- sum(as.vector(post1)) + sum(as.vector(post2))
      smooth[1,t,,] <- post1 / fac
      smooth[2,t,,] <- post2 / fac 
      post1 <- post1 * f$phi[1,t-1,,]
      post2 <- post2 * f$phi[2,t-1,,]
      fac <- sum(as.vector(post1)) + sum(as.vector(post2))
      smooth[1,t-1,,] <- post1 / fac
      smooth[2,t-1,,] <- post2 / fac
    }else{
      post1 <- post1 * f$phi[1,t-1,,]
      post2 <- post2 * f$phi[2,t-1,,]
      fac <- sum(as.vector(post1)) + sum(as.vector(post2))
      smooth[1,t-1,,] <- post1 / fac
      smooth[2,t-1,,] <- post2 / fac
    }
  }
  
  return(smooth)
  
}



