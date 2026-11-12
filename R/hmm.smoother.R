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

hmm.smoother <- function(f, L, K, P){
  
  m <- dim(f$phi)[1]
  
  if (m == 1){
    
    s <- hmm.smoother1(f=f, K=K, L=L, P=NULL)
    
  } else if (m == 2){
    
    s <- hmm.smoother2(f=f, K1=K[[1]], K2=K[[2]], L=L, P=P)
    
  } 
  
  return(s)
  
}

hmm.smoother1 <- function(f, L, K, P=NULL){
  ## Smoothing the filtered estimates
  ## The equations for smoothing are presented in Pedersen et al. 2011, Oikos, Appendix
  m <- dim(f$phi)[1]
  T <- dim(f$phi)[2]
  row <- dim(f$phi)[3]
  col <- dim(f$phi)[4]
  
  for (ii in 1:m){
    # convert movement kernels from matrix to cimg for convolution
    K[[ii]] <- imager::as.cimg(K[[ii]])
  }

  smooth <- array(0, dim = dim(f$phi))
  sum_phi_T <- sum(f$phi[,T,,], na.rm = TRUE)
  if (is.finite(sum_phi_T) && sum_phi_T > 0) {
    smooth[,T,,] <- f$phi[,T,,] / sum_phi_T
  } else {
    sum_LT <- max(sum(L[T,,], na.rm = TRUE), 1e-300)
    for (ii in 1:m) smooth[ii,T,,] <- L[T,,] / (m * sum_LT)
  }
  
  ## create empty lists for upcoming vars
  p <- list(); Rp <- list(); post <- list()
  
  for(t in T:2){
    
    RAT <- smooth[,t,,] / pmax(f$pred[,t,,], 1e-300)
    ## preserve the same array dims even if m == 1
    if (m == 1) RAT <- abind::abind(RAT, RAT, along = 0)[1,,, drop=F]

    for (ii in 1:m){
      
      # convolve today's smoother prediction with movement kernel
      p[[ii]] = imager::as.cimg(t(RAT[ii,,]))
      Rp[[ii]] <- t(as.matrix(imager::convolve(p[[ii]], K[[ii]])))
          
    }

    # multiply by transition probability 
    if (m == 1){
      post[[1]] <- matrix(Rp[[1]], row, col)
      
    } else if (m == 2){
      post[[1]] <- matrix(P[1,1] * Rp[[1]] + P[1,2] * Rp[[2]], row, col)
      post[[2]] <- matrix(P[2,1] * Rp[[1]] + P[2,2] * Rp[[2]], row, col)
       
    } else if (m > 2){
      stop('More than 2 behavior states is not currently supported.')
    } else{
      print('Skipping mult trans prob')
    }

    for (ii in 1:m){
      post[[ii]] <- post[[ii]] * f$phi[ii,t-1,,]
    }
    fac <- sum(vapply(post, sum, numeric(1), na.rm = TRUE))
    if (is.finite(fac) && fac > 0) {
      for (ii in 1:m){
        smooth[ii,t-1,,] <- post[[ii]] / fac
      }
    } else {
      sum_phi_prev <- sum(f$phi[,t-1,,], na.rm = TRUE)
      if (is.finite(sum_phi_prev) && sum_phi_prev > 0) {
        smooth[,t-1,,] <- f$phi[,t-1,,] / sum_phi_prev
      }
    }
  }
  
  return(smooth)
  
}


hmm.smoother2 <- function(f, K1, K2, L, P){
  ## Smoothing the filtered estimates
  ## The equations for smoothing are presented in Pedersen et al. 2011, Oikos, Appendix
  T <- dim(f$phi)[2]
  row <- dim(f$phi)[3]
  col <- dim(f$phi)[4]
  
  # convert movement kernel from matrix to cimg for convolution
  K1 <- imager::as.cimg(K1)
  K2 <- imager::as.cimg(K2)
  
  smooth <- array(0, dim = dim(f$phi))
  sum_phi_T <- sum(f$phi[,T,,], na.rm = TRUE)
  if (is.finite(sum_phi_T) && sum_phi_T > 0) {
    smooth[,T,,] <- f$phi[,T,,] / sum_phi_T
  } else {
    sum_LT <- max(sum(L[T,,], na.rm = TRUE), 1e-300)
    smooth[1,T,,] <- L[T,,] / (2 * sum_LT)
    smooth[2,T,,] <- L[T,,] / (2 * sum_LT)
  }
  
  #smooth <- f$phi  #default; fill in as the prediction step.
  
  for(t in T:2){
    RAT <- smooth[,t,,] / pmax(f$pred[,t,,], 1e-300)
    
    # convolve today's smoother prediction with movement kernel
    p1 = imager::as.cimg(t(RAT[1,,]))
    Rp1 <- imager::convolve(p1, K1)
    p2 = imager::as.cimg(t(RAT[2,,]))
    Rp2 <- imager::convolve(p2, K2)
    Rp1 = t(as.matrix(Rp1))
    Rp2 = t(as.matrix(Rp2))
    
    post1 <- matrix(P[1,1] * Rp1 + P[1,2] * Rp2, row, col) * f$phi[1,t-1,,]
    post2 <- matrix(P[2,1] * Rp1 + P[2,2] * Rp2, row, col) * f$phi[2,t-1,,]
    
    fac <- sum(post1, na.rm = TRUE) + sum(post2, na.rm = TRUE)
    if (is.finite(fac) && fac > 0) {
      smooth[1,t-1,,] <- post1 / fac
      smooth[2,t-1,,] <- post2 / fac
    } else {
      sum_phi_prev <- sum(f$phi[,t-1,,], na.rm = TRUE)
      if (is.finite(sum_phi_prev) && sum_phi_prev > 0) {
        smooth[,t-1,,] <- f$phi[,t-1,,] / sum_phi_prev
      }
    }
  }
  
  return(smooth)
  
}



