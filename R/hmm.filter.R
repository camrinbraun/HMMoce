#' HMM filter functions
#' 
#' 
#' @param g grid from \code{\link{setup.grid}}
#' @param L is likelihood array output from \code{make.L}
#' @param K is list of length m where each element represents the movement (diffusion) kernel for that behavior state \code{\link{gausskern}}
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

hmm.filter <- function(g, L, K, P, m, maskL = FALSE, bound.thr = 0.1, minBounds = 10){
  
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
  
  # Build separable 1D Toeplitz convolution operators when K[[ii]] is low-rank (rank 1 or 2)
  build_toeplitz_1d <- function(n, kvec) {
    lk <- length(kvec)
    ctr <- floor(lk / 2) + 1L
    idx <- outer(seq_len(n), seq_len(n), "-") + ctr
    valid <- idx >= 1L & idx <= lk
    Tmat <- matrix(0, n, n)
    Tmat[valid] <- kvec[idx[valid]]
    Tmat
  }
  
  sep_ops <- vector("list", m)
  for (ii in 1:m){
    K_mat <- as.matrix(K[[ii]])
    sv <- svd(K_mat)
    # Check if rank-1 or rank-2 decomposition captures >= 99.99% of kernel energy
    r_use <- if (length(sv$d) == 1L || sv$d[2] < 1e-4 * sv$d[1]) {
      1L
    } else if (length(sv$d) == 2L || sv$d[3] < 1e-3 * sv$d[1]) {
      2L
    } else {
      0L
    }
    if (r_use > 0L) {
      ops_r <- vector("list", r_use)
      for (rk in seq_len(r_use)) {
        u_vec <- sqrt(sv$d[rk]) * sv$u[, rk]
        v_vec <- sqrt(sv$d[rk]) * sv$v[, rk]
        ops_r[[rk]] <- list(
          Tv = build_toeplitz_1d(col, v_vec),
          Tu_t = t(build_toeplitz_1d(row, u_vec))
        )
      }
      sep_ops[[ii]] <- ops_r
    } else {
      K[[ii]] <- imager::as.cimg(K_mat)
      sep_ops[[ii]] <- NULL
    }
  }

  psi <- rep(0, T - 1) # sum of the probability of states at each step

  ## create empty lists for upcoming vars
  q <- vector("list", m); post <- vector("list", m)
  
  # Start filter iterations
  for(t in 2:T){
    
    for (ii in 1:m){
      # convolve previous day's likelihood with movement kernels
      M_prev <- phi[ii, t-1, , ]
      if (!is.null(sep_ops[[ii]])) {
        ops_r <- sep_ops[[ii]]
        q_mat <- ops_r[[1]]$Tv %*% M_prev %*% ops_r[[1]]$Tu_t
        if (length(ops_r) > 1L) {
          for (rk in 2:length(ops_r)) {
            q_mat <- q_mat + (ops_r[[rk]]$Tv %*% M_prev %*% ops_r[[rk]]$Tu_t)
          }
        }
        q_mat[q_mat < 0] <- 0
        q[[ii]] <- q_mat
      } else {
        p_img <- imager::as.cimg(t(M_prev))
        q[[ii]] <- t(as.matrix(imager::convolve(p_img, K[[ii]])))
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
      
      denom <- if (is.finite(psi[t-1]) && psi[t-1] > 0) psi[t-1] else 1
      for (ii in 1:m){
        phi[ii,t,,] <- post[[ii]] / denom
      }
      #phi[1,t,,] <- post1 / (psi[t-1] + 1e-15)
      #phi[2,t,,] <- post2 / (psi[t-1] + 1e-15)
      
  } # t loop
  
  list(phi = phi, pred = pred, psi = psi)
  
}