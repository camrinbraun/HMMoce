#' Negative Log likelihood of parameters
#' 
#' @param parvec vector of length 6 containing * kernel 1 * kernel 2 * diagonal
#'   of 2x2 matrix
#' @param g grid from \code{\link{setup.grid}}
#' @param L final likelihood (2D)
#'   
#' @return parameter values
#' @export
#' @references Pedersen MW, Patterson TA, Thygesen UH, Madsen H (2011)
#'   Estimating animal behavior and residency from movement data. Oikos
#'   120:1281-1290. doi: 10.1111/j.1600-0706.2011.19044.x
#'   

get.nll.fun <- function(parvec = c(10, 30, 5, 2, .707, .8), g, L){
  
  K1 <- gausskern(parvec[1], parvec[2], muadv = 0)
  K2 <- gausskern(parvec[3], parvec[4], muadv = 0)
  P <- matrix(c(parvec[5], 1-parvec[5], 1-parvec[6], parvec[6]), 2, 2, byrow = TRUE)

  f = hmm.filter(g, L, K1, K2, P)
  nllf <- -sum(log(f$psi[f$psi>0]))
  print(paste0("HMM -log(L): ", nllf))
  #flush.console()
  nllf
  
}

