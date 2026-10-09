#' Calculate negative log likelihood
#' 
#' \code{neg.log.lik.fun} calculates the negative log likelihood
#' 
#' @param pars kernel parameters. Usually output from \code{est_params}
#' @param g see \code{hmm.filter}
#' @param L see \code{hmm.filter}
#' @param maskL see \code{hmm.filter}
#' @param bound.thr see \code{hmm.filter}
#' @param minBounds see \code{hmm.filter}
#' @param kernel character either "normal" or "isotropic", depending on desired kernel characteristics
#' 
#' @return a negative log likelihood
#' @author Paul Gatti
#' 

neg.log.lik.fun <- function(pars, g, L, maskL = TRUE, bound.thr = 0.1, minBounds = 10, kernel = 'normal'){
  
  if (length(pars) == 4){
    sigmas = pars[1:2]
    sizes = rep(ceiling(sigmas[1]*4),2)
    pb = pars[3:4]
    muadvs = c(0,0)
  } else if (length(pars) == 1){
    sigmas = pars[1]
    sizes = rep(ceiling(sigmas[1]*4),2)
    pb = NULL
    muadvs = c(0)
  }
  
  if (kernel == 'normal'){
    # behav 1
    K1 = gausskern.pg(sizes[1], sigmas[1], muadv = muadvs[1])
    
    # behav 2
    if (!is.null(pb)) K2 = gausskern.pg(sizes[2], sigmas[2], muadv = muadvs[2])
    
  } else if (kernel == 'isotropic'){
    
    # behav 1
    K1 = gausskern.isotrop(sizes[1], sigmas[1], muadv = muadvs[1], ratio.xy = 1)
    
    # behav 2
    if (!is.null(pb)) K2 = gausskern.isotrop(sizes[2], sigmas[2], muadv = muadvs[2], ratio.xy = 1)
    
  }

  if (!is.null(pb)) {
    P <- matrix(c(pb[1], 1 - pb[1], 1 - pb[2], pb[2]), 2, 2, byrow = TRUE)
  } else {
    P <- NULL
  }
  
  if (!is.null(pb)) {
    f. <- hmm.filter(g=g, L=L, K = list(K1, K2), P=P, m=2, maskL = maskL, bound.thr = bound.thr, minBounds = minBounds)
  } else{
    f. <- hmm.filter(g=g, L=L, K = list(K1), P=P, m=1, maskL = maskL, bound.thr = bound.thr, minBounds = minBounds)
  }
  
  nllf. <- -sum(log(f.$psi))
  return(nllf.)
}