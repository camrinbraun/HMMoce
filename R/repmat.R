#' Repeat your matrix?
#' 
#' @param X is matrix you want to repeat
#' @param m is output row dimension
#' @param n is output col dimension
#' @return repeated matrix
#' @references Pedersen MW, Patterson TA, Thygesen UH, Madsen H (2011)
#'   Estimating animal behavior and residency from movement data. Oikos
#'   120:1281-1290.
#'   

repmat <- function(X, m, n){
  # function from Pedersen et al 2011
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
}