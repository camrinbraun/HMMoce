#' Creates grid from matrices
#' 
#' @param x vector, usually of longitude data
#' @param y vector, usually of latitude data
#' @return list of 2 matrices for lon/lat values
#'   
#' @examples
#' \dontrun{
#' x <- c(1, 2, 3)
#' y <- c(10, 9, 8)
#' meshgrid(x, y)
#' }
#' @references Pedersen MW, Patterson TA, Thygesen UH, Madsen H (2011)
#'   Estimating animal behavior and residency from movement data. Oikos
#'   120:1281-1290. doi: 10.1111/j.1600-0706.2011.19044.x


meshgrid <- function(x, y){
  Y <- repmat(as.matrix(y), 1, length(x))
  X <- repmat(t(as.matrix(x)), length(y), 1)
  list(X = X, Y = Y)
}
