#' Create Gaussian Kernel
#' 
#' \code{gausskern.nostd} calculates 2D Gaussian kernel based on kernel size, 
#' deviation, and advection. Differs from \code{gausskern} in that it should result in a more circular, and perhaps realistic, kernel. This *should* avoid the issue with square kernels allowing more diffusion in the diagonal than horizontal or vertical.
#' 
#' @param siz size of the kernel, siz x siz. Must be a positive integer.
#' @param sigma standard deviation of the kernel. Unit is cell width. Must be a 
#'   positive number.
#' @param muadv advection of the kernel. Unit of the input is cell width. 
#'   Defaults to 0.
#' @return Gaussian kernel as a 2D matrix of size (siz x siz)
#' @export
#' 
#' @examples
#' kern = gausskern.nostd(3, 0.5)
#' 
#' @author Paul Gatti

gausskern.nostd <- function (siz, sigma, muadv = 0) {
  if (round(siz) < 1) siz = 1
  x = 1:round(siz)
  mu = c(mean(x), mean(x)) + muadv

  fx = exp(-0.5 * ((x - mu[1])/sigma)^2)/sqrt((2 * pi) * sigma**2)
  fy = exp(-0.5 * ((x - mu[2])/sigma)^2)/sqrt((2 * pi) * sigma**2)
  fx[!is.finite(fx)] = 0
  fy[!is.finite(fy)] = 0
  kern = (fx %*% t(fy))
  kern[is.nan(kern)] = 0
  kern
}
