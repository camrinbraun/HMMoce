#' Create isotropic Gaussian Kernel
#' 
#' \code{gausskern.isotrop} calculates 2D Gaussian kernel based on kernel size, 
#' deviation, and advection
#' 
#' @param siz size of the kernel, siz x siz. Must be a positive integer.
#' @param sigma standard deviation of the kernel. Unit is cell width. Must be a 
#'   positive number.
#' @param muadv advection of the kernel. Unit of the input is cell width. 
#'   Defaults to 0.
#' @param ratio.xy ratio in x:y
#' @return Gaussian kernel as a 2D matrix of size (siz x siz)
#' @export
#' 
#' @examples
#' kern = gausskern.isotrop(3, 0.5)
#' @note No final standardization by the sum. This function is used to design a pseudo-optimal size for the kernel based on sigma while maintaining 99% of the allowable diffusion.
#' 
#' @author Paul Gatti

gausskern.isotrop <- function (siz, sigma, muadv = 0, ratio.xy = 1) {
  
  # x longitude
  sizx = siz / ratio.xy
  sigmax = siz / ratio.xy
  if (round(sizx) < 1) sizx = 1
  x = 1:round(sizx)
  mux = c(mean(x), mean(x)) + muadv/ratio.xy
  fx = exp(-0.5 * ((x - mux[1])/sigmax)^2)/sqrt((2 * pi) * sigmax**2) 
  fx[!is.finite(fx)] = 0
  
  # y latitude
  sizy = siz
  sigmay = siz
  if (round(sizy) < 1) sizy = 1
  y = 1:round(sizy)
  muy = c(mean(y), mean(y)) + muadv
  fy = exp(-0.5 * ((y - muy[2])/sigmay)^2)/sqrt((2 * pi) * sigmay**2) 
  fy[!is.finite(fy)] = 0
  
  # kernel
  kern = (fx %*% t(fy))
  kern = kern/(sum(kern, na.rm = T))#;kern
  kern[is.nan(kern)] = 0
  kern
}
