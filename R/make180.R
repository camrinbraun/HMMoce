#' convert  0, 360 --> -180, 180
#' @param lon is numeric or integer longitude to be converted
#' @return converted lon

make180 <- function(lon) {
  isnot180 <- max(lon) > 180
  if (isnot180) {
    ind <- which(lon > 180)
    lon[ind] <- lon[ind] - 360
  }
  return(lon)
}
