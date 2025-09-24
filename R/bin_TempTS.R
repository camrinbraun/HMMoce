#' Convert time series to PDT-like data for use in HMMoce.
#' 
#' \code{bin_TempTS} compares tag SST to remotely sensed SST and calculates 
#' likelihoods
#' 
#' @param ts is data frame of time series data. Currently this is built for WC 
#'   tags but can handle others if you know what youre doing (e.g. proper input 
#'   format). Required inputs cols are just date, Depth and Temperature.
#' @param res is numeric indicating desired depth bin resolution. Default is 8 
#'   to match native PDT format on WC tags; however, theres no reason why this 
#'   couldnt be higher resolution, especially if consider we later match this to
#'   a 3D product such as HYCOM that is higher depth resolution (usually). **NOT WORKING**
#' @param tsDates is vector of desired dates as output. This is designed to be
#'   used for creating higher temporal resolution PDT datasets from time series
#'   data.
#' @return a dataframe that matches the required format for a PDT data product 
#'   for further use in HMMoce
#'   
#' @export
#' 
#' @examples
#' \dontrun{
#' }

bin_TempTS <- function (ts, tsDates=NULL, res = 8) 
{
	sm0 <- plyr::ddply(ts, c("date"), function(x) c(nrecTempTS = nrow(x[which(!is.na(x$Temperature) & !is.na(x$Depth)), ])))

	if (any(sm0$nrecTempTS == 0)) 
		warning("no Temperature/DepthTS available for dates:\n", 
						paste(sm0$date[which(sm0$nrecTempTS == 0)], collapse = ", "))
	if (!is.null(ts$Corrected.Depth)) {
		warning("Corrected.Depth column found. This data will be used instead!")
		ts$Depth <- ts$Corrected.Depth
	}
	ts <- ts[which(!is.na(ts$Temperature) & !is.na(ts$Depth)), ]
	sm <- plyr::ddply(ts, c("date"), function(x) c(nrecTempTS = nrow(x)))
	d.interval <- findInterval(sm$date, dateVec)
	fdates <- c()
	pdt.rec <- c()
	ud <- unique(d.interval)
	
	for (d in ud) {
		idx <- which(d.interval %in% d)
	  print(ts$date[idx][1])
		#i <- which(as.character(ts$date) == d)
		x <- ts[idx, ]
		depth.range <- range(x$Depth)
		depth <- x$Depth/res
		ii <- which((x$Depth%%res)/res >= 0.5)
		depth[ii] <- ceiling(depth[ii]) * res
		ii <- which((x$Depth%%res)/res < 0.5)
		depth[ii] <- trunc(depth[ii]) * res
		x$Depth <- depth
		unique_depths <- unique(depth)
		if (length(unique_depths) > 2) {
			h <- hist(depth, breaks = unique_depths, plot = F)
			identifiers <- c("DeployID", "Serial", "Ptt")
			identifiers <- identifiers[which(identifiers %in% names(ts))]
			add <- plyr::ddply(x[, which(names(x) %in% c(identifiers, "Depth", "Temperature"))], c(identifiers, "Depth"), function(x) c(nrecs = nrow(x), MeanTemp = mean(round(x$Temp, 2)), MinTemp = min(round(x$Temp, 2)), MaxTemp = max(round(x$Temp, 2))))
			add$date <- x$date[1]
			add$bin <- 1:nrow(add)
			pdt.rec <- rbind(pdt.rec, add)
		}
	}
	pdt.rec$MeanPDT <- rowMeans(cbind(pdt.rec$MaxTemp, pdt.rec$MinTemp))
	#pdt.rec$date <- as.Date(pdt.rec$date)
	return(pdt.rec)
}