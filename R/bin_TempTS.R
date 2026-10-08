#' Convert time series to PDT-like data for use in HMMoce.
#' 
#' \code{bin_TempTS} converts time series data (such as from MT tags) to PDT-like data for use in HMMoce
#' 
#' @author Ben Galuardi
#' @param ts Data frame of time series data (Date, Depth, Temperature).
#' @param out_dates POSIX vector of desired output dates.
#' @param bin_res Numeric indicating desired depth bin resolution (default = 8).
#' @return A dataframe of summarized depth-temperature information.
#' @export
#' 
#' @examples
#' \dontrun{
#' tsFile <- system.file("extdata", "141259-Series.csv", package = "HMMoce")
#' ts <- read.table(tsFile, sep=',', header=T)
#' ts$Date <- as.POSIXct(paste(ts$Day, ts$Time), format='%d-%b-%Y %H:%M:%S', tz='UTC')
#' ts <- ts[,c('Date','Depth','Temperature')]
#' 
#' ## generate example daily summary of depth-temp profiles
#' pdt_dates <- seq.POSIXt(tag, pop, by = 'day') 
#' 
#' ## generate depth-temp summary from time series
#' pdt <- bin_TempTS(ts, out_dates = pdt_dates, bin_res = 25)
#' pdt <- pdt[,c('Date','Depth','MinTemp','MaxTemp')]
#' }
bin_TempTS <- function(ts, out_dates, bin_res = 8) {
  # Drop missing depth or temperature records immediately
  ts <- ts[!is.na(ts$Temperature) & !is.na(ts$Depth), ]
  
  # Assign time indices and vectorize depth binning
  ts <- ts %>%
    dplyr::mutate(
      time_idx = findInterval(Date, out_dates),
      Depth = round(Depth / bin_res) * bin_res
    )
  
  # Group, summarize, and filter in a single pipeline
  pdt_rec <- ts %>%
    dplyr::group_by(time_idx, Depth) %>%
    dplyr::summarize(
      nrecs = dplyr::n(),
      MeanTemp = mean(round(Temperature, 2), na.rm = TRUE),
      MinTemp = min(round(Temperature, 2), na.rm = TRUE),
      MaxTemp = max(round(Temperature, 2), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(time_idx) %>%
    # Mimics original logic: require >2 unique depth bins per time step
    dplyr::filter(dplyr::n_distinct(Depth) > 2) %>%
    dplyr::mutate(
      Date = out_dates[time_idx],
      bin = dplyr::row_number(),
      MeanPDT = (MaxTemp + MinTemp) / 2
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Depth, nrecs, MeanTemp, MinTemp, MaxTemp, Date, bin, MeanPDT)
  
  return(as.data.frame(pdt_rec))
}