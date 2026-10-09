#' Read and format tag data
#' 
#' \code{read.wc} reads and formats tag data output from Wildlife Computers Data
#' Portal
#' 
#' @param filename is path to the file where your data lives
#' @param tag is POSIXct object of the tagging date
#' @param pop is POSIXct object of the pop-up date
#' @param type is character indicating which type of data to read. Choices are 
#'   'sst', 'pdt', 'light' corresponding to those data files output from WC Data
#'   Portal
#' @param dateFormat is input format argument for \code{as.POSIXct} or NULL.
#'   Default is NULL and if NULL the code will attempt to detect date format but
#'   weird things can happen so its best to input your own.
#' @param verbose is logical indicating whether a verbose output with more 
#'   details on the loaded files is desired. Default is FALSE.
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr n
#'   
#' @return a the data read as a data.frame
#' 
#' @examples
#' \dontrun{
#' # example data in the package
#' sstFile <- system.file("extdata", "141259-SST.csv", package = "HMMoce")
#' ptt <- 141259
#' 
#' # set temporal and spatial bounds
#' iniloc <- data.frame(matrix(c(13, 10, 2015, 41.3, -69.27, 10, 4, 2016, 40.251, -36.061),
#'  nrow = 2, ncol = 5, byrow = TRUE))
#'  names(iniloc) <- list('day','month','year','lat','lon')
#'  tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), 
#'  format = '%d/%m/%Y', tz='UTC')
#'  pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), 
#'  format = '%d/%m/%Y', tz='UTC')
#'
#' # read and format the example data
#' tag.sst <- read.wc(ptt, sstFile, type = 'sst', tag=tag, pop=pop)
#' }
#' @export

read.wc <- function(filename, tag, pop, type = 'sst', dateFormat=NULL, verbose=FALSE){
  
  #if(substr(wd, nchar(wd), nchar(wd)) == '/'){
  #} else{
  #  wd <- paste(wd, '/', sep='')
  #}
  
  if(type == 'pdt'){
    # READ IN PDT DATA FROM WC FILES
    #data <- utils::read.table(paste(wd, ptt,'-PDTs.csv', sep=''), sep=',',header=T,blank.lines.skip=F, skip = 0)
    #data <- utils::read.csv(filename)# sep=',', header=T, blank.lines.skip=F, skip = 0)
    data <- utils::read.table(filename, sep=',', header=T, blank.lines.skip=F, skip = 0)
    #if (length(grep('Discont16', names(data))) == 0 & ncol(data) > 89) names(data)[90:94] <- c('Depth16','MinTemp16','MaxTemp16','X.Ox16','Discont16')
    if (verbose) print(paste('If read.wc() fails for type=pdt, check the number of column headers in the PDTs.csv file.'))
    data <- extract.pdt(data)
    if (is.null(dateFormat)){
      dts <- as.POSIXct(data$Date, format = findDateFormat(data$Date))
      
    } else{
      dts <- as.POSIXct(data$Date, format = dateFormat)
      
    }
    data$Date <- dts
    d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
    #didx <- dts >= (tag + d1) & dts <= (pop - d1)
    didx <- dts >= tag & dts <= pop
    data <- data[didx,]; dts <- dts[didx]
    
    # check for days with not enough data for locfit
    data1 <- data
    data1$dts <- as.Date(dts)
    dt.cut <- data.frame(group_by(data1, dts) %>% summarise(n=n()))
    #dt.cut <- dt.cut[which(dt.cut[,2] < 3),1]
    #if (length(dt.cut) == 0){
    #  
    #} else{
    #  data <- data1[-which(data1$dts %in% dt.cut),]
    #}
    
    # get unique dates
    udates <- unique(as.Date(data$Date))
    
    # get data gaps
    #print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have PDT data...', sep=''))
    #gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
    #print(paste('Data gaps are ', paste(gaps[gaps > 1], collapse=', '), ' days in PDT...'))
    gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
    if(verbose){
      print(utils::head(data))
      print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have PDT data...', sep=''))
      print(paste('Data gaps are ', paste(gaps[gaps > 1], collapse=', '), ' days in PDT...'))
    }

  } else if(type == 'sst'){
    # READ IN TAG SST FROM WC FILES
    #data <- utils::read.table(paste(wd, ptt, '-SST.csv', sep=''), sep=',',header=T, blank.lines.skip=F)
    data <- utils::read.table(filename, sep=',',header=T, blank.lines.skip=F)
    if (is.null(dateFormat)){
      dts <- as.POSIXct(data$Date, format = findDateFormat(data$Date))
      
    } else{
      dts <- as.POSIXct(data$Date, format = dateFormat)
      
    }
    data$Date <- dts
    d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
    #didx <- dts >= (tag + d1) & dts <= (pop - d1)
    didx <- dts >= tag & dts <= pop
    data <- data[didx,]
    if (length(data[,1]) < 1){
      stop('Something wrong with reading and formatting of tags SST data. Check date format.')
    }
    dts <- as.POSIXct(data$Date, format = findDateFormat(data$Date))
    udates <- unique(as.Date(dts))
    data <- data[,c('Date','Depth','Temperature')]
    
    # get data gaps
    gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
    if (verbose){
      print(utils::head(data))
      print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have SST data...', sep=''))
      print(paste('Data gaps are ', paste(gaps[gaps > 1], collapse=', '), ' days...'))
    }
    
  } else if(type == 'light'){
    # READ IN LIGHT DATA FROM WC FILES
    #data <- utils::read.table(paste(wd,'/', ptt, '-LightLoc.csv', sep=''), sep=',',header=T, blank.lines.skip=F,skip=2)
    data <- try(utils::read.table(filename, sep=',',header=T, blank.lines.skip=F, skip=2), TRUE)
    
    if (class(data) == 'try-error'){
      data <- try(utils::read.table(filename, sep=',',header=T, blank.lines.skip=F, skip=0), TRUE)
      if (class(data) == 'try-error') stop('Tried reading light data with skip=2 (old WC format) and skip=0 (new WC format) but both failed. Check source light data file and try again.')
    }
    
    if(!any(grep('depth', names(data), ignore.case=T))) data <- utils::read.table(filename, sep=',',header=T, blank.lines.skip = F, skip = 1)
    if(!any(grep('depth', names(data), ignore.case=T))) data <- utils::read.table(filename, sep=',',header=T, blank.lines.skip = F, skip = 0)
    data <- data[which(!is.na(data[,1])),]
    
    #dts <- as.POSIXct(data$Day, format = findDateFormat(data$Day))
    dts <- lubridate::parse_date_time(data$Day, orders=c('dby', 'dbY'), tz='UTC')
    
    if(as.Date(dts[1]) > as.Date(Sys.Date()) | as.Date(dts[1]) < '1990-01-01'){
      stop('Error: dates are in the future or before 1990 and thus likely did not parse correctly.')    
    }
    
    d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
    didx <- (dts > (tag + d1)) & (dts < (pop - d1))
    data <- data[didx,]
    data$dts <- as.Date(dts[didx])
    udates <- unique(as.Date(dts))
    
    # get data gaps
    gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
    if (verbose){
      print(utils::head(data))
      print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have light data...', sep=''))
      print(paste('Data gaps are ', paste(gaps[gaps > 1], collapse=', '), ' days...'))
    }
    
  }
  
  #return(list(data = data, udates = udates))
  return(data)
  
}

