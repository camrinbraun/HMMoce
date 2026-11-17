## ----init, echo=FALSE, message=FALSE-------------------------------------
#library(sp)
#library(HMMoce)
devtools::load_all('../../HMMoce')
#library(raster)
#library(fields)
#library(tidyverse)
options(tidy=TRUE)
knitr::opts_chunk$set(
  comment = '', fig.width = 6, fig.height = 6, tidy = TRUE)
dir <- getwd()
  
## vignette figure generation script is at github.com/camrinbraun/HMMoce_run/vignette/


## ----iniloc,size='small', tidy=TRUE--------------------------------------

# SET START/END LOCATIONS
## iniloc is dataframe containing cols: day, month, year, lat, lon and rows: start, end
iniloc <- data.frame(matrix(c(13, 10, 2015, 41.3, -69.27,
                              10, 4, 2016, 40.251, -36.061), nrow = 2, ncol = 5, byrow = T))
names(iniloc) <- list('day','month','year','lat','lon')
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')

# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- seq.POSIXt(tag, pop, by = '24 hours')

## ----load_sst,size='small'-----------------------------------------------
sstFile <- system.file("extdata", "141259-SST.csv", package = "HMMoce")
tag.sst <- read.wc(sstFile, type = 'sst', tag=tag, pop=pop, verbose=F)
tag.sst <- tag.sst[,c('Date','Depth','Temperature')]
head(tag.sst)


## ----load_pdt,size='small'-----------------------------------------------
# DEPTH-TEMPERATURE PROFILE DATA
## example is output from Wildlife Computer Portal
## pdt needs to contain at least:
##    - Date (POSIXct)
##    - Depth
##    - MinTemp
##    - MaxTemp
##    - MeanTemp (optional): if meantemp doesn't exist for whatever reason, HMMoce will calculate the midpoint between min/max temps and use that

pdtFile <- system.file("extdata", "141259-PDTs.csv", package = "HMMoce")
pdt <- read.wc(pdtFile, type = 'pdt', tag = tag, pop = pop, verbose = F)
pdt <- pdt[,c('Date','Depth','MinTemp','MaxTemp')]
head(pdt)

## ----load_series,size='small'--------------------------------------------
# DEPTH-TEMPERATURE TIME SERIES DATA
## exampling showing how to coerce depth-temp time series to a PDT-like summarized product

tsFile <- system.file("extdata", "141259-Series.csv", package = "HMMoce")
## bandaid to access the series data without installing new HMMoce
tsFile <- '~/work/RCode/HMMoce/inst/extdata/141259-Series.csv'
ts <- read.table(tsFile, sep=',', header=T)
ts$Date <- as.POSIXct(paste(ts$Day, ts$Time), format='%d-%b-%Y %H:%M:%S', tz='UTC')
ts <- ts[,c('Date','Depth','Temperature')]

## generate depth-temp summary from time series
pdt <- bin_TempTS(ts, out_dates = dateVec, bin_res = 25)
pdt <- pdt[,c('Date','Depth','MinTemp','MaxTemp')]
head(pdt)

## ----light_raw,size='small', cache=F-------------------------------------

lightFile <- system.file("extdata", "141259-LightLoc.csv", package = "HMMoce")
light <- read.wc(lightFile, type = 'light', tag=tag, pop=pop, verbose=F)
## combine character vectors "Day" and "Time" to generate POSIXct object
light$Date <- lubridate::dmy_hms(paste(light$Day, light$Time, sep = ' '))
light <- light[,c('Date','Type')]

## ----light_raw_ex,size='small'-------------------------------------------
head(light)

## ----light_est,size='small'----------------------------------------------

# LIGHT BASED POSITIONS FROM GPE2 (INSTEAD OF RAW LIGHTLOCS FROM PREVIOUS)
llFile <- system.file("extdata", "141259-Locations-GPE2.csv", package = "HMMoce")
lightloc <- read.table(llFile, sep = ',', header = T, blank.lines.skip = F)
lightloc <- lightloc[which(lightloc$Type != 'Argos'),]
lightloc <- lightloc[,c('Date','Longitude','Error.Semi.minor.axis','Latitude','Error.Semi.major.axis','Offset','Offset.orientation')]
lightloc$Date <- as.POSIXct(lightloc$Date, format = findDateFormat(lightloc$Date))
head(lightloc)

## ----light_ell,size='small', eval=FALSE, echo=FALSE, message=FALSE-------
#  
#  plot(c(-75,-65), c(25,60), type='n', xlab='Longitude', ylab='Latitude', main='Example light-based error ellipse')
#  points(lightloc$Longitude[1], lightloc$Latitude[1], pch=16)
#  plotrix::draw.ellipse(lightloc$Longitude[1], lightloc$Latitude[1], lightloc$Error.Semi.minor.axis[1] / 1000 / 111, lightloc$Error.Semi.major.axis[1] / 1000 / 111, col=alpha("red", 0.1))
#  text(-72, 58, 'Error semi minor = 71 km')
#  text(-72, 56, 'Error semi major = 1,500 km')
#  text(-72, 54, 'Offset and orientation = 0')
#  world(add=T)
#  

## ----mwtlight_err,size='small', eval=FALSE-------------------------------
#  
#  MWTdata ## from Lat&Long sheet from x-tag
#  names(MWTdata) = c('Date','Latitude','Longitude')
#  ## set fixed longitude error estimate in METERS
#  MWTdata$Error.Semi.minor.axis = .7 * 1000 * 111
#  

## ----mmd,size='small'----------------------------------------------------
mmdFile <- system.file("extdata", "141259-MinMaxDepth.csv", package = "HMMoce")
mmd <- read.table(mmdFile, sep = ',', header = T, blank.lines.skip = F)[,c('Date','MinDepth','MaxDepth')]
mmd$Date <- as.POSIXct(mmd$Date, format = findDateFormat(mmd$Date))
head(mmd)


## ----mmd_series,size='small'---------------------------------------------
seriesFile <- system.file("extdata", "141259-Series.csv", package = "HMMoce")
series <- read.table(seriesFile, sep = ',', header = T, blank.lines.skip = F)[,c('Day','Time','Depth','Temperature')]
series$Date <- as.Date(as.POSIXct(paste(series$Day, series$Time), format = '%d-%b-%Y %H:%M:%S', tz='UTC'))
mmd <- series %>% group_by(Date) %>% dplyr::summarise(n=n(), MinDepth = min(Depth, na.rm=T), MaxDepth = max(Depth, na.rm=T), .groups='keep')
head(mmd)


## ----splim,size='small'--------------------------------------------------
# SET SPATIAL LIMITS
# these are the lat/lon bounds of your study area (e.g. where you think the animal went)
sp.lim <- list(lonmin = -80,
               lonmax = -25,
               latmin = 15,
               latmax = 50)

## setup the spatial grid to base likelihoods on
locs.grid <- setup.locs.grid(sp.lim, res='quarter')


## ----env-sst,size='small', cache=FALSE, eval=FALSE-----------------------
#  ## here we get use a ridculously small spatial extent for this species, just to show some example environmental data
#  udates <- seq.Date(as.Date(tag), as.Date(pop), by = 'day')
#  sst.dir <- paste0(dir, '/EnvData/sst/')
#  if (!dir.exists(sst.dir)) dir.create(sst.dir, recursive = TRUE)
#  if (!file.exists(paste0(sst.dir, 'oisst_', udates[1], '.nc'))) get.env(udates[1], filename='oisst', type = 'sst', sst.type='oi', spatLim = sp.lim, save.dir = sst.dir)
#  
#  sst <- terra::rast(paste0(sst.dir, 'oisst_', udates[1], '.nc'))
#  plot(sst, main='Example SST field')
#  world(add=T)

## ----env-hycom,size='small', eval=FALSE, cache=FALSE---------------------
#  ## you need some representation of environmental depth-temperature
#  ## in this case we're using hycom
#  
#  dir <- getwd()
#  hycom.dir <- paste0(dir,'/EnvData/hycom/')
#  if (!dir.exists(hycom.dir)) dir.create(hycom.dir, recursive = TRUE)
#  if (!file.exists(paste0(hycom.dir, 'hycom_', udates[1], '.nc'))) get.env(udates[1], filename='hycom', type = 'hycom', spatLim = sp.lim_small, save.dir = hycom.dir)
#  
#  b <- terra::rast(paste0(hycom.dir, 'hycom_', udates[1], '.nc'))
#  
#  ## plot top 3 depth levels from HYCOM
#  #par(mfrow=c(3,1))
#  for (i in c(1,10,20)){
#    plot(b[[i]])
#    world(add=T)
#  }

## ----env-bathy, size='small', eval=FALSE---------------------------------
#  
#  bathy.dir <- paste0(dir, '/EnvData/bathy/')
#  if (!dir.exists(bathy.dir)) dir.create(bathy.dir, recursive = TRUE)
#  
#  if (!file.exists(paste0(bathy.dir, 'bathy.nc'))){
#    bathy <- HMMoce::get.bath.data(sp.lim_small$lonmin, sp.lim_small$lonmax, sp.lim_small$latmin, sp.lim_small$latmax, folder = bathy.dir, res=1)
#  } else{ ## OR (once downloaded and reading the .nc later)
#    bathy <- irregular_ncToRaster(paste0(bathy.dir, 'bathy.nc'), varid = 'topo')
#  }
#  
#  ## example bathymetry data
#  plot(bathy, main='Example bathymetry field')
#  world(add=T)

## ----lik-light, size='small', eval=FALSE---------------------------------
#  
#    # LIGHT-BASED LIKELIHOODS
#    L.light.srss <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec, res = 1, focalDim = 15) # if using raw light data
#    L.light.lonlat <- calc.lightloc(lightloc, locs.grid = locs.grid, dateVec = dateVec, errEll = TRUE) # if using light-based positions (w/ Lat & Lon)
#  
#      L.light.lon <- calc.lightloc(lightloc, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE) # if using light-based positions (w/ Lon only)
#  

## ----lik-sst, size='small', eval=FALSE-----------------------------------
#  
#    # SST LIKELIHOODS
#    L.sst <- calc.sst(tag.sst, filename='oisst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
#  
#    ## identical to above but in parallel
#    L.sst <- calc.sst.par(tag.sst, filename='oisst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1, ncores=2)
#  

## ----lik-ohc, size='small', eval=FALSE-----------------------------------
#  
#    # OCEAN HEAT CONTENT (INTEGRATED PDTS)
#    L.ohc <- calc.ohc(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
#    ## the parallel version:
#    #L.ohc <- calc.ohc.par(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)

## ----lik-hycom, size='small', eval=FALSE---------------------------------
#  
#    # HYCOM PROFILE BASED LIKELIHOODS
#    L.hycom <- calc.hycom(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = as.POSIXct(dateVec), use.se = T)
#    ## the parallel version:
#    #L.hycom <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = as.POSIXct(dateVec), use.se = T)

## ----lik-bathy, size='small', eval=FALSE---------------------------------
#  
#    ## bathymetry based likelihood
#    ## resample bathy to a more reasonable (coarse) grid for likelihood calculations
#    ## hi-res bathy grid will work but will take longer
#    bathy_resamp <- terra::resample(bathy, L.sst) # or whatever grid makes sense to resample to
#    L.bathy <- calc.bathy(mmd, bathy_resamp, dateVec, focalDim = 3, sens.err = 5, lik.type = 'dnorm')
#    L.bathy <- calc.bathy.par(mmd, bathy_resamp, dateVec, focalDim = 3, sens.err = 5, lik.type = 'max', ncores = 4)
#  

## ----lik-bt, size='small', eval=FALSE------------------------------------
#  
#    ## bottom temperature based likelihood
#    L.bt <- calc.bottomTemp(tag.bt, dateVec, focalDim = 3, sens.err = 1, bt.dir = bt.dir, filename = 'bottomT', varName = 'Temperature')
#    #L.bt <- calc.bottomTemp.par(tag.bt, dateVec, focalDim = 3, sens.err = 1, bt.dir = bt.dir, filename = 'bottomT', varName = 'Temperature', ncores = 4)

## ----makeL, size='small', eval=FALSE-------------------------------------
#  
#  # COMBINE LIKELIHOOD MATRICES
#  # make list of rasters
#  L.rasters <- list(L.light, L.sst, L.ohc)
#  
#  ## typically these will need to be resampled to have matching resolution and extent
#  resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) terra::res(x)[1]))
#  L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
#  
#  ## good idea to save these likelihood rasters at some point to keep from having to re-calculate them if (when) R crashes
#  
#  ## finally, combine the likelihoods into a master observation likelihood that will be used for the modeling
#  L <- make.L(ras.list,
#              iniloc,
#              dateVec,
#              known.locs = NULL)
#  
#  image.plot(res$g$lon[1,], res$g$lat[,1], L[12,,],
#             main='Example overall likelihood at single time step',
#             zlim=c(.01, 1),
#             xlab='Longitude', ylab='Latitude')

## ----param, size='small', eval=FALSE-------------------------------------
#  
#  ## if you want to try coarse grids for parameter estimation
#  ## use coarse.L() and supply the outputs in place of L and g below
#  L.mle <- coarse.L(L, L.res$L.rasters)$L.mle
#  g.mle <- coarse.L(L, L.res$L.rasters)$g.mle
#  
#  ## opt.params is a wrapper for the various optimization routines in HMMoce
#  pars.optim <- opt.params(pars.init = c(2,.2,.6,.8),
#                             lower.bounds = c(0.1, 0.001, .1, .1),
#                             upper.bounds = c(6, .6, .9, .9),
#                             g = L.res$g,
#                             L = L,
#                             alg.opt = 'optim',
#                             write.results = FALSE)
#  ## about 22 mins on example blue shark 141259 with full grid
#  
#  pars.optim.mle <- opt.params(pars.init = c(2,.2,.6,.8),
#                             lower.bounds = c(0.1, 0.001, .1, .1),
#                             upper.bounds = c(6, .6, .9, .9),
#                             g = g.mle,
#                             L = L.mle,
#                             alg.opt = 'optim',
#                             write.results = FALSE)
#  
#  ## about 1.5 mins on example blue shark 141259 with coarse grid
#  
#  ## nlminb is also supported in HMMoce but testing suggests
#  ## this is rarely the best choice due to lack of convergence and other issues
#  pars.nlminb <- opt.params(pars.init = c(2,.2,.6,.8),
#                              lower.bounds = c(0.1, 0.001, .1, .1),
#                              upper.bounds = c(5, .5, .9, .9),
#                              g = L.res$g,
#                              L = L,
#                              alg.opt = 'nlminb',
#                              write.results = FALSE)
#  ## about 30 mins on blue shark 141259
#  
#  pars.ga <- opt.params(pars.init = c(2,.2,.6,.8),
#                            lower.bounds = c(0.1, 0.001, .1, .1),
#                              upper.bounds = c(6, .6, .9, .9),
#                              g = L.res$g,
#                              L = L,
#                              alg.opt = 'ga',
#                              write.results = FALSE,
#                              ncores = ceiling(parallel::detectCores() * .9))
#    ## about 1.6 hrs on blue shark 141259 w 15 cores
#  
#    pars.ga.mle <- opt.params(pars.init = c(2,.2,.6,.8),
#                          lower.bounds = c(0.1, 0.001, .1, .1),
#                          upper.bounds = c(6, .6, .9, .9),
#                          g = g.mle,
#                          L = L.mle,
#                          alg.opt = 'ga',
#                          write.results = FALSE,
#                          ncores = ceiling(parallel::detectCores() * .9))
#     ## about 2 mins with MLE grid but results way different
#  
#    ## example using the genetic algorithm and only one behavior state
#    pars.ga.one <- opt.params(pars.init = c(2),
#                          lower.bounds = c(1),
#                          upper.bounds = c(8),
#                          g = L.res$g,
#                          L = L,
#                          alg.opt = 'ga',
#                          write.results = FALSE,
#                          ncores = 4)
#  

## ----pars1, size='small', eval=FALSE-------------------------------------
#  
#  ## as an example, grab the pars from the GA using both states
#  #pars <- pars.ga$par
#  
#  ## or use values from a previous run as done here:
#  pars <- c(4.964, 0.217, 0.367, 0.484)
#  
#  ## if only one state is used in the estimation routines (above), this will catch that here and set the parameters accordingly
#  if (length(pars) == 4){
#    sigmas = pars[1:2]
#    sizes = rep(ceiling(sigmas[1]*4),2)
#    pb = pars[3:4]
#    muadvs = c(0,0)
#  } else if (length(pars) == 1){
#    sigmas = pars[1]
#    sizes = rep(ceiling(sigmas[1]*4),2)
#    pb = NULL
#    muadvs = c(0)
#  }
#  
#  K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
#  if (!is.null(pb)) K2 <- gausskern.pg(sizes[2], sigmas[2], muadv=muadvs[2])
#  
#  image.plot(K1, main='Diffusion kernel for "migratory" state')
#  image.plot(K2, main='Diffusion kernel for "resident" state')
#  
#  # ** show kernel and how it works **
#  ## convolution?
#  

## ----pars2, size='small', eval=FALSE-------------------------------------
#  
#  ## set transition matrix, if applicable
#  if (!is.null(pb)){
#   P <- matrix(c(pb[1], 1- pb[1], 1 - pb[2], pb[2]), 2, 2, byrow=TRUE)
#  } else{
#   P <- NULL
#  }
#  

## ----hmm, size='small', eval=FALSE---------------------------------------
#  
#  # RUN THE FILTER STEP
#  if (!is.null(pb)){
#    K <- list(K1,K2)
#    f <- hmm.filter(g = L.res$g, L = L, K = K, P = P, m = 2)
#  } else{
#    K <- list(K1)
#    f <- hmm.filter(g = L.res$g, L = L, K = K, P = P, m = 1)
#  }
#  nllf <- -sum(log(f$psi[f$psi>0])) # negative log-likelihood
#  AIC <- 2 * nllf + 2 * length(which(!is.na(pars)))
#  
#  # RUN THE SMOOTHING STEP
#  s <- hmm.smoother(f, K = K, L = L, P = P)
#  
#  # GET THE MOST PROBABLE TRACK AS MEAN OF POSTERIOR DISTRIBUTION OF STATE
#  tr <- calc.track(s, g = L.res$g, dateVec, iniloc, method='mean')
#  

## ----plot-hmm-F, size='small', eval=FALSE--------------------------------
#  
#  plotHMM(s, tr, dateVec, ptt='141259_example', save.plot = F)
#  
#  plotRD(s, tr, ptt='141259_example', g = g, makePlot = TRUE, save.plot=FALSE)
#  

