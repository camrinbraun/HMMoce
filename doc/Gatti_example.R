## ----setup, eval=FALSE---------------------------------------------------
#  #---------------------------------------------------------------------------------------------
#  # pgatti
#  # creation 15/10/19
#  # last update 21/10/19
#  #  R version 3.4.4
#  #=============================================================================================
#  # HMMoce
#  # shark example
#  #=============================================================================================
#  # 1. load tag & env data + compute likelihoods
#  # 2. Bathymetry and bottom temp based likelihoods
#  # 3. model fitting
#  #---------------------------------------------------------------------------------------------
#  
#  rm(list=ls())
#  graphics.off()
#  #
#  
#  # shortcut functions------
#  ac=function(x){return(as.character(x))}
#  an=function(x){return(as.numeric(ac(x)))}
#  an.=function(x){return(as.numeric(x))}
#  av=function(x){return(as.vector(x))}
#  gf=function(){graphics.off()}
#  #-------------------------
#  
#  dir='C:/Users/pgatti/Documents/Halibut/R/functions/whoi_hmmoce_pkg/'
#  source(paste0(dir,'pgatti_potential_functions.r'))
#  
#  setwd(dir)
#  library(HMMoce)
#  library(GA)
#  
#  #----------------------------------------------
#  # 1. load tag & env data + compute likelihoods
#  #----------------------------------------------
#  preprocessing=F
#  if(preprocessing){
#    print('not runned yet')
#    stop()
#  
#    # following copy/pasted from Camrin example
#  
#    #------------
#    # LOAD THE TAG DATA
#    #------------
#    # setwd()
#  
#    # SET INITIAL LOCATIONS (TAG AND POP-UP)
#    iniloc <- data.frame(matrix(c(13, 10, 2015, 41.3, -69.27,
#                                  10, 4, 2016, 40.251, -36.061), nrow = 2, ncol = 5, byrow = T))
#    names(iniloc) <- list('day','month','year','lat','lon')
#    tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
#    pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')
#  
#    # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
#    dateVec <- as.Date(seq(tag, pop, by = 'day'))
#  
#    # READ IN DATA AS OUTPUT FROM WC PORTAL
#    # SST DATA
#    sstFile <- system.file("extdata", "141259-SST.csv", package = "HMMoce")
#    tag.sst <- read.wc(sstFile, type = 'sst', tag=tag, pop=pop, verbose=T)
#    sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data
#  
#    # DEPTH-TEMPERATURE PROFILE DATA
#    pdtFile <- system.file("extdata", "141259-PDTs.csv", package = "HMMoce")
#    pdt <- read.wc(pdtFile, type = 'pdt', tag=tag, pop=pop, verbose=T)
#    pdt.udates <- pdt$udates; pdt <- pdt$data
#  
#    # RAW LIGHT DATA
#    #lightFile <- system.file("extdata", "141259-LightLoc.csv", package = "HMMoce")
#    #light <- read.wc(ptt, lightFile, type = 'light', tag=tag, pop=pop);
#    #light.udates <- light$udates; light <- light$data
#  
#    # LIGHT BASED POSITIONS FROM GPE2 (INSTEAD OF RAW LIGHTLOCS FROM PREVIOUS)
#    locsFile <- system.file("extdata", "141259-Locations-GPE2.csv", package = "HMMoce")
#    locs <- read.table(locsFile, sep = ',', header = T, blank.lines.skip = F)
#    locDates <- as.Date(as.POSIXct(locs$Date, format=findDateFormat(locs$Date)))
#  
#    # SET SPATIAL LIMITS
#    # these are the lat/lon bounds of your study area (e.g. where you think the animal went)
#    sp.lim <- list(lonmin = -82,
#                   lonmax = -25,
#                   latmin = 15,
#                   latmax = 50)
#  
#    #------------
#    # GET ENVIRONMENTAL DATA
#    #------------
#    # env data downloads can be
#    #large, depending on application for 180 days of data spanning the NW Atlantic
#    #(the example application), the downloads will take ~10mins on Amazon EC2.
#    #Personal computers will likely be slower.
#  
#    # DOWNLOAD SST DATA
#    #sst.dir <- paste(tempdir(), '/sst/', sep='')
#  #  sst.dir <- paste('~/work/EnvData/hmmoce', '/sst/', sep='')
#    dir.create(sst.dir, recursive = TRUE)
#    get.env(sst.udates, filename='oisst', type = 'sst', sst.type='oi', spatLim = sp.lim, save.dir = sst.dir)
#  
#    # YOU NEED SOME REPRESENTATION OF ENVIRONMENTAL DEPTH-TEMPERATURE
#    # HYCOM DATA
#    #hycom.dir <- paste0(dir,'EnvData/hycom/', sep='')
#    hycom.dir <- paste('~/work/EnvData/hmmoce', '/hycom/', sep='')
#    dir.create(hycom.dir, recursive = TRUE)
#    get.env(pdt.udates[90:length(pdt.udates)],
#            filename='hycom', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)
#  
#    # OR WORLD OCEAN ATLAS DATA
#    #woa.dir <- paste(tempdir(), '/woa/', sep='')
#    #dir.create(woa.dir, recursive = TRUE)
#    #get.env(type = 'woa', resol = 'quarter', save.dir = woa.dir)
#    # THEN LOAD AND CHECK THE DOWNLOADED RDA FILE FOR WOA
#    #load(paste(woa.dir,'woa.quarter.rda',sep=''))
#    #str(woa.quarter)
#    #List of 4
#    #$ watertemp: num [1:44, 1:46, 1:57, 1:12] 26.5 26.5 26.4 26.3 26.2 ...
#    #$ lon      : num [1:44(1d)] -95.5 -94.5 -93.5 -92.5 -91.5 -90.5 -89.5 -88.5 -87.5 -86.5 ...
#    #$ lat      : num [1:46(1d)] 9.5 10.5 11.5 12.5 13.5 14.5 15.5 16.5 17.5 18.5 ...
#    #$ depth    : num [1:57(1d)] 0 5 10 15 20 25 30 35 40 45 ...
#  
#    # BATHYMETRY
#    bathy.dir <- paste('~/work/EnvData/hmmoce', '/bathy/', sep='')
#    dir.create(bathy.dir, recursive = TRUE)
#    bathy <- get.bath.data(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax, folder = bathy.dir)
#    #library(raster); plot(bathy)
#    # OR READ IT FROM NETCDF
#    #bathy.nc <- RNetCDF::open.nc(paste(bathy.dir, 'bathy.nc', sep=''))
#  
#    #------------
#    # CALCULATE LIKELIHOODS
#    #------------
#    # .par functions are the same calculations as those lacking .par, except they have been parallelized to leverage multiple CPUs
#    locs.grid <- setup.locs.grid(sp.lim)
#  
#    # vector indicating which likelihoods to run (e.g. 1=light, 2=sst, 5=hycom)
#    # can be combined with if() statements around calc functions: if (any(likVec == 5) & !exists('L.5')){calc.hycom(...)}
#    likVec <- c(1,2,5)
#  
#    # LIGHT-BASED LIKELIHOODS
#    #L.1 <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec, res=0.25) # if trying to use raw light levels, not currently recommended (v0.2)
#    L.1 <- calc.lightloc(locs, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE)
#  
#    #library(fields);library(raster)
#    #plot(L.1[[12]]); world(add=T)
#  
#    # SST LIKELIHOODS
#    #L.2 <- calc.sst(tag.sst, filename='oisst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
#    L.2 <- calc.sst(tag.sst, filename='oisst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
#    # save.image() # good idea to save after these larger calculations in case the next one causes problems
#    # gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running
#  
#    # PDT LIKELIHOODS
#    # OCEAN HEAT CONTENT (INTEGRATED PDTS)
#    L.3 <- calc.ohc.par(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
#    # save.image() # good idea to save after these larger calculations in case the next one causes problems
#    # gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running
#  
#    # WORLD OCEAN ATLAS-BASED LIKELIHOODS
#    L.4 <- calc.woa.par(pdt, ptt=ptt, woa.data = woa.quarter, sp.lim=sp.lim, focalDim = 9, dateVec = dateVec, use.se = T)
#    # save.image() # good idea to save after these larger calculations in case the next one causes problems
#    # gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running
#  
#    # HYCOM PROFILE BASED LIKELIHOODS
#    L.5 <- calc.hycom(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = as.POSIXct(dateVec), use.se = T)
#    # save.image() # good idea to save after these larger calculations in case the next one causes problems
#    # gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running
#  
#  }else{
#    load(paste0(dir,"likelihoods_example_20191014.rda"))
#  }
#  
#  #----------------------------------------------
#  # 2. Bathymetry and bottom temp based likelihoods
#  # (likely non sense here (for a pelagic) but just to illlustrate the computation & check if this work)
#  #----------------------------------------------
#  calc.lik.bottom=FALSE
#  if(calc.lik.bottom){
#    rm(L.5)
#    L.rasters <- mget(ls(pattern = 'L\\.')) # use with caution as all workspace items containing 'L.' will be listed. We only want the likelihood outputs calculated above
#    resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) terra::res(x)[1]))
#    L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
#  
#    dim(bathy);res(bathy)
#    # aggregate etopo bathy grid 0.08-->.25
#    b2=terra::resample(bathy, L.res[[1]][[1]])
#    dim(b2);res (b2)
#    # aggregate etopo bathy grid 0.08-->.25
#    b3=terra::resample(bathy, L.res$L.mle.res)
#    dim(b3);res(b3)
#  
#    # there is probably a way to avoid the loop here
#    tag.pdt=NULL
#    for(dd in unique(ac(pdt$Date))){
#      # dd=unique(ac(pdt$Date))[1]
#      xx=pdt[ac(pdt$Date)==dd,]
#      tag.pdt=rbind(tag.pdt,xx[which.max(xx$Depth),c('Date','Depth','MinTemp', 'MaxTemp')])
#    }
#  
#    # bathymetry based likelihood
#    L.6<-calc.bathy(tag.pdt,b3,dateVec, focalDim = 3, sens.err = 1)
#    L.62<-calc.bathy.par(tag.pdt,b3,dateVec,  focalDim = 3, sens.err = 1,ncores=4)
#  
#    # par(mfrow=c(1,2));plot(L.6[[2]],main='linear loop');plot(L.62[[2]],main='parallel')
#  
#    # bottom temperature based likelihood
#    L.7<-calc.bottomTemp(tag.pdt[1:4,],b3,dateVec[1:4], focalDim = 3, sens.err = 1,hycom.dir='C:/Users/pgatti/Documents/Halibut/R/functions/whoi_hmmoce_pkg/EnvData/hycom/',filename='hycom')
#    L.72<-calc.bottomTemp.par(tag.pdt[1:4,],b3,dateVec[1:4], focalDim = 3, sens.err = 1,hycom.dir='C:/Users/pgatti/Documents/Halibut/R/functions/whoi_hmmoce_pkg/EnvData/hycom/',filename='hycom')
#  
#    # par(mfrow=c(1,2));plot(L.7[[2]],main='linear loop');plot(L.72[[2]],main='parallel')
#  
#  }
#  #----------------------------------------------
#  
#  #----------------------------------------------
#  # 3. model fitting
#  #------------
#  # PREPARE TO RUN THE MODEL
#  #------------
#  
#  # some issue with L.5...
#  print("L.5 remove --> no way to display it")
#  rm(L.5)
#  
#  L.rasters <- mget(ls(pattern = 'L\\.')) # use with caution as all workspace items containing 'L.' will be listed. We only want the likelihood outputs calculated above
#  resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) terra::res(x)[1]))
#  L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
#  
#  # Figure out appropriate L combinations
#  # use this if you have a vector (likVec) indicating which likelihoods you are calculating
#  # for example, likVec <- c(1,2,5) for light, sst, and hycom likelihoods
#  if (length(likVec) > 2){
#    L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
#  } else{
#    L.idx <- utils::combn(likVec, 2, simplify=F)
#  }
#  
#  run.idx <- c(1,2,4)
#  
#  # vector of appropriate migr kernel speed. see ?makePar for more info.
#  parVec <- c(2, 4)
#  
#  # GOOD IDEA TO CLEAN THINGS UP AND SAVE
#  #rm(list=c('L.1','L.2','L.3','L.4','L.5', 'woa.quarter'))
#  # setwd(); base::save.image('.rda')
#  
#  #------------
#  # RUN THE MODEL
#  #------------
#  # CAN BE PARALLELIZED...
#  #require(foreach)
#  #print('Processing in parallel... ')
#  #ncores <- ceiling(parallel::detectCores() * .25)
#  #cl = parallel::makeCluster(ncores)
#  #doParallel::registerDoParallel(cl, cores = ncores)
#  #ans = foreach::foreach(tt = run.idx) %dopar%{
#  
#  #for (tt in run.idx){
#  #  for (bnd in bndVec){
#  #    for (i in parVec){
#  
#  tt=1
#  i=parVec[1]
#  
#  ptt <- 141259
#  #runName <- paste(ptt,'_idx',tt,'_bnd',bnd,'_par',i,sep='')
#  
#  # COMBINE LIKELIHOOD MATRICES
#  # L.idx combination indicates likelihood surfaces to consider
#  L <- make.L(L.res[[1]][L.idx[[tt]]],
#              iniloc = iniloc, dateVec = dateVec)
#  
#  
#  L <- L$L
#  g <- L.res$g
#  lon <- g$lon[1,]
#  lat <- g$lat[,1]
#  
#  rm(L.1,L.2,L.3,L.4,L.5,L.6)
#  #rm(L.res,L.rasters,locs.grid,bathy)
#  
#  #-----------------------------------------------
#  # Parameter estimation
#  #-----------------------------------------------
#  # a. nlminb or optim (gradient based)
#  # b. GA genetic algorithm
#  #-----------------------------------------------
#  
#  #-----------------------------------------------
#  # a. nlminb or optim
#  #-----------------------------------------------
#  maskL=FALSE
#  
#  # pars sigma1 (resident), sigma2 (migrant), p11 switch from behavior 1 to 1 (remain resident), p22
#  # c(1,2)/1000 * 3600 * 24/111/g$dla/(.5*pi)
#  pars.init=c(2,.2,.6,.8)
#  lower.bounds=c(0.1,0.001,.1,.1)
#  upper.bounds=c(5,.5,.9,.9)
#  
#  alg.opt='optim'
#  if(alg.opt=='optim'){
#    t0=Sys.time()
#    res.optim=optim(par=pars.init, fn=neg.log.lik.fun,g=g.mle,L=L.mle,maskL=maskL, gr=NULL,
#                    lower=lower.bounds, upper=upper.bounds,
#                    method='L-BFGS-B',hessian=FALSE)
#    # hessian=T --> compute hessian matrix for estimate of parameter SD --> usually slow down the computations
#    t1=Sys.time()
#    print(t1-t0)
#    save(res.optim,file=paste0('optim.',ptt,'_idx',tt,'_bnd',bnd,'_maskL',maskL,'.mle.RData'))
#  
#    t0=Sys.time()
#    res.optim=optim(par=pars.init, fn=neg.log.lik.fun,g=g,L=L,maskL=maskL, gr=NULL,
#                    lower=lower.bounds, upper=upper.bounds,
#                    method='L-BFGS-B',hessian=FALSE)
#    t1=Sys.time()
#    print(t1-t0)
#    save(res.optim,file=paste0('optim.',ptt,'_idx',tt,'_bnd',bnd,'_maskL',maskL,'.RData'))
#  }
#  
#  alg.opt='nlminb'
#  if(alg.opt=='nlminb'){
#    t0=Sys.time()
#    #res.nlminb=nlminb(pars.init, objective=neg.log.lik.fun, gradient = NULL, hessian = NULL,g=g,L=L,
#    #                  scale = 1, lower = 0, upper = max.sigma)
#    res.nlminb=tryCatch({nlminb(pars.init, objective=neg.log.lik.fun,maskL=maskL, gradient = NULL, hessian = NULL,
#                                g=g.mle,L=L.mle,
#                                scale = 1, lower=lower.bounds, upper=upper.bounds,
#                                control=list(iter.max=100))},error=function(error_message){return(error_message)})
#    t1=Sys.time()
#    print(t1-t0)
#    save(res.optim,file=paste0('nlminb.',ptt,'_idx',tt,'_bnd',bnd,'_maskL',maskL,'.mle.RData'))
#  
#    # print('nlmib has failed with the acurate grid ?!')
#    t0=Sys.time()
#    #res.nlminb=nlminb(pars.init, objective=neg.log.lik.fun, gradient = NULL, hessian = NULL,g=g,L=L,
#    #                  scale = 1, lower = 0, upper = max.sigma)
#    res.nlminb=tryCatch({nlminb(pars.init, objective=neg.log.lik.fun,maskL=maskL, gradient = NULL, hessian = NULL,
#                                g=g,L=L,
#                                scale = 1, lower=lower.bounds, upper=upper.bounds,
#                                control=list(iter.max=100))},error=function(error_message){return(error_message)})
#    t1=Sys.time()
#    print(t1-t0)
#    save(res.optim,file=paste0('nlminb.',ptt,'_idx',tt,'_bnd',bnd,'_maskL',maskL,'.RData'))
#  }
#  
#  #-----------------------------------------------
#  # b. GA genetic algorithm
#  #-----------------------------------------------
#  # alg.opt='GA'
#  # evolutionnary algorithm, completely different approach from gradient based such as optim or nlminb
#  # do not need initialization value,
#  # although you can still provide a parameter set that might work,
#  # and the algorithm will integrate it in its initial population of parameter set
#  # IT IS MUCH MORE time/computing consumming BUT also more resilient to local minima
#  # has a common feature with SEM, i.e. random sampling (mutations)
#  
#  alg.opt='GA'
#  if(alg.opt='GA'){
#  
#    pm=.1 # probability of mutation # default value of GA
#    mit=100 # maximum number of generations (i.e. iterations)
#    Run=100 # number of iteration without improvement befor ethe algorithm is stopped
#    pSz=100 # population size (each individual is a parameter set)
#    Ncores=4 # number of core to use for the parallelization, usually its good to use half of the cores (and let some for the background processes)
#  
#   # pm=.1 # probability of mutation # default value of GA
#  #  mit=1 # maximum number of generations (i.e. iterations)
#  #  Run=20 # number of iteration without improvement befor ethe algorithm is stopped
#  #  pSz=50 # population size (each individual is a parameter set)
#  #  Ncores=4 # number of core to use for the parallelization, usually its good to use half of the cores (and let some for the background processes)
#  
#    print('GA with the heavy grid might be long to run (~7 hours here)')
#    t0=Sys.time()
#    res.ga=ga(type='real-valued',fitness=pos.log.lik.fun,
#              g=g.mle,L=L.mle,maskL=maskL,
#              lower=lower.bounds,
#              upper=upper.bounds,
#              popSize=pSz,#crossover=gareal_blxCrossover,
#              maxiter=mit,run=Run,#pmutation=pm,
#              names=c('sigma1.ncell','sigma2.ncell','pswitch11','pswitch22'),#suggestions=2,
#              keepBest=T,
#              parallel=Ncores#'snow'
#    )
#    t1=Sys.time()
#    print(t1-t0)
#    save(res.ga,file=paste0('GA.',ptt,'_idx',tt,'_bnd',bnd,'_maskL',maskL,'.mle.RData'))
#  
#  
#    pm=.1 # probability of mutation # default value of GA
#    mit=100 # maximum number of generations (i.e. iterations)
#    Run=30 # number of iteration without improvement befor ethe algorithm is stopped
#    pSz=100 # population size (each individual is a parameter set)
#    Ncores=4
#  
#    t0=Sys.time()
#    res.ga=ga(type='real-valued',fitness=pos.log.lik.fun,
#              g=g,L=L,maskL=maskL,
#              lower=lower.bounds,
#              upper=upper.bounds,
#              popSize=pSz,#crossover=gareal_blxCrossover,
#              maxiter=mit,run=Run,#pmutation=pm,
#              names=c('sigma1.ncell','sigma2.ncell','pswitch11','pswitch22'),#suggestions=2,
#              keepBest=T,
#              parallel=Ncores#'snow'
#    )
#    t1=Sys.time()
#    print(t1-t0)
#    save(res.ga,file=paste0('GA.',ptt,'_idx',tt,'_bnd',bnd,'_maskL',maskL,'.RData'))
#  }
#  
#  optim.post.ga=FALSE
#  if(optim.post.ga){
#    pars.init=an(res.ga@solution)
#    t0=Sys.time()
#    res.optim=optim(par=pars.init, fn=neg.log.lik.fun,g=g,L=L,maskL=maskL, gr=NULL,
#                    lower=lower.bounds, upper=upper.bounds,
#                    method='L-BFGS-B',hessian=FALSE)
#    t1=Sys.time()
#    print(t1-t0)
#    save(res.optim,file=paste0('postGA.',ptt,'_idx',tt,'_bnd',bnd,'_maskL',maskL,'.RData'))
#  
#  }
#  
#  # duration filter with and without maskL
#  test=FALSE
#  if(test){
#    t0=Sys.time()
#    #ll=pos.log.lik.fun(pars,g.mle,L.mle)
#    f. <- hmm.filter(g.mle, L.mle, K1,K2, P)
#    t1=Sys.time()
#    print(t1-t0)
#  
#    t0=Sys.time()
#    #ll=pos.log.lik.fun(pars,g.mle,L.mle)
#    f.2 <- hmm.filter(g.mle, L.mle, K1,K2, P,maskL=F)
#    t1=Sys.time()
#    print(t1-t0)
#  }
#  
#  #----------------------------------
#  # filter / smooth
#  #----------------------------------
#  alg.opt='optim'
#  
#  # parameter set
#  if(alg.opt='GA'){
#    load(paste0('GA.',ptt,'_idx',tt,'_bnd',bnd,'_maskL',maskL,'.RData'))
#    plot(res.ga)
#    pars=an(res.ga@solution)
#    pars
#  }
#  
#  if(alg.opt=='nlminb'){
#    load(paste0('nlminb.',ptt,'_idx',tt,'_bnd',bnd,'_mask',maskL,'.RData'))
#    pars=res.nlminb$par
#    pars
#  }
#  
#  if(alg.opt=='optim'){
#    load(paste0('optim.',ptt,'_idx',tt,'_bnd',bnd,'_maskL',maskL,'.RData'))
#    pars=res.optim$par
#    pars
#  }
#  # pars=pars.init
#  
#  sigmas=pars[1:2]
#  sizes=ceiling(sigmas*4)
#  pb=pars[3:4]
#  muadvs=c(0,0)
#  
#  # if(!is.na(pars[5])) muadvs[1]<-pars[5]
#  #  if(!is.na(pars[6])) muadvs[2]<-pars[6]
#  
#  gausskern.PG=F
#  if(gausskern.PG){
#    # behav 1
#    if(sizes[1]%%2==0){sizes[1]=sizes[1]+1}
#    ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
#    while(ss<.999){
#      sizes[1]=sizes[1]+2
#      ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
#    }
#    K1 = gausskern.pg(sizes[1],sigmas[1],muadv=muadvs[1])
#    rm(ss)
#    K1=mask.K(K1)
#  
#    # behav 2
#    if(sizes[2]%%2==0){sizes[2]=sizes[2]+1}
#    ss=sum(gausskern.nostd(sizes[2],sigmas[2],muadv=muadvs[2]))
#    while(ss<.999){
#      sizes[2]=sizes[2]+2
#      ss=sum(gausskern.nostd(sizes[2],sigmas[2],muadv=muadvs[2]))
#    }
#    K2 = gausskern.pg(sizes[2],sigmas[2],muadv=muadvs[2])
#    rm(ss)
#    K2=mask.K(K2)
#  
#  }else{
#    sizes=rep(ceiling(sigmas[1]*4),2)
#    K1 = gausskern.pg(sizes[1],sigmas[1],muadv=muadvs[1])
#    K2 = gausskern.pg(sizes[2],sigmas[2],muadv=muadvs[2])
#  }
#  
#  P <- matrix(c(pb[1],1-pb[1],1-pb[2],pb[2]),2,2,byrow=TRUE)
#  
#  # RUN THE FILTER STEP
#  t0=Sys.time()
#  f <- hmm.filter(g, L, K=list(K1, K2), P = P, m = 2)
#  
#  #if(!is.na(bnd)){
#  #  f <- hmm.filter(g, L, K1, K2, maskL=T, P, minBounds = bnd)
#  #  maskL.logical <- TRUE
#  #} else{
#  # f <- hmm.filter(g, L, K1, K2, P, maskL=F)
#  # # f2 <- hmm.filter(g.mle, L.mle, K1, K2, P, maskL=F)
#  #  maskL.logical <- FALSE
#  #}
#  t1=Sys.time()
#  t1-t0
#  
#  # RUN THE SMOOTHING STEP
#  t0=Sys.time()
#  s <- hmm.smoother(f, K = list(K1, K2), L, P)
#  t1=Sys.time()
#  t1-t0
#  
#  # GET THE MOST PROBABLE TRACK
#  tr <- calc.track(s, g, dateVec, iniloc)
#  #setwd(myDir);
#  plotHMM(s, tr, dateVec, ptt=runName, save.plot = F)
#  
#  # WRITE OUT RESULTS
#  #outVec <- matrix(c(ptt=ptt, minBounds = bnd, migr.spd = i,
#  #                   Lidx = paste(L.idx[[tt]],collapse=''), P1 = P[1,1], P2 = P[2,2],
#  #                   spLims = sp.lim[1:4], resol = terra::res(L.rasters[[resamp.idx]]),
#  #                   maskL = maskL, NLL = nllf, name = runName), ncol=15)
#  #write.table(outVec,paste(dataDir, 'outVec_results.csv', sep=''), sep=',', col.names=F, append=T)
#  #names(outVec) <- c('ptt','bnd','migr.spd','Lidx','P1','P2','spLims','resol','maskL','nll','name')
#  #res <- list(outVec = outVec, s = s, g = g, tr = tr, dateVec = dateVec, iniloc = iniloc, grid = terra::res(L.res[[1]]$L.5)[1])
#  #setwd(myDir); save(res, file=dir,paste(runName, '-HMMoce_res.rda', sep=''))
#  #save.image(file=paste(ptt, '-HMMoce.RData', sep=''))
#  #source('~/HMMoce/R/hmm.diagnose.r') # not yet functional
#  #hmm.diagnose(res, L.idx, L.res, dateVec, locs.grid, iniloc, bathy, pdt, plot=T)
#  
#  #write.table(outVec, file='HMMoce_results_outVec.csv', sep=',', append=T)
#  
#  #    } # parVec loop
#  #  } # bndVec loop
#  #} # L.idx loop
#  
#  
#  #parallel::stopCluster(cl)
#  #closeAllConnections()

