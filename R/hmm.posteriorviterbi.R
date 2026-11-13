#########################
#Viterbi most probable track
##########################

###########################################################
#This code finds the most probable sequence of grid cells through the posteriors of HMMoce
###########################################################
#' Collaborative development by Drs. Julie Nielsen, Martin Arostegui, Camrin Braun, and Ben Galuardi
#' Cite: 
#' Nielsen JK, Bryan DR, Rand KM, Arostegui MC, Braun CD, Galuardi B, McDermott SF (2023) Geolocation of a demersal fish (Pacific cod) in a high-latitude island chain (Aleutian Islands, Alaska). Animal Biotelemetry 11: 29.
#'
#' Parameters
#' @param pars HMM parameters output from opt.params()
#' @param s HMM backward smoothed output of hmm.smoother()
#' @param g original grid on which HMM was run, assigned from L.res$g
#' @param iniloc dataframe of tag and popup locations/dates (2 rows, 5 columns)
#' @param dateVec timestep sequence used in HMM
#' @param lims.lat vector of c(min,max) desired latitude for reduced spatial extent in which to run algorithm
#' @param lims.lon vector of c(min,max) desired longitude for reduced spatial extent in which to run algorithm
#' @param percentile uncertainty percentile of the daily posterior to include in admissible paths (default is 99.9% of daily posterior)
#' @param approach for two-state outputs, use a priori most likely behavior state sequence with separate posteriors and kernels ('separate'), or joint posteriors and kernels weighted by relative behavior state probability ('joint')
#' 
hmm.posteriorviterbi <- function(pars, s, g, iniloc, dateVec, lims.lat, lims.lon, percentile = 99.9, approach = "separate"){

#### For two-state, join behavior state-specific posteriors and kernels weighted by daily relative probability ####
if (approach == "joint") {
  # calculate the estimated behavior state
  if(dim(s)[1] > 1){
    p.resid <- apply(s, c(1,2), sum)[2,]
    p.migrant <- 1-p.resid
  } else{
    p.resid <- rep(NA, length.out = length(dateVec))
  }
  
  if (dim(s)[1] == 2){
    p.1 <- apply(s[1,,,], 1, sum)
    p.2 <- apply(s[2,,,], 1, sum)
    sv <- -(apply(s[1,,,], 1, sum) > apply(s[2,,,], 1, sum)) + 1
    #sv[sv == 0] <- NA
    
    # normalize for each behav at each time point
    s_norm <- array(NA, dim=dim(s)[c(3,4,2,1)])
    for (i in 1:dim(s_norm)[3]){
      s_norm[,,i,1] <- (s[1,i,,] / max(s[1,i,,], na.rm=T)) * p.1[i]
      s_norm[,,i,2] <- (s[2,i,,] / max(s[2,i,,], na.rm=T)) * p.2[i]
    }
    
    
    ## sum across behavior states
    s_norm <- rowSums(s_norm, dims=3, na.rm=TRUE)
    s_norm <- aperm(s_norm,c(3,1,2))
    
    ## set up the dims for the s output needed in viterbi
    T <- dim(s)[2] # dimension of time
    col <- dim(s)[3] # ncols
    row <- dim(s)[4] # nrows
    
    s_viterbi <- array(0, dim = c(1, T, col, row))
    for (tt in 1:T){
      s_viterbi[1,tt,,] <- s_norm[tt,,]
    }
    
    # replace original 2-state posteriors with across-state, weighted posterior
    s <- s_viterbi
  }
  
  # normalize the posteriors from 0-1 for each day
  s_idx = s
  for(i in 1:(dim(s)[2])){
    s_idx[1, i , ,] = imager::renorm(s[1, i, , ], 0, 1)
  }
  
  # set < X% posterior probability to 0
  s[s_idx < (1-percentile/100)] = 0
  
  #Crop s to outer spatial extent of overall UD of focal deployment (degrees lat/lon) [manually determined with plotRD]
  min.lat.crop<- lims.lat[1]
  max.lat.crop<- lims.lat[2]
  min.lon.crop<- lims.lon[1]
  max.lon.crop<- lims.lon[2]
  
  # Latitudes to remove
  crop.lat.index <- which(g$lat[, 1] < min.lat.crop | g$lat[, 1] > max.lat.crop)
  # Longitudes to remove
  crop.lon.index <- which(g$lon[1, ] < min.lon.crop | g$lon[1, ] > max.lon.crop)
  
  # Crop s (spatially)
  if (length(crop.lat.index)==0 & length(crop.lon.index)==0) {
    L.cropped  <- s[,,,,drop=F]
  } else if (length(crop.lat.index)==0 & !length(crop.lon.index)==0) {
    L.cropped  <- s[,,-crop.lon.index,,drop=F]
  } else if (!length(crop.lat.index)==0 & length(crop.lon.index)==0) {
    L.cropped  <- s[,,,-crop.lat.index,drop=F]
  } else {
    L.cropped  <- s[,,-crop.lon.index,-crop.lat.index,drop=F]
  }
  
  # Match HMMoce L format to Nielsen L format (ready for Viterbi)
  L.vit <- aperm(L.cropped, c(4,3,2,1))
  L.vit <- aperm(L.vit,c(2,1,3,4))
  L.vit <- array(reverse(reverse(reverse(apply(L.vit,2,c)))), dim = c(dim(L.vit)[2],dim(L.vit)[1],dim(L.vit)[3],dim(L.vit)[4]))
  
  # Generate other Viterbi-specific objects
  if (length(crop.lon.index)==0) {
    longs.vit <- g$lon[1, ]
  } else {
    longs.vit <- g$lon[1, -crop.lon.index]
  }
  
  if (length(crop.lat.index)==0) {
    lats.vit <- g$lat[ ,1]
  } else {
    lats.vit <- g$lat[-crop.lat.index, 1]
  }
  rows.vit <- dim(L.vit)[1]
  cols.vit <- dim(L.vit)[2]
  timesteps.vit <- dim(L.vit)[3]
  behavs.vit <- dim(L.vit)[4]
  # Identify starting cell
  ilo.vit <- which.min(abs(iniloc[1,"lon"]-longs.vit))
  ila.vit <- which.min(abs(iniloc[1,"lat"]-rev(lats.vit)))
  
  # Generate kernels normally
  if (length(pars) == 4){ # if two-state model
    sigmas = pars[1:2]
    sizes = rep(ceiling(sigmas[1]*6),2)
    pb = pars[3:4]
    muadvs = c(0,0)
    K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
    K2 <- gausskern.pg(sizes[2], sigmas[2], muadv=muadvs[2])
    K <- list(K1,K2)
  } else if (length(pars) == 1){ # if one-state model
    sigmas = pars[1]
    sizes = rep(ceiling(sigmas[1]*6),2)
    pb = NULL
    muadvs = c(0)
    K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
    K <- list(K1)
  }
  
  # Ensure diffusion kernel has a central cell (requires odd dimensions)
  if (is.na(p.resid[1])) {
    # One-state model
    if (sizes[1] %% 2 == 0) { # if kernel has even dimensions
      sizes <- sizes + 1
      pb = NULL
      muadvs = c(0)
      K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
      K <- list(K1)
    } else {}
  } else {
    # Two-state model
    if (sizes[1] %% 2 == 0 | sizes[2] %% 2 == 0) { # if kernel(s) have even dimensions
      if (sizes[1] %% 2 == 0) {sizes[1] <- sizes[1] + 1} # yields odd dimensions
      if (sizes[2] %% 2 == 0) {sizes[2] <- sizes[2] + 1} # yields odd dimensions
      pb = pars[3:4]
      muadvs = c(0,0)
      K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
      K2 <- gausskern.pg(sizes[2], sigmas[2], muadv=muadvs[2])
    } else {}
  }
  
  #mpt.calc: 4-dimensional array to hold pathway and probability of that pathway [x,y,time,index# of most probable previous cell]
  # and [x,y,time,cumulative probability]
  #Keeps track of the highest probability pathway between each grid cell and grid cells from previous day
  
  mpt.calc<- array(NA, dim=c(rows.vit,cols.vit,timesteps.vit,2)) #X,Y,days, pathway (which index of matrix highest), and probability at node
  mpt.calc[ila.vit,ilo.vit,1,1]<- which(row(mpt.calc[,,1,1])==ila.vit & col(mpt.calc[,,1,1])==ilo.vit) #index # of release location
  mpt.calc[ila.vit,ilo.vit,1,2]<- 1 #Probability at release location
  mpt.calc[,,1,2]<- log(mpt.calc[,,1,2])
  
  #Day 1: initial probability
  #image.plot(reverse(mpt.calc[,,1,1]))#index# of release location on first day (this is just a single pixel)
  #image.plot(reverse(mpt.calc[,,1,2]))#log probability of release location (this is just a single pixel)
  
  #i=2 #Recommend running this through a day at a time, starting at i=2, to understand how it works, etc.
  #Iterate through each day
  start.time<- Sys.time()
  for(i in 2:timesteps.vit){ #Start at 2
    # One-state model
    if (is.na(p.resid[1])) {
      kern <- K[[1]]
      z <- 1 # behavior state used in testing path admissibility for percentile uncertainty
    } else {
      # Two-state model transformed in single joint kernel, weighted by state probabilities
      kern <- K1*p.migrant[i] + K2*p.resid[i]
    }
    
    #Define how different kernels will fit into viterbi study area matrix
    kern.dim<- dim(kern)[1] #symmetrical
    kern.y.lower<- ceiling(kern.dim/2) #lowest row where kernel can be fully placed into larger matrix
    kern.y.upper<- rows.vit - floor(kern.dim/2) #largest row where kernel can be fully placed into larger matrix
    kern.x.lower<- ceiling(kern.dim/2) #lowest col where kernel can be fully placed into larger matrix
    kern.x.upper<- cols.vit - floor(kern.dim/2) #largest col where kernel can be fully placed into larger matrix
    kern.offset<- floor(kern.dim/2) #Value subtracted/added from each x and y when kernel fully fits into larger matrix
    
    #For each grid cell, find which previous day's grid cell has highest probability of transition
    L_day <- L.vit[,,i,1]
    valid_prev <- which(!is.na(mpt.calc[,,i-1,1]), arr.ind = TRUE)
    if (nrow(valid_prev) > 0) {
      for (k_idx in seq_len(nrow(valid_prev))) {
        y <- valid_prev[k_idx, 1]
        x <- valid_prev[k_idx, 2]
        if (y > kern.y.lower && y < kern.y.upper) {
          ymin <- y - kern.offset; ymax <- y + kern.offset; kern.y <- 1:kern.dim
        } else if (y <= kern.y.lower) {
          ymin <- 1; y.offset <- kern.y.lower - y; kern.y <- (1:kern.dim)[-(1:y.offset)]; ymax <- y + kern.offset
        } else {
          ymax <- rows.vit; y.offset <- y - kern.y.upper; kern.y <- 1:(kern.dim - y.offset); ymin <- y - kern.offset
        }
        if (x > kern.x.lower && x < kern.x.upper) {
          xmin <- x - kern.offset; xmax <- x + kern.offset; kern.x <- 1:kern.dim
        } else if (x <= kern.x.lower) {
          xmin <- 1; x.offset <- kern.x.lower - x; kern.x <- (1:kern.dim)[-(1:x.offset)]; xmax <- x + kern.offset
        } else {
          xmax <- cols.vit; x.offset <- x - kern.x.upper; kern.x <- 1:(kern.dim - x.offset); xmin <- x - kern.offset
        }
        if (length(ymin:ymax) == kern.dim) kern.y <- 1:kern.dim
        if (length(xmin:xmax) == kern.dim) kern.x <- 1:kern.dim
        
        sub_upd <- L_day[ymin:ymax, xmin:xmax, drop = FALSE] * kern[kern.y, kern.x, drop = FALSE]
        pos_mask <- !is.na(sub_upd) & (sub_upd > 0)
        if (any(pos_mask)) {
          update.index <- y + (x - 1L) * rows.vit
          cand_prob <- log(sub_upd) + mpt.calc[y, x, i-1, 2]
          cur_idx  <- mpt.calc[ymin:ymax, xmin:xmax, i, 1]
          cur_prob <- mpt.calc[ymin:ymax, xmin:xmax, i, 2]
          better <- pos_mask & (is.na(cur_idx) | (cand_prob > cur_prob))
          if (any(better)) {
            cur_idx[better]  <- update.index
            cur_prob[better] <- cand_prob[better]
            mpt.calc[ymin:ymax, xmin:xmax, i, 1] <- cur_idx
            mpt.calc[ymin:ymax, xmin:xmax, i, 2] <- cur_prob
          }
        }
      }
    }
    # Fallback if bounded kernel support did not reach any positive cell in L_day
    if (all(is.na(mpt.calc[,,i,1]))) {
      best_prev_lin <- which.max(mpt.calc[,,i-1,2])
      if (length(best_prev_lin) > 0) {
        bp <- arrayInd(best_prev_lin[1], c(rows.vit, cols.vit))
        pos_l <- which(!is.na(L_day) & L_day > 0, arr.ind = TRUE)
        if (nrow(pos_l) == 0) pos_l <- bp
        d2 <- (pos_l[, 1] - bp[1, 1]) ^ 2 + (pos_l[, 2] - bp[1, 2]) ^ 2
        sig_eff <- max(sigmas[1], 1)
        lp <- log(pmax(L_day[pos_l], 1e-300)) - 0.5 * (d2 / (sig_eff ^ 2)) + mpt.calc[bp[1, 1], bp[1, 2], i-1, 2]
        mpt.calc[,,i,1][pos_l] <- best_prev_lin[1]
        mpt.calc[,,i,2][pos_l] <- lp
      }
    }
    print (i)
  }
  Sys.time()-start.time
  
  #Reconstruct MPT from grid cell at last time step highest cumulative probability
  mpt<- matrix(NA,nrow=timesteps.vit,ncol=2) #1 = row, 2=col
  end.loc<- arrayInd(which.max(mpt.calc[,,timesteps.vit,2]),c(rows.vit,cols.vit)) #Grid cell with most cumulative probability at the end
  mpt[timesteps.vit,]<- end.loc
  
  prev.loc<- mpt.calc[end.loc[,1],end.loc[,2],timesteps.vit,1] #Index # of grid cell from previous step that had highest probability moving to end loc
  for(i in (timesteps.vit-1):1){
    mpt[i,]<- arrayInd(prev.loc,c(rows.vit,cols.vit)) #index number of grid cell with highest probability
    prev.loc<-  mpt.calc[mpt[i,1],mpt[i,2],i,1] #re-set prev loc to next time step
  }
  
  mpt.latlong<- matrix(NA,nrow=timesteps.vit,ncol=2) #Convert row/col to lat/long
  for(i in 1:timesteps.vit){
    mpt.latlong[i,1]<- lats.vit[((rows.vit+1)-mpt[i,1])]
    mpt.latlong[i,2]<- longs.vit[mpt[i,2]]
  }
  
  # Generate dataframe of most probable track by date, lat/lon, probability of resident state
  vit.mpt<- data.frame("date"=dateVec,"lat"=mpt.latlong[,1],"lon"=mpt.latlong[,2],"p.resid"=p.resid)
  return(vit.mpt)
}


#### For two-state, use a priori behavior state sequence with separate posteriors and kernels based on which state is most likely on each day ####
if (approach == "separate") {
  if (length(pars) == 4){ # if two-state model
    sigmas = pars[1:2]
    sizes = rep(ceiling(sigmas[1]*4),2)
    pb = pars[3:4]
    muadvs = c(0,0)
    K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
    if (!is.null(pb)) K2 <- gausskern.pg(sizes[2], sigmas[2], muadv=muadvs[2])
    K <- list(K1,K2)
  } else if (length(pars) == 1){ # if one-state model
    sigmas = pars[1]
    sizes = rep(ceiling(sigmas[1]*4),2)
    pb = NULL
    muadvs = c(0)
    K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
    K <- list(K1)
  }
  
  
  # normalize the posteriors from 0-1 for each day
  s_idx = s
  for(i in 1:(dim(s)[2])){
    if (dim(s)[1]==2) {
      s_idx[1, i , ,] = imager::renorm(s[1, i, , ], 0, 1)
      s_idx[2, i , ,] = imager::renorm(s[2, i, , ], 0, 1)
    } else {
      s_idx[1, i , ,] = imager::renorm(s[1, i, , ], 0, 1)
    }
    
  }
  
  # set < X% posterior probability to 0
  s[s_idx < (1-percentile/100)] = 0
  
  #Crop s to outer spatial extent of overall UD of focal deployment (degrees lat/lon) [manually determined with plotRD]
  min.lat.crop<- lims.lat[1]
  max.lat.crop<- lims.lat[2]
  min.lon.crop<- lims.lon[1]
  max.lon.crop<- lims.lon[2]
  
  # Latitudes to remove
  crop.lat.index <- which(g$lat[, 1] < min.lat.crop | g$lat[, 1] > max.lat.crop)
  # Longitudes to remove
  crop.lon.index <- which(g$lon[1, ] < min.lon.crop | g$lon[1, ] > max.lon.crop)
  
  # Crop s (spatially)
  if (length(crop.lat.index)==0 & length(crop.lon.index)==0) {
    L.cropped  <- s[,,,,drop=F]
  } else if (length(crop.lat.index)==0 & !length(crop.lon.index)==0) {
    L.cropped  <- s[,,-crop.lon.index,,drop=F]
  } else if (!length(crop.lat.index)==0 & length(crop.lon.index)==0) {
    L.cropped  <- s[,,,-crop.lat.index,drop=F]
  } else {
    L.cropped  <- s[,,-crop.lon.index,-crop.lat.index,drop=F]
  }
  
  # Match HMMoce L format to Nielsen L format (ready for Viterbi)
  L.vit <- aperm(L.cropped, c(4,3,2,1))
  L.vit <- aperm(L.vit,c(2,1,3,4))
  L.vit <- array(reverse(reverse(reverse(apply(L.vit,2,c)))), dim = c(dim(L.vit)[2],dim(L.vit)[1],dim(L.vit)[3],dim(L.vit)[4]))
  
  # Generate other Viterbi-specific objects
  if (length(crop.lon.index)==0) {
    longs.vit <- g$lon[1, ]
  } else {
    longs.vit <- g$lon[1, -crop.lon.index]
  }
  
  if (length(crop.lat.index)==0) {
    lats.vit <- g$lat[ ,1]
  } else {
    lats.vit <- g$lat[-crop.lat.index, 1]
  }
  rows.vit <- dim(L.vit)[1]
  cols.vit <- dim(L.vit)[2]
  timesteps.vit <- dim(L.vit)[3]
  behavs.vit <- dim(L.vit)[4]
  # Identify starting cell
  ilo.vit <- which.min(abs(iniloc[1,"lon"]-longs.vit))
  ila.vit <- which.min(abs(iniloc[1,"lat"]-rev(lats.vit)))
  
  # Determine behavior state sequence if two-state model used
  distr <- s
  if (dim(distr)[1]==2){
    sv <- -(apply(distr[1,,,], 1, sum) > apply(distr[2,,,], 1, sum)) + 1
    sv <- sv + 1 # two-state model, where 1 is migratory and 2 is resident
  } else {
    sv <- NA # one-state model
  }
  
  # Ensure diffusion kernel has a central cell (requires odd dimensions)
  if (is.na(sv[1])==TRUE) {
    # One-state model
    if (sizes[1] %% 2 == 0) { # if kernel has even dimensions
      sizes <- sizes + 1
      pb = NULL
      muadvs = c(0)
      K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
      K <- list(K1)
    } else {}
  } else {
    # Two-state model
    if (sizes[1] %% 2 == 0 | sizes[2] %% 2 == 0) { # if kernel(s) have even dimensions
      if (sizes[1] %% 2 == 0) {sizes[1] <- sizes[1] + 1} # yields odd dimensions
      if (sizes[2] %% 2 == 0) {sizes[2] <- sizes[2] + 1} # yields odd dimensions
      pb = pars[3:4]
      muadvs = c(0,0)
      K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
      if (!is.null(pb)) K2 <- gausskern.pg(sizes[2], sigmas[2], muadv=muadvs[2])
      K <- list(K1,K2)
    } else {}
  }
  
  #mpt.calc: 4-dimensional array to hold pathway and probability of that pathway [x,y,time,index# of most probable previous cell]
  # and [x,y,time,cumulative probability]
  #Keeps track of the highest probability pathway between each grid cell and grid cells from previous day
  
  mpt.calc<- array(NA, dim=c(rows.vit,cols.vit,timesteps.vit,2)) #X,Y,days, pathway (which index of matrix highest), and probability at node
  mpt.calc[ila.vit,ilo.vit,1,1]<- which(row(mpt.calc[,,1,1])==ila.vit & col(mpt.calc[,,1,1])==ilo.vit) #index # of release location
  mpt.calc[ila.vit,ilo.vit,1,2]<- 1 #Probability at release location
  mpt.calc[,,1,2]<- log(mpt.calc[,,1,2]) 
  
  #Day 1: initial probability
  #image.plot(reverse(mpt.calc[,,1,1]))#index# of release location on first day (this is just a single pixel)
  #image.plot(reverse(mpt.calc[,,1,2]))#log probability of release location (this is just a single pixel)
  
  #i=2 #Recommend running this through a day at a time, starting at i=2, to understand how it works, etc.
  #Iterate through each day
  start.time<- Sys.time()
  for(i in 2:timesteps.vit){ #Start at 2
    # One-state model
    if (is.na(sv[1])==TRUE) {
      kern <- K[[1]]
      z <- 1 # behavior state used in testing path admissibility for percentile uncertainty
    } else {
      # Two-state model
      if (sv[i]==1){
        kern <- K[[1]] #first movement state (HMMoce: migratory), D1, defined during forward filter/reverse smoothing
        z <- 1 # behavior state used in testing path admissibility for percentile uncertainty
      } else if (sv[i]==2) {
        kern <- K[[2]] #second movement state (HMMoce: resident), D2, defined during forward filter/reverse smoothing
        z <- 2 # behavior state used in testing path admissibility for percentile uncertainty
      }
    }
    
    #Define how different kernels will fit into viterbi study area matrix
    kern.dim<- dim(kern)[1] #symmetrical
    kern.y.lower<- ceiling(kern.dim/2) #lowest row where kernel can be fully placed into larger matrix
    kern.y.upper<- rows.vit - floor(kern.dim/2) #largest row where kernel can be fully placed into larger matrix
    kern.x.lower<- ceiling(kern.dim/2) #lowest col where kernel can be fully placed into larger matrix
    kern.x.upper<- cols.vit - floor(kern.dim/2) #largest col where kernel can be fully placed into larger matrix
    kern.offset<- floor(kern.dim/2) #Value subtracted/added from each x and y when kernel fully fits into larger matrix
    
    #For each grid cell, find which previous day's grid cell has highest probability of transition
    L_day <- if (behavs.vit == 1 || sv[i] == 1) L.vit[,,i,1] else L.vit[,,i,2]
    valid_prev <- which(!is.na(mpt.calc[,,i-1,1]), arr.ind = TRUE)
    if (nrow(valid_prev) > 0) {
      for (k_idx in seq_len(nrow(valid_prev))) {
        y <- valid_prev[k_idx, 1]
        x <- valid_prev[k_idx, 2]
        if (y > kern.y.lower && y < kern.y.upper) {
          ymin <- y - kern.offset; ymax <- y + kern.offset; kern.y <- 1:kern.dim
        } else if (y <= kern.y.lower) {
          ymin <- 1; y.offset <- kern.y.lower - y; kern.y <- (1:kern.dim)[-(1:y.offset)]; ymax <- y + kern.offset
        } else {
          ymax <- rows.vit; y.offset <- y - kern.y.upper; kern.y <- 1:(kern.dim - y.offset); ymin <- y - kern.offset
        }
        if (x > kern.x.lower && x < kern.x.upper) {
          xmin <- x - kern.offset; xmax <- x + kern.offset; kern.x <- 1:kern.dim
        } else if (x <= kern.x.lower) {
          xmin <- 1; x.offset <- kern.x.lower - x; kern.x <- (1:kern.dim)[-(1:x.offset)]; xmax <- x + kern.offset
        } else {
          xmax <- cols.vit; x.offset <- x - kern.x.upper; kern.x <- 1:(kern.dim - x.offset); xmin <- x - kern.offset
        }
        if (length(ymin:ymax) == kern.dim) kern.y <- 1:kern.dim
        if (length(xmin:xmax) == kern.dim) kern.x <- 1:kern.dim
        
        sub_upd <- L_day[ymin:ymax, xmin:xmax, drop = FALSE] * kern[kern.y, kern.x, drop = FALSE]
        pos_mask <- !is.na(sub_upd) & (sub_upd > 0)
        if (any(pos_mask)) {
          update.index <- y + (x - 1L) * rows.vit
          cand_prob <- log(sub_upd) + mpt.calc[y, x, i-1, 2]
          cur_idx  <- mpt.calc[ymin:ymax, xmin:xmax, i, 1]
          cur_prob <- mpt.calc[ymin:ymax, xmin:xmax, i, 2]
          better <- pos_mask & (is.na(cur_idx) | (cand_prob > cur_prob))
          if (any(better)) {
            cur_idx[better]  <- update.index
            cur_prob[better] <- cand_prob[better]
            mpt.calc[ymin:ymax, xmin:xmax, i, 1] <- cur_idx
            mpt.calc[ymin:ymax, xmin:xmax, i, 2] <- cur_prob
          }
        }
      }
    }
    # Fallback if bounded kernel support did not reach any positive cell in L_day
    if (all(is.na(mpt.calc[,,i,1]))) {
      best_prev_lin <- which.max(mpt.calc[,,i-1,2])
      if (length(best_prev_lin) > 0) {
        bp <- arrayInd(best_prev_lin[1], c(rows.vit, cols.vit))
        pos_l <- which(!is.na(L_day) & L_day > 0, arr.ind = TRUE)
        if (nrow(pos_l) == 0) pos_l <- bp
        d2 <- (pos_l[, 1] - bp[1, 1]) ^ 2 + (pos_l[, 2] - bp[1, 2]) ^ 2
        sig_eff <- max(sigmas[1], 1)
        lp <- log(pmax(L_day[pos_l], 1e-300)) - 0.5 * (d2 / (sig_eff ^ 2)) + mpt.calc[bp[1, 1], bp[1, 2], i-1, 2]
        mpt.calc[,,i,1][pos_l] <- best_prev_lin[1]
        mpt.calc[,,i,2][pos_l] <- lp
      }
    }
    print (i)
  }
  Sys.time()-start.time
  
  #Reconstruct MPT from grid cell at last time step highest cumulative probability
  mpt<- matrix(NA,nrow=timesteps.vit,ncol=2) #1 = row, 2=col
  end.loc<- arrayInd(which.max(mpt.calc[,,timesteps.vit,2]),c(rows.vit,cols.vit)) #Grid cell with most cumulative probability at the end
  mpt[timesteps.vit,]<- end.loc
  
  prev.loc<- mpt.calc[end.loc[,1],end.loc[,2],timesteps.vit,1] #Index # of grid cell from previous step that had highest probability moving to end loc 
  for(i in (timesteps.vit-1):1){
    mpt[i,]<- arrayInd(prev.loc,c(rows.vit,cols.vit)) #index number of grid cell with highest probability
    prev.loc<-  mpt.calc[mpt[i,1],mpt[i,2],i,1] #re-set prev loc to next time step
  }
  
  mpt.latlong<- matrix(NA,nrow=timesteps.vit,ncol=2) #Convert row/col to lat/long
  for(i in 1:timesteps.vit){
    mpt.latlong[i,1]<- lats.vit[((rows.vit+1)-mpt[i,1])]
    mpt.latlong[i,2]<- longs.vit[mpt[i,2]]
  }
  
  
  # calculate the estimated behavior state
  if(dim(s)[1] > 1){
    p.resid <- apply(s, c(1,2), sum)[2,]
  } else{
    p.resid <- rep(NA, length.out = length(dateVec))
  }
  
  # Generate dataframe of most probable track by date, lat/lon, probability of resident state
  vit.mpt<- data.frame("date"=dateVec,"lat"=mpt.latlong[,1],"lon"=mpt.latlong[,2],"p.resid"=p.resid)
  return(vit.mpt)
}

}
