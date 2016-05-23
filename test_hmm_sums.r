# fix the last day's weird probability sum...

load("~/GitHub/CAM_DATA/sword795_example.RData")

# use better rparameters

par0= c(12,36,12,6,.707,.866)
D1 <- par0[1:2] # parameters for kernel 1. this is behavior mode transit
D2 <- par0[3:4] # parameters for kernel 2. resident behavior mode
p <- par0[5:6]
#----------------------------------------------------------------------------------#
# GENERATE MOVEMENT KERNELS. D VALUES ARE MEAN AND SD PIXELS
K1 = as.cimg(gausskern(D1[1], D1[2], muadv = 0))
K2 = as.cimg(gausskern(D2[1], D2[2], muadv = 0))
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

#----------------------------------------------------------------------------------#
# RUN THE FILTER STEP
f2 = hmm.filter_test(g,L,K1,K2,P)
# res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
# fields::image.plot(lon, lat, res/max(res), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# RUN THE SMOOTHING STEP
s2 = hmm.smoother_test(f2, K1, K2, P, plot = F)
apply(s2, c(2), sum)

# PLOT IT IF YOU WANT TO SEE LIMITS (CI)
sres = apply(s2[1,,,], 2:3, sum, na.rm=T)
fields::image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# GET THE MOST PROBABLE TRACK
#----------------------------------------------------------------------------------#
distr = s2
meanlat <- apply(apply(distr, c(2, 4), sum) * repmat(t(as.matrix(g$lat[,1])), T, 1), 1, sum)
meanlon <- apply(apply(distr, c(2, 3), sum) * repmat(t(as.matrix(g$lon[1,])), T, 1), 1, sum)

lines(meanlon, meanlat, col = 2)


# Things to look at... 
# L
# L.locs
# 2:(T-1); (T-1:2) in filter/smooth.. need to leave the first/last day alone while updating all others. but, we want to use first/last days info going forward and backward...

