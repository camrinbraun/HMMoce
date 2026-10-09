mask.K=function(K){
  #K=matrix(data=NA)
  coord.Kx=matrix(data=rep(1:ncol(K),each=nrow(K)),ncol=ncol(K),nrow=nrow(K))-(trunc(ncol(K)/2)+1)
  coord.Ky=matrix(data=rep(1:ncol(K),nrow(K)),ncol=ncol(K),nrow=nrow(K))-(trunc(ncol(K)/2)+1)
  coord=cbind(as.vector(coord.Kx),as.vector(coord.Ky))
  # euclidean distance with the point at the center of the kernel (coordinates 0 0)
  coord=apply(coord,1,function(x){return(sqrt(x[1]**2+x[2]**2))})
  coord=matrix(data=coord,ncol=ncol(K),nrow=nrow(K))
  coord[round(coord)>(trunc(ncol(K)/2))]<-NA
  mask=which(is.na(coord))
  K[mask]<-0
  K=K/sum(K,na.rm=T)
  return(K)
}