### trace of a matrix
tr<-function(x)
{
  return(sum(diag(x)))
}

### calculate the RMSE of beta
RMSE<-function(x)
{
  sqrt(mean((x)^2))
}

### calculate the F-norm of matrix
F.norm<-function(x)
{
  sqrt(tr(t(x)%*%x))
}

### adjust the sign of the fac_hat in order to measure the RMSE
adjust.fac<-function(fac_hat, fac)
{
  for (i in 1:nrow(fac))
  {
    if ((cor(fac_hat[i,], fac[i, ]))<0)
      fac_hat[i,] = -fac_hat[i,]
  }
  return(fac_hat)
}


###performance of group estimation
ranindex=function(glocN,TglocN){
Tlabel=rep(0,length(TglocN[[1]]))
Elabel=rep(0,length(glocN[[1]]))
#true
nTgroup=length(TglocN)
for(i in 1:nTgroup){
  id=which(as.vector(TglocN[[i]]))
  Tlabel[id]=i 
}
#estimate 
nEgroup=length(glocN)
for(i in 1:nEgroup){
  id=which(as.vector(glocN[[i]]))
  Elabel[id]=i 
}
return(adjustedRandIndex(Tlabel, Elabel))
}
