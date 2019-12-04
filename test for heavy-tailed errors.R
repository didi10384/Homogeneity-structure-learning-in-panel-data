# X (Time*N*p)
# Y (Time*N)
library(e1071)
Test.HeavyTail=function(X,Y){
  N=dim(X)[2]
  p=dim(X)[3]
  X.EK=matrix(0,N,p)
  Y.EK=numeric(N)
  for(i in 1:N){
     Y.EK[i]=kurtosis(Y[,i])  
     for(j in 1:p){
     X.EK[i,j]=kurtosis(X[,i,j])  
   }
  }
  return (list(X.EK=X.EK,Y.EK=Y.EK))
}