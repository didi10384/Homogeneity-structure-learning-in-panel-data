# X (Time*N*p)
# Y (Time*N)
library(tseries)
Test.Stationarity=function(X,Y){
  N=dim(X)[2]
  p=dim(X)[3]
  X.pvalue=matrix(0,N,p)
  Y.pvalue=numeric(N)
  for(i in 1:N){
      Y.pvalue[i]=adf.test(Y[,i])$p.value
    for(j in 1:p){
      X.pvalue[i,j]=adf.test(X[,i,j])$p.value
    }
  }
  return (list(X.pvalue=X.pvalue, Y.pvalue=Y.pvalue))
}
  