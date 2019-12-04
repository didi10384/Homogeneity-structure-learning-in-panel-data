############################ Input ############################
# X (Time*N*p)
#alpha (N*1)
#beta (N*p)
#factor(q*Time)
#lambda (N*q)
############################ Output ############################
# Y (Time*N)
library(pracma) 
prediction=function(X,alpha,beta,factor,lambda){
  Time=dim(X)[1]
  N=dim(X)[2]
  Y=matrix(0,nrow=Time,ncol=N)
  for (i in 1:Time){
    for (j in 1:N){
      Y[i,j]=alpha[j]+dot(X[i,j,],beta[j,])+dot(factor[,i],lambda[j,])
    }
  }
  return (Y)
}