############################ Input ############################
# X (Time*N*p)
# Y (Time*N)
# covariance ("sample" or "robust")
# regression ("ols" or "robust")
# factor (default=NULL or a q*Time matrix)
# glocNT (default=NULL or a list of true structures)
############################ Output ############################
# beta_hat1(estimated beta before homogeneity detection)(N*p)
# beta_hat2(estimated beta after homogeneity detection)(N*p)
# q (estimated number of factors)
# B_hat (estimated factor loadings)
# fac (estimated factors if not provided)
# covZ(estimated covariance matrix for averaged Zt)
# glocN(grouping location: list of true/false for each group) 
# ng (number of groups)
# lambda (factor coefficients) (N*q)
# alpha (intercept)(N*1)
############################ source the estimtion functions
source("estimator/huber_loss.R")
source("estimator/huber_regression.R")
source("estimator/search_tau_lm_huber.R")
source("estimator/cv_lm_huber.R")
source("estimator/auxiliary.R")
source("estimator/cca.R")
source("estimator/solve_tau_eq.R")
source("estimator/beta_estimation.R")
source("estimator/group_estimation.R")
source("estimator/factor_estimation.R")

############################ load libraries
library(MASS)
library(wbs)
library(plyr)
library(mclust)

RobustPanel=function(X,Y,covariance,regression,factor=NULL,glocNT=NULL){
##### Step 1: get common parameters
Time=dim(X)[1]
N=dim(X)[2]
p=dim(X)[3]

##### Step 2: estimate beta for each individual 
#factor is unknown(default)
if (is.null(factor)){
  ### get Zt (averaged X over all individuals for each time)(Time*p)
  Zt = apply(X, c(1,3),mean) 
  
  ### estimate covariance matrix of Zt
  if (covariance == "sample")
  {
    covZ=cov(Zt)
  }
  if (covariance == "robust")
  {
    covZ=cov.robust(Zt, verbose=F)$covz
  }
  
  ### estimate number of factors and factor loading B
  estB=est.B(covZ)
  B_hat = estB$B_hat
  q=estB$q

  ### estimate factor 
  if (covariance == "sample")
  {
    fac=est.fac.ols(Zt, B_hat) # ols method
  }
  if (covariance == "robust")
  {
    fac_tau = est.fac.robust(Zt, B_hat, n_tau = Time, verbose = F) # robust method
    fac = fac_tau$fac_hat
  } 
  
  ### estimate beta before homogeneity detection
  if (regression == "ols")
  { 
    Mbefore=panel.factor.ols(X, Y, fac,Time,N, p, q)
    beta_hat1=Mbefore$beta_ols
    lambda=Mbefore$lambda_ols
    alpha=Mbefore$alpha_ols
    Y_defac=Mbefore$Y_defac
  }
  if (regression == "robust")
  {
    Mbefore=panel.factor.robust(X, Y, fac, Time,N, p, q, n_tau = Time)
    beta_hat1=Mbefore$beta_robust
    lambda=Mbefore$lambda_robust
    alpha=Mbefore$alpha_robust
    Y_defac=Mbefore$Y_defac
  }
}
# factor is known (oracal)
else
{ 
  fac=factor
  q=nrow(fac)
  B_hat=NULL
  covZ=NULL
  if (regression == "ols")
  { 
    Mbefore=panel.factor.ols(X, Y, fac,Time,N, p, q)
    beta_hat1=Mbefore$beta_ols
    lambda=Mbefore$lambda_ols
    alpha=Mbefore$alpha_ols
    Y_defac=Mbefore$Y_defac
  }
  if (regression == "robust")
  {
    Mbefore=panel.factor.robust(X, Y, fac, Time,N, p, q, n_tau = Time)
    beta_hat1=Mbefore$beta_robust
    lambda=Mbefore$lambda_robust
    alpha=Mbefore$alpha_robust
    Y_defac=Mbefore$Y_defac
  }
}
  
##### Step 3: homogeneity detection
#version 1: grouping beta together 
if (is.null(glocNT)){
G=6
betaGroup=beta.group(beta_hat1,G=G)
glocN=betaGroup$glocN
gmat=betaGroup$g_mat
ng=betaGroup$ng
}
else{
glocN=glocNT
ng=length(glocN)
}

##### Step 4: estimate beta after homogeneity detection
if (regression == "ols")
{ 
  Mafter=panel.factor.ols.group(X, Y_defac, glocN, ng, N, Time)
  beta_hat2=Mafter$beta_ols
}
if (regression == "robust")
{
  Mafter=panel.factor.robust.group(X, Y_defac, glocN, ng, N, Time)
  beta_hat2=Mafter$beta_robust
}

return(list(beta_hat1=beta_hat1,beta_hat2=beta_hat2,
            lambda=lambda,alpha=alpha,
            #lambda2=lambda2,alpha2=alpha2,
            q=q,B_hat=B_hat,fac=fac,
            covZ=covZ,glocN=glocN,ng=ng))
}