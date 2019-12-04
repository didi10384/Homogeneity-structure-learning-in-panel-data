### estimate the model before homogeneity detection 
panel.factor.robust<-function(X, Y, fac, Time,N, p, q, n_tau = 10)
{
  para_robust = matrix(0, nrow = N, ncol = p+q+1)
  taus = rep(0, n_tau)
  Y_defac = matrix(0, Time, N)
  ### estimate the parameters for i = 1, ..., N
  for (i in 1:N)
  {
    Yi = Y[,i]
    Xi = X[,i,]
    
    W = cbind(Xi, t(fac))
    
    ### Huber estimator of the paramters
    if (i <=n_tau)
    {
      taus[i] = search.tau.lm.huber(W, Yi, fit.intcp = T, opt=list(method='default'))
    }
    tau = mean(taus[1:min(i, n_tau)])
    para_robust[i,] = lm.huber.fixTau(W, Yi, tau)$Coefficients
    Y_defac[,i] = Yi - t(fac)%*%para_robust[i,(p+2):(p+q+1)] - para_robust[i,1]
  }
  return(list(lambda_robust = para_robust[,(p+2):(p+q+1)], 
              beta_robust = para_robust[,2:(p+1)],
              alpha_robust = para_robust[,1],
              tau = taus, Y_defac = Y_defac))
}

### estimte the model by ols estimtor before homogeneity detection 
panel.factor.ols<-function(X, Y, fac,Time,N, p, q)
{
  para_ols = matrix(0, nrow = N, ncol = p+q+1)
  Y_defac = matrix(0, Time, N)
  for (i in 1:N)
  {
    Yi = Y[,i]
    Xi = X[,i,]
    
    W = cbind(Xi, t(fac))
    #W1 = cbind(1, W)
    ### OLS
    #para_ols[i,] = solve(crossprod(W1))%*%crossprod(W1, Yi)
    para_ols[i,] = lm(Yi~W)$coefficients
    Y_defac[,i] = Yi - t(fac)%*%para_ols[i,(p+2):(p+q+1)] - para_ols[i,1]
  }
  return(list(lambda_ols = para_ols[,(p+2):(p+q+1)], 
              beta_ols = para_ols[,2:(p+1)],
              alpha_ols = para_ols[,1],
              Y_defac = Y_defac))
}

## estimate the model when the group information is known
panel.factor.robust.group<-function(X, Y_defac, glocN, ng, N, Time)
{
  ### obtain sum of all subjects at the same time point t
  Y=apply(Y_defac,1,sum)
  X_new=matrix(0,Time,ng)
  for(t in 1:Time)
  {
    Xt=X[t,,]
    X_new[t,]=sapply(glocN, function(x){sum(x*Xt)})
  }
  tau = search.tau.lm.huber(X_new, Y, fit.intcp = F, opt=list(method='default'))
  beta_robust = lm.huber.fixTau(X_new, Y, tau, fit.intcp = F)$Coefficients[-1]
  beta_l = lapply(1:length(beta_robust), function(i) beta_robust[i]*glocN[[i]])
  beta_N = Reduce("+", beta_l)

  return(list(beta_robust = beta_N, beta_robustk = beta_robust, tau = tau))
}


### estimate the model by ols when the group information is known
panel.factor.ols.group<-function(X, Y_defac, glocN, ng, N, Time)
{
  ### obtain sum of all subjects at the same time point t
  Y=apply(Y_defac,1,sum)
  X_new=matrix(0,Time,ng)
  for(t in 1:Time)
  {
    Xt=X[t,,]
    X_new[t,]=sapply(glocN, function(x){sum(x*Xt)})
  }

  #beta_ols = solve(crossprod(X_new))%*%crossprod(X_new, Y)
  beta_ols = lm(Y~X_new-1)$coefficients
  beta_l = lapply(1:length(beta_ols), function(i) beta_ols[i]*glocN[[i]])
  beta_N = Reduce("+", beta_l)

  return(list(beta_ols = beta_N, beta_olsk = beta_ols))
}
