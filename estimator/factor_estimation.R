#robust covariance estimation 
huber.est<-function(x, tau)
{
  ltau = rep(0, length(x))
  ind_small = abs(x)<=tau
  ind_large = !ind_small
  ltau[ind_small] = x[ind_small]
  ltau[ind_large] = tau*sign(x[ind_large])
  return(mean(ltau))
}

cov.robust<-function(Zt, verbose = F)
{
  p = ncol(Zt)
  ### define the outer function on Zt[,i]
  zz<-function(z){
    z2 = outer(z, z, "-")
    zz_vec = z2[upper.tri(z2)]
    return(zz_vec)
  }
  ### caculate the outer product for (j1,j2) column
  ZZt = apply(Zt, 2, zz)
  covz = matrix(0, p, p)
  for (i in 1:p)
  {
    for (j in i:p)
    {
      if (verbose)
        cat(i, j, "\r")
      
      ### run a huber estimator to estimate out the intercept, which is the estimator
      y = ZZt[,i]*ZZt[,j]/2
      tau = get.tau(tau_range = c(1, 1.7), U = y, p = p, Time = nrow(Zt))

      covz[i,j] = huber.est(y, tau)
    }
  }
  covz[lower.tri(covz)] = t(covz)[lower.tri(covz)]
  return(list(covz = covz, tau = tau))
}

### estimate B matrix 
est.B<-function(covZ)
{
  eig_Z = eigen(covZ)
  
  #estimate number of factors
  eigvalues=eig_Z$values
  p=sum(eigvalues>0)
  eigratios=eigvalues[1:(p-1)]/eigvalues[2:p]
  q=which.max(eigratios)
  q=min(q,1)

  ### estimate B & factor
  if (q>1)
  B_hat = as.matrix(eig_Z$vectors[,1:q])%*%diag(sqrt(eig_Z$values[1:q]))
  if (q==1)
  B_hat = as.matrix(eig_Z$vectors[,1])*sqrt(eig_Z$values[1])  
  
  return(list(B_hat=B_hat,q=q))
}

### estimate factor by robust covariance estimation
est.fac.robust<-function(Zt, B_hat, n_tau = 10, verbose = F)
{
  q = ncol(B_hat)
  Time = nrow(Zt)
  n_tau = min(n_tau, Time)
  p = ncol(Zt)
  fac_hat = matrix(0, nrow = q, ncol = Time) 
  Zt = as.matrix(Zt)
  taus = rep(0, n_tau)
  for (t in 1:Time)
  {
    if (verbose)
      cat(t, "\r")
    if (t <= n_tau)
    {
      taus[t] = search.tau.lm.huber(B_hat, as.vector(Zt[t,]), fit.intcp = F, opt=list(method='default'))
    }
    tau = mean(taus[1:min(t, n_tau)])
    hat_huber = lm.huber.fixTau(B_hat, Zt[t,], tau, fit.intcp = F)
    fac_hat[,t] = hat_huber$Coefficients[-1]
  }
  return(list(fac_hat = fac_hat, taus = taus))
}

### estimate the factor by the ols estimator
est.fac.ols<-function(Zt, B_hat)
{
  solve(crossprod(B_hat))%*%crossprod(B_hat, t(Zt))
}