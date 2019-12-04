################ search tau for best huber regression #######################
# search for best tau producing best huber mean for a list of input
# input: x (vector), search method, options
search.tau.lm.huber = function(X,Y,fit.intcp=TRUE,
                               opt=list(method='default')){
  if (opt$method == 'default') {
    min_tau = mad(lm(Y~X)$residuals)/2
    tau.set = seq(log(min_tau),log(max(Y)-min(Y)),by=0.1)
    tau.set = exp(tau.set)
    return( search.tau.lm.huber(X,Y,fit.intcp=fit.intcp, 
                                opt=list(method='set',tau_set=tau.set)))
  }
  
  if(opt$method == 'set') {
    tau.set = opt$tau_set
    scores = vector(mode='numeric',length=length(tau.set)) 
    
    for (i in 1:length(tau.set))
      scores[i] = cv.lm.huber(X,Y,tau.set[i],fit.intcp=fit.intcp, 
                              nfold=5, nrepeat=1)
    
    return(max(tau.set[min(scores)==scores]))
  }
}
#############################################################################


