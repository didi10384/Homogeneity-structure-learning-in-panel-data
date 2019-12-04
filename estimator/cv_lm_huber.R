############ cross-validation given tau for huber regression ################
cv.lm.huber = function(X,Y,tau,fit.intcp=TRUE,loss.validation='mse', 
                       seed=42, nfold=5, nrepeat=1) {
  set.seed(seed)
  
  # set loss function for validation
  if (loss.validation == 'mse') # mean square error
    loss.func = function(x) { return( mean(x^2) ) }
  else if (loss.validation == 'mae') # mean absolute deviation
    loss.func = function(x) { return( mean(abs(x)) ) }
  else if (loss.validation == 'huber') # mean huber loss 
    loss.func = function(x) {return(loss.huber(x,tau))}
  
  loss.set = c()
  for (rep in 1:nrepeat) {
    folds = cut(sample(1:nrow(X)),breaks=nfold,labels=FALSE, 
                ordered_result = FALSE)
    
    for (i in 1:nfold) {
      index.test = folds==i
      mdl = lm.huber.fixTau(X[!index.test,,drop = F], Y[!index.test], 
                            fit.intcp=fit.intcp, tau=tau)
      resids = Y[index.test]-mdl$predict(X[index.test,,drop = F])
      loss.set = append(loss.set,loss.func(resids) )
    }
  }
  return( mean(loss.set) )
}
#############################################################################