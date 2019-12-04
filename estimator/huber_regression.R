################ huber regression given fixed tau #############################
# fit linear model given tau value
# data: X=n*p, y=n*1 
# min: sum_i huber_tau( Y_i - X_i*beta - miu )
lm.huber.fixTau = function(X,Y,tau=NA,fit.intcp=TRUE,tol=1e-4,
                           iter.max=100,batch_size=1000, verbose=FALSE) {
  batch_size = min(batch_size,nrow(X))
  
  # Initialize: use OLS weights as initializer
  if(fit.intcp) {
    mdl = lm(Y~X)
    intcp = mdl$coefficients[1] # intercept
    wX = mdl$coefficients[2:length(mdl$coefficients)] # vector, len=p
  } else {
    mdl = lm(Y~X-1)
    intcp = 0 # intercept
    wX = mdl$coefficients # vector, len=p
  }
  resids = mdl$residuals
  if( is.na(tau) )
    tau = median(abs(resids))
  rm(mdl)
  
  iter = 0
  # use gradient descent on huber loss until converge or exceed max.iter
  repeat {
    # check if exceed max iteration
    iter = iter+1
    if( iter > iter.max ) {
      if( verbose )
        print(c('stopped: exceed max iteration, current gradient size:', 
                grad.size))
      # print(c('stopped: exceed max iteration, current gradient size:',
      #         grad.size))
      break
    }
    
    # specify data for gradient and calculate gradient
    # spread gradient among coordinates and update weights & residuals
    idx = sample(1:nrow(X), batch_size)
    resids = Y[idx] - X[idx,, drop = F]%*%wX - intcp
    grad.resids = 2*pmin(tau, pmax(-tau,resids)) # 1*batch_size
    if(fit.intcp) grad.intcp = sum(grad.resids)/batch_size 
    else grad.intcp = 0
    grad.wX = as.vector(grad.resids%*%X[idx,])/batch_size
    grad.size = sum(grad.intcp^2+grad.wX^2)
    
    # in the first iteration, adjust the scale of 
    # gradient to the right order
    if (iter == 1) {
      step.scale = 1024
      base = loss.huber(resids,tau) 
      # should decrease huber loss, or at least not increase too much
      step = X[idx,, drop = F]%*%grad.wX + grad.intcp
      loss = base + 1
      while( loss > base ) {
        step.scale = step.scale/2
        loss = loss.huber(resids-step.scale*step,tau)
      }
      if(verbose)
        print(c('loss',loss,'step.scale',step.scale))
      rm(base,step,loss)
    }
    
    # check if converged
    if( grad.size < tol )
      break
    
    # else update wX, intcp
    wX = wX + grad.wX*step.scale/sqrt(iter)
    intcp = intcp + grad.intcp*step.scale/sqrt(iter)
    
    if(verbose)
      print(c('Iter: ',iter,'grad.size',grad.size, 
              'huber_loss', loss.huber(resids,tau)))
  }
  
  # return parameters, intercept, loss, n_iter, iter_max
  pred.func = function(X) {return(as.vector(X%*%wX + intcp))}
  result = list(Coefficients = c(intcp, wX), 
                tau = tau,
                loss=loss.huber(resids,tau), 
                nIter=iter, 
                maxIter=iter.max, 
                tolerance=tol, 
                predict=pred.func)
  if(verbose) 
    print(result)
  return(result)
}
###############################################################################