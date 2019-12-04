###############################################################################
# return mean huber loss of a vector 
loss.huber = function(x,tau) {
  x = abs(x)
  tmp = sum(x[x<tau]^2)+sum(2*tau*x[x>=tau])-tau^2*length(x[x>=tau])
  return(tmp/length(x))
}
###############################################################################
