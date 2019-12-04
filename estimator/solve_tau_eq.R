tau.equation<-function(tau, U, p, Time)
{
  m = floor(Time/2)
  T1 = Time*(Time - 1)/2
  eq_diff = sum(pmin(U^2, tau^2))/tau^2/T1 - 3*log(p)/m
  return(eq_diff)
}

get.tau<-function(tau_range, U, p, Time, verbose = F)
{
  tau_min = tau_range[1]
  tau_max = tau_range[2]
  tau_min_eq = tau.equation(tau_min, U, p, Time)
  while(tau_min_eq<0)
  {
    if (verbose)
      cat("change tau min with eq = ", tau_min_eq, "\n")
    tau_min = tau_min - 0.1
    tau_min_eq = tau.equation(tau_min, U, p, Time)
  }
  tau_max_eq = tau.equation(tau_max, U, p, Time)
  while(tau_max_eq>0)
  {
    if (verbose)
      cat("change tau max with eq = ", tau_max_eq, "\n")
    tau_max = tau_max + 0.1
    tau_max_eq = tau.equation(tau_max, U, p, Time)
  }
  del = tau_min_eq - tau_max_eq
  while(del>10^{-4})
  {
    tau_seq = seq(tau_min, tau_max, length.out = 10)
    eq_value = sapply(tau_seq, tau.equation, U, p = p, Time)
    tau_min_ind = max(which(eq_value>0))
    tau_max_ind = min(which(eq_value<0))
    del = eq_value[tau_min_ind] - eq_value[tau_max_ind]
    tau_min = tau_seq[tau_min_ind]
    tau_max = tau_seq[tau_max_ind]
    if (verbose)
      cat("shrinkage tau_min and tau_max", tau_min, tau_max, del, "\n")
  }
  return(mean(c(tau_min, tau_max)))
}

#get.tau(tau_range = c(1, 1.7), U = y, p = p, Time = nrow(Zt))
