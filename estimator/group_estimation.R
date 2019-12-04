# ### estimate the groups
beta.group<-function(beta_est,G)
{
  N = nrow(beta_est)
  vb = as.vector(beta_est)
  ord = order(vb)
  vec = sort(vb)

  ### estimate groups
  w = wbs(vec)
  #w_cpt = changepoints(w, penalty="ssic.penalty", alpha = 3,Kmax=G)
  w_cpt = changepoints(w, Kmax=G)

  ind = sort(w_cpt$cpt.ic$ssic.penalty)
  if (length(vec) - max(ind)< 10)
  {
    ind[length(ind)] = length(vec)
  }
  po = unique(c(ind, length(vec)))
  ng = length(po)
  gg = rep(1:ng, c(po[1], diff(po)))
  gvb = rep(0, length(vb))
  gvb[ord] = gg
  g_mat = matrix(gvb, nrow = N)
  glocN = lapply(1:ng, function(i) g_mat==i)
  return(list(glocN = glocN, po = po, ng = ng, g_mat = g_mat))
}
