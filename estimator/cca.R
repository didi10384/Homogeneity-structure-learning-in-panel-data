cc <- function (X, Y) 
{
  Xnames = dimnames(X)[[2]]
  Ynames = dimnames(Y)[[2]]
  ind.names = dimnames(X)[[1]]
  res = rcc(X, Y, 0, 0)
  return(res)
}

rcc <- function (X, Y, lambda1, lambda2) 
{
  Xnames <- dimnames(X)[[2]]
  Ynames <- dimnames(Y)[[2]]
  ind.names <- dimnames(X)[[1]]
  Cxx <- var(X, na.rm = TRUE, use = "pairwise") + diag(lambda1, 
                                                       ncol(X))
  Cyy <- var(Y, na.rm = TRUE, use = "pairwise") + diag(lambda2, 
                                                       ncol(Y))
  Cxy <- cov(X, Y, use = "pairwise")
  res <- geigen(Cxy, Cxx, Cyy)
  names(res) <- c("cor", "xcoef", "ycoef")
  scores <- comput(X, Y, res)
  return(list(cor = res$cor, names = list(Xnames = Xnames, 
                                          Ynames = Ynames, ind.names = ind.names), xcoef = res$xcoef, 
              ycoef = res$ycoef, scores = scores))
}

geigen <- function (Amat, Bmat, Cmat) 
{
  Bdim <- dim(Bmat)
  Cdim <- dim(Cmat)
  if (Bdim[1] != Bdim[2]) 
    stop("BMAT is not square")
  if (Cdim[1] != Cdim[2]) 
    stop("CMAT is not square")
  p <- Bdim[1]
  q <- Cdim[1]
  s <- min(c(p, q))
  if (max(abs(Bmat - t(Bmat)))/max(abs(Bmat)) > 1e-10) 
    stop("BMAT not symmetric.")
  if (max(abs(Cmat - t(Cmat)))/max(abs(Cmat)) > 1e-10) 
    stop("CMAT not symmetric.")
  Bmat <- (Bmat + t(Bmat))/2
  Cmat <- (Cmat + t(Cmat))/2
  Bfac <- chol(Bmat)
  Cfac <- chol(Cmat)
  Bfacinv <- solve(Bfac)
  Cfacinv <- solve(Cfac)
  Dmat <- t(Bfacinv) %*% Amat %*% Cfacinv
  if (p >= q) {
    result <- svd2(Dmat)
    values <- result$d
    Lmat <- Bfacinv %*% result$u
    Mmat <- Cfacinv %*% result$v
  }
  else {
    result <- svd2(t(Dmat))
    values <- result$d
    Lmat <- Bfacinv %*% result$v
    Mmat <- Cfacinv %*% result$u
  }
  geigenlist <- list(values, Lmat, Mmat)
  names(geigenlist) <- c("values", "Lmat", "Mmat")
  return(geigenlist)
}

comput <- function (X, Y, res) 
{
  X.aux = scale(X, center = TRUE, scale = FALSE)
  Y.aux = scale(Y, center = TRUE, scale = FALSE)
  X.aux[is.na(X.aux)] = 0
  Y.aux[is.na(Y.aux)] = 0
  xscores = X.aux %*% res$xcoef
  yscores = Y.aux %*% res$ycoef
  corr.X.xscores = cor(X, xscores, use = "pairwise")
  corr.Y.xscores = cor(Y, xscores, use = "pairwise")
  corr.X.yscores = cor(X, yscores, use = "pairwise")
  corr.Y.yscores = cor(Y, yscores, use = "pairwise")
  return(list(xscores = xscores, yscores = yscores, corr.X.xscores = corr.X.xscores, 
              corr.Y.xscores = corr.Y.xscores, corr.X.yscores = corr.X.yscores, 
              corr.Y.yscores = corr.Y.yscores))
}


svd2 <- function (x, nu = min(n, p), nv = min(n, p), LINPACK = FALSE) 
{
  dx <- dim(x)
  n <- dx[1]
  p <- dx[2]
  svd.x <- try(svd(x, nu, nv, LINPACK))
  if (class(svd.x) == "try-error") {
    nNA <- sum(is.na(x))
    nInf <- sum(abs(x) == Inf)
    if ((nNA > 0) || (nInf > 0)) {
      msg <- paste("sum(is.na(x)) = ", nNA, "; sum(abs(x)==Inf) = ", 
                   nInf, ".  'x stored in .svd.x.NA.Inf'", sep = "")
      stop(msg)
    }
    tf <- tempfile("svd.LINPACK.error.matrix", tmpdir = getwd(), 
                   fileext = ".rda")
    save(x, nu, nv, LINPACK, svd.x, file = tf)
    msg <- paste("svd failed using LINPACK = ", LINPACK, 
                 " with n = ", n, " and p = ", p, ";  x stored in '", 
                 tf, "'", sep = "")
    warning(msg)
    svd.x <- try(svd(x, nu, nv, !LINPACK))
    if (class(svd.x) == "try-error") {
      stop("svd also failed using LINPACK = ", !LINPACK, 
           ";  x stored in '", tf, "'")
    }
  }
  svd.x
}
