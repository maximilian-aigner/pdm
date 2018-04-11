hsvmpathc <-
function(x, y, KK ,  nlam, flmin, ulam, isd, 
                      eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, delta, nobs, nvars, 
                      vnames ) {
  #################################################################################
  #data setup
  y <- as.factor(y)
  y <- c(-1, 1)[as.numeric(y)]
  if (!all(y %in% c(-1, 1))) 
    stop("y should be a factor with two levels")
  if (delta < 0) 
    stop("delta must be non-negative")
  delta <- as.double(delta)
 # KK<-  as.integer(KK)

  #################################################################################
  # call Fortran core
  fit1 <- .Fortran("hsvmlassoNETc", delta, lam2, nobs, nvars, 
                   as.double(x), as.double(y) ,KK = as.integer(KK),  jd, pf, pf2, dfmax, pmax, nlam, 
                   flmin, ulam, eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
                   beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                   alam = double(nlam), npass = integer(1), jerr = integer(1) ,cluster=double(nvars*nlam))
  #################################################################################
  # output
  
  outlist <- getoutput(fit1, maxit, pmax, nvars, vnames)
  outlist <- c(outlist, list(npasses = fit1$npass, jerr = fit1$jerr))
  class(outlist) <- c("hsvmpathc")
  outlist
  
}
