# sample weight functions from the conditional distributions 
# (mixing distribution parameters included)
cw_gig0 <- function(di, ri, a=1, b=1, p= - 0.5) {
  return(rgig(n=1, chi=b, psi=a + ri, lambda=di / 2 + p))
}

cw_gamma0 <- function(di, ri, a=2, b=2) {
  return(rgamma(n=1, shape=di / 2 + a, rate=ri / 2 + b))
}
cw_mass0 <- function(di, ri, gamma=0.6, lambda=0.8)
{ 
  # contaminated normal
  t1 <- lambda ^ {di / 2} * exp( - ri * lambda / 2) * gamma
  t2 <- exp( - ri / 2) * (1 - gamma)
  if (runif(1) < (t1 / (t1 + t2))) 
    return(lambda)
  else 
    return(1)
}


##
# sample weight functions from the conditional distributions
cw_gig <- function(di, ri,...) {
  # Sample weight using GIG(a=1, b=1, p=-0.5 default) mixing distribution (I step in DA; I1 step in DAI)
  #
  # Args:
  #   di: Integer. The number of components observed in the response Y_i
  #   ri: Non-negative real number. r_{i,\bm{k}} in the paper
  #   ...: Mixing distribution parameters. a,b, and p in GIG(a=1, b=1, p=-0.5 default)
  #
  # Returns:
  #   Positive real number. Sampled weight from the conditional distribution.  
  return(cw_gig0(di=di, ri=ri, ...))
}
cw_gamma <- function(di, ri, ...) {
  # Sample weight using gamma(a=2, b=2) mixing distribution 
  #   ...: Mixing distribution parameters. a and b in gamma(a=2, b=2)
  
  return(cw_gamma0(di=di, ri=ri, ...))
}
cw_constant <- function(di, ri){ 
  # Constant weight w=1
  return(1)
}
cw_mass <- function(di, ri, ...) {
  # Sample weight using mixing distribution w = lambda with prob gamma; w = 1 with prob 1-gamma 
  #   ...: Mixing distribution parameters. w and lambda in w = lambda with prob gamma; w = 1 with prob 1-gamma 
  return(cw_mass0(di, ri, ...))
}

