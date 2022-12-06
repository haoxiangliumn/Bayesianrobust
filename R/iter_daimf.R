# DAI (monotone to full) each iteration (impute a monotone pattern to fully observed) (I step (I1 and I2 if contain missing values) and P step).

# 
iter_daimf <- function(beta, sigma, yobs, X, cw, A, m, I, ni, d, n, p, yfull=FALSE){
  # DAI (monotone to full) each iteration (impute a monotone pattern to fully observed) (I step (I1 and I2 if contain missing values) and P step).
  # Args:
  #   beta: Vector. Previous linear regression parameters \tilde{\beta}_j, j\in\{1,...,d\}.
  #   sigma: Matrix. Previous covariance matrix \tilde{\Sigma}_j.
  #   yobs: Matrix. Y observed matrix with Y_obs and the missing values NA.
  #   X: Matrix. Feature matrix X.
  #   cw: Function. Function to generate samples from conditional weight distribution. cw(d, r_i); d: dim of y_i.
  #   A: Matrix. Positive semidefinite matrix in prior.
  #   m: parameter in prior.
  #   I: !is.na(yobs) missing index matrix. Observed TRUE;  missing FALSE.
  #   ni: Positive integer. colSums(I,na.rm = T). Number of observations in each pattern.
  #   d: Positive integer. Number of response components.
  #   n: Positive integer. Number of observations.
  #   p: Positive integer. Number of features.
  #   yfull: Logical, default FALSE. TRUE: impute all missing components in the response and output them. FALSE: no output (missing componnents in the response).
  # 
  # Returns:
  #   List. Matrix sigma = sigma(t+1); Matrix beta = beta(t+1)
  #   beta: linear regression parameters
  #   sigma: covariance matrix
  
  #I <- ! is.na(yobs)
  #ni <- colSums(I, na.rm=T)
  #if(all(ni == cummax(ni))){
  #  stop("yobs is not monotone")
  #} 
  #d <- ncol(I)
  #n <- nrow(I)
  #p <- ncol(X)
  # S1 I step
  # impute weight
  w <- rep(0, n) # w_i
  nobs0 <- d # number of obs flag
  for (i in 1:n) {
    obs <- I[i, ]
    d_i <- sum(obs)
    ri <- t(t(yobs[i, obs, drop=FALSE]) - t(beta[, obs, drop=FALSE]) %*% t(X[i, , drop=FALSE])) %*% 
          mat_inv(sigma[obs, obs, drop=FALSE]) %*% (t(yobs[i, obs, drop=FALSE]) - t(beta[, obs, drop=FALSE]) %*% t(X[i, , drop=FALSE]))
    w[i] <- cw(d_i, ri)
    ## impute missing value
    nobs <- length(obs[obs == TRUE])
    if (nobs == d){
      next # no missing no need to impute
    }
    mu1 <- t(beta[, ! obs, drop=FALSE]) %*% t(X[i, , drop=FALSE])
    diff2 <- t(yobs[i, obs, drop=FALSE]) - t(beta[, obs, drop=FALSE]) %*% t(X[i, , drop=FALSE]) # y_obs-beta_obs^T x  vertical vector
    if (nobs < nobs0){
    # no need to update if in the same pattern
    sigma11 <- sigma[ ! obs, ! obs, drop=FALSE]
    sigma21 <- sigma[obs, ! obs, drop=FALSE]
    sigma12 <- t(sigma21)
    sigma22inv <- mat_inv(sigma[obs, obs, drop=FALSE])
    nobs0 <- nobs
    }
    yobs[i, ! obs] <- mvrnorm(n=1, mu=mu1 + sigma12 %*% sigma22inv %*% diff2, Sigma=(sigma11 - sigma12 %*% sigma22inv %*% sigma21) / w[i])
  }

  # p step
  ## sigma
  Lambda <- diag(w)
  Omega <- mat_inv(t(X) %*% Lambda %*% X)
  hatbeta <- Omega %*% t(X) %*% Lambda %*% yobs
  S <- t(yobs - X %*% hatbeta) %*% Lambda %*% (yobs - X %*% hatbeta)
  sigma <- rinvwishart(nu=n - p + m - d, S=S + A)
  ## beta
  beta <- t(rmatnorm(s=1, M=t(hatbeta), U=sigma, V=Omega))
  if(yfull) return(list(sigma = sigma, beta = beta, yimp = yobs[!I]))
  else return(list(sigma = sigma, beta = beta))
}






