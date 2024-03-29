# DA (monotone) each iteration. 

iter_da <- function(beta, sigma, yobs, X, cw, A, m, I, ni, d, n, p, yfull=FALSE, ...){
  # DA each iteration. (yobs monotone)
  # Args:
  #   beta: Vector. Previous linear regression parameters \tilde{\beta}_j, j\in\{1,...,d\}.
  #   sigma: Matrix. Previous covariance matrix \tilde{\Sigma}_j.
  #   yobs: Monotone matrix. Y observed matrix with Y_obs and the missing values NA.
  #   X: Matrix. Feature matrix X.
  #   cw: Function. Function to generate samples from conditional weight distribution. cw(d, r_i); d: dim of y_i.  
  #   A: Matrix. Positive semidefinite matrix in prior.
  #   m: parameter in prior.
  #   I: !is.na(yobs) missing index matrix. Observed TRUE;  missing FALSE.
  #   ni: Positive integer. colSums(I,na.rm = T). Number of observations in each pattern
  #   d: Positive integer. Number of response components
  #   n: number of observations
  #   p: number of features
  #   yfull: Logical, default FALSE. TRUE: impute all missing components in the response and output them. FALSE: no output (missing componnents in the response).
  #   ... Parameters in the mixing distribution. 
  # Returns:
  #   List. Matrix sigma = sigma(t+1); Matrix beta = beta(t+1)
  #   beta: linear regression parameters
  #   sigma: covariance matrix
  
  #I <- !is.na(yobs) # missing index : observed TRUE  missing FALSE
  #ni <- colSums(I, na.rm = T) # \sum_{i=1}^l n_i
  #if(all(ni == cummax(ni))){
  #stop("yobs is not monotone")
  #} 
  #d <- ncol(I)
  #n <- nrow(I)
  #p <- ncol(X)
  
  # S1 I step 
  ## impute weight
  
  w <- rep(0, n) # w_i

  for (i in 1:n) {
    obs <- I[i, ]
    d_i <- sum(obs)
    ri <- t(t(yobs[i, obs, drop=FALSE]) - t(beta[, obs, drop=FALSE]) %*% t(X[i, , drop=FALSE])) %*% 
          mat_inv(sigma[obs, obs, drop=FALSE]) %*% (t(yobs[i, obs, drop=FALSE]) - t(beta[, obs, drop=FALSE]) %*% t(X[i, , drop=FALSE]))
    w[i] <- cw(d_i, ri, ...)
  }
  
  # p step
  #G0 <- NULL
  #hatbeta0 <- NULL
  H <- matrix(0, nrow=d, ncol=d)
  GZ <- matrix(0, nrow=p, ncol=d) #(G1Z1,...,Gdzd)
  hatbetah <- matrix(0, nrow=p, ncol=d) #(hatbeta1h1,...,hatbetadhd)
  num <- 0
  for (k in 1:d) {
    ## sigma
    index <- c(k:d)
    if(num < ni[k])
    {
      num <- ni[k]
      numobs <- c(1:num)
      Lambda <- diag(c(w[1:num]))
      X1 <- X[numobs, , drop=FALSE] # x(\bm(k),k)
      yobs1 <- yobs[numobs, index, drop=FALSE] # y(\bm(k),k)
      Omega <- mat_inv(t(X1) %*% Lambda %*% X1)
      hatbeta <- Omega %*% t(X1) %*% Lambda %*% yobs1
      G <- t(chol(Omega)) # u
      #Omega0 <- append(G0, list(list(G)))
      #hatbeta0 <- append(hatbeta0, list(list(hatbeta)))
      S <- t(yobs1- X1 %*% hatbeta) %*% Lambda %*% (yobs1 - X1 %*% hatbeta)
      B <- A[index, index, drop=FALSE] + S
      L <- t(chol(mat_inv(B))) # e
    }
    else {
      hatbeta <- hatbeta[, -1, drop=FALSE]
      L <- L[-1, -1, drop=FALSE]
    }
    u <- as.matrix(c(sqrt(rchisq(n=1, df=num - k + (m - p - d + 1))), rnorm(n=d - k))) # F
    h <- L %*% u # h_l = e_lF_l
    H[, k] <- c(rep(0, k - 1), h)
    Z <- as.matrix(mvrnorm(n=1, mu=rep(0, p), Sigma=diag(p)))
    GZ[, k] <- G %*% Z # uz
    hatbetah[, k] <- hatbeta %*% h # (\hatbeta_lg_l)_{l=1}^d
  }
  invh <- mat_inv(H)
  sigma <- t(invh) %*% invh #\Sigma = H^{-1}^T H^{-1}
  # beta
  beta <- (GZ + hatbetah) %*% invh # (u_1Z_1,...,u_dZ_d)h^{-1} + (\hatbeta_1g_1,...,\hatbeta_dg_d)h^{-1}
  if(!yfull) return(list(sigma = sigma, beta = beta))
  else{
    nobs0 <- d # number of obs flag
    for (i in 1:n) {
      obs <- I[i, ]
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
    return(list(sigma = sigma, beta = beta, yimp = yobs[!I]))
  }
}






