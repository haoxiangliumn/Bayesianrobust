# monotone data augmentation gibbs DA each iteration
da <- function(X, yobs, m=ncol(yobs), A=diag(0, ncol(yobs)), cw=cw_gamma, iter=100, yfull=FALSE)
{
  # DA each iteration (yobs has to be monotone).
  # Args:
  #   X: matrix. Feature matrix.
  #   yobs: Monotone Matrix (need to reshape into monotone structure). Y observed matrix with Y_obs and the missing values NA.
  #   m: Real number. Prior parameter.
  #   A: Matrix. Positive semidefinite matrix in prior.
  #   cw: Function. Function to generate samples from conditional weight distribution. cw(d, r_i); d: dim of y_i.
  #   iter: Positive integer. Number of iterations
  #   yfull: Logical, default FALSE. TRUE: impute all missing components in the response and output them. FALSE: no output (missing componnents in the response).
  # 
  # Returns:
  #   List. First iter+1 elements: Matrix sigma and Matrix beta t=0 to t=iter; iter+2 element: time to implement the algorithm.
  #         beta: linear regression parameters
  #         sigma: covariance matrix  
  #
  # TODO: add I1/I2 y imputation to monotone step for non-monotone responses.
  
  I <- ! is.na(yobs) # missing index : observed TRUE  missing FALSE
  ni <- colSums(I, na.rm = T) # \sum_{i=1}^l n_i
  if ( ! all(ni == cummax(ni))) {
    stop("yobs is not monotone")
  } 
  d <- ncol(I)
  n <- nrow(I)
  p <- ncol(X)
  for (i in 1:d) {
    if (sum(I[c(1:ni[i]), i]) != ni[i]) {
      stop("yobs is not monotone")
      break
    } 
  }
  
  if(yfull) samples <- matrix(0, nrow=iter + 1, ncol=p * d + d * d + sum(!I))
  else samples <- matrix(0, nrow=iter + 1, ncol=p * d + d * d)
  # initial (use OLS from the pattern 1 data)
  first <- c(1:ni[1])
  beta <- mat_inv(t(X[first, , drop=FALSE]) %*% X[first, , drop=FALSE]) %*% t(X[first, , drop=FALSE]) %*% yobs[first, , drop=FALSE]
  sigma <- t(yobs[first, , drop=FALSE] - X[first, , drop=FALSE] %*% beta) %*% (yobs[first, , drop=FALSE] - X[first, , drop=FALSE] %*% beta)
  if(!yfull) samples[1, ] <- c(as.vector(beta), as.vector(sigma))
  else samples[1, ] <- c(as.vector(beta), as.vector(sigma), matrix(0,nrow=1, ncol=sum(!I)))
  start.time <- Sys.time()
  for (i in 1:iter) {
    new <- iter_da(beta=beta, sigma=sigma, yobs=yobs, X=X, cw=cw, A=A, m=m, I=I, ni=ni, d=d, n=n, p=p, yfull=yfull)
    beta <- new$beta
    sigma <- new$sigma
    if(!yfull) samples[i + 1, ] <- c(as.vector(beta), as.vector(sigma))
    else {
      yimp <- new$yimp
      samples[i + 1, ] <- c(as.vector(beta), as.vector(sigma), as.vector(yimp))
    }
  }
  end.time <- Sys.time()
  result <- list(samples=samples, time=end.time-start.time)
  return(result)
}

#result_mda = mda_gibbs(X, yobs, m, A, cw=cw_gamma, iter=iter)
