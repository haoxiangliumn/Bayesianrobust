# DAI (impute y_{(\bm{k})} any missing structure to fully observed if y has missing values)
daif <- function(X, yobs, m=ncol(yobs), A=diag(0, ncol(yobs)), cw=cw_gamma, iter=100, yfull=FALSE, ...)
{
  # DAI (impute y_{(\bm{k})} any missing structure to fully observed if y has missing values)
  # Args:
  #   X: matrix. Feature matrix.
  #   yobs: Matrix (with missing as NA). Y observed matrix with Y_obs and the missing values NA.
  #   m: Real number. Prior parameter.
  #   A: Matrix. Positive semidefinite matrix in prior.
  #   cw: Function. Function to generate samples from conditional weight distribution. cw(d, r_i); d: dim of y_i.
  #   iter: Positive integer. Number of iterations
  #   yfull: Logical, default FALSE. TRUE: impute all missing components in the response and output them. FALSE: no output (missing componnents in the response).
  #   ... Parameters in the mixing distribution. 
  # Returns:
  #   List. First iter+1 elements: Matrix sigma and Matrix beta t=0 to t=iter; iter+2 element: time to implement the algorithm.
  #         beta: linear regression parameters
  #         sigma: covariance matrix
  
  
  I <- ! is.na(yobs) # missing index : observed TRUE  missing FALSE
  d <- ncol(I)
  n <- nrow(I)
  p <- ncol(X)
  
  if(yfull) samples <- matrix(0, nrow=iter + 1, ncol=p * d + d * d + sum(!I))
  else samples <- matrix(0, nrow=iter + 1, ncol=p * d + d * d)
  # initial values: 0
  beta <- matrix(0, nrow=p, ncol=d)
  sigma <- diag(d)
  if(!yfull) samples[1, ] <- c(as.vector(beta), as.vector(sigma))
  else samples[1, ] <- c(as.vector(beta), as.vector(sigma), matrix(0,nrow=1, ncol=sum(!I)))
  
  start.time <- Sys.time()
  for (i in 1:iter) {
    new <- iter_daif(beta=beta, sigma=sigma, yobs=yobs, X=X, cw=cw, A=A, m=m, I=I, d=d, n=n, p=p, yfull=yfull, ...)
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
