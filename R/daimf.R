# DAI (impute a monotone pattern y_{(\bm{k})} to fully observed if y has missing values)
daimf <- function(X, yobs, m=ncol(yobs), A=diag(0, ncol(yobs)), cw=cw_gamma, iter=100, yfull=FALSE, ...)
{
  # DAI (impute a monotone pattern y_{(\bm{k})} to fully observed if y has missing values)
  # Args:
  #   X: matrix. Feature matrix.
  #   yobs: Monotone Matrix (need to reshape into monotone structure). Y observed matrix with Y_obs and the missing values NA.
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
  
  
  if(dim(t(yobs))[1] == 1) # vector
  {
    yobs <- as.matrix(yobs) 
  }
  I <- ! is.na(yobs) # missing index : observed TRUE  missing FALSE
  ni <- colSums(I, na.rm = T) # \sum_{i=1}^l n_i
  if (ni[1] == 0){
    stop("An observation should contain at least one observed entry")
  }
  # check NA monotone
  if ( all(ni == cummax(ni))) {
    for (i in 1:length(ni)){
      if (sum(I[c(1:ni[i]), i]) != ni[i]){
        stop("yobs is not monotone")
      }
    } 
  }
  
  d <- ncol(I)
  n <- nrow(I)
  p <- ncol(X)
  
  # initial value OLS 
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
    new <- iter_daimf(beta=beta, sigma=sigma, yobs=yobs, X=X, cw=cw, A=A, m=m, I=I, ni=ni, d=d, n=n, p=p, yfull=yfull, ...)
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