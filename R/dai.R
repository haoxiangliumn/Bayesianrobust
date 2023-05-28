# Gibbs sampler implementing DAI algorithms (any missing structure to fully observed/monotone)
# DAI (DAI) (impute y_{(I)} to y_{(Ik)} if y has missing values)
dai <- function(X, yobs, m=ncol(yobs), A=diag(0, ncol(yobs)), cw=cw_gamma, Ik=matrix(TRUE, nrow=nrow(yobs), ncol=ncol(yobs)), iter=100, yfull=FALSE, ...)
{
  # DAI (impute y_{(I)} to y_{(Ik)}) (any missing structure to fully observed/monotone)
  # Args:
  #   X: matrix. Feature matrix.
  #   yobs: Matrix. Y observed matrix with Y_obs and the missing values NA.
  #   m: Real number. Prior parameter.
  #   A: Matrix. Positive semidefinite matrix in prior.
  #   cw: Function. Function to generate samples from conditional weight distribution. cw(d, r_i); d: dim of y_i.
  #   Ik: Imputation matrix structure (\bm{k}'). Ik_{ij} = TRUE: y_{ij} has value; Ik_{ij} = FALSE: y_{ij} no value.
  #   iter: Positive integer. Number of iterations.
  #   yfull: Logical, default FALSE. TRUE: impute all missing components in the response and output them. FALSE: no output (missing componnents in the response).
  #   ... Parameters in the mixing distribution. 
  # Returns:
  #   List. First iter+1 elements: Matrix sigma and Matrix beta t=0 to t=iter; iter+2 element: time to implement the algorithm.
  #         beta: linear regression parameters
  #         sigma: covariance matrix

  I <- ! is.na(yobs) # missing index : observed TRUE  missing FALSE
  ni <- colSums(Ik, na.rm = T) # \sum_{i=1}^l n_i
  if ( ! all(ni == cummax(ni))) {
    stop("Ik does not represent a monotone pattern")
  }
  d <- ncol(I)
  n <- nrow(I)
  p <- ncol(X)
  for (i in 1:d) {
    if (sum(Ik[c(1:ni[i]), i]) != ni[i]) {
      stop("Ik does not represent a monotone pattern")
    }
    for (j in 1:n){
      if(I[j, i]){
        if(!Ik[j, i]){
          stop("Ik should be a larger missing pattern than yobs with missing pattern I")
        }
      }
    }
  }
  
  if(!is.logical(Ik)) {
    stop("All elements in Ik should be logical: TRUE or FALSE")
  }
  if(-1 %in% (Ik - I)) {
    stop("Ik should be larger than I")
  }

  if(yfull) samples <- matrix(0, nrow=iter + 1, ncol=p * d + d * d + sum(!I))
  else samples <- matrix(0, nrow=iter + 1, ncol=p * d + d * d)
  # initial values: 0
  beta <- matrix(0, nrow=p, ncol=d)
  sigma <- diag(d)
  if(!yfull) samples[1, ] <- c(as.vector(beta), as.vector(sigma))
  else samples[1, ] <- c(as.vector(beta), as.vector(sigma), matrix(0,nrow=1, ncol=sum(!I)))
  
  start.time <- Sys.time()
  for (i in 1:iter) {
    new <- iter_dai(beta=beta, sigma=sigma, yobs=yobs, X=X, cw=cw, A=A, m=m, I=I, Ik=Ik, ni=ni, d=d, n=n, p=p, yfull=yfull, ...)
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

