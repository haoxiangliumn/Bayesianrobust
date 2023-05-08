# DA(I) algorithm
dageneral <- function(X, yobs, m=ncol(yobs), A=diag(0, ncol(yobs)), cw=cw_gamma, Ik, iter=100, yfull=FALSE)
{
  # General function. DA(I) sampler.
  # Args:
  #   X: matrix. Feature matrix.
  #   yobs: Monotone Matrix (need to reshape into monotone structure). Y observed matrix with Y_obs and the missing values NA.
  #   m: Real number. Prior parameter.
  #   A: Matrix. Positive semidefinite matrix in prior.
  #   cw: Function. Function to generate samples from conditional weight distribution. cw(d, r_i); d: dim of y_i.
  #   Ik: Matrix. All elements should be logical values. Imputation matrix structure (k'). Ik should represent a monotone missing pattern. If Ik is missing, for monotone yobs, the algorithm will implement DA(monotone), and for non-monotone yobs, the algorithm will set Ik to TRUE for all the elements. Ik_{ij} = TRUE: y_{ij} has value; Ik_{ij} = FALSE: y_{ij} NA.
  #   iter: Positive integer. Number of iterations
  #   yfull: Logical, default FALSE. TRUE: impute all missing components in the response and output them. FALSE: no output (missing componnents in the response).
  #
  # Returns:
  #   List. First iter+1 elements: Matrix sigma and Matrix beta t=0 to t=iter; iter+2 element: time to implement the algorithm.
  #         beta: linear regression parameters
  #         sigma: covariance matrix
  
  # A: positive semidefinite
  if( ! isSymmetric(A)){
    stop("A should be Positive semidefinite")
  }
  if(min(eigen(A, symmetric=TRUE)$values) < 0){
    stop("A should be Positive semidefinite")
  } 
  
  yobs <- as.matrix(yobs) 
  X <- as.matrix(X) 
  
  if( ! missing(Ik)) {
    if( ! is.logical(Ik)){
      stop("Ik elements should be TRUE/FALSE logical values")
    }
  }
  I <- ! is.na(yobs) # missing index : observed TRUE  missing FALSE
  ni <- colSums(I, na.rm=T) # \sum_{i=1}^l n_i
  monotone <- TRUE
  if (ni[1] == 0){
    stop("An observation should contain at least one observed entry")
  }
  # check NA monotone
  if ( all(ni == cummax(ni))) {
    for (i in 1:length(ni)){
      if (sum(I[c(1:ni[i]), i]) != ni[i]){
        monotone <- FALSE
        break
      }
    } 
  }
  else{
    monotone <- FALSE
  }
  if(!missing(Ik)){
    nik <- colSums(I, na.rm=T)
    # check NA monotone
    if ( all(nik == cummax(nik))) {
      for (i in 1:length(nik)){
        if (sum(I[c(1:nik[i]), i]) != nik[i]){
          stop("Ik should be monotone")
        }
      } 
    }
    else{
      stop("Ik should be monotone")
    }
  }
  
  
  d <- ncol(I)
  n <- nrow(I)
  p <- ncol(X)
  for (i in 1:d) {
    if (sum(I[c(1:ni[i]), i]) != ni[i]) {
      monotone <- 0
      break
    } 
  }
  if (monotone) { # monotone
    if(missing(Ik)){
      result <- da(X, yobs, m=m, A=A, cw=cw, iter=iter, yfull=yfull) # DA
      print("da")
      return(result)
    }
    else if(sum(Ik) == nrow(yobs) * ncol(yobs)) {
      result <- daimf(X, yobs, m=m, A=A, cw=cw, iter=iter, yfull=yfull) # monotone and impute to fully observed
      print("daimf")
      return(result)
      }
    else {
      result <- dai(X, yobs, m, A, cw=cw, Ik=Ik, iter=iter, yfull=yfull) # dai general
      print("dai")
      return(result)
    } 
  } 
  else {
    if(missing(Ik)){
      result <- daif(X, yobs, m=m, A=A, cw=cw, iter=iter, yfull=yfull) # impute to full
      print("daif")
      return(result)
    }
    else if(sum(Ik) == nrow(yobs) * ncol(yobs)){
      result <- daif(X, yobs, m=m, A=A, cw=cw, iter=iter, yfull=yfull) # impute to full
      print("daif")
      return(result)
    }
    else{
      result <- dai(X, yobs, m, A, cw=cw, Ik=Ik, iter=iter, yfull=yfull) # impute to Ik \bm{k'}
      print("dai")
      return(result)
    }
  }
} 