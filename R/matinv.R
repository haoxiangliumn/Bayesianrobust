# matrix inverse
mat_inv <- function(X) {
  # Calculate matrix (n*n or 1*1) inverse.
  # 
  # Args: 
  #   X: matrix (n*n or 1*1).
  # 
  # Returns:
  #   Inverse of the matrix X.
  if (all(dim(X) == c(1,1))) {
    return(1 / X)
  }
  else { 
    return(inv(X))
  }
}
