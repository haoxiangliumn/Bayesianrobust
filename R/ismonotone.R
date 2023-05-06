# monotone check
ismonotone <- function(yobs)
{
  # Check whether yobs is monotone or not, and permute yobs into the standard monotone format after permuting the rows and columns.
  # Args:
  #   yobs: matrix. Monotone Matrix 
  # Returns:
  #   List. monotone: logic. TRUE if yobs is monotone; FALSE otherwise.
  #         ynew: matrix. yobs in standard monotone format/ clothest to monotone.
  #         ynew_column: vector. ynew (permuted) original column location in yobs.
  #         ynew_row: vector. ynew (permuted) original row location in yobs.

  if(dim(t(yobs))[1] == 1) # vector
  {
    yobs <- as.matrix(yobs) 
  }
  I <- ! is.na(yobs) # missing index : observed TRUE  missing FALSE
  row_obs <- rowSums(I, na.rm = T) 
  row_sort <- order(row_obs, decreasing = TRUE) 
  col_obs <- colSums(I, na.rm = T)
  col_sort <- order(col_obs, decreasing = FALSE) 
  I_new <- I[row_sort, col_sort, drop = FALSE]
  ynew <- yobs[row_sort, col_sort, drop = FALSE]
  ni <- colSums(I_new, na.rm = T)
  monotone <- TRUE
  if (ni[1] == 0){
    error("An observation should contain at least one observed entry")
  }
  # check NA monotone
  if ( all(ni == cummax(ni))) {
   for (i in 1:length(ni)){
     if (sum(I_new[c(1:ni[i]), i]) != ni[i]){
       monotone <- FALSE
       break
     }
   } 
  }
  else{
    monotone <- FALSE
  }
  names(col_sort) = names(yobs)
  return(list(monotone = monotone, ynew = ynew, ynew_column = col_sort, ynew_row = row_sort))
}

