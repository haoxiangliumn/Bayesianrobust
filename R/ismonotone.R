# monotone check
ismonotone <- function(yobs, x)
{
  # Check whether yobs is monotone or not, and permute yobs into the standard monotone format after permuting the rows and columns.
  # Args:
  #   yobs: matrix. Monotone matrix.
  #   x: matrix. Feature matrix.
  # Returns:
  #   List. monotone: logic. TRUE if yobs is monotone; FALSE otherwise.
  #         ynew: matrix. yobs in standard monotone format/ closest to monotone.
  #         ynew_column: vector. ynew (permuted) original column location in yobs.
  #         ynew_row: vector. ynew (permuted) original row location in yobs.
  #         xnew: matrix. the corresponding feature matrix after rearranging the rows. 
  #         permute: logic. TRUE if permutation is needed (to be monotone/ closest to monotone); FALSE if no need to permute.

  yobs <- as.matrix(yobs) 
  x <- as.matrix(x)
  I <- ! is.na(yobs) # missing index : observed TRUE  missing FALSE
  row_obs <- rowSums(I, na.rm = T) 
  row_sort <- order(row_obs, decreasing = TRUE) 
  col_obs <- colSums(I, na.rm = T)
  col_sort <- order(col_obs, decreasing = FALSE) 
  I_new <- I[row_sort, col_sort, drop = FALSE]
  ynew <- yobs[row_sort, col_sort, drop = FALSE]
  xnew <- x[row_sort, , drop = FALSE]
  ni <- colSums(I_new, na.rm = T)
  monotone <- TRUE
  permute = !all(row_sort == c(1:length(row_sort))) * all(col_sort == c(1:length(col_sort))) 
  if (ni[1] == 0){
    stop("An observation should contain at least one observed entry")
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
  ynew = as.matrix(ynew)
  return(list(monotone = monotone, ynew = ynew, ynew_column = col_sort, ynew_row = row_sort, xnew = xnew, permute = permute))
}

