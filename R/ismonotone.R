# monotone check
ismonotone <- function(yobs)
{
  # Check whether yobs is monotone or not, and permute yobs into the standard monotone format after permuting the rows and columns.
  # Args:
  #   yobs: Monotone Matrix 
  # Returns:
  #   List. monotone: TRUE if yobs is monotone; FALSE otherwise.
  #         ynew: yobs in standard monotone format/ clothest to monotone.
  #         ynew_column: permuted yobs original column location in yobs.
  #         ynew_row: permuted yobs original row location in yobs.
  I <- ! is.na(yobs) # missing index : observed TRUE  missing FALSE
  row_obs <- rowSums(I, na.rm = T) 
  row_sort <- order(row_obs, decreasing = TRUE) 
  col_obs <- colSums(I, na.rm = T)
  col_sort <- order(col_obs, decreasing = FALSE) 
  I_new <- I[row_sort, col_sort]
  ynew <- yobs[row_sort, col_sort]
  ni <- colSums(I_new, na.rm = T)
  monotone <- TRUE
  if ( ! all(ni == cummax(ni))) {
    monotone <- FALSE
  }
  names(col_sort) = names(yobs)
  return(list(monotone = monotone, ynew = ynew, col_sort = col_sort, row_sort = row_sort))
}

