#' Return the positive part of a value
#' @param x the vector
#' @return the maximum of x and 0
#' @export
pos <- function(x) pmax(x, 0)

#' Create folds for time series cross validation. Folds are continguous and sequential.
#' @param n number of rows of a matrix/data_frame to which ids are assigned
#' @param fold_len length of a fold
#' @return a vector of integers assigning values in 1 - nfold to each row.
create_folds <- function(n, fold_len){
  foldid  <- integer(n)
  nfold <- ceiling(n / fold_len)
  for(i in seq.int(0, nfold - 1L)){
    if(i < nfold - 1L){
      foldid[(i*fold_len + 1):((i+1)*fold_len)] = i+1
    }else{
      foldid[(i*fold_len + 1):n] = i+1
    }
  }
  foldid
}

