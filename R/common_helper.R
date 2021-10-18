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

#' Moving average of value for specific day of week over last x weeks
#' @param v a vector of values (assume they are ordered chronologically and contiguous)
#' @param x integer number of weeks
#' @return a vector average of the values for the last x weeks, scaled appropriately to avoid NAs
#' @export
ma_week <- function(v, x = 6L) {
  if (length(v) < 1L) stop("v must have at least length 1.")
  r <- v
  week_num <- sort(rep(1:ceiling(length(v) / 7), 7))[1:length(v)]
  number_of_lags <- week_num
  ind = which(week_num >= x)
  number_of_lags[ind] <- x

  for (i in seq(x - 1))
    r[week_num > i] <- r[week_num > i] + dplyr::lag(v, 7L * i)[week_num > i]
  r / number_of_lags
}

