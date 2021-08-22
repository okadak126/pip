#' Return the positive part of a value
#' @param x the vector
#' @return the maximum of x and 0
#' @export
pos <- function(x) pmax(x, 0)

## expandWeek <- function(days){
##   weeks0 = days%%7
##   weeks = matrix(0, ncol = 7, nrow = length(weeks0))
##   for(i in 0:6){
##     if(i != 0){
##       weeks[weeks0 == i, i] = 1
##     } else {
##       weeks[weeks0 == i, 7] = 1
##     }
##   }
##   weeks
## }
## ma <- function(z, n=5) {
##   res <- numeric(length(z)-n)
##   for(i in 1:(length(res))){
##     res[i] <- mean(z[i:(i+n-1)])}
##   return(res)
## }

# moving average of past n days, padding the front with NAs
## mamiss <- function(z, n=7){
##   c(rep(NA, n), ma(z, n))
## }

## Does cross validation

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
      foldid[(i*fold_len + 1):n] = i+1 #
    }
  }
  foldid
}

