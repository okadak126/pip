#' Compute the prediction statistics
#'
#' Compute three-day predicted usage, waste, remaining and shortage
#' using the actual usage and predicted usage.
#'
#' @param y is the number of units used at the current day i
#' @param t_pred is the sum of the predicted number of units for days
#'     i+1, i+2, i+3
#' @param initial_expiry_data the number of units that can be used the
#'     current day and the next (2-length vector)
#' @param initial_collection_data is the number of units to collect
#'     for days i, i+1, i+2
#' @param start is when we start the clock. So for example, in this
#'     default invocation the shelf storage (initial_collection_data)
#'     was 60, 60, 60, on days 10, 11, and 12. So the evaluation of
#'     the model begins at day 11 but the decision to order collection
#'     starts on day 13.
#' @return a list with four components, \code{x} is the number of
#'     units to collect, \code{r} is a matrix of two columns
#'     indicating remaining units usable on day i+1 and day i+2,
#'     \code{w} is waste, and \code{s} is shortage.
#' @export
compute_prediction_statistics <- function(y, #actual usage arranged according to date
                                          t_pred,
                                          initial_expiry_data = c(0, 0),
                                          initial_collection_data = c(60, 60, 60),
                                          start = 10) {
  N <- length(y)## remove day column
  r <- matrix(0, nrow = N , ncol = 2)
  x <- numeric(N + 3)
  w <- s <- numeric(N)
  x[seq.int(start, start + 2)] <- initial_collection_data

  for (i in seq.int(start, N)) {
    if(i == start){
      w[i] <- pos(initial_expiry_data[1] - y[i])
      r[i, 1] <- pos(initial_expiry_data[1] + initial_expiry_data[2] - y[i] - w[i])
      s[i] <- pos(y[i] - initial_expiry_data[1] - initial_expiry_data[ 2] - x[i])
      r[i, 2] <- pos(x[i] - pos(y[i] - initial_expiry_data[1] - initial_expiry_data[2]))
    }else{
      w[i] <- pos(r[i - 1 , 1] - y[i])
      r[i, 1] <- pos(r[i - 1, 1] + r[i - 1, 2] - y[i] - w[i])
      s[i] <- pos(y[i] - r[i - 1, 1] - r[i - 1, 2] - x[i])
      r[i, 2] <- pos(x[i] - pos(y[i] - r[i - 1, 1] - r[i - 1, 2]))
    }
    x[i + 3] <- max(floor(pos(t_pred[i] - x[i + 1] - x[i + 2] - r[i, 1] - r[i, 2] + 1)))
  }
  return(list(x = x, r = r, w = w, s = s))
}

#' Compute the Cross-Validation Loss for each hyperparameter
#' @param w_cv a vector of waste determined by pip::compute_prediction_statistics
#' @param r_cv a vector of remaining inventory determined by pip::compute_prediction_statistics
#' @param s_cv a vector of inventory shortage as determined by pip::compute_prediction_statistics
#' @param n number of days to consider (length of w_cv, s_cv)
#' @param lo_inv_limit threshold for insufficient inventory below which to penalize (shortage risk)
#' @param hi_inv_limit threshold for surfeit inventory above which to penalize (wastage risk)
compute_cv_loss <- function(w_cv, r_cv, s_cv, n,
                            lo_inv_limit, hi_inv_limit) {

  waste_loss <- apply(w_cv, 2, function(x) sum(x)) # average waste
  invmax_loss <- apply(r_cv, 2, function(x) sum(pos(x - hi_inv_limit))) # penalty on surfeit
  invmin_loss <- apply(r_cv, 2, function(x) sum(pos(lo_inv_limit - x))) # penalty on dearth
  short_loss <- apply(s_cv, 2, function(x) sum(x)) # average shortage

  #pos_rmse <- apply(predictions_cv, 2, function(x) sqrt(sum((pos(t_true - x)^2)[-i_ignore]) / max(sum(t_true - x > 0), 1) ))
  #neg_rmse <- apply(predictions_cv, 2, function(x) sqrt(sum((pos(t_true - x)^2)[-i_ignore]) / max(sum(t_true - x <= 0), 1) ))

  cv_loss <- ((waste_loss + invmax_loss) + (invmin_loss  + short_loss)) / n

  # smooth the loss to avoid random dips
  cv_loss[-c(1, length(cv_loss))] <- stats::filter(cv_loss, rep(1/3, 3))[-c(1, length(cv_loss))]
  cv_loss
}

#' Run cross-validation method
#' @param data the full dataset (data.frame) used as input to the model
#' @param pred_var_indices the column indices corresponding to the features in "data"
#' @param resp_var_index the (single) column index corresponding to the response in "data
#' @param l1_bounds vector containing the possible values of the l1 bound on coefs aside from
#'     days of the week and seven day moving average
#' @param lag_bounds vector containing possible values of the bound on the seven day moving average
#'     (lag) coefficient
#' @param c0 the lower bound on remaining "fresh" inventory (used in optimization)
#' @param start the first day in the dataset that the model's predictions are evaluated
#' @param seed_prop the proportion of seed data (treated as fold 1 for training)
#' @param fold_size the desired size of each fold (aside from the seed fold 1)
#' @param cv_type the variety of cross_validation used. The model can either train on all previous
#'     folds (exp for exanding) or on all other folds besides the left-out fold (loo for leave-one-out, default)
#' @param show_progress a TRUE/FALSE flag for showing progress, default FALSE
cross_validate <- function(data, pred_var_indices, resp_var_index, l1_bounds, lag_bounds,
                           c0, start, seed_prop, fold_size, cv_type = "loo", show_progress = FALSE) {

  XX <- cbind(1, as.matrix(data[, pred_var_indices])) # this is the input - add intercept
  y <- data[, resp_var_index]                         # true usage values

  nhypers <- length(l1_bounds) * length(lag_bounds)

  ## Set up cross-validation scheme
  N <- floor(nrow(data) * (1 - seed_prop))
  seed_length <- nrow(data) - N
  nfolds <- ceiling(N / fold_size)

  foldid <- create_folds(N, fold_size) # assign CV fold ids to each row
  foldid <- c(rep(0, seed_length), foldid) + 1 # Tack on the seed length (add 1 to prevent 0 ind)
  print(foldid)

  preds_cv <- w_cv <- r_cv <- s_cv <- matrix(0,
                                             nrow = nrow(data),
                                             ncol = nhypers)

  #the optimization always assume that we make prediction by the end of the day, say 23:59:59
  pb <- if (show_progress) utils::txtProgressBar(min = 0, max = nfolds, style = 3) else NULL

  for (k in seq_len(nfolds)) {

    # Solve for the coefficients and predict
    indices <- if (cv_type == "exp") { which(foldid < k + 1) } else { which(foldid != k + 1) }

    for (l in seq_len(nhypers)) {
      coefs <- single_lpSolve(data,
                              l1_bounds = l1_bounds[(l - 1) / length(lag_bounds) + 1],
                              lag_bounds = lag_bounds[(l - 1) %% length(lag_bounds) + 1],
                              num_vars = length(pred_var_indices), # exclude date and response
                              ind = indices,
                              start = start,
                              c = c0)

      t_pred <- XX %*% coefs # usage for next 3 days (for each L and lag bound)

      # There are predictions for each l1_bound.
      # Compute waste, shortage, and remaining inventory for each prediction

      rr <- compute_prediction_statistics(y,
                                          t_pred = t_pred,
                                          initial_expiry_data = c(0, 0),
                                          initial_collection_data = y[seq.int(start + 1L, start + 3L)],
                                          start = start)

      w_cv[foldid == k + 1, l] <- rr$w[foldid == k + 1]
      r_cv[foldid == k + 1, l] <- rr$r[foldid == k + 1, 1] + rr$r[foldid == k + 1, 2]
      s_cv[foldid == k + 1, l] <- rr$s[foldid == k + 1]
      preds_cv[foldid == k + 1, l] <- t_pred[foldid == k + 1]
    }

    if (show_progress) utils::setTxtProgressBar(pb, k)
  }
  if (show_progress) close(pb)

  i_ignore <- c(seq_len(seed_length), (nrow(data) - 2):nrow(data))

  list(w_cv = w_cv[-i_ignore, ],
       r_cv = r_cv[-i_ignore, ],
       s_cv = s_cv[-i_ignore, ],
       predictions_cv = preds_cv[-i_ignore, ])

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#' Run a fit on a given dataset
#'
#' The dataset should have one more column than the predictors and
#' response. This additional column, which is typically a date or day
#' number, has to be named "date" or "day"
#' @param data the dataset
#' @param c0 the lower bound on remaining "fresh" inventory (used in optimization)
#' @param history_window Number of days to look back
#' @param penalty_factor penalty for shortage specified by doctors
#' @param start the first day in the dataset that the model's predictions are evaluated
#' @param l1_bounds vector containing the possible values of the l1 bound on coefs aside from
#'     days of the week and seven day moving average
#' @param lag_bounds vector containing possible values of the bound on the seven day moving average
#'     (lag) coefficient
#' @param date_column the name of the date or day number column as a
#'     regex, default is "day|date" i.e. day or date
#' @param response_column the name of the response column, default is
#'     "plt_used"
#' @param show_progress a TRUE/FALSE flag for showing progress,
#'     default TRUE
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' build_model(c0 = 30,
#'            history_window = 200,
#'            penalty_factor = 15,
#'            start = 10,
#'            response_column = "plt.used",
#'            data = day3_data,
#'            show_progress = FALSE)
#'
#' @export
#'
build_model <- function(data, ## The data set
                        c0 = 15, ## the minimum number of units to keep on shelf (c)
                        history_window = 200, ## Number of days to look back
                        penalty_factor = 15, ## the penalty factor for shortage, specified by doctor
                        start = 10,   ## the day you start evaluating
                        l1_bounds = seq(from = 200, to = 0, by = -2), ## allowed values of the l1 bound
                        lag_bounds = c(NA), ## vector of bounds on coefficient for seven day moving average of daily use (NA = no bound)
                        date_column = "day|date",  ## we assume column is named date or day by default
                        response_column = "plt_used",
                        show_progress = TRUE) {

  ## Check that we have enough history!
  n <- nrow(data)
  if (history_window > n) {
    stop(sprintf("build_model: history_window(%d) must be less than data rows(%d)",
                 history_window, n))
  }

  ## Preprocess data
  d <- data[seq.int(n - history_window + 1, n), ]

  data_column_names <- names(data)
  number_of_columns <- ncol(data)
  print(paste0("Number of features:", number_of_columns - 2)) # This is the number of features used

  resp_var_index <- grep(response_column, data_column_names)
  date_var_index <- grep(date_column, data_column_names)
  if (length(resp_var_index) < 1 || length(date_var_index) < 1) {
    stop("Response and/or date columns not found!")
  }
  pred_var_indices <- setdiff(seq.int(2L, number_of_columns), c(resp_var_index, date_var_index))
  predictor_names <- data_column_names[pred_var_indices]

  cv_result <- cross_validate(d, pred_var_indices, resp_var_index,
                              l1_bounds, lag_bounds, c0, start, seed_prop = 0.00,
                              fold_size = 7L, cv_type = "loo", show_progress = TRUE)
  #cv_result <- cross_validate(d, pred_var_indices, resp_var_index,
  #                            l1_bounds, lag_bounds, c0, start, seed_prop = 0.44,
  #                            fold_size = 7L, cv_type = "exp", show_progress = TRUE)

  cv_loss <- compute_cv_loss(cv_result$w_cv,
                             cv_result$r_cv,
                             cv_result$s_cv,
                             n, penalty_factor, 60)

  # ignore the first batch - there seems to be some instability in the first few solves.
  index <- which.min(cv_loss)
  l1_bound_best <- l1_bounds[(index - 1) / length(lag_bounds) + 1]
  lag_bound_best <- lag_bounds[(index - 1) %% length(lag_bounds) + 1]
  print(l1_bound_best)
  print(lag_bound_best)

  coefs <- single_lpSolve(d,
                          l1_bounds = l1_bound_best,
                          lag_bounds = lag_bound_best,
                          num_vars = length(pred_var_indices),
                          start = start,
                          c = c0)

  names(coefs) <- c("intercept", predictor_names)
  list(
    l1_bound = l1_bound_best,
    lag_bound = lag_bound_best,
    w = cv_result$w_cv[ , index],
    r = cv_result$r_cv[ , index],
    s = cv_result$s_cv[ , index],
    cv_loss_order = order(order(cv_loss)), # return list of cv values for each hyperparameter to compare with eval
    coefs = coefs
  )
}


#'
#' Predict the next three day sum and units to order if available
#' @param model the trained model
#' @param new_data a new data frame with the same predictor names as in
#'     the training data, even if one row
#' @return the predicted three day total of how many units to collect
#' @export
#'
predict_three_day_sum <- function(model, ## the trained model
                                  new_data) { ## the new dataset with specified number of features
  predictor_names <- names(model$coefs)[-1]
  common_names <- intersect(predictor_names, names(new_data))
  if (length(common_names) != length(predictor_names)) {
    stop("New data set column names don't match training data set!")
  }
  new_data <- as.matrix(new_data[, predictor_names])
  ceiling(sum(new_data %*% model$coefs[-1] ) + model$coefs[1])
}

#' Run a fit on a given dataset
#'
#' The dataset should have one more column than the predictors and
#' response. This additional column, which is typically a date or day
#' number, has to be named "date" or "day"
#' @param data the dataset
#' @param c0 the c0 value
#' @param train_window Number of days to look back
#' @param test_window Number of days on which to evaluate model
#' @param penalty_factor penalty for shortage specified by doctors
#' @param start the day you start evaluating the model (for CV)
#' @param l1_bounds vector containing the possible values of the l1 bound on coefs aside from
#'     days of the week and seven day moving average
#' @param lag_bounds vector containing possible values of the bound on the seven day moving average
#'     (lag) coefficient
#' @param date_column the name of the date or day number column as a
#'     regex, default is "day|date" i.e. day or date
#' @param response_column the name of the response column, default is
#'     "plt_used"
#' @examples
#'
#' @export
#'
evaluate_model <- function(data, ## The data set
                        c0 = 15, ## the minimum number of units to keep on shelf (c)
                        train_window = 200, ## Number of days to look back
                        test_window = 7, ## Number of days on which to evaluate model
                        penalty_factor = 15, ## the penalty factor for shortage, specified by doctor
                        start = 10, ## the day you start evaluating
                        l1_bounds = seq(from = 200, to = 0, by = -2), ## allowed values of the l1 bound
                        lag_bounds = c(NA),
                        date_column = "day|date",  ## we assume column is named date or day by default
                        response_column = "plt_used") {

  # Validate input dataset
  n <- nrow(data)
  if (train_window + test_window > n) {
    stop(sprintf("build_model: history_window(%d) must be less than data rows(%d)",
                 train_window + test_window, n))
  }
  resp_var_index <- grep(response_column, names(data))
  date_var_index <- grep(date_column, names(data))
  if (length(resp_var_index) < 1 || length(date_var_index) < 1) {
    stop("Response and/or date columns not found!")
  }
  pred_var_indices <- setdiff(seq.int(2L, ncol(data)), c(resp_var_index, date_var_index))

  # Partition data into train and test sets
  data_train <- data[seq.int(n - train_window - test_window + 1, n - test_window), ]
  data_test <- data[seq.int(n - test_window + 1, n), ]
  y_test <- data_test[, resp_var_index]

  train_prop <- nrow(data_train) / (nrow(data_test) + nrow(data_train))

  # Run CV procedure to build model on training data and obtain the "best" hyperparams
  model <- build_model(data = data_train,
                      c0 = c0,
                      history_window = train_window,
                      penalty_factor = penalty_factor,
                      start = start,
                      l1_bounds = l1_bounds,
                      lag_bounds = lag_bounds)

  plot(x = l1_bounds, y = model$cv_loss_order[seq(from = 1, to = length(l1_bounds), by = length(lag_bounds))])
  lines(x = l1_bounds, y = model$cv_loss_order[seq(from = 1, to = length(l1_bounds), by = length(lag_bounds))], type="l")

  # Train the model on training window for all hyperparams and compute stats in test window
  eval_result <- cross_validate(rbind(data_train, data_test), pred_var_indices, resp_var_index,
                                l1_bounds, lag_bounds, c0, start, seed_prop = train_prop,
                                fold_size = test_window, cv_type = "exp")


  # Compute the true 3 day usage
  t_true <- dplyr::lead(y_test, n = 1L, default = 0) +
    dplyr::lead(y_test, n = 2L, default = 0) +
    dplyr::lead(y_test, n = 3L, default = 0)

  i_ignore <- c((length(y_test) - 2):length(y_test))

  avg_waste    <- apply(eval_result$w_cv, 2, function(x) sum(x))
  avg_shortage <- apply(eval_result$s_cv, 2, function(x) sum(x))
  rmse <- apply(eval_result$predictions_cv, 2, function(x) {
    sqrt(sum(((t_true[-i_ignore] - x)^2)) / (test_window - 3))
    })

  # Compute loss for each hyperparameter set for the test set
  cv_loss <- compute_cv_loss(eval_result$w_cv,
                             eval_result$r_cv,
                             eval_result$s_cv,
                             test_window,
                             penalty_factor, 60)

  plot(x = l1_bounds, y = order(order(cv_loss))[seq(from = 1, to = length(l1_bounds), by = length(lag_bounds))])
  lines(x = l1_bounds, y = order(order(cv_loss))[seq(from = 1, to = length(l1_bounds), by = length(lag_bounds))], type="l")

  min_cvloss_idx <- which.min(cv_loss)

  top_idx <- 1 # best hyperparameter combo
  min_short_idx <- which(avg_shortage == min(avg_shortage))
  print(min_short_idx)
  min_waste_idx <- intersect(which(avg_waste == min(avg_waste[min_short_idx])), min_short_idx)
  print(min_waste_idx)
  min_rmse_idx <- max(intersect(which(abs(rmse - min(rmse[min_waste_idx])) <= 1.0), min_waste_idx))

  if (length(min_short_idx) < 2L) {
    top_idx <- min_short_idx
  } else if (length(min_waste_idx) < 2L) {
    top_idx <- min_waste_idx
  } else {
    top_idx <- min_rmse_idx
  }

  l1_bound_train <- model$l1_bound
  lag_bound_train <- model$lag_bound
  l1_bound_test <- l1_bounds[(top_idx - 1) / length(lag_bounds) + 1]
  lag_bound_test <- lag_bounds[(top_idx - 1) %% length(lag_bounds) + 1]
  lag_offset <- if (length(which(lag_bounds == lag_bound_train)) > 0L) which(lag_bounds == lag_bound_train) else 1
  model_idx <- (which(l1_bounds == l1_bound_train) - 1) * length(lag_bounds) + lag_offset

  l1_bound_cvloss <- l1_bounds[(min_cvloss_idx - 1) / length(lag_bounds) + 1]
  lag_bound_cvloss <- lag_bounds[(min_cvloss_idx - 1) %% length(lag_bounds) + 1]

  print(sprintf("Distance between train and test loss vectors: %f", sqrt(sum((model$cv_loss_order - order(order(cv_loss)))^2) / length(cv_loss)) ))
  print(sprintf("Model found best l1 bound was %d, and the best was %d", l1_bound_train, l1_bound_test))
  print(sprintf("According to our CV loss function, the best l1_bound is %d with loss %f", l1_bound_cvloss, cv_loss[min_cvloss_idx]))
  print(sprintf("Model's Waste: %f, Shortage: %f, RMSE: %f", avg_waste[model_idx], avg_shortage[model_idx], rmse[model_idx]))

  print(sprintf("Best Waste: %f, Shortage: %f, RMSE: %f", avg_waste[top_idx], avg_shortage[top_idx], rmse[top_idx]))
  #warning(sprintf("Model found best lag bound was %d, and it was actually %d", lag_bound_train, lag_bound_test))
}
