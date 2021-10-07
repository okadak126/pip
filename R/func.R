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
#' @param c the ideal minimum number of fresh units remaining at the end
#'     of the day. this is also the minimum number that should be collected
#' @return a list with four components, \code{x} is the number of
#'     units to collect, \code{r} is a matrix of two columns
#'     indicating remaining units usable on day i+1 and day i+2,
#'     \code{w} is waste, and \code{s} is shortage.
#' @export
compute_prediction_statistics <- function(y, #actual usage arranged according to date
                                          t_pred,
                                          initial_expiry_data = c(0, 0),
                                          initial_collection_data = c(60, 60, 60),
                                          start = 10,
                                          c = 30) {
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
    x[i + 3] <- max(floor(pos(t_pred[i] - x[i + 1] - x[i + 2] - r[i, 1] - r[i, 2] + 1)), c)
  }
  return(list(x = x, r = r, w = w, s = s))
}

#' Compute the Cross-Validation Loss for each hyperparameter
#' @param preds a vector of 3 day usage predictions
#' @param y a vector of actual daily usage
#' @param w a vector of waste determined by pip::compute_prediction_statistics
#' @param r1 a vector of remaining inventory determined by pip::compute_prediction_statistics
#' @param r2 a vector of remaining inventory determined by pip::compute_prediction_statistics
#' @param s a vector of inventory shortage as determined by pip::compute_prediction_statistics
#' @param penalty_factor factor to additionally penalize shortage terms over waste
#' @param lo_inv_limit threshold for insufficient inventory below which to penalize (shortage risk)
#' @param hi_inv_limit threshold for surfeit inventory above which to penalize (wastage risk)
#' @export
compute_loss <- function(preds, y, w, r1, r2, s, penalty_factor,
                         lo_inv_limit, hi_inv_limit) {

  if (nrow(w) != nrow(s) || nrow(w) != nrow(r1) || nrow(w) != nrow(r2)) {
    stop("Waste, inventory, and/or shortage vectors of different lengths.")
  }

  waste_loss <- apply(w, 2, function(x) sum(x^2)) # average squared waste
  short_loss <- apply(s, 2, function(x) sum(x^2)) # average squared shortage

  ## true next-three-day usage for eachd day i
  t_true <- (dplyr::lead(y, n = 1L, default = 0) +
    dplyr::lead(y, n = 2L, default = 0) +
    dplyr::lead(y, n = 3L, default = 0))[-c((length(y) - 2):length(y))]

  ## Add bias
  pos_rss <- apply(preds, 2, function(x) sum((pos(x - t_true - lo_inv_limit)^2) ))
  neg_rss <- apply(preds, 2, function(x) sum((pos(t_true + lo_inv_limit - x)^2) ))

  ## Loss = average(RSS_pos_bias + RSS_neg_bias + squared_waste + pi^2 * squared_shortage)
  ## This takes into account how much we prefer waste over shortage situations,
  ## but biased RSS makes the function more sensitive (waste and shortage are relatively rare)
  loss <- (pos_rss + neg_rss + waste_loss + penalty_factor^2 * short_loss) / nrow(w)
  loss
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
#' @param penalty_factor how much we penalize shortage over waste
#' @param seed_prop the proportion of seed data (treated as fold 1 for training)
#' @param fold_size the desired size of each fold (aside from the seed fold 1)
#' @param cv_type the variety of cross_validation used. The model can either train on all previous
#'     folds (exp for exanding) or on all other folds besides the left-out fold (loo for leave-one-out, default)
#' @param show_progress a TRUE/FALSE flag for showing progress, default FALSE
cross_validate <- function(data, pred_var_indices, resp_var_index, l1_bounds, lag_bounds,
                           c0, start, penalty_factor, seed_prop, fold_size, cv_type = "loo", show_progress = FALSE) {

  XX <- cbind(1, as.matrix(data[, pred_var_indices])) # this is the input - add intercept
  y <- data[, resp_var_index]                         # true usage values

  nhypers <- length(l1_bounds) * length(lag_bounds)

  ## Set up cross-validation scheme
  N <- floor(nrow(data) * (1 - seed_prop))
  seed_length <- nrow(data) - N
  nfolds <- ceiling(N / fold_size)

  foldid <- create_folds(N, fold_size) # assign CV fold ids to each row
  foldid <- c(rep(0, seed_length), foldid) + 1 # Tack on the seed length (add 1 to prevent 0 ind)

  preds_cv <- w_cv <- r2_cv <- r1_cv <- s_cv <- matrix(0,
                                             nrow = nrow(data),
                                             ncol = nhypers)

  #the optimization always assume that we make prediction by the end of the day, say 23:59:59
  pb <- if (show_progress) utils::txtProgressBar(min = 0, max = nfolds, style = 3) else NULL

  for (k in seq_len(nfolds)) {


    # Solve for the coefficients and predict
    indices <- if (cv_type == "exp") { which(foldid < k + 1) } else { which(foldid != k + 1) }

    # Solve all at once.
    coefs <- single_lpSolve(data,
                            l1_bounds = l1_bounds,
                            lag_bounds = lag_bounds,
                            num_vars = length(pred_var_indices), # exclude date and response
                            ind = indices,
                            start = start,
                            c = c0,
                            shortage_factor = penalty_factor)

    t_pred <- ceiling(XX %*% coefs) # usage for next 3 days (for each L and lag bound)

    for (l in seq_len(nhypers)) {
      rr <- compute_prediction_statistics(y,
                                          t_pred = t_pred[, l],
                                          initial_expiry_data = c(0, 0),
                                          initial_collection_data = y[seq.int(start, start + 2L)],
                                          start = start,
                                          c = c0)

      w_cv[foldid == k + 1, l] <- rr$w[foldid == k + 1]
      r1_cv[foldid == k + 1, l] <- rr$r[foldid == k + 1, 1]
      r2_cv[foldid == k + 1, l] <- rr$r[foldid == k + 1, 2]
      s_cv[foldid == k + 1, l] <- rr$s[foldid == k + 1]
      preds_cv[foldid == k + 1, l] <- t_pred[foldid == k + 1, l]
    }

    if (show_progress) utils::setTxtProgressBar(pb, k)
  }
  if (show_progress) close(pb)

  initial_mask <- start + 5L # mask the first 6 days including the "start" day. This is so that
                            # the initial values do not affect the waste or shortage results

  # Also drop the last 3 days, since the 3 day predictions cannot include the last 3 days in dataset.
  i_ignore <- c(union(seq_len(initial_mask), seq_len(seed_length)), (nrow(data) - 2):nrow(data))

  list(w = w_cv[-i_ignore, , drop = FALSE],
       r1 = r1_cv[-i_ignore, , drop = FALSE],
       r2 = r2_cv[-i_ignore, , drop = FALSE],
       s = s_cv[-i_ignore, , drop = FALSE],
       preds = preds_cv[-i_ignore, , drop = FALSE])

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
#' @param lo_inv_limit inventory count below which we penalize result (potential short)
#' @param hi_inv_limit inventory count above which we penalize result (potential waste)
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
                        penalty_factor = 5, ## the penalty factor for shortage, specified by doctor
                        lo_inv_limit = 30,  ## inventory count below which we penalize result
                        hi_inv_limit = 60,  ## inventory count above which we penalize result
                        start = 10,   ## the day you start evaluating
                        l1_bounds = seq(from = 200, to = 0, by = -2), ## allowed values of the l1 bound
                        lag_bounds = -1, ## vector of bounds on coefficient for seven day moving average of daily use (NA = no bound)
                        date_column = "day|date",  ## we assume column is named date or day by default
                        response_column = "plt_used",
                        show_progress = TRUE,
                        plot_losses = TRUE) {

  ## Check that we have enough history
  n <- nrow(data)
  if (history_window > n) {
    stop(sprintf("build_model: history_window(%d) must be less than data rows(%d)",
                 history_window, n))
  }

  ## Preprocess data
  d <- data[seq.int(n - history_window + 1, n), , drop = FALSE]
  data_column_names <- names(data)
  print(paste0("Number of features:", ncol(data) - 2)) # This is the number of features used

  ## Find and exclude the date and response columns from the dataset
  resp_var_index <- grep(response_column, data_column_names)
  date_var_index <- grep(date_column, data_column_names)
  if (length(resp_var_index) < 1 || length(date_var_index) < 1) {
    stop("Response and/or date columns not found!")
  }
  pred_var_indices <- setdiff(seq.int(2L, ncol(data)), c(resp_var_index, date_var_index))
  predictor_names <- data_column_names[pred_var_indices]

  ## Obtain predictions, waste, loss, remaining inventory, and shortage from
  ## best hyperparameter via cross-validation
  cv_result <- cross_validate(d, pred_var_indices, resp_var_index,
                              l1_bounds, lag_bounds, c0, start, penalty_factor, seed_prop = 0.00,
                              fold_size = 14L, cv_type = "loo", show_progress = TRUE)

  ## Calculate losses for each of the CV results and obtain coefficients from best loss
  cv_loss <- compute_loss(cv_result$preds,
                          d[1:(nrow(d) - start - 5), resp_var_index],
                          cv_result$w,
                          cv_result$r1,
                          cv_result$r2,
                          cv_result$s,
                          penalty_factor, lo_inv_limit, hi_inv_limit)


  index <- which.min(cv_loss)
  l1_bound_min <- l1_bounds[(index - 1) / length(lag_bounds) + 1]
  lag_bound_min <- lag_bounds[(index - 1) %% length(lag_bounds) + 1]

  if (plot_losses) {
    plot(x = l1_bounds, y = cv_loss[seq(from = 1, to = length(cv_loss), by = 2)])
    print(l1_bound_min)
    print(lag_bound_min)
  }

  ### Solve using hyperparameter value with least loss
  coefs <- single_lpSolve(d,
                          l1_bounds = l1_bound_min,
                          lag_bounds = lag_bound_min,
                          num_vars = length(pred_var_indices),
                          start = start,
                          c = c0,
                          shortage_factor = penalty_factor)


  # Return coefficients and other features of the CV result
  names(coefs) <- c("intercept", predictor_names)
  list(
    l1_bound = l1_bound_min,
    lag_bound = lag_bound_min,
    w = cv_result$w[ , index],
    r1 = cv_result$r1[ , index],
    r2 = cv_result$r2[ , index],
    s = cv_result$s[ , index],
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
  ceiling(new_data %*% model$coefs[-1] + model$coefs[1])
}

#' Evaluate the model on a held out validation set within the original training set
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
                        penalty_factor = 5, ## the penalty factor for shortage, specified by doctor
                        lo_inv_limit = 30,  ## inventory count below which we penalize result
                        hi_inv_limit = 60,  ## inventory count above which we penalize result
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
  data_train <- data[seq.int(n - train_window - test_window + 1, n - test_window), , drop = FALSE]
  data_test <- data[seq.int(n - test_window + 1, n), , drop = FALSE]
  data_full <- data[seq.int(n - train_window - test_window + 1, n), , drop = FALSE]

  y_test <- data_test[, resp_var_index]
  y_full <- data_full[, resp_var_index]

  train_prop <- nrow(data_train) / (nrow(data_full))

  # Run CV procedure to build model on training data and obtain the "best" hyperparams
  model <- build_model(data = data_train,
                      c0 = c0,
                      history_window = train_window,
                      penalty_factor = penalty_factor,
                      start = start,
                      l1_bounds = l1_bounds,
                      lag_bounds = lag_bounds)


  t_pred_cv <- predict_three_day_sum(model, data_full)

  # Compute the prediction statistics for the best param
  model_result <- compute_prediction_statistics(y = y_full,
                                                t_pred = t_pred_cv,
                                                initial_expiry_data = c(0, 0),
                                                initial_collection_data = y_full[seq.int(start, start + 2L)],
                                                start = start,
                                                c = c0)


  min_cv_loss <- compute_loss(matrix(t_pred_cv[seq.int(n - test_window + 1L, n - 3L)], ncol = 1),
                          y_test,
                          matrix(model_result$w[seq.int(n - test_window + 1L, n - 3L)], ncol = 1),
                          matrix(model_result$r[seq.int(n - test_window + 1L, n - 3L), 1L], ncol = 1),
                          matrix(model_result$r[seq.int(n - test_window + 1L, n - 3L), 2L], ncol = 1),
                          matrix(model_result$s[seq.int(n - test_window + 1L, n - 3L)], ncol = 1),
                          penalty_factor,
                          lo_inv_limit,
                          hi_inv_limit)

  nhypers <- length(l1_bounds) * length(lag_bounds)
  eval_losses <- rep(0, nhypers)
  X_full <- cbind(1, as.matrix(data_full[, pred_var_indices]))

  for (l in seq_len(nhypers)) {
      coefs <- single_lpSolve(data_train,
                              l1_bounds = l1_bounds[(l - 1) / length(lag_bounds) + 1],
                              lag_bounds = lag_bounds[(l - 1) %% length(lag_bounds) + 1],
                              num_vars = length(pred_var_indices), # exclude date and response
                              start = start,
                              c = c0,
                              shortage_factor = penalty_factor)

      t_pred_test <- ceiling(X_full %*% coefs) # usage for next 3 days (for each L and lag bound)

      # There are predictions for each l1_bound.
      # Compute waste, shortage, and remaining inventory for each prediction
      rr <- compute_prediction_statistics(y_full,
                                          t_pred = t_pred_test,
                                          initial_expiry_data = c(0, 0),
                                          initial_collection_data = y_full[seq.int(start, start + 2L)],
                                          start = start,
                                          c = c0)

      eval_loss <- compute_loss(matrix(t_pred_test[seq.int(n - test_window + 1L, n - 3L)], ncol = 1L),
                                y_test,
                                matrix(rr$w[seq.int(n - test_window + 1L, n - 3L)], ncol = 1L),
                                matrix(rr$r[seq.int(n - test_window + 1L, n - 3L),1L], ncol = 1L),
                                matrix(rr$r[seq.int(n - test_window + 1L, n - 3L),2L], ncol = 1L),
                                matrix(rr$s[seq.int(n - test_window + 1L, n - 3L)], ncol = 1L),
                                penalty_factor,
                                lo_inv_limit,
                                hi_inv_limit)

      eval_losses[l] <- eval_loss
  }

  min_eval_loss <- min(eval_losses)
  min_eval_loss_idx <- which.min(eval_losses)

  l1_bound_best_eval <- l1_bounds[(min_eval_loss_idx - 1) / length(lag_bounds) + 1]
  lag_bound_best_eval <- lag_bounds[(min_eval_loss_idx - 1) %% length(lag_bounds) + 1]

  print(sprintf("CV Found best L was %d, loss was %f", model$l1_bound, min_cv_loss))
  print(sprintf("Test Set found best L was %d, loss was %f", l1_bound_best_eval, min_eval_loss))

}
