#This is the function tranining a model with fixed lambda(1)
#1:p: betas
#p+1:p+N:w
#p+N+1:p+2N:r1
#p+2N+1:p+3N:r2
#p+3N+1:p+4N:x
#p+4N+1:p+5N: t
#p+5N+1:intercept
#p+5N+2:p+5N+1+p:l1 bounds
#########


#' Solve the LP problem for platelet usage prediction
#' @param d the data frame. We assume that the first column contains a date and the response is the
#'   last column.
#' @param l1_bounds the upper bound on the l1 norm of coefficient vector (aside from lag and days of week)
#' @param lag_bounds the upper bound on the seven day moving average of usage coefficient (lag)
#' @param num_vars the number of features/covariates in the dataset.
#' @param start the 0th index of the date when we begin evaluating the model
#' @param c the value for c_0
#' @param shortage_factor the amount by which we multiply shortage over waste (ratio)
#' @param buffer a value to be added to account for training in the initial stages, default 10
#' @param ind the indices for cross validation used only if non-null. (Should be sorted!!)
#' @param reset_basis TRUE or FALSE for whether we reset the LP basis before each solve (for non-uniqueness)
#' @importFrom lpSolveAPI make.lp add.constraint set.bounds lp.control delete.constraint
#' @importFrom lpSolveAPI set.objfn get.variables set.rhs get.basis set.basis
#' @examples
#' data_file <- system.file("extdata", "sample_data.Rdata", package = "pip")
#' load(data_file)
#'
#' #' # Note that SBCpip::create_dataset() produces data of the form lpdata from source files
#' sol <- single_lpSolve(lpdata, l1_bounds, lag_bounds, num_vars)
#'
#' # Add the coefficient names to the resulting coefficients. Each
#' # column represents a pairing of hyperparameter values from l1_bounds
#' # and lag_bounds. The effects of the L1 penalty should be apparent.
#' rownames(sol) <- c("Intercept", colnames(lpdata[2:(ncol(lpdata) - 1)]))
#' @export
#'
single_lpSolve <- function(d, l1_bounds, lag_bounds, num_vars,
                           start = 10, c = 10, shortage_factor = 15,
                           buffer = 10, ind = NULL, reset_basis = FALSE) {
  doing_cv <- !is.null(ind)

  ## ind and l1_bounds must be sorted (ascending and descending, respectively)

  # Assumes that date is the first column and the output (usage) is the last column.
  XX <- d[, 2:(num_vars + 1)]
  y <- d[, (num_vars + 2)]

  N <- nrow(XX)
  p <- ncol(XX)

  if (!doing_cv) ind <- 1:N

  ## y is usage, so y[i] is known usage on day i
  ## p betas, N ws, N r1s, N r2s, N xs, N ts, N s's, 1 intercept, p bounds
  N_var <- p + 6*N + 1 + p

  # renaming 0th index for each set of decision variables for clarity
  w_idx0 <- p         # waste
  r1_idx0 <- p + N    # remaining inventory expiring in 1 day
  r2_idx0 <- p + 2*N  # remaining inventory expiring in 2 days
  x_idx0 <- p + 3*N   # collection
  t_idx0 <- p + 4*N   # predicted 3 day usage (i+1, i+2, i+3)
  s_idx0 <- p + 5*N   # shortage (KO Added)
  bb_idx0 <- p + 6*N + 1 # absolute value bounds on beta (coefficients) - we skip the intercept

  # Altered to penalize r1
  pind_w  <- w_idx0 + ind
  pind_s <- s_idx0 + ind # indices of waste for CV fold (minimize sum of these in obj func)
  pind_r1 <- r1_idx0 + ind
  pind_r2 <- r2_idx0 + ind

  my.lp <- lpSolveAPI::make.lp(0,N_var)
  lp.control(my.lp,
             timeout=60) # Should be included in config

  obj_coefficients <- rep(0,N_var)

  ## Set coefs of w to be 1s and those of and s to be 5s (penalty factor) (eq 9)
  obj_coefficients[pind_w] <- 1
  obj_coefficients[pind_s] <- shortage_factor
  obj_coefficients[pind_r1] <- 0.001 # These add stability / uniqueness to the solution (break ties)
  obj_coefficients[pind_r2] <- -0.001
  lpSolveAPI::set.objfn(my.lp,obj_coefficients)

  ## Constrain sum of day of week coefs of beta to be zero
  ## to prevent non-identifiability and improve interpretability
  constraint_coefficients <- rep(0,N_var)
  constraint_coefficients[1:7] <- 1
  lpSolveAPI::add.constraint(my.lp, constraint_coefficients, type = "=", 0)

  ## w >=0, s >= 0, r1 >= 0, r2 >= 0, x >= 0, t >= 0
  for (i in 1:N){
    ## w_{i} >= 0
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[w_idx0+i] <- 1
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, type = ">=", 0)

    ## s_{i} >= 0
    constraint_coefficients <- rep(0, N_var)
    constraint_coefficients[s_idx0+i] <- 1
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, type = ">=", 0)

    ## r1_{i} >= 0
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[r1_idx0+i] <- 1 # r_i(1) >= 0
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    ## r2_{i} >= 0 (further constrained >= c later)
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[r2_idx0+i] <- 1 ## r_i(2) >= 0
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    ## x_{i} >= 0
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[x_idx0+i] <- 1 ## x_i >= 0
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    # t_{i} >= 0
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[t_idx0+i] <- 1 ## t_i >= 0
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

  }

  ## The initial values for x_{start+1}, x_{start+2}, x_{start+3} ensure
  ## the feasibility of this problem.
  rhs_offset <- c(c + buffer, 0, 0)
  for (i in seq_len(3L)) {
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[x_idx0+start+i] <- 1 ## x_{start + i}
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", y[start+i] + rhs_offset[i])
  }

  ## Set r1 to zero for start day
  constraint_coefficients <- rep(0,N_var)
  constraint_coefficients[r1_idx0+start] <- 1 ## r_{start}(1)
  lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=",0)

  ## Setting r2 to zero for start day
  constraint_coefficients <- rep(0,N_var)
  constraint_coefficients[r2_idx0+start] <- 1 ## r_{start}(2)
  lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=",0)

  ## Constraints on three-day prediction (according to a linear model)
  ## t_i = XX*beta + intercept [Eqn 11]
  ## Note that t_start is completely determined by our initial collection and inventory values
  for(i in (start+1):N) {
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[t_idx0+i] <- 1 ## t_i
    constraint_coefficients[1:p] <- -as.numeric(XX[i,]) ## p features on day i
    constraint_coefficients[p+6*N+1] <- -1 ## negate intercept
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", 0)
  }

  ## Constraints on remaining inventories, collection, waste, and shortage
  ## r1_{start} and r2_{start} are already fixed
  for (i in (start+1):N){

    ## Inequalities constraining r1_{i} [3]
    ## r1_{i} - r2_{i-1} - r1_{i-1} + w_{i} >= -y_{i} (Eqn 14)
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[r1_idx0+i] <- 1 #r_i(1)
    constraint_coefficients[w_idx0+i] <- 1 #w_i
    constraint_coefficients[r1_idx0+i-1] <- -1 #r_{i-1}(1)
    constraint_coefficients[r2_idx0+i-1] <- -1 # r_{i-1}(2)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", -y[i])

    ## r1_{i} - r2_{i-1} <= 0
    constraint_coefficients <- rep(0, N_var)
    constraint_coefficients[r1_idx0+i] <- 1 #r_i(1)
    constraint_coefficients[r2_idx0+i-1] <- -1  #r_{i-1}(2)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", 0)

    ## Inequalities constraining r2_{i}
    ## Note: r2_{i} >= x_{i} + r1_{i-1} + r2_{i-1} - r1_{i} - w_{i} - y_{i}
    ## is implied by the equation below

    ## r2_{i} - x_{i} <= 0
    constraint_coefficients <- rep(0, N_var)
    constraint_coefficients[r2_idx0+i] <- 1
    constraint_coefficients[x_idx0+i] <- -1
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", 0)

    ## Inequalities constraining waste [2]
    ##w_{i} >= r1_{i-1} - y_{i} (Eqn 12)
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[w_idx0+i] <- 1 ## w_i
    constraint_coefficients[r1_idx0+i-1] <- -1 ## r_{i-1}(1)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", -y[i])

    ## w_i <= r1_{i-1}
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", 0)

    ## Inequalities constraining shortage [2]

    ## s_{i} + r1_{i-1} + r2_{i-1} + x_{i} - w_{i} >= y_{i}
    ## Implied by the equation below.
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[w_idx0+i] <- -1 ## w_i
    constraint_coefficients[x_idx0+i] <- 1 ## x_i
    constraint_coefficients[r1_idx0+i-1] <- 1 ## r_{i-1}(1)
    constraint_coefficients[r2_idx0+i-1] <- 1 ## r_{i-1}(2)
    constraint_coefficients[s_idx0+i] <- 1 ## s_i
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", y[i])

    # s_{i} <= y_{i} (we can only be short as much as we need that day)
    constraint_coefficients <- rep(0, N_var)
    constraint_coefficients[s_idx0+i] <- 1
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", y[i])

    ## Primary Equation [Eq 15] - This makes the problem solvable (KO added shortage)
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[r1_idx0+i] <- 1 ## r_i(1)
    constraint_coefficients[r2_idx0+i] <- 1 ## r_i(2)
    constraint_coefficients[w_idx0+i] <- 1 ## w_i
    constraint_coefficients[x_idx0+i] <- -1 ## x_i
    constraint_coefficients[r1_idx0+i-1] <- -1 ## r_{i-1}(1)
    constraint_coefficients[r2_idx0+i-1] <- -1 ## r_{i-1}(2)
    #if (c == 0)
    constraint_coefficients[s_idx0+i] <- -1 ## s_i
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", -y[i])

  }

  ## Constraint all r2 after start to be greater than c (Equation 17)
  for (i in ind) {
    if(i >= start+4){
      constraint_coefficients <- rep(0,N_var)
      constraint_coefficients[r2_idx0+i] = 1 ## r_i(2)
      lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", c)
    }
  }

  ## collection (note that this is an optimistic strategy, anticipating no waste)
  ## x_{i+3}+x_{i+2}+x_{i+1} + r2_{i} + r1_{i} -t_i== 0:start+N
  ## We have x_{start + 1}, x_{start + 2}, x_{start + 3}, r1_{start}, and r2_{start} fixed,
  ## so we begin evaluating from start + 4
  for (i in seq.int(start+4, N)) {

    ## Implements equation 8
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[x_idx0+i] <- 1 ## x_{i}
    if (i %in% ind) {
      constraint_coefficients[x_idx0+i-1] <- 1 ## x_{i-1}
      constraint_coefficients[x_idx0+i-2] <- 1 ## x_{i-2}
      constraint_coefficients[r1_idx0+i-3] <- 1 ## r_{i-3}(1)
      constraint_coefficients[r2_idx0+i-3] <- 1 ## r_{i-3}(2)
      constraint_coefficients[t_idx0+i-3] <- -1 ## t_{i-3}

      # Take into account potential for waste and shortage (partial)
      lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", 0)
    } else {
      lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", y[i])
    }
  }

  ## Coefficient Constraints
  for (j in 1:p) {
    ## constrain bounds (l1 constraints) KO added this
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[bb_idx0+j] <- 1 ## l1 constraint
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    ## penalty term constraints (constrain below 0)
    constraint_coefficients[j] <- 1 ## beta_j
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    ## penalty term constraints (constrain above 0)
    constraint_coefficients[j] <- -1  ## beta_j
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)
  }

  lpSolveAPI::set.bounds(my.lp, lower = rep(-1000,ncol(my.lp)), upper = rep(1000,ncol(my.lp)))
  lpSolveAPI::lp.control(my.lp,sense='min')

  coefficients_matrix <- matrix(0, ncol = length(l1_bounds) * length(lag_bounds), nrow = p+1)

  # Create a constant basis (this is in case there are multiple optimal solutions)
  pre_rhs <- lpSolveAPI::get.rhs(my.lp) # grab current RHS

  lpSolveAPI::set.rhs(my.lp, rep(0, length(pre_rhs))) # set RHS to 0
  lpSolveAPI::set.bounds(my.lp, lower = rep(0,ncol(my.lp)), upper = rep(0,ncol(my.lp))) # set bounds to 0
  lpSolveAPI::set.objfn(my.lp, rep(1, N_var)) # reset objfn to sum of all vars (unique solution)
  status_init <- solve(my.lp)

  obj_coefficients[pind_w] <- 1
  obj_coefficients[pind_s] <- shortage_factor
  obj_coefficients[pind_r1] <- 0.001
  obj_coefficients[pind_r2] <- -0.001
  lpSolveAPI::set.objfn(my.lp, obj_coefficients)
  pre_basis <- lpSolveAPI::get.basis(my.lp)
  lpSolveAPI::set.rhs(my.lp, pre_rhs)
  lpSolveAPI::set.bounds(my.lp, lower = rep(-1000,ncol(my.lp)), upper = rep(1000,ncol(my.lp)))

  for (i in seq_along(l1_bounds)) {

    # Reset basis
    if (reset_basis)
      lpSolveAPI::set.basis(my.lp, pre_basis)

    ## Vars 1-8 are not constrained by the l1_bound
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[(bb_idx0 + 9):(bb_idx0 + p)] <- 1
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", l1_bounds[i])

    nConstraints <- nrow(my.lp) # index refers to l1 constraint

    for (j in seq_along(lag_bounds)) {

      # Apply a common constraint to the first and second lags
      if (lag_bounds[j] > 0) {
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[bb_idx0 + 8] <- 1
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", lag_bounds[j])
      }

      status <- solve(my.lp)

      # If a timeout occurs, try to solve the problem with looser constraints
      if (status == 7 && !doing_cv) {
        while (status == 7) {
          warning("A timeout occurred. Loosening L1 constraints", call. = FALSE)
          lpSolveAPI::delete.constraint(my.lp, nConstraints) # delete the L1 constraint

          # relax the L1 constraints aside from day of week and seven day moving average
          constraint_coefficients <- rep(0,N_var)
          constraint_coefficients[(bb_idx0 + 9):(bb_idx0 + p)] <- 1
          lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", l1_bounds[i] + 8)
          status <- solve(my.lp)
        }
      }

      if (doing_cv && (status > 0))  {
        # If we are running cross validation, we can accept some failures
        warning(sprintf("lpSolveStatus is %d", status))
        coefficients_matrix[ , (i - 1) * length(lag_bounds) + j] <- 0.0 # set all coefs to 0

      }
      else if (status > 0) {
        stop(sprintf("lpSolveStatus is %d", status))
      }
      else {
        coeffs <- c(lpSolveAPI::get.variables(my.lp)[6*N+p+1], lpSolveAPI::get.variables(my.lp)[1:p])
        coeffs[abs(coeffs) < 1e-10] <- 0.0 # set coeffs below tolerance to 0 (may fix some numeric issues)
        coefficients_matrix[ , (i - 1) * length(lag_bounds) + j] <- coeffs  ## all the beta_j's
      }

      # If a constraint was added, clean up the lag bound
      while (nConstraints < nrow(my.lp)) {
        lpSolveAPI::delete.constraint(my.lp, nrow(my.lp)) # delete the lag bound constraint
      }

    }
    # delete the L1 bound constraint
    lpSolveAPI::delete.constraint(my.lp, nConstraints)

  }

  rm(my.lp)

  coefficients_matrix

}

#' Solve the LP problem for platelet usage prediction (ignoring shortage in the optimization)
#' @param d the data frame
#' @param l1_bounds the upper bound on the l1 norm of coefficient vector (aside from lag and days of week)
#' @param lag_bounds the upper bound on the seven day moving average of usage coefficient (lag)
#' @param num_vars the number of features/covariates in the dataset.
#' @param start the 0th index of the date when we begin evaluating the model
#' @param c the value for c_0
#' @param shortage_factor the amount by which we multiply shortage over waste (ratio)
#' @param buffer a value to be added to account for training in the initial stages, default 10
#' @param ind the indices for cross validation used only if non-null. (Should be sorted!!)
#' @param reset_basis TRUE or FALSE for whether we reset the LP basis before each solve (for non-uniqueness)
#' @importFrom lpSolveAPI make.lp add.constraint set.bounds lp.control delete.constraint
#' @importFrom lpSolveAPI set.objfn get.variables set.rhs get.basis set.basis
#'
single_lpSolve_wasteOnly <- function(d, l1_bounds, lag_bounds, num_vars,
                           start = 10, c = 10, shortage_factor = 15,
                           buffer = 10, ind = NULL, reset_basis = FALSE) {
  doing_cv <- !is.null(ind)

  ## ind and l1_bounds must be sorted (ascending and descending, respectively)

  # Assumes that date is the first column and the output (usage) is the last column.
  XX <- d[, 2:(num_vars + 1)]
  y <- d[, (num_vars + 2)]

  N <- nrow(XX)
  p <- ncol(XX)

  if (!doing_cv) ind <- 1:N

  ## y is usage, so y[i] is known usage on day i
  ## p betas, N ws, N r1s, N r2s, N xs, N ts, N s's, 1 intercept, p bounds
  N_var <- p + 5*N + 1 + p

  # renaming 0th index for each set of decision variables for clarity
  w_idx0 <- p         # waste
  r1_idx0 <- p + N    # remaining inventory expiring in 1 day
  r2_idx0 <- p + 2*N  # remaining inventory expiring in 2 days
  x_idx0 <- p + 3*N   # collection
  t_idx0 <- p + 4*N   # predicted 3 day usage (i+1, i+2, i+3)
  bb_idx0 <- p + 5*N + 1 # absolute value bounds on beta (coefficients) - we skip the intercept

  # Altered to penalize r1
  pind_w  <- w_idx0 + ind
  pind_r1 <- r1_idx0 + ind
  pind_r2 <- r2_idx0 + ind

  my.lp <- lpSolveAPI::make.lp(0,N_var)
  lp.control(my.lp,
             timeout=60) # Should be included in config

  obj_coefficients <- rep(0,N_var)

  ## Set coefs of w to be 1s and those of and s to be 5s (penalty factor) (eq 9)
  obj_coefficients[pind_w] <- 1
  obj_coefficients[pind_r1] <- 0.001 # These add stability / uniqueness to the solution (break ties)
  obj_coefficients[pind_r2] <- -0.001
  lpSolveAPI::set.objfn(my.lp,obj_coefficients)

  ## Constrain sum of day of week coefs of beta to be zero
  ## to prevent non-identifiability and improve interpretability
  constraint_coefficients <- rep(0,N_var)
  constraint_coefficients[1:7] <- 1
  lpSolveAPI::add.constraint(my.lp, constraint_coefficients, type = "=", 0)

  ## w >=0, s >= 0, r1 >= 0, r2 >= 0, x >= 0, t >= 0
  for (i in 1:N){
    ## w_{i} >= 0
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[w_idx0+i] <- 1
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, type = ">=", 0)

    ## r1_{i} >= 0
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[r1_idx0+i] <- 1 # r_i(1) >= 0
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    ## r2_{i} >= 0 (further constrained >= c later)
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[r2_idx0+i] <- 1 ## r_i(2) >= 0
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    ## x_{i} >= 0
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[x_idx0+i] <- 1 ## x_i >= 0
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    # t_{i} >= 0
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[t_idx0+i] <- 1 ## t_i >= 0
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

  }

  ## The initial values for x_{start+1}, x_{start+2}, x_{start+3} ensure
  ## the feasibility of this problem.
  rhs_offset <- c(c + buffer, 0, 0)
  for (i in seq_len(3L)) {
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[x_idx0+start+i] <- 1 ## x_{start + i}
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", y[start+i] + rhs_offset[i])
  }

  ## Set r1 to zero for start day
  constraint_coefficients <- rep(0,N_var)
  constraint_coefficients[r1_idx0+start] <- 1 ## r_{start}(1)
  lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=",0)

  ## Setting r2 to zero for start day
  constraint_coefficients <- rep(0,N_var)
  constraint_coefficients[r2_idx0+start] <- 1 ## r_{start}(2)
  lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=",0)

  ## Constraints on three-day prediction (according to a linear model)
  ## t_i = XX*beta + intercept [Eqn 11]
  ## Note that t_start is completely determined by our initial collection and inventory values
  for(i in (start+1):N) {
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[t_idx0+i] <- 1 ## t_i
    constraint_coefficients[1:p] <- -as.numeric(XX[i,]) ## p features on day i
    constraint_coefficients[p+5*N+1] <- -1 ## negate intercept
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", 0)
  }

  ## Constraints on remaining inventories, collection, and waste
  ## r1_{start} and r2_{start} are already fixed
  for (i in (start+1):N){

    ## Inequalities constraining r1_{i}
    ## r1_{i} - r2_{i-1} - r1_{i-1} + w_{i} >= -y_{i} (Eqn 14)
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[r1_idx0+i] <- 1 #r_i(1)
    constraint_coefficients[w_idx0+i] <- 1 #w_i
    constraint_coefficients[r1_idx0+i-1] <- -1 #r_{i-1}(1)
    constraint_coefficients[r2_idx0+i-1] <- -1 # r_{i-1}(2)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", -y[i])

    ## r1_{i} - r2_{i-1} <= 0
    constraint_coefficients <- rep(0, N_var)
    constraint_coefficients[r1_idx0+i] <- 1 #r_i(1)
    constraint_coefficients[r2_idx0+i-1] <- -1  #r_{i-1}(2)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", 0)

    ## Inequalities constraining r2_{i} [1]

    ## r2_{i} - x_{i} <= 0
    constraint_coefficients <- rep(0, N_var)
    constraint_coefficients[r2_idx0+i] <- 1
    constraint_coefficients[x_idx0+i] <- -1
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", 0)

    ## Inequalities constraining waste [2]
    ##w_{i} >= r1_{i-1} - y_{i} (Eqn 12)
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[w_idx0+i] <- 1 ## w_i
    constraint_coefficients[r1_idx0+i-1] <- -1 ## r_{i-1}(1)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", -y[i])

    ## w_i <= r1_{i-1}
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", 0)

    ## Primary Equation [Eq 15] - This makes the problem solvable
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[r1_idx0+i] <- 1 ## r_i(1)
    constraint_coefficients[r2_idx0+i] <- 1 ## r_i(2)
    constraint_coefficients[w_idx0+i] <- 1 ## w_i
    constraint_coefficients[x_idx0+i] <- -1 ## x_i
    constraint_coefficients[r1_idx0+i-1] <- -1 ## r_{i-1}(1)
    constraint_coefficients[r2_idx0+i-1] <- -1 ## r_{i-1}(2)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", -y[i])

  }

  ## Constraint all r2 after start to be greater than c (Equation 17)
  for (i in ind) {
    if(i >= start+4){
      constraint_coefficients <- rep(0,N_var)
      constraint_coefficients[r2_idx0+i] = 1 ## r_i(2)
      lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", c)
    }
  }

  ## collection (note that this is an optimistic strategy, anticipating no waste)
  ## x_{i+3}+x_{i+2}+x_{i+1} + r2_{i} + r1_{i} -t_i== 0:start+N
  ## We have x_{start + 1}, x_{start + 2}, x_{start + 3}, r1_{start}, and r2_{start} fixed,
  ## so we begin evaluating from start + 4
  for (i in seq.int(start+4, N)) {

    ## Implements equation 8
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[x_idx0+i] <- 1 ## x_{i}
    if (i %in% ind) {
      constraint_coefficients[x_idx0+i-1] <- 1 ## x_{i-1}
      constraint_coefficients[x_idx0+i-2] <- 1 ## x_{i-2}
      constraint_coefficients[r1_idx0+i-3] <- 1 ## r_{i-3}(1)
      constraint_coefficients[r2_idx0+i-3] <- 1 ## r_{i-3}(2)
      constraint_coefficients[t_idx0+i-3] <- -1 ## t_{i-3}

      lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", 0)
    } else {
      lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", y[i])
    }
  }

  ## Coefficient Constraints
  for (j in 1:p) {
    ## constrain bounds (l1 constraints) KO added this
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[bb_idx0+j] <- 1 ## l1 constraint
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    ## penalty term constraints (constrain below 0)
    constraint_coefficients[j] <- 1 ## beta_j
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    ## penalty term constraints (constrain above 0)
    constraint_coefficients[j] <- -1  ## beta_j
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)
  }

  lpSolveAPI::set.bounds(my.lp, lower = rep(-1000,ncol(my.lp)), upper = rep(1000,ncol(my.lp)))
  lpSolveAPI::lp.control(my.lp,sense='min')

  coefficients_matrix <- matrix(0, ncol = length(l1_bounds) * length(lag_bounds), nrow = p+1)

  # Create a constant basis (this is in case there are multiple optimal solutions)
  pre_rhs <- lpSolveAPI::get.rhs(my.lp) # grab current RHS

  lpSolveAPI::set.rhs(my.lp, rep(0, length(pre_rhs))) # set RHS to 0
  lpSolveAPI::set.bounds(my.lp, lower = rep(0,ncol(my.lp)), upper = rep(0,ncol(my.lp))) # set bounds to 0
  lpSolveAPI::set.objfn(my.lp, rep(1, N_var)) # reset objfn to sum of all vars (unique solution)
  status_init <- solve(my.lp)

  obj_coefficients[pind_w] <- 1
  obj_coefficients[pind_r1] <- 0.001
  obj_coefficients[pind_r2] <- -0.001
  lpSolveAPI::set.objfn(my.lp, obj_coefficients)
  pre_basis <- lpSolveAPI::get.basis(my.lp)
  lpSolveAPI::set.rhs(my.lp, pre_rhs)
  lpSolveAPI::set.bounds(my.lp, lower = rep(-1000,ncol(my.lp)), upper = rep(1000,ncol(my.lp)))

  for (i in seq_along(l1_bounds)) {

    # Reset basis
    if (reset_basis)
      lpSolveAPI::set.basis(my.lp, pre_basis)

    ## Vars 1-8 are not constrained by the l1_bound
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[(bb_idx0 + 9):(bb_idx0 + p)] <- 1
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", l1_bounds[i])

    nConstraints <- nrow(my.lp) # index refers to l1 constraint

    for (j in seq_along(lag_bounds)) {

      # Apply a common constraint to the first lags
      if (lag_bounds[j] > 0) {
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[bb_idx0 + 8] <- 1
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", lag_bounds[j])
      }

      status <- solve(my.lp)

      # If a timeout occurs, try to solve the problem with looser constraints
      if (status == 7 && !doing_cv) {
        while (status == 7) {
          warning("A timeout occurred. Loosening L1 constraints", call. = FALSE)
          lpSolveAPI::delete.constraint(my.lp, nConstraints) # delete the L1 constraint

          # relax the L1 constraints aside from day of week and seven day moving average
          constraint_coefficients <- rep(0,N_var)
          constraint_coefficients[(bb_idx0 + 9):(bb_idx0 + p)] <- 1
          lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", l1_bounds[i] + 8)
          status <- solve(my.lp)
        }
      }

      if (doing_cv && (status > 0))  {
        # If we are running cross validation, we can accept some failures
        warning(sprintf("lpSolveStatus is %d", status))
        coefficients_matrix[ , (i - 1) * length(lag_bounds) + j] <- 0.0 # set all coefs to 0

      }
      else if (status > 0) {
        stop(sprintf("lpSolveStatus is %d", status))
      }
      else {
        coeffs <- c(lpSolveAPI::get.variables(my.lp)[5*N+p+1], lpSolveAPI::get.variables(my.lp)[1:p])
        coeffs[abs(coeffs) < 1e-10] <- 0.0 # set coeffs below tolerance to 0 (may fix some numeric issues)
        coefficients_matrix[ , (i - 1) * length(lag_bounds) + j] <- coeffs  ## all the beta_j's
      }

      # If a constraint was added, clean up the lag bound
      while (nConstraints < nrow(my.lp)) {
        lpSolveAPI::delete.constraint(my.lp, nrow(my.lp)) # delete the lag bound constraint
      }

    }
    # delete the L1 bound constraint
    lpSolveAPI::delete.constraint(my.lp, nConstraints)

  }

  rm(my.lp)

  coefficients_matrix

}


