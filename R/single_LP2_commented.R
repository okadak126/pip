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


#' Solve an LP problem
#' @param d the data frame
#' @param l1_bounds the upper bound on the l1 norm of coefficient vector (aside from lag and days of week)
#' @param lag_bounds the upper bound on the seven day moving average of usage coefficient (lag)
#' @param num_vars the number of features/covariates in the dataset.
#' @param start the start date
#' @param c the value for c_0
#' @param buffer a value to be added to account for training in the initial stages, default 10
#' @param ind the indices for cross validation used only if non-null. (Should be sorted!!)
#' @importFrom lpSolveAPI make.lp add.constraint set.bounds lp.control delete.constraint
#' @importFrom lpSolveAPI set.objfn get.variables
#' @export
#'
single_lpSolve <- function(d, l1_bounds, lag_bounds, num_vars, start = 10, c = 30, buffer = 10, ind = NULL) {
  doing_cv <- !is.null(ind)

  ## ind and l1_bounds must be sorted (ascending and descending, respectively)

  # Assumes that date is the first column and the output (usage) is the last column.
  XX <- d[, 2:(num_vars + 1)]
  y <- d[, (num_vars + 2)]

  N <- nrow(XX)
  p <- ncol(XX)

  if (!doing_cv) ind <- 1:N

  # Altered to penalize r1
  pind_r1 <- (p + N) + ind # indices of waste for CV fold (minimize sum of these in obj func)
  pind_w  <- p + ind

  ## y is usage, so y[i] is known usage on day i
  ## p betas, N ws, N r1s, N r2s, N xs, N ts, 1 intercept, p bounds
  N_var <- p + 5*N + 1 + p

  my.lp <- lpSolveAPI::make.lp(0,N_var)
  lp.control(my.lp,
             timeout=45) # Should be included in config

  obj_coefficients <- rep(0,N_var)

  ## Set coefs of w and r1 to be 1's (eq 9)
  obj_coefficients[pind_r1] <- 1
  obj_coefficients[pind_w] <- 40
  lpSolveAPI::set.objfn(my.lp,obj_coefficients)

  ## Constrain sum of day of week coefs of beta to be zero??  Why?
  ## Reply: There is non-identifiability issue with the day of week
  ## information. We can add the same value its all coefficients of
  ## these seven features, and shift the intercept correspondingly,
  ## the final model will be equivalent. Wihout loss of generality,
  ## we constrain the week sum of coefs to be 0, it leads to better
  ## intepretability and solves the identifiability issue. That's
  ## why I do not need to remove the first feature. If you have
  ## removed the first feature, then this constraint is no longer
  ## needed. However, later, when we apply L1 penalty, I have left
  ## out the day of week features(7)+average PLT consumption
  ## feature(1)(because we know that then are important and a linear
  ## model with them will not lead to overfitting), we need to
  ## change that accordingly as a result. I have removed the first
  ## two lines.
  constraint_coefficients <- rep(0,N_var)
  constraint_coefficients[1:7] <- 1
  lpSolveAPI::add.constraint(my.lp, constraint_coefficients, type = "=", 0)

  ##w >=0
  for (i in 1:N){
    ## Setting each coeff of w to be non-negative: (eq 13)
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[i+p] <- 1
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, type = ">=", 0)
  }

  ##start
  ## model evaluation begins on day start + 1 but collection
  ## begins start + 3
  ##
  rhs_offset <- c(c + buffer, 0, 0)
  for (i in seq_len(3L)) {
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[start+i+p+3*N] <- 1 ## x_{start + 1}
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", y[start+i] + rhs_offset[i])
  }

  ## Set r1 to zero for start day
  constraint_coefficients <- rep(0,N_var)
  constraint_coefficients[start+p+N] <- 1 ## r_{start}(1)
  lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=",0)

  ## Setting r2 to zero for start day
  ## The initial values for x_{start+1}, x_{start+2}, x_{start+3} have made sure
  ## the feasibility of this problem. The model training part here
  ## does not involve any initial information we input as a
  ## practitioner, including remaining bags, new bags will be arriving
  ## for the next three days.
  constraint_coefficients <- rep(0,N_var)
  constraint_coefficients[start+p+N*2] <- 1 ## r_{start}(2)
  lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=",0)

  ##prediction
  ## Eqn 11
  ## Total need t_i. Z is the covariate matrix (p columns)
  for(i in (start+1):N) {
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[i+p+4*N] <- 1 ## t_i
    constraint_coefficients[1:p] <- -as.numeric(XX[i,]) ## p features on day i
    constraint_coefficients[p+1+5*N] <- -1 ## negate intercept
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", 0)
  }

  for (i in (start+1):N){

    ##remaining
    ##r1_{i} - r2_{i-1}- r_{i-1} + w_{i} >= -y_{i}
    ##w_{i} + r2_{i} + r1_{i}- r1_{i-1}-r2_{i-1} - x_{i} == -y_{i}
    ## r2_{i} - x_{i} <= 0
    # r1_{i} - r2_{i-1} <= 0 (KO added this)

    ## Eqn 14
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[i+p+N] <- 1 #r_i(1)
    constraint_coefficients[i+p] <- 1 #w_i
    constraint_coefficients[i+p+N-1] <- -1 #r_{i-1}(1)
    constraint_coefficients[i+p+2*N-1] <- -1 # r_{i-1}(2)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", -y[i])

    ## Eqn 16
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[i+p+N] <- 1 # r_i(1) >= 0
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    ## KO Added this
    constraint_coefficients <- rep(0, N_var)
    constraint_coefficients[i+p+N] <- 1 #r_i(1)
    constraint_coefficients[i+p+2*N-1] <- -1  #r_{i-1}(2)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", 0)

    ## Eqn 15 (KO corrected)
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[i+p+N] <- 1 ## r_i(1)
    constraint_coefficients[i+p+N*2] <- 1 ## r_i(2)
    constraint_coefficients[i+p] = 1 ## w_i
    constraint_coefficients[i+p+3*N] <- -1 ## x_i
    constraint_coefficients[i+p+N-1] <- -1 ## r_{i-1}(1)
    constraint_coefficients[i+p+N*2-1] <- -1 ## r_{i-1}(2)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", -y[i])

    ## Eqn 17. Why is it missing c_0? Maybe because of start day??
    ## Reply: Yes! Later from day start+4->last day, we have added
    ## the c there. You have noticed it there. Here we only need
    ## to initialize start+1->start+3, but I was being lazy and
    ## had not written them in a seperate small loop.

    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[i+p+N*2] <- 1 ## r_i(2)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    ##wasted
    ##w_{i} >= r1_{i-1} - y_{i}
    ## Eqn 12
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[i+p] <- 1 ## w_i
    constraint_coefficients[i+p+N-1] <- -1 ## r_{i-1}(1)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", -y[i])
    ##constraint_coefficients[i+p+N-1] = -1 ## repeated from 2 lines above

    # w_i <= r_{i-1}(1)
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", 0)
  }

  ## (Equation 17)
  for (i in ind) {
    if(i >= start+4){
      constraint_coefficients <- rep(0,N_var)
      constraint_coefficients[i+p+N*2] = 1 ## r_i(2)
      lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", c)
    }
  }

  for (i in seq.int(start+4, N)) {

    ##collecting
    ##x_{i+3} +  x_{i+2}+x_{i+1}  + r2_{i} + r1_{i} -t_i== 0:start+N

    ## Implements equation 8 (equation 18 in the slide deck is incorrect - should be 8)
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[i+3*N+p] <- 1 ## x_{i}
    if (i %in% ind) {
      constraint_coefficients[i+3*N-1+p] <- 1 ## x_{i-1}
      constraint_coefficients[i+3*N-2+p] <- 1 ## x_{i-2}
      constraint_coefficients[i+p+N-3] <- 1 ## r_{i-3}(1)
      constraint_coefficients[i+p+N*2-3] <- 1 ## r_{i-3}(2)
      constraint_coefficients[i+p+N*4-3] <- -1 ## t_{i-3}
      lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", 0)
    } else {
      lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "=", y[i])

    }
  }

  for (j in 1:p) {
    ## constrain bounds (l1 constraints) KO added this
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[j+N*5+p+1] <- 1 ## l1 constraint
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    ## penalty term constraints (constrain below 0)
    constraint_coefficients[j] <- 1 ## beta_j
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)

    ## penalty term constraints (constrain above 0)
    constraint_coefficients[j] <- -1  ## beta_j
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, ">=", 0)
  }

  lpSolveAPI::set.bounds(my.lp, lower = rep(-1000,ncol(my.lp)), upper = rep(1000,ncol(my.lp)))
  #lpSolveAPI::set.type(my.lp, columns = (p+1):(p+5*N), type = c("integer"))
  lpSolveAPI::lp.control(my.lp,sense='min')

  coefficients_matrix <- matrix(0, ncol = length(l1_bounds) * length(lag_bounds), nrow = p+1)


  pre_rhs <- lpSolveAPI::get.rhs(my.lp) # grab current RHS
  # Obtain blank state of LP by having it solve all 0
  lpSolveAPI::set.rhs(my.lp, rep(0, length(pre_rhs))) # set RHS to 0
  lpSolveAPI::set.bounds(my.lp, lower = rep(0,ncol(my.lp)), upper = rep(0,ncol(my.lp))) # set bounds to 0
  lpSolveAPI::set.objfn(my.lp, rep(1, N_var)) # reset objfn to sum of all vars (unique solution)
  status_init <- solve(my.lp)
  #print(sprintf("Initial Solve Status = %d", status_init))

  obj_coefficients[pind_r1] <- 1
  obj_coefficients[pind_w] <- 40
  lpSolveAPI::set.objfn(my.lp, obj_coefficients)
  pre_basis <- lpSolveAPI::get.basis(my.lp)
  lpSolveAPI::set.rhs(my.lp, pre_rhs)
  lpSolveAPI::set.bounds(my.lp, lower = rep(-1000,ncol(my.lp)), upper = rep(1000,ncol(my.lp)))
  #print("Ready for evaluation")

  for (i in seq_along(l1_bounds)) {
    #print(paste0("solving for ", l1_bounds[i]))

    # Reset basis
    lpSolveAPI::set.basis(my.lp, pre_basis)

    ## Vars 1-8 are not constrained by the l1_bound
    constraint_coefficients <- rep(0,N_var)
    constraint_coefficients[(p + 5*N + 1 + 9):(p + 5*N + 1 + p)] <- 1
    lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", l1_bounds[i])

    nConstraints <- nrow(my.lp) # index refers to l1 constraint

    for (j in seq_along(lag_bounds)) {

      # Apply a constraint to the seven day moving average (lag) parameter if applicable
      if (lag_bounds[j] > 0) {
        constraint_coefficients <- rep(0,N_var)
        constraint_coefficients[p + 5*N + 1 + 8] <- 1
        lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", lag_bounds[j])
      }

      # Reset LP before it is solved
      status <- solve(my.lp)
      #print(paste0("Solution Count = ",  lpSolveAPI::get.solutioncount(my.lp)))

      # If a timeout occurs, try to solve the problem with looser constraints
      if (status == 7 && !doing_cv) {
        while (status == 7) {
          warning("A timeout occurred. Loosening L1 constraints", call. = FALSE)
          lpSolveAPI::delete.constraint(my.lp, nConstraints) # delete the L1 constraint
          # relax the L1 constraints aside from day of week and seven day moving average
          constraint_coefficients <- rep(0,N_var)
          constraint_coefficients[(p + 5*N + 1 + 9):(p + 5*N + 1 + p)] <- 1
          lpSolveAPI::add.constraint(my.lp, constraint_coefficients, "<=", l1_bounds[i] + 8)
          status <- solve(my.lp)
        }
      }

      if (doing_cv && (status == 5 || status == 7))  {
        # If we are running cross validation, we can accept some failures
        coefficients_matrix[ , (i - 1) * length(lag_bounds) + j] <- 0.0 # set all coefs to 0
      }
      else if (status > 0) {
        stop(sprintf("lpSolveStatus is %d", status))
      }
      else {
        coeffs <- c(lpSolveAPI::get.variables(my.lp)[5*N+p+1], lpSolveAPI::get.variables(my.lp)[1:p])
        coeffs[abs(coeffs) < 1e-10] <- 0.0 # set coeffs below tolerance to 0 (may fix some numeric issues)
        coefficients_matrix[ , (i - 1) * length(lag_bounds) + j] <- coeffs  ## all the beta_j's
        #print(sprintf("Value of objective (waste) is %f", lpSolveAPI::get.objective(my.lp)))
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

