#' @importFrom stats p.adjust runif
NULL

#' #' Fast Survival Probability Extraction at Specific Times
#'
#' Computes survival probabilities at specific time points using a fitted survival model,
#' optionally interpolating if the target time points are not in the model's time grid.
#'
#' @param model A fitted survival model (must implement `$predict()` method).
#' @param data A data frame of covariates for prediction.
#' @param target_times A numeric vector of times at which to extract survival probabilities.
#' @param time_grid Optional numeric vector of times over which to predict survival curves.
#'                  Defaults to `sort(unique(target_times))`.
#' @param transform A function applied to the survival probabilities (default: `identity`).
#'
#' @return A numeric vector of survival probabilities at the specified target times.
#' @keywords internal
fast_predict_at_times <- function(model, data, target_times, time_grid = NULL, transform = identity) {
  # Use unique target times as time grid if none provided
  if (is.null(time_grid)) {
    time_grid <- sort(unique(target_times))
  }
  
  # Predict survival curves over minimal time grid
  pred <- model$predict(data, time_grid)
  surv_curves <- pred$predictions  # rows: patients, cols: time points
  
  # If all target_times are in time_grid, avoid interpolation
  if (all(target_times %in% time_grid)) {
    idx_time <- match(target_times, time_grid)
    surv_vals <- surv_curves[cbind(seq_along(target_times), idx_time)]
  } else {
    # Otherwise interpolate
    surv_vals <- mapply(function(i) {
      approx(time_grid, surv_curves[i, ], xout = target_times[i], rule = 2)$y
    }, seq_along(target_times))
  }
  
  return(transform(surv_vals))
}

#' Compute Conformal P-values
#'
#' Computes one-sided conformal p-values (right-tailed or left-tailed) for survival predictions
#' at a grid of time points, using inverse probability of censoring weights (IPCW).
#'
#' @param data.test Test data (patients to predict).
#' @param data.cal Calibration data (used for p-value computation).
#' @param surv_model A fitted survival model wrapper.
#' @param cens_model A fitted censoring model wrapper.
#' @param time_points Time grid for predictions. If `NULL`, uses equally spaced times.
#' @param num_time_points Number of time points (used if `time_points = NULL`).
#' @param alternative `"greater"` or `"smaller"` — direction of the one-sided test.
#' @param break_ties Whether to add random noise to break score ties (default: `FALSE`).
#' @param fast Whether to use vectorized approximations (recommended: `TRUE`).
#'
#' @return A matrix of p-values (rows: patients, columns: time points).
#' @keywords internal
compute_cp <- function(data.test, data.cal, surv_model, cens_model, time_points=NULL, num_time_points=100, alternative="greater",
                       break_ties=FALSE, fast=TRUE) {
  if(is.null(time_points)) {
    time_points <- seq(0, max(data.cal$time), length.out=num_time_points)
  }
  n <- nrow(data.cal)
  if(alternative=="greater") {
    scoring_fun <- function(x, t) { surv_model$predict(x, t)$predictions}
  } else if (alternative=="smaller") {
    scoring_fun <- function(x, t) { 1 - surv_model$predict(x, t)$predictions }
  } else {
    stop("Error: unknown alternative!")
  }
  ## Compute calibration scores
  scores.cal <- rep(NA, n)
  idx.event <- which(data.cal$status==1)
  if(fast) {
    if (length(idx.event) > 0) {
      new_data <- data.cal[idx.event, , drop = FALSE]
      event_times <- data.cal$time[idx.event]
      transform_fn <- switch(alternative,
                             "greater" = identity,
                             "smaller" = function(x) 1 - x,
                             stop("Unknown alternative")
      )
      scores.cal[idx.event] <- fast_predict_at_times(surv_model, new_data, event_times, transform=transform_fn)
    }
  } else {
    scores.cal[idx.event] <- sapply(idx.event, function(i) { scoring_fun(data.cal[i,], data.cal$time[i])} )
  }
  
  scores.test <- scoring_fun(data.test, time_points)
  if(break_ties) {
    scores.cal[idx.event] <- scores.cal[idx.event] + runif(length(scores.cal[idx.event]), min = -1e-6, max = 1e-6)
    scores.test <- scores.test + runif(length(scores.test), min = -1e-6, max = 1e-6)
  }
  comparison.operator <- ">="
  compare_fun <- match.fun(comparison.operator)
  compare_fun_equal <- match.fun("==")
  
  ## Compute IPC weights
  weights <- rep(NA, n)
  if (fast) {
    weights <- rep(NA, n)
    if (length(idx.event) > 0) {
      new_data <- data.cal[idx.event, , drop = FALSE]
      event_times <- data.cal$time[idx.event]            
      weights[idx.event] <- 1 / pmax(fast_predict_at_times(cens_model, new_data, event_times, transform=identity), 1e-6)
    }
  } else {
    weights[idx.event] <- sapply(idx.event, function(i) {
      prob.cens <- cens_model$predict(data.cal[i,], data.cal$time[i])$predictions
      return(1/prob.cens)
    } )
  }
  ## For numerical stability, do not allow extremely large weights
  weights[weights>n] <- n
  
  if(length(idx.event)>0) {
    den <- 1+sum(weights[idx.event])
  } else {
    den <- 1
  }
  ##print("weights")
  ##print(weights[idx.event])
  
  ## Initialize matrix for p-values
  pvals_matrix <- matrix(NA, nrow = nrow(data.test), ncol = length(time_points))
  
  ## Loop over different horizon values, but vectorized inside
  ##print(alternative)
  ##print(summary(scores.cal[idx.event]))
  for (h in seq_along(time_points)) {
    ##print(time_points[h])
    
    ## Compute score comparisons efficiently using outer()
    score_comparison <- outer(as.numeric(scores.cal[idx.event]), as.numeric(scores.test[, h]), FUN = compare_fun)
    ##score_comparison_equal <- outer(as.numeric(scores.cal[idx.event]), as.numeric(scores.test[, h]), FUN = compare_fun_equal)
    ##cat(sprintf("Number of ties for t = %f: %d\n", time_points[h], sum(score_comparison_equal)))
    ##print(summary(as.numeric(scores.test[, h])))
    
    ## Compute the numerator using matrix operations
    num_values <- 1 + colSums(weights[idx.event] * score_comparison, na.rm = TRUE)
    
    ## Compute and store p-values for the current horizon
    pvals_matrix[,h] <- pmin(1, num_values / den)
    ##print(pmin(1, num_values / den))
  }
  
  ## Enforce mononicity with respect to horizong by computing the running max
  if(FALSE) {
    if(alternative=="greater") {
      pvals_matrix <- t(apply(pvals_matrix, 1, function(row) cummax(row)))
    } else {
      pvals_matrix <- t(apply(pvals_matrix, 1, function(row) rev(cummax(rev(row)))))
    }
  }
  
  pvals_matrix <- matrix(pvals_matrix, nrow(data.test), ncol = length(time_points))
  colnames(pvals_matrix) <- time_points
  return(pvals_matrix)
}

#' Construct Conformal Survival Bands
#'
#' Computes marginal conformal prediction bands for survival probabilities at specified time points
#' using right-tailed and left-tailed p-values. Optionally applies a doubly robust correction
#' that clips the band using the fitted survival model's predictions.
#'
#' @param data.test A data frame of test data. Must include `time` and `status` columns.
#' @param data.cal A data frame of calibration data (also with `time` and `status` columns).
#' @param surv_model A fitted survival model object (e.g., created by `GRF_SurvivalForestWrapper$new()`).
#' @param cens_model A fitted censoring model object (same interface as `surv_model`, but fit to censored data).
#' @param time_points Optional numeric vector of time points at which to compute the bands. If `NULL`, a regular grid is used.
#' @param num_time_points Integer specifying number of time points to use if `time_points` is `NULL`. Default is 100.
#' @param doubly_robust Logical. If `TRUE`, the bands are extended to include the fitted model predictions (default: `TRUE`).
#' @param fast Logical. Whether to use fast, vectorized approximations for score and weight computations. Recommended. Default is `TRUE`.
#' @param use_bh Logical. Whether to use the Benjamini–Hochberg adjustment to the conformal p-values. Recommended. Default is `TRUE`.
#' @param estimate_pi0 Logical. Whether to estimate the null proportion, to increase power of BH. Default is `TRUE`.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{lower}{A matrix of lower survival bounds (patients × time points).}
#'   \item{upper}{A matrix of upper survival bounds.}
#'   \item{time_points}{The time grid used.}
#'   \item{model_pred}{Matrix of fitted survival probabilities from the black-box model.}
#' }
#'
#' @details
#' This function implements the conformal survival band (CSB) construction method described in the associated paper.
#' For each time point, it computes one-sided conformal p-values and applies Benjamini-Hochberg (BH) adjustment
#' to construct valid marginal prediction intervals. If \code{doubly_robust = TRUE}, the lower and upper bounds
#' are clipped to ensure they lie within the model-predicted curve, improving practical performance.
#'
#' @export
conformal_survival_band <- function(data.test, data.cal, surv_model, cens_model, time_points=NULL, num_time_points=100,
                                    doubly_robust=TRUE, fast=TRUE, use_bh=TRUE, estimate_pi0=FALSE) {
    n.test <- nrow(data.test)
    if(is.null(time_points)) {
        time_points <- seq(0, max(data.cal$time), length.out=num_time_points)
    }
    ## Calculate conformal pvalues
    pvals_rt <- compute_cp(data.test, data.cal, surv_model, cens_model, time_points, alternative="greater", fast=fast)
    pvals_lt <- compute_cp(data.test, data.cal, surv_model, cens_model, time_points, alternative="smaller", fast=fast)
    ## ## Define Storey's estimator
    ## estimate_pi0 <- function(pvals, lambda = 0.5) {
    ##     m <- length(pvals)
    ##     pi0_hat <- (mean(pvals > lambda)+1/m) / (1-lambda)
    ##     pi0_hat <- min(pi0_hat, 1)  # Ensure pi0 is not greater than 1
    ##     return(pi0_hat)
    ## }
    ## Calculate lower and upper bounds using BH
    if(use_bh) {
        if(estimate_pi0) {
            ##cat("NOTE: using Storey's method.\n")
            #pi0_lt_storey <- apply(pvals_lt, 2, estimate_pi0, lambda = 0.5)
            ## Estimate an upper bound for the LT null proportion (proportion alive at time t)
            pi0_lt <- sapply(time_points, function(t) {            
                ## Number of individuals with observed time >= t
                n1 <- sum(data.cal$time >= t)
                ## Number of censored individuals with censoring time <= t (may or may not be alive, count them as alive)
                n2 <- sum(data.cal$event[data.cal$time <= t] == 0)
                ## Estimate of proportion alive
                (1 + n1 + n2) / (1 + nrow(data.cal))
            })
            #print("Estimated pi0_lt:")
            #df_lt = tibble(time=time_points, storey=pi0_lt_storey, new =pi0_lt)
            #print(df_lt, n=200)
            #pi0_rt_storey <- apply(pvals_rt, 2, estimate_pi0, lambda = 0.5)
            ## Estimate an upper bound for the RT null proportion (proportion dead at time t)
            pi0_rt <- sapply(time_points, function(t) {            
                ## Number of individuals with observed time <= t (count them all as dead, even if censored)
                n1 <- sum(data.cal$time <= t)
                ## Estimate of proportion alive
                (1 + n1) / (1 + nrow(data.cal))
            })
            #print("Estimated pi0_rt:")
            #df_rt = tibble(time=time_points, storey=pi0_rt_storey, new =pi0_rt)
            #print(df_rt, n=200)
            ## Repeat each value across all rows
            pi0_lt <- matrix(rep(pi0_lt, each = nrow(pvals_lt)), nrow = nrow(pvals_lt), byrow = FALSE)
            pi0_rt <- matrix(rep(pi0_rt, each = nrow(pvals_rt)), nrow = nrow(pvals_rt), byrow = FALSE)
            upper <- pmin(1, pi0_lt*apply(pvals_lt, 2, p.adjust, method = "BH"))
            lower <- 1 - pmin(1, pi0_rt*apply(pvals_rt, 2, p.adjust, method = "BH"))
        } else {
            upper <- apply(pvals_lt, 2, p.adjust, method = "BH")
            lower <- 1 - apply(pvals_rt, 2, p.adjust, method = "BH")
        }
    } else {
        upper <- pvals_lt
        lower <- 1 - pvals_rt
    }
    upper <- matrix(upper, nrow(data.test), ncol = length(time_points))
    lower <- matrix(lower, nrow(data.test), ncol = length(time_points))
    colnames(lower) <- time_points
    colnames(upper) <- time_points
    ## Compute model predictions
    model_pred <- matrix(surv_model$predict(data.test, time_points)$predictions, nrow(data.test), ncol = length(time_points))
    colnames(model_pred) <- time_points
    if(doubly_robust) {
        ## Make results doubly robust (not used)
        lower_dr <- pmin(lower, model_pred)
        lower_dr <- matrix(lower_dr, nrow(data.test), ncol = length(time_points))
        upper_dr <- pmax(upper, model_pred)
        upper_dr <- matrix(upper_dr, nrow(data.test), ncol = length(time_points))
        colnames(lower_dr) <- time_points
        colnames(upper_dr) <- time_points
        lower <- lower_dr
        upper <- upper_dr
    }
    out <- list(lower=lower, upper=upper, time_points=time_points, model_pred=model_pred)
    return(out)
}
