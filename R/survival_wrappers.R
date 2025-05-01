#' @importFrom stats approx approxfun splinefun model.response model.frame model.matrix predict pnorm qnorm runif
#' @importFrom survival survreg survreg.control coxph basehaz
#' @importFrom randomForestSRC rfsrc
#' @importFrom grf survival_forest
#' @importFrom R6 R6Class
NULL

#' #' Initialize a Survival Model
#'
#' Instantiates a survival model wrapper object based on the selected model type.
#'
#' @param model_type A character string specifying the model type. One of:
#'   \itemize{
#'     \item \code{"grf"} — Generalized Random Forest via the \code{grf} package.
#'     \item \code{"survreg"} — Parametric survival model via \code{survreg}.
#'     \item \code{"rf"} — Random survival forest via \code{randomForestSRC}.
#'     \item \code{"cox"} — Cox proportional hazards model.
#'   }
#'
#' @return An object of class inheriting from \code{SurvivalModelWrapper}.
#' @export
init_surv_model <- function(model_type) {
  surv_model <- switch(model_type,
    "grf" = GRF_SurvivalForestWrapper$new(),
    "survreg" = SurvregModelWrapper$new(dist = "lognormal"),
    "rf" = randomForestSRC_SurvivalWrapper$new(),
    "cox" = CoxphModelWrapper$new(),
    stop("Unknown model type!")
  )
  return(surv_model)
}

#' Ensure Input is Matrix
#'
#' Utility function to safely convert vectors to matrices for consistent input handling.
#'
#' @param new_data A vector or matrix.
#' @return A matrix with one row (if input was a vector) or the original matrix.
#' @keywords internal
ensure_matrix <- function(new_data) {
  ## Check if the input is a vector
  if (is.vector(new_data)) {
    ## Convert the vector to a matrix with one row
    new_data <- matrix(new_data, nrow = 1)
  ## Check if the input is already a matrix
  } else if (is.matrix(new_data)) {
    ## If it's already a matrix, return it unchanged
    return(new_data)
  ## Raise an error if the input is neither a vector nor a matrix
  } else {
    stop("Input must be either a vector or a matrix")
  }
  ## Return the resulting matrix
  return(new_data)
}

#' Interpolate Survival Probabilities at Given Time Points
#'
#' Applies monotonic interpolation (PCHIP-style) to estimate survival probabilities
#' at arbitrary time points from a survival object produced by \code{randomForestSRC::rfsrc()}.
#'
#' @param pred A prediction object containing survival curves and time points.
#' @param time_points A numeric vector of time points at which to interpolate survival probabilities.
#'
#' @return A numeric matrix of interpolated survival probabilities.
#' @keywords internal
get_survival_prob_at_time <- function(pred, time_points) {
  ## Extract the survival times and survival probabilities from the prediction object
  times <- pred$time.interest
  survival_probs <- pred$survival

  ## Initialize a matrix to hold the interpolated probabilities
  num_individuals <- nrow(survival_probs)
  num_time_points <- length(time_points)
  interpolated_probs <- matrix(NA, nrow = num_individuals, ncol = num_time_points)

  ## Perform monotone interpolation using PCHIP for each individual and each time point
  for (i in 1:num_individuals) {
    ## Create the monotone interpolation function for the current individual
    interp_function <- splinefun(x = times, y = survival_probs[i, ], method = "monoH.FC")

    ## Interpolate at each time point and ensure the probabilities are within [0, 1]
    interpolated_probs[i, ] <- pmin(pmax(interp_function(time_points), 0), 1)
  }

  return(interpolated_probs)
}

#' @title Abstract Survival Model Wrapper
#'
#' @description
#' This R6 class provides a common interface for fitting and predicting survival models.
#' Subclasses should implement the \code{fit()} method, and optionally override \code{predict()} and \code{predict_quantiles()}.
#'
#' @field model Fitted model object.
#' @field formula Formula used for model fitting.
#' @field time.points Default failure times used for interpolation.
#' @field use_covariates Optional vector of covariate names to include.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{fit(formula, data)}}{Abstract method to fit a survival model.}
#'   \item{\code{predict(new_data, time.points = NULL)}}{Predict survival curves.}
#'   \item{\code{predict_quantiles(new_data, probs)}}{Return survival quantiles for each individual.}
#'   \item{\code{predict_interpolate(survival_probs, original_failure_times, time.points)}}{Linear interpolation utility.}
#'   \item{\code{parse_formula(formula, data)}}{Extracts time, status, and covariates from formula.}
#'   \item{\code{select_columns(data)}}{Restricts columns to specified covariates.}
#' }
#'
#' @export
SurvivalModelWrapper <- R6Class("SurvivalModelWrapper",
  public = list(
    model = NULL,                ## Holds the trained model object.
    formula = NULL,              ## Stores the formula used to fit the model.
    time.points = NULL,          ## A sequence of time points for which survival probabilities are calculated.
    use_covariates = NULL,       ## List of relevant covariates (e.g., c("X1", "X3")) to use when fitting the censoring model


    #' @description Constructor for SurvivalModelWrapper.
    #' @param use_covariates Optional character vector of covariates to be used.
    initialize = function(use_covariates = NULL) {
      self$use_covariates <- use_covariates
    },

    #' @description Abstract method for fitting a survival model. Should be implemented in subclass.
    #' @param formula A survival formula (e.g., Surv(time, status) ~ predictors).
    #' @param data A data.frame with columns for time, status, and covariates.
    fit = function(formula, data) {
      stop("This method should be implemented in the subclass.")
    },


    ## ' @description Predicts survival probabilities for individuals in new data.
    #' @param new_data A data.frame of new individuals with covariates.
    #' @param time.points Optional vector of time points. Defaults to model's time.points.
    #' @return A list with `predictions` (matrix of survival probabilities) and `time.points`.
    predict = function(new_data, time.points = NULL) {
      ## If time.points is not provided, use the default values from the model
      if (is.null(time.points)) {
        time.points <- self$time.points
      }

      ## Ensure time.points is provided either by the user or from the model
      if (is.null(time.points)) {
        stop("Error: time.points are not provided, and there is no default value set.")
      }

      ## Step 1: Predict quantiles for each probability in `probs`
      probs <- seq(0.01, 0.99, by = 0.01)
      survival_times <- self$predict_quantiles(new_data, probs = probs)

      ## Step 2: Initialize a matrix to hold survival probabilities
      survival_probs <- matrix(NA, nrow = nrow(survival_times), ncol = length(time.points))

      ## Step 3: Loop through each individual and interpolate the survival probabilities
      for (i in 1:nrow(survival_times)) {
        if(!all(is.na(survival_times[i, ]))) {
          ## Create an interpolation function for each individual
          interp_fun <- approxfun(rev(as.numeric(survival_times[i, ])), rev(probs), rule = 2, method = "linear", ties = "ordered")
          ## Apply the interpolation function to the provided time.points
          survival_probs[i, ] <- interp_fun(time.points)
        } else {
            survival_probs[i, ] <- 1
        }
      }

      ## Step 4: Return the predicted survival curves and failure times
      list(predictions = survival_probs, time.points = time.points)
    },

    #' @description Predicts survival quantiles (e.g. median, quartiles).
    #' @param new_data New data with covariates.
    #' @param probs Numeric vector of quantile probabilities. Defaults to c(0.25, 0.5, 0.75).
    #' @return A data.frame of predicted quantiles for each individual.
    predict_quantiles = function(new_data, probs = c(0.25, 0.5, 0.75)) {
        ## Predict survival curves
        predictions <- self$predict(new_data)
        survival_curves <- predictions$predictions
        time_points <- predictions$time.points  ## Time points associated with the survival curves

        ## Add padding to ensure interpolation will work
        survival_curves <- cbind(1,survival_curves,0)
        time_points <- c(0,time_points,max(time_points)+1)

        ## Function to find the survival time corresponding to a given survival percentile using built-in interpolation
        find_quantile <- function(survival_probs, time_points, percentile) {
            target_prob <- 1 - percentile  ## Convert percentile to survival probability threshold
            ## Use linear interpolation
            interpolated_time <- approx(x = rev(survival_probs), y = rev(time_points), xout = target_prob, rule = 2, ties="ordered")$y
            return(interpolated_time)
        }

        ## Initialize a list to store quantiles for each individual
        quantiles_list <- list()

        ## Loop over each individual
        for (i in 1:nrow(survival_curves)) {
            ## For each individual, find the survival times at the specified percentiles
            quantiles <- sapply(probs, function(p) find_quantile(survival_curves[i, ], time_points, p))
            quantiles_list[[i]] <- quantiles
        }

        ## Convert the list of quantiles to a data frame
        quantiles_df <- do.call(rbind, quantiles_list)
        colnames(quantiles_df) <- paste0("Q", probs * 100, "%")
        rownames(quantiles_df) <- paste0("Individual_", 1:nrow(quantiles_df))

        return(as.data.frame(quantiles_df))
    },

    #' @description Linearly interpolates survival probabilities at new time points.
    #' @param survival_probs Matrix of survival probabilities.
    #' @param original_failure_times Original time points corresponding to survival_probs.
    #' @param time.points New time points to interpolate to.
    #' @return Matrix of interpolated probabilities.
    predict_interpolate = function(survival_probs, original_failure_times, time.points) {
      ## Initialize a matrix to store interpolated survival probabilities
      survival_probs_interp <- matrix(NA, nrow = nrow(survival_probs), ncol = length(time.points))

      ## Loop through each individual and interpolate the survival probabilities at the custom times
      for (i in 1:nrow(survival_probs)) {
        ## Create an interpolation function for each individual
        interp_fun <- approxfun(original_failure_times, survival_probs[i, ], rule = 2, ties = "ordered")

        ## Apply the interpolation function to the custom failure times
        survival_probs_interp[i, ] <- interp_fun(time.points)
      }

      ## Return the interpolated survival probabilities
      return(survival_probs_interp)
    },

    #' @description Parses a survival formula and extracts components.
    #' @param formula A survival formula, e.g., Surv(time, status) ~ X1 + X2.
    #' @param data A data.frame containing survival and covariate columns.
    #' @return A list with `time`, `status`, and `covariates`.
    parse_formula = function(formula, data) {
      self$formula <- formula
      response <- model.response(model.frame(formula, data))
      time <- response[, 1]
      status <- response[, 2]
      covariates <- model.matrix(formula, data)[, -1, drop = FALSE]  ## Remove intercept
      list(time = time, status = status, covariates = covariates)
    },

    #' @description Subsets the data to only include selected covariates.
    #' @param new_data A data.frame to filter.
    #' @return A filtered data.frame with only the relevant covariates.
    select_columns = function(new_data) {
      if(!is.null(self$use_covariates)) {
          new_data_sel <- new_data %>% select(time, status, self$use_covariates)
      } else {
          new_data_sel <- new_data
      }
      return(new_data_sel)
  }

  ),

)

#' @title GRF Survival Forest Model
#'
#' @description
#' Implements a wrapper for generalized random forest survival models using the \pkg{grf} package.
#'
#' @export
GRF_SurvivalForestWrapper <- R6Class("GRF_SurvivalForestWrapper",
  inherit = SurvivalModelWrapper,
  public = list(
 
    #' @description
    #' Fit a generalized random forest survival model using the \pkg{grf} package.
    #'
    #' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
    #' @param data A data.frame containing columns for time, status, and covariates.
    fit = function(formula, data) {
      data <- self$select_columns(data)
      parsed_data <- self$parse_formula(formula, data)

      ## Fit the survival forest model
      self$model <- grf::survival_forest(parsed_data$covariates,
                                         Y = parsed_data$time,
                                         D = parsed_data$status,
                                         num.trees = 100)
      self$time.points <- self$model$failure.times  ## Extract the default failure times
      if(length(self$time.points)==0) {
          self$time.points = seq(min(data$time), max(data$time), length.out=100)
      }
    },

    #' @description
    #' Predicts survival probabilities at specified time points for new data.
    #'
    #' @param new_data A data.frame of new observations with the same structure as the training data.
    #' @param time.points Optional numeric vector of time points at which to compute predictions.
    #'
    #' @return A list with:
    #' \describe{
    #'   \item{\code{predictions}}{Matrix of survival probabilities (rows = individuals, cols = time points).}
    #'   \item{\code{time.points}}{Vector of time points used in prediction.}
    #' }
    predict = function(new_data, time.points = NULL) {
      new_data <- self$select_columns(new_data)

      ## Generate the design matrix from the new data using the stored formula
      covariates_new <- model.matrix(self$formula, new_data)[, -1, drop = FALSE]

      ## Predict survival curves
      predictions <- predict(self$model, newdata = covariates_new)

      ## Use default failure times if custom ones are not provided
      original_failure_times <- self$time.points
      survival_probs <- predictions$predictions

      ## If predictions fail, replace survival_probs with ones
      if(length(survival_probs)==0){
          survival_probs = matrix(1, nrow = nrow(new_data), ncol = length(self$time.points))
      }

      ## If custom failure times are provided, use the default interpolation method
      if (!is.null(time.points)) {
        survival_probs_interp <- self$predict_interpolate(survival_probs, original_failure_times, time.points)
        return(list(predictions = survival_probs_interp, time.points = time.points))
      }

      ## If no custom failure times are provided, return the original predictions
      return(list(predictions = survival_probs, time.points = original_failure_times))
    }
  )
)

#' @title RandomForestSRC Survival Wrapper
#'
#' @description
#' Wrapper for the \code{randomForestSRC::rfsrc()} survival model.
#'
#' @export
randomForestSRC_SurvivalWrapper <- R6Class("randomForestSRC_SurvivalWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    #' @description
    #' Fits a survival forest model using the \pkg{randomForestSRC} package.
    #' This method parses the survival formula and fits a random survival forest to the data.
    #'
    #' @param formula A survival formula, typically \code{Surv(time, status) ~ predictors}.
    #' @param data A data.frame containing time, status, and covariates.
    #' @param ntree Number of trees to grow in the forest (default = 100).
    #' @param ... Additional arguments passed to \code{randomForestSRC::rfsrc()}.
    fit = function(formula, data, ntree = 100, ...) {
      data <- self$select_columns(data)

      ## Extract the model frame based on the formula and data
      mf <- model.frame(formula, data)

      ## Ensure the response is a survival object (Surv)
      response <- model.response(mf)
      if (!inherits(response, "Surv")) {
        stop("The left-hand side of the formula must be a survival object (Surv(time, status)).")
      }
      ## Fit the survival forest model using randomForestSRC
      self$model <- randomForestSRC::rfsrc(formula, data = data, ntree = ntree, save.memory=TRUE, ...)

      ## Store the failure times from the model
      self$time.points <- self$model$time.interest

      if(length(self$time.points)==0) {
          self$time.points = seq(min(data$time), max(data$time), length.out=100)
      }
    },

    #' @description
    #' Predicts survival probabilities at specified time points for new data.
    #'
    #' @param new_data A data.frame of new observations with the same structure as the training data.
    #' @param time.points Optional numeric vector of time points at which to compute predictions.
    #'
    #' @return A list with:
    #' \describe{
    #'   \item{\code{predictions}}{Matrix of survival probabilities (rows = individuals, cols = time points).}
    #'   \item{\code{time.points}}{Vector of time points used in prediction.}
    #' }
    predict = function(new_data, time.points = NULL) {
      new_data <- self$select_columns(new_data)
      ## Ensure that new_data is correctly formatted
      if (!is.data.frame(new_data)) {
        stop("new_data must be a data frame.")
      }

      ## Use the model's failure times if custom failure times are not provided
      if (is.null(time.points)) {
        time.points <- self$time.points
      }

      ## Predict survival curves using the trained random forest model
      predictions <- predict(self$model, newdata = new_data)
      ## Extract survival probabilities for each individual at each failure time
      survival_probs <- get_survival_prob_at_time(predictions, time.points)

      ## If predictions fail, replace survival_probs with ones
      if(length(survival_probs)==0){
          survival_probs = matrix(1, nrow = nrow(new_data), ncol = length(self$time.points))
      }

        ## Return predictions and failure times
        list(predictions = survival_probs, time.points = time.points)
    }
  )
)

#' @title Parametric Survival Model Wrapper
#'
#' @description
#' Wrapper for parametric survival models via \code{survival::survreg()}.
#'
#' @field dist A character string indicating the distribution used by the model (e.g., "weibull", "lognormal").
#'
#' @export
SurvregModelWrapper <- R6Class("SurvregModelWrapper",
  inherit = SurvivalModelWrapper,
  public = list(
    dist = NULL,  ## Distribution parameter

    #' @description Constructor for SurvregModelWrapper.
    #' @param use_covariates Optional character vector of covariate names to use.
    #' @param dist Distribution to use in the survreg model (e.g., \"weibull\", \"lognormal\").
    initialize = function(use_covariates = NULL, dist = "weibull") {
      self$use_covariates <- use_covariates
      self$dist <- dist
    },
 
    #' @description Fit a parametric survival model using survreg.
    #' @param formula A formula of the form \code{Surv(time, status) ~ predictors}.
    #' @param data A data.frame with survival outcome and covariates.
    fit = function(formula, data) {
      data <- self$select_columns(data)
      parsed_data <- self$parse_formula(formula, data)
      self$model <- survival::survreg(formula, data = data, dist = self$dist, control = survival::survreg.control(maxiter = 1000))
      self$time.points <- seq(min(parsed_data$time), max(parsed_data$time), length.out = 100)
    },

    #' @description Predict survival time quantiles.
    #' @param new_data A data.frame of covariates.
    #' @param probs A numeric vector of quantile probabilities.
    #' @return A data.frame of quantile predictions.
    predict_quantiles = function(new_data, probs = c(0.25, 0.5, 0.75)) {

      ## Predict quantiles for each probability in probs
      quantiles_matrix <- sapply(probs, function(p) {
        predict(self$model, newdata = new_data, type = "quantile", p = p)
      })

      ## Ensure quantiles_matrix is a matrix
      quantiles_matrix <- ensure_matrix(quantiles_matrix)

      ## Set column names to represent the quantiles
      colnames(quantiles_matrix) <- paste0("Q", probs * 100, "%")
      rownames(quantiles_matrix) <- paste0("Individual_", 1:nrow(quantiles_matrix))

      ## Return the quantile estimates as a data frame
      as.data.frame(quantiles_matrix)
    }
  )
)

#' @title Cox Proportional Hazards Wrapper
#'
#' @description
#' Wrapper for Cox proportional hazards models via \code{survival::coxph()}.
#'
#' @field model The fitted Cox model object from \code{survival::coxph()}.
#' @field formula_env Environment where the formula is evaluated (optional).
#'
#' @export
CoxphModelWrapper <- R6Class("CoxphModelWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    model = NULL,             ## To store the fitted coxph model
    formula_env = NULL,        ## To store the formula environment

    #' @description Fit a Cox proportional hazards model using coxph().
    #' @param formula A formula like \code{Surv(time, status) ~ predictors}.
    #' @param data A data.frame with survival times, censoring indicator, and predictors.
    fit = function(formula, data) {
        data <- self$select_columns(data)
        ## Fit the coxph model and store it
        parsed_data <- self$parse_formula(formula, data)
        self$model <- survival::coxph(formula, data = data, x = TRUE, y = TRUE)
        self$time.points <- seq(min(parsed_data$time), max(parsed_data$time), length.out = 100)

        ## Capture and store the formula's environment
        #self$formula_env <- environment(formula)

        ## Ensure the environment is set correctly for terms evaluation
        #environment(self$model$terms) <- self$formula_env
    },
 
    #' @description Predict survival probabilities for new individuals using Cox model.
    #' @param new_data A data.frame with covariates.
    #' @param time.points Optional numeric vector of time points to predict at.
    #' @return A list with `predictions` matrix and `time.points`.
    predict = function(new_data, time.points = NULL) {
        new_data <- self$select_columns(new_data)

        ## Use the model's failure times if custom failure times are not provided
        if (is.null(time.points)) {
            time.points <- self$time.points
        }


        ## Step 1: Calculate linear predictors for the new data
        linear_predictors <- predict(self$model, newdata = new_data, type = "lp")

        ## Step 2: Extract the baseline cumulative hazard
        baseline_hazard <- survival::basehaz(self$model, centered = TRUE)


        ## Step 3: Interpolate baseline hazard to match the requested time.points
        extended_hazard_times <- c(0, baseline_hazard$time, max(baseline_hazard$time) * 10)
        extended_hazard_probs <- c(0, baseline_hazard$hazard, 1)
        baseline_hazard_at_times <- approx(extended_hazard_times, extended_hazard_probs, xout = time.points)$y

        ## Step 4: Compute survival probabilities for all individuals
        ## Survival probability is calculated as: exp(-baseline hazard * exp(linear predictor))
        survival_probs <- t(exp(-outer(baseline_hazard_at_times, exp(linear_predictors), "*")))

        return(list(predictions = survival_probs, time.points = time.points))
    }


  )
)
