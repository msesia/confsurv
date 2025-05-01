#' @title SurvivalDataGenerator: Synthetic Survival Data Generator
#'
#' @description
#' An R6 class for generating synthetic survival data using user-defined covariate generators,
#' survival time models, and censoring mechanisms. The class composes these components to simulate
#' realistic right-censored survival datasets.
#'
#' @docType class
#' @name SurvivalDataGenerator
#'
#' @field covariate_generator Function taking an integer `num_samples` and returning a matrix or data.frame.
#' @field survival An object with a `sample(X)` method to simulate event times based on covariates.
#' @field censoring An object with a `sample(X)` method to generate censoring times based on covariates.
#'
#' @return An object of class `SurvivalDataGenerator`, with a `$sample()` method that returns
#'   a data.frame including: `event_time`, `censoring_time`, `time`, `status`, and covariates `X1, X2, ..., Xp`.
#'
#' @examples
#' cov_gen <- function(n) matrix(rnorm(n * 2), nrow = n)
#' surv_gen <- ExponentialDistribution$new(function(X) rep(0.1, nrow(X)))
#' cens_gen <- ExponentialDistribution$new(function(X) rep(0.05, nrow(X)))
#' gen <- SurvivalDataGenerator$new(cov_gen, surv_gen, cens_gen)
#' data <- gen$sample(100)
#' head(data)
#'
#' @export
SurvivalDataGenerator <- R6::R6Class("SurvivalDataGenerator",
  public = list(
    covariate_generator = NULL,
    survival = NULL,
    censoring = NULL,

    #' @description Constructor for SurvivalDataGenerator.
    #' @param covariate_generator Function to generate covariates.
    #' @param survival_generator Object for generating event times.
    #' @param censoring_generator Object for generating censoring times.
    initialize = function(covariate_generator, survival_generator, censoring_generator) {
      self$covariate_generator <- covariate_generator
      self$survival <- survival_generator
      self$censoring <- censoring_generator
    },

    #' @description Generate a dataset of survival data with censoring.
    #' @param num_samples Number of observations to simulate.
    sample = function(num_samples) {
      X <- self$covariate_generator(num_samples)
      colnames(X) <- paste("X", 1:ncol(X), sep = "")
      T <- self$survival$sample(X)
      C <- self$censoring$sample(X)
      time <- pmin(T, C)
      status <- as.integer(time == T)
      data.frame(event_time = T, censoring_time = C, time = time, status = status, X)
    }
  )
)

#' @title SurvivalDistribution: Abstract Base Class for Survival Models
#'
#' @description
#' An abstract R6 base class for defining survival distributions used in synthetic data generation.
#' Subclasses must implement the core methods: \code{sample}, \code{predict}, and \code{predict_quantiles}.
#'
#' @docType class
#' @name SurvivalDistribution
#'
#' @return An abstract R6 class to be extended by concrete survival distribution implementations.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{initialize(...)}}{
#'     Constructor. Should be implemented in subclasses.
#'   }
#'   \item{\code{sample(X, T = NULL)}}{
#'     Abstract method to sample survival or censoring times based on covariates.
#'   }
#'   \item{\code{predict(X, time.points)}}{
#'     Abstract method to return survival probabilities at given time points.
#'   }
#'   \item{\code{predict_quantiles(X, probs = c(0.25, 0.5, 0.75))}}{
#'     Abstract method to return quantile survival times for given probabilities.
#'   }
#' }
#'
#' @export
SurvivalDistribution <- R6::R6Class("SurvivalDistribution",
  public = list(

    #' @description
    #' Abstract constructor. Should be implemented by subclass.
    #' @param ... Subclass-specific parameters.
    initialize = function(...) {
      stop("This method should be implemented by the derived class.")
    },

    #' @description
    #' Abstract method to sample survival times given covariates.
    #' @param X Covariate matrix or data.frame.
    #' @param T Optional vector of lower bounds (e.g., event times); returned values must exceed T.
    #' @return A numeric vector of sampled survival times.
    sample = function(X, T = NULL) {
      stop("This method should be implemented by the derived class.")
    },

    #' @description
    #' Abstract method to predict survival probabilities at specific time points.
    #' @param X Covariate matrix or data.frame.
    #' @param time.points Numeric vector of time values to evaluate survival.
    #' @return A list with two components:
    #'   \itemize{
    #'     \item \code{predictions}: Matrix of survival probabilities.
    #'     \item \code{time.points}: Vector of corresponding time points.
    #'   }
    predict = function(X, time.points) {
      stop("This method should be implemented by the derived class.")
    },

    #' @description
    #' Abstract method to compute quantiles of survival time distribution.
    #' @param X Covariate matrix or data.frame.
    #' @param probs Numeric vector of quantile probabilities.
    #' @return A data.frame of survival time quantiles per individual.
    predict_quantiles = function(X, probs = c(0.25, 0.5, 0.75)) {
      stop("This method should be implemented by the derived class.")
    }
  )
)

#' @title LogNormalDistribution: Covariate-dependent Log-normal Survival Model
#'
#' @description
#' An R6 class for modeling survival times using a log-normal distribution with covariate-dependent
#' location (\code{mu}) and scale (\code{sigma}). Supports sampling, survival prediction, and quantiles.
#'
#' @docType class
#' @name LogNormalDistribution
#'
#' @field mu_fun A function that maps covariates \code{X} to the log-mean parameter \code{mu}.
#' @field sigma_fun A function that maps covariates \code{X} to the log-scale parameter \code{sigma}.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{initialize(mu_fun, sigma_fun)}}{
#'     Constructor method.
#'     \itemize{
#'       \item \code{mu_fun}: Function returning log-mean values from covariates.
#'       \item \code{sigma_fun}: Function returning log-scale values from covariates.
#'     }
#'   }
#'
#'   \item{\code{sample(X, T = NULL, max_reps = 100)}}{
#'     Samples times from the log-normal distribution.
#'     \itemize{
#'       \item \code{X}: A matrix or data.frame of covariates.
#'       \item \code{T}: Optional lower bound vector for rejection sampling.
#'       \item \code{max_reps}: Max attempts to satisfy constraint (default 100).
#'     }
#'     Returns a numeric vector of sampled times.
#'   }
#'
#'   \item{\code{predict(X, time.points)}}{
#'     Computes survival probabilities at given times.
#'     \itemize{
#'       \item \code{X}: Covariates.
#'       \item \code{time.points}: Vector of time values.
#'     }
#'     Returns a list with:
#'     \itemize{
#'       \item \code{predictions}: Matrix of survival probabilities.
#'       \item \code{time.points}: The input time values.
#'     }
#'   }
#'
#'   \item{\code{predict_quantiles(X, probs = c(0.25, 0.5, 0.75))}}{
#'     Computes quantiles for each observation.
#'     \itemize{
#'       \item \code{X}: Covariates.
#'       \item \code{probs}: Vector of desired quantiles.
#'     }
#'     Returns a data.frame with quantile values per individual.
#'   }
#'
#'   \item{\code{fit(formula, data)}}{
#'     Placeholder method for compatibility; does nothing.
#'   }
#' }
#'
#' @examples
#' mu_fn <- function(X) 0.5 * X[,1]
#' sigma_fn <- function(X) rep(0.5, nrow(X))
#' dist <- LogNormalDistribution$new(mu_fun = mu_fn, sigma_fun = sigma_fn)
#' X <- data.frame(X1 = rnorm(5), X2 = rnorm(5))
#' dist$sample(X)
#' dist$predict(X, time.points = c(1, 2, 3))
#' dist$predict_quantiles(X, probs = c(0.25, 0.5, 0.75))
#'
#' @export
LogNormalDistribution <- R6::R6Class("LogNormalDistribution",
  inherit = SurvivalDistribution,
  public = list(
      
    mu_fun = NULL,                 ## Function to calculate mean based on covariates X
    sigma_fun = NULL,              ## Function to calculate standard deviation based on covariates X

    #' @description Initialize a new LogNormalDistribution object.
    #' @param mu_fun A function mapping covariates X to log-mean values.
    #' @param sigma_fun A function mapping covariates X to log-scale values.
    initialize = function(mu_fun, sigma_fun) {
        self$mu_fun <- mu_fun
        self$sigma_fun <- sigma_fun
    },

    #' @description Sample survival times.
    #' @param X Covariate matrix or data.frame.
    #' @param T Optional lower bounds; if provided, samples must exceed `T[i]`.
    #' @param max_reps Max rejection sampling attempts per sample (default = 100).
    #' @return A numeric vector of sampled survival times.
    sample = function(X, T=NULL, max_reps=100) {
        ## Ensure that we are working with a data frame excluding time and status if present
        if (is.data.frame(X)) {
            X <- X[, !(colnames(X) %in% c("time", "status")), drop = FALSE]
        }
        num_samples <- nrow(X)
        mu_x <- self$mu_fun(X)
        sigma_x <- self$sigma_fun(X)
        if(is.null(T)) {
            log_out <- rnorm(num_samples, mean = mu_x, sd = sigma_x)
            out <- exp(log_out)
        } else {
            ## Initialize output vector for conditional samples
            out <- rep(NA, num_samples)

            ## Loop through each individual and sample conditionally
            for (i in 1:num_samples) {
                iter <- 0
                while (iter < max_reps) {
                    log_C <- rnorm(1, mean = mu_x[i], sd = sigma_x[i])
                    C <- exp(log_C)

                    ## Check if the sampled value satisfies the condition C > `T[i]`
                    if (is.null(T) || C > T[i]) {
                        out[i] <- C
                        break  ## Exit the loop for this individual when C > `T[i]`
                    }
                    iter <- iter + 1
                }

                ## If maximum iterations exceeded, set the sample equal to `T[i]`
                if (is.na(out[i]) && !is.null(T)) {
                    out[i] <- T[i]
                    warning(paste("Max iterations exceeded for individual", i, "; sample set to `T[i]`."))
                }
            }
        }
        return(out)
    },

    #' @description Predict survival probabilities at specified time points.
    #' @param X Covariate matrix or data.frame.
    #' @param time.points Numeric vector of time points.
    #' @return A list with `predictions` (matrix) and `time.points`.
    predict = function(X, time.points) {
        ## Ensure that we are working with a data frame excluding time and status if present
        if (is.data.frame(X)) {
            X <- X[, !(colnames(X) %in% c("time", "status")), drop = FALSE]
        }

        ## Compute the mean (mu) and standard deviation (sigma) for each individual
        mu_x <- self$mu_fun(X)
        sigma_x <- self$sigma_fun(X)

        ## Initialize a matrix to store survival probabilities
        log_t <- log(time.points)
        log_t[time.points <= 0] <- -Inf  ## Handle non-positive time points

        ## Vectorized calculation of survival probabilities
        log_t_matrix <- matrix(log_t, nrow = nrow(X), ncol = length(time.points), byrow = TRUE)
        mu_matrix <- matrix(mu_x, nrow = nrow(X), ncol = length(time.points))
        sigma_matrix <- matrix(sigma_x, nrow = nrow(X), ncol = length(time.points))

        ## Compute the survival probabilities for all individuals at all time points
        z_scores <- (log_t_matrix - mu_matrix) / sigma_matrix
        predictions <- 1 - pnorm(z_scores)
        predictions[log_t_matrix == -Inf] <- 1  ## Survival probability is 1 at time 0 or negative times

        return(list(predictions = predictions, time.points = time.points))
    },

    #' @description Predict survival quantiles for each individual.
    #' @param X Covariate matrix or data.frame.
    #' @param probs Numeric vector of quantiles to compute (default: c(0.25, 0.5, 0.75)).
    #' @return A data.frame with one row per individual and one column per quantile.
    predict_quantiles = function(X, probs = c(0.25, 0.5, 0.75)) {
        ## Ensure that we are working with a data frame excluding time and status if present
        if (is.data.frame(X)) {
            X <- X[, !(colnames(X) %in% c("time", "status")), drop = FALSE]
        }

        ## Compute the mean (mu) and standard deviation (sigma) for each individual
        mu_x <- self$mu_fun(X)
        sigma_x <- self$sigma_fun(X)

        ## Vectorized calculation of quantiles
        qnorm_matrix <- matrix(qnorm(probs), nrow = nrow(X), ncol = length(probs), byrow = TRUE)
        mu_matrix <- matrix(mu_x, nrow = nrow(X), ncol = length(probs))
        sigma_matrix <- matrix(sigma_x, nrow = nrow(X), ncol = length(probs))

        ## Calculate quantiles using the inverse survival function formula
        quantiles_matrix <- exp(mu_matrix + sigma_matrix * qnorm_matrix)

        ## Convert the quantiles matrix to a data frame
        quantiles_df <- as.data.frame(quantiles_matrix)
        colnames(quantiles_df) <- paste0("Q", probs * 100, "%")
        rownames(quantiles_df) <- paste0("Individual_", 1:nrow(quantiles_df))

        return(quantiles_df)
    },

    #' @description Dummy fit method (does nothing).
    #' @param formula Not used.
    #' @param data Not used.
    fit = function(formula, data) {
    }

    )
  )

#' ExponentialDistribution: Parametric Exponential Survival Model
#'
#' This class defines a parametric survival model where times follow an exponential distribution.
#' The hazard rate (lambda) is modeled as a function of covariates through a user-provided function.
#' The class supports sampling, survival probability prediction, and quantile estimation.
#'
#' @section Inherits:
#' \code{\link{SurvivalDistribution}}: A generic base class for survival distributions.
#'
#' @section Fields:
#' \describe{
#'   \item{\code{rate_fun}}{Function that maps covariates to the exponential rate parameter (λ).}
#' }
#'
#' @section Methods:
#' \describe{
#'
#'   \item{\code{initialize(rate_fun)}}{
#'     Constructor that initializes the exponential distribution model.
#'     \itemize{
#'       \item \code{rate_fun}: A function of the covariate matrix \code{X} that returns a numeric vector of rates.
#'     }
#'   }
#'
#'   \item{\code{sample(X, T = NULL)}}{
#'     Draws samples from the exponential distribution using covariate-specific rates.
#'     If \code{T} is provided, the sampled values are conditionally shifted so that \code{C > T[i]} for each sample.
#'     \itemize{
#'       \item \code{X}: A matrix or data.frame of covariates.
#'       \item \code{T}: Optional vector of true survival times to enforce right-censoring.
#'     }
#'     Returns a numeric vector of sampled times.
#'   }
#'
#'   \item{\code{predict(X, time.points)}}{
#'     Computes survival probabilities \code{P(T > t)} for each individual at given time points.
#'     \itemize{
#'       \item \code{X}: Covariates.
#'       \item \code{time.points}: Numeric vector of time points.
#'     }
#'     Returns a list:
#'     \itemize{
#'       \item \code{predictions}: Matrix of survival probabilities (n_individuals × n_time_points).
#'       \item \code{time.points}: The vector of time points used.
#'     }
#'   }
#'
#'   \item{\code{predict_quantiles(X, probs = c(0.25, 0.5, 0.75))}}{
#'     Computes the survival quantiles (inverse CDF) for each individual.
#'     \itemize{
#'       \item \code{X}: Covariates.
#'       \item \code{probs}: Numeric vector of probabilities for quantile computation.
#'     }
#'     Returns a data.frame with quantile estimates per individual.
#'   }
#'
#'   \item{\code{fit(formula, data)}}{
#'     No-op method included for compatibility. Does not perform model fitting.
#'   }
#'
#' }
#'
#' @examples
#' # Define a rate function depending on the first covariate
#' rate_fun <- function(X) { exp(0.2 * X[, 1]) }
#'
#' # Instantiate the exponential distribution
#' exp_model <- ExponentialDistribution$new(rate_fun = rate_fun)
#'
#' # Generate some covariate data
#' X <- data.frame(X1 = rnorm(5), X2 = rnorm(5))
#'
#' # Sample times
#' sampled_times <- exp_model$sample(X)
#'
#' # Predict survival probabilities at t = 1, 2, 3
#' exp_model$predict(X, c(1, 2, 3))
#'
#' # Predict quantiles (e.g., median survival time)
#' exp_model$predict_quantiles(X, probs = 0.5)
#'
#' @export
ExponentialDistribution <- R6::R6Class("ExponentialDistribution",
  inherit = SurvivalDistribution,
  public = list(
      #' @field rate_fun A function that returns a numeric vector of rates given covariates.
      rate_fun = NULL,

    #' @description Constructor for ExponentialDistribution.
    #' @param rate_fun A function mapping covariates \code{X} to exponential rates.      
    initialize = function(rate_fun) {
        self$rate_fun <- rate_fun
    },

    #' @description Sample survival or censoring times from the exponential distribution.
    #' @param X A data.frame or matrix of covariates.
    #' @param T Optional numeric vector of lower bounds; sampled values will satisfy \code{sample > T[i]}.
    #' @return A numeric vector of sampled times.
    sample = function(X, T = NULL) {
        ## Ensure that we are working with a data frame excluding time and status if present
        if (is.data.frame(X)) {
            X <- X[, !(colnames(X) %in% c("time", "status")), drop = FALSE]
        }
        num_samples <- nrow(X)
        ## Compute censoring rates
        censoring_rate <- self$rate_fun(X)

        ## Sample times from an exponential distribution
        out <- rexp(num_samples, rate = censoring_rate)

        ## If conditioning values T are provided, shift samples so that C > T
        if (!is.null(T)) {
            out <- out + T
        }

        return(out)
    },

    #' @description Compute survival probabilities \code{P(T > t)} for each row of X at each time point.
    #' @param X Covariate matrix or data.frame.
    #' @param time.points Vector of time values.
    #' @return A list with:
    #' \itemize{
    #'   \item \code{predictions}: Matrix of survival probabilities.
    #'   \item \code{time.points}: The time grid used.
    #' }
    predict = function(X, time.points) {
        ## Ensure that X is a matrix (this should be the case for survival models)
        if (is.data.frame(X)) {
            X <- X[, !(colnames(X) %in% c("time", "status")), drop = FALSE]
        }

        ## Compute censoring rates using the rate function
        censoring_rate <- self$rate_fun(X)  ## Vector of length n_individuals

        ## Compute survival probabilities in a vectorized way
        ## pexp applied element-wise to time.points across censoring rates
        ## Expand censoring_rate to match the dimensions of time.points
        censoring_matrix <- matrix(censoring_rate, nrow = length(censoring_rate), ncol = length(time.points), byrow = FALSE)
        time_matrix <- matrix(time.points, nrow = length(censoring_rate), ncol = length(time.points), byrow = TRUE)

        predictions <- pexp(time_matrix, rate = censoring_matrix, lower.tail = FALSE)

        ## Return predictions and time points
        return(list(predictions = predictions, time.points = time.points))
    },


    #' @description Estimate survival time quantiles using the inverse CDF.
    #' @param X Covariate matrix or data.frame.
    #' @param probs Vector of survival probabilities (e.g., 0.5 for median).
    #' @return A data.frame of quantile estimates (rows = individuals, cols = quantiles).
    predict_quantiles = function(X, probs = c(0.25, 0.5, 0.75)) {
        ## Ensure that X is a matrix (this should be the case for survival models)
        if (is.data.frame(X)) {
            X <- X[, !(colnames(X) %in% c("time", "status")), drop = FALSE]
        }

        ## Compute the censoring rates (λ) for each individual based on covariates
        censoring_rate <- self$rate_fun(X)  ## λ: a vector of length n_individuals

        ## Compute the exponential quantile function for each probability
        log_term <- -log(1 - probs)  ## Vectorized for probs: a vector of length n_probs

        ## Calculate quantiles for all individuals and all probabilities in a vectorized manner
        quantiles_matrix <- outer(1 / censoring_rate, log_term, "*")  ## n_individuals x n_probs matrix

        ## Convert the quantiles matrix to a data frame for easier handling
        quantiles_df <- as.data.frame(quantiles_matrix)

        ## Set appropriate column names for quantiles
        colnames(quantiles_df) <- paste0("Q", probs * 100, "%")

        ## Set appropriate row names for individuals
        rownames(quantiles_df) <- paste0("Individual_", 1:nrow(quantiles_df))

        return(quantiles_df)
    },

 
    #' @description Dummy method for API compatibility. Does not do anything.
    #' @param formula Not used.
    #' @param data Not used.
    fit = function(formula, data) {
    }

    )
  )
