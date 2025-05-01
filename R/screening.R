#' Select Patients Based on Survival Bands
#'
#' Applies a screening rule to select patients whose conformal survival bands meet a specified risk threshold.
#'
#' @param time.points A numeric vector of time points corresponding to the survival curves.
#' @param survival_lb A matrix of lower survival bounds (patients × time points).
#' @param survival_ub A matrix of upper survival bounds (patients × time points).
#' @param screening_time A scalar time at which to apply the screening rule.
#' @param screening_prob A scalar survival probability threshold for screening.
#' @param screening_crit A string indicating the screening criterion. One of `"low risk"` (select if survival > threshold)
#'   or `"high risk"` (select if survival < threshold).
#'
#' @return A list with:
#' \describe{
#'   \item{time}{The screening time used.}
#'   \item{selected}{An integer vector of patient indices who meet the screening criterion.}
#' }
#'
#' @details
#' This function selects patients whose **conformal survival band** at a given time satisfies a risk threshold.
#' Specifically, it uses the lower bound of the band if screening for "low risk", and the upper bound if screening for "high risk".
#'
#' Internally, it delegates to [select_patients_point()] for interpolation and comparison.
#'
#' @seealso [select_patients_point()], [conformal_survival_band()]
#' 
#' @export
select_patients_band <- function(time.points, survival_lb, survival_ub, screening_time, screening_prob, screening_crit) {
    stopifnot(screening_crit %in% c("low risk", "high risk"))
    if(screening_crit == "low risk") {
        survival_dist <- survival_lb
    } else {
        survival_dist <- survival_ub
    }
    selected <- select_patients_point(time.points, survival_dist, screening_time, screening_prob, screening_crit)
    return(selected)
}

#' Select Patients Based on Predicted Survival Probabilities
#'
#' Given a survival matrix and a screening rule, identifies patients who meet the rule.
#'
#' @param time.points Time grid (must match columns of `survival_dist`).
#' @param survival_dist Matrix of survival probabilities (rows: patients).
#' @param screening_time Time at which to apply the rule.
#' @param screening_prob Survival threshold.
#' @param screening_crit `"low risk"` (survival >= threshold) or `"high risk"` (<=).
#'
#' @return A list with selected indices and the screening time.
#' @export
select_patients_point <- function(time.points, survival_dist, screening_time, screening_prob, screening_crit) {
  stopifnot(screening_crit %in% c("low risk", "high risk"))
  
  comparison.operator <- ifelse(screening_crit == "low risk", ">=", "<=")
  compare_fun <- match.fun(comparison.operator)
  
  n <- nrow(survival_dist)
  
  is.selected <- sapply(1:n, function(i) {
    surv_values <- as.numeric(survival_dist[i, ])
    
    if (length(time.points) == 1) {
      if (!isTRUE(all.equal(screening_time, time.points))) {
        stop("screening_time does not match time.points when length == 1")
      }
      estimated_prob <- surv_values[1]
    } else {
      interp_fun <- approxfun(time.points, surv_values, rule = 2, method = "linear", ties = "ordered")
      estimated_prob <- interp_fun(screening_time)
    }   
    return(compare_fun(estimated_prob, screening_prob))
  })  
  idx.selected <- which(is.selected)  
  out <- list(time = screening_time, selected = idx.selected)
  return(out)
}


#' Evaluate Screening Selections Without Oracle Event Times
#'
#' Computes conservative bounds on the survival rate of selected patients
#' without access to true event times. This function evaluates selection performance
#' under two extreme assumptions for censored data:
#' - lower bound assumes censored patients failed at the censoring time
#' - upper bound assumes censored patients never fail (event time = ∞)
#'
#' @param data.test A data frame containing the test dataset. Must include columns `time` and `status`.
#' @param idx.selected A vector of indices for selected patients (row indices into `data.test`).
#' @param screening_time Numeric, the time horizon for screening (e.g., 5 months).
#' @param screening_prob Numeric between 0 and 1; the survival threshold used for classification.
#' @param screening_crit A string, either `"low risk"` or `"high risk"`, indicating the direction of the screening rule.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Screening time}{The time horizon used for screening}
#'   \item{Screening rule}{A textual description of the rule applied}
#'   \item{Screened}{Number of selected patients}
#'   \item{Proportion survived (lower)}{Conservative estimate assuming censored events occurred at censoring time}
#'   \item{Proportion survived (upper)}{Liberal estimate assuming censored events never occurred}
#' }
#'
#' @seealso \code{\link{evaluate_selections}}, \code{\link{select_patients_band}}
#' @export
evaluate_selections_without_oracle <- function(data.test, idx.selected, screening_time, screening_prob, screening_crit) {
  stopifnot(screening_crit %in% c("low risk", "high risk"))
  
  # Create oracle versions of the data
  data.test.oracle.lower <- data.test
  data.test.oracle.lower$event_time <- data.test$time  # always use time, regardless of censoring
  
  data.test.oracle.upper <- data.test
  data.test.oracle.upper$event_time <- ifelse(data.test$status == 1, data.test$time, Inf)

  # Evaluate selections using lower oracle
  res.lower <- evaluate_selections(data.test.oracle.lower, idx.selected, screening_time, screening_prob, screening_crit)
  res.lower$`Proportion survived (lower)` <- res.lower$`Proportion survived`
  res.lower$`Proportion survived` <- NULL
  res.lower$Valid <- NULL

  # Evaluate selections using upper oracle
  res.upper <- evaluate_selections(data.test.oracle.upper, idx.selected, screening_time, screening_prob, screening_crit)
  res.upper$`Proportion survived (upper)` <- res.upper$`Proportion survived`
  res.upper$`Proportion survived` <- NULL
  res.upper$Valid <- NULL

  # Merge the two results by common keys
  merge(res.lower, res.upper, by = c("Screening time", "Screening rule", "Screened"))
}


#' Evaluate Screening Rule Against True or Imputed Event Times
#'
#' Given a set of selected patients, compute the proportion who survive beyond a given time,
#' and determine whether this satisfies the conservative screening criterion.
#'
#' @param data.test.oracle A data frame that must contain a column `event_time`, representing the known or imputed event time.
#' @param idx.selected Integer vector of row indices indicating selected patients.
#' @param screening_time Numeric value representing the time horizon for evaluation.
#' @param screening_prob Numeric value between 0 and 1, used as a threshold for the screening rule.
#' @param screening_crit Character string; either `"low risk"` (expecting high survival) or `"high risk"` (expecting low survival).
#'
#' @return A named list with the screening time, rule description, number of patients screened,
#'         observed survival proportion, and whether the result satisfies the conservative screening rule.
#' @export
evaluate_selections <- function(data.test.oracle, idx.selected, screening_time, screening_prob, screening_crit) {
  stopifnot(screening_crit %in% c("low risk", "high risk"))

  is_alive <- data.test.oracle$event_time >= screening_time
  n <- nrow(data.test.oracle)

  num_selections <- length(idx.selected)
  if (num_selections > 0) {
    selected_alive <- sum(is_alive[idx.selected])
    selected_dead <- num_selections - selected_alive
    prop_selected_alive <- selected_alive / num_selections
  } else {
    selected_alive <- 0
    selected_dead <- 0
    prop_selected_alive <- if (screening_crit == "low risk") 1 else 0
  }

  if (screening_crit == "low risk") {
    conservative <- prop_selected_alive > screening_prob
  } else {
    conservative <- prop_selected_alive < screening_prob
  }

  # Fallback in degenerate case
  if (is.na(conservative)) conservative <- TRUE

  rule_desc <- interpret_screening_rule(screening_time, screening_prob, screening_crit)

  result <- list(
    `Screening time` = screening_time,
    `Screening rule` = rule_desc,
    Screened = num_selections,
    `Proportion survived` = prop_selected_alive,
    Valid = conservative
  )

  return(result)
}


#' Format Screening Rule as String
#'
#' @param time Screening time.
#' @param prob Survival probability threshold.
#' @param crit Screening direction (`"low risk"` or `"high risk"`).
#'
#' @return A formatted string describing the rule.
#' @export
interpret_screening_rule <- function(time, prob, crit) {
    stopifnot(crit %in% c("low risk", "high risk"))
    stopifnot(time>=0)
    stopifnot(prob>=0 & prob<=1)
    if(crit=="high risk") {
        interpretation <- sprintf("high-risk patients with P(survival at time %.2f) < %.2f", time, prob)
    } else {
        interpretation <- sprintf("low-risk patients with P(survival at time %.2f) > %.2f", time, prob)
    }
    return(interpretation)
}
