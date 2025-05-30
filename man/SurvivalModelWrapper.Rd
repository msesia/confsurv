% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_wrappers.R
\name{SurvivalModelWrapper}
\alias{SurvivalModelWrapper}
\title{Abstract Survival Model Wrapper}
\description{
This R6 class provides a common interface for fitting and predicting survival models.
Subclasses should implement the \code{fit()} method, and optionally override \code{predict()} and \code{predict_quantiles()}.
}
\section{Methods}{

\describe{
\item{\code{fit(formula, data)}}{Abstract method to fit a survival model.}
\item{\code{predict(new_data, time.points = NULL)}}{Predict survival curves.}
\item{\code{predict_quantiles(new_data, probs)}}{Return survival quantiles for each individual.}
\item{\code{predict_interpolate(survival_probs, original_failure_times, time.points)}}{Linear interpolation utility.}
\item{\code{parse_formula(formula, data)}}{Extracts time, status, and covariates from formula.}
\item{\code{select_columns(data)}}{Restricts columns to specified covariates.}
}
}

\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{model}}{Fitted model object.}

\item{\code{formula}}{Formula used for model fitting.}

\item{\code{time.points}}{Default failure times used for interpolation.}

\item{\code{use_covariates}}{Optional vector of covariate names to include.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-SurvivalModelWrapper-new}{\code{SurvivalModelWrapper$new()}}
\item \href{#method-SurvivalModelWrapper-fit}{\code{SurvivalModelWrapper$fit()}}
\item \href{#method-SurvivalModelWrapper-predict}{\code{SurvivalModelWrapper$predict()}}
\item \href{#method-SurvivalModelWrapper-predict_quantiles}{\code{SurvivalModelWrapper$predict_quantiles()}}
\item \href{#method-SurvivalModelWrapper-predict_interpolate}{\code{SurvivalModelWrapper$predict_interpolate()}}
\item \href{#method-SurvivalModelWrapper-parse_formula}{\code{SurvivalModelWrapper$parse_formula()}}
\item \href{#method-SurvivalModelWrapper-select_columns}{\code{SurvivalModelWrapper$select_columns()}}
\item \href{#method-SurvivalModelWrapper-clone}{\code{SurvivalModelWrapper$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SurvivalModelWrapper-new"></a>}}
\if{latex}{\out{\hypertarget{method-SurvivalModelWrapper-new}{}}}
\subsection{Method \code{new()}}{
Constructor for SurvivalModelWrapper.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SurvivalModelWrapper$new(use_covariates = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{use_covariates}}{Optional character vector of covariates to be used.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SurvivalModelWrapper-fit"></a>}}
\if{latex}{\out{\hypertarget{method-SurvivalModelWrapper-fit}{}}}
\subsection{Method \code{fit()}}{
Abstract method for fitting a survival model. Should be implemented in subclass.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SurvivalModelWrapper$fit(formula, data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{formula}}{A survival formula (e.g., Surv(time, status) ~ predictors).}

\item{\code{data}}{A data.frame with columns for time, status, and covariates.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SurvivalModelWrapper-predict"></a>}}
\if{latex}{\out{\hypertarget{method-SurvivalModelWrapper-predict}{}}}
\subsection{Method \code{predict()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SurvivalModelWrapper$predict(new_data, time.points = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{new_data}}{A data.frame of new individuals with covariates.}

\item{\code{time.points}}{Optional vector of time points. Defaults to model's time.points.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list with \code{predictions} (matrix of survival probabilities) and \code{time.points}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SurvivalModelWrapper-predict_quantiles"></a>}}
\if{latex}{\out{\hypertarget{method-SurvivalModelWrapper-predict_quantiles}{}}}
\subsection{Method \code{predict_quantiles()}}{
Predicts survival quantiles (e.g. median, quartiles).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SurvivalModelWrapper$predict_quantiles(new_data, probs = c(0.25, 0.5, 0.75))}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{new_data}}{New data with covariates.}

\item{\code{probs}}{Numeric vector of quantile probabilities. Defaults to c(0.25, 0.5, 0.75).}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A data.frame of predicted quantiles for each individual.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SurvivalModelWrapper-predict_interpolate"></a>}}
\if{latex}{\out{\hypertarget{method-SurvivalModelWrapper-predict_interpolate}{}}}
\subsection{Method \code{predict_interpolate()}}{
Linearly interpolates survival probabilities at new time points.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SurvivalModelWrapper$predict_interpolate(
  survival_probs,
  original_failure_times,
  time.points
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{survival_probs}}{Matrix of survival probabilities.}

\item{\code{original_failure_times}}{Original time points corresponding to survival_probs.}

\item{\code{time.points}}{New time points to interpolate to.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Matrix of interpolated probabilities.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SurvivalModelWrapper-parse_formula"></a>}}
\if{latex}{\out{\hypertarget{method-SurvivalModelWrapper-parse_formula}{}}}
\subsection{Method \code{parse_formula()}}{
Parses a survival formula and extracts components.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SurvivalModelWrapper$parse_formula(formula, data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{formula}}{A survival formula, e.g., Surv(time, status) ~ X1 + X2.}

\item{\code{data}}{A data.frame containing survival and covariate columns.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list with \code{time}, \code{status}, and \code{covariates}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SurvivalModelWrapper-select_columns"></a>}}
\if{latex}{\out{\hypertarget{method-SurvivalModelWrapper-select_columns}{}}}
\subsection{Method \code{select_columns()}}{
Subsets the data to only include selected covariates.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SurvivalModelWrapper$select_columns(new_data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{new_data}}{A data.frame to filter.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A filtered data.frame with only the relevant covariates.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SurvivalModelWrapper-clone"></a>}}
\if{latex}{\out{\hypertarget{method-SurvivalModelWrapper-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SurvivalModelWrapper$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
