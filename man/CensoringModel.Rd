% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/censoring_wrappers.R
\name{CensoringModel}
\alias{CensoringModel}
\title{CensoringModel Class}
\description{
This class inherits from \code{SurvivalModelWrapper} and wraps any compatible survival model.
It fits a censoring distribution (i.e., the survival function of the censoring time)
by inverting the event indicator: \code{status := 1 - status}.
}
\details{
Wraps a \code{SurvivalModelWrapper} to model censoring distributions by flipping event indicators.
}
\section{Methods}{

\describe{
\item{\code{new(model)}}{Creates a new \code{CensoringModel} with a specified survival model wrapper.}
\item{\code{fit(data)}}{Fits the censoring distribution using the flipped event indicators.}
\item{\code{predict(new_data, time.points = NULL)}}{Predicts the censoring survival curve.}
}
}

\section{Super class}{
\code{\link[confsurv:SurvivalModelWrapper]{confsurv::SurvivalModelWrapper}} -> \code{CensoringModel}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-CensoringModel-new}{\code{CensoringModel$new()}}
\item \href{#method-CensoringModel-fit}{\code{CensoringModel$fit()}}
\item \href{#method-CensoringModel-predict}{\code{CensoringModel$predict()}}
\item \href{#method-CensoringModel-clone}{\code{CensoringModel$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="confsurv" data-topic="SurvivalModelWrapper" data-id="parse_formula"><a href='../../confsurv/html/SurvivalModelWrapper.html#method-SurvivalModelWrapper-parse_formula'><code>confsurv::SurvivalModelWrapper$parse_formula()</code></a></span></li>
<li><span class="pkg-link" data-pkg="confsurv" data-topic="SurvivalModelWrapper" data-id="predict_interpolate"><a href='../../confsurv/html/SurvivalModelWrapper.html#method-SurvivalModelWrapper-predict_interpolate'><code>confsurv::SurvivalModelWrapper$predict_interpolate()</code></a></span></li>
<li><span class="pkg-link" data-pkg="confsurv" data-topic="SurvivalModelWrapper" data-id="predict_quantiles"><a href='../../confsurv/html/SurvivalModelWrapper.html#method-SurvivalModelWrapper-predict_quantiles'><code>confsurv::SurvivalModelWrapper$predict_quantiles()</code></a></span></li>
<li><span class="pkg-link" data-pkg="confsurv" data-topic="SurvivalModelWrapper" data-id="select_columns"><a href='../../confsurv/html/SurvivalModelWrapper.html#method-SurvivalModelWrapper-select_columns'><code>confsurv::SurvivalModelWrapper$select_columns()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-CensoringModel-new"></a>}}
\if{latex}{\out{\hypertarget{method-CensoringModel-new}{}}}
\subsection{Method \code{new()}}{
Constructor for CensoringModel.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CensoringModel$new(model)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{model}}{A survival model wrapper (e.g., output of \code{init_censoring_model()}).}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-CensoringModel-fit"></a>}}
\if{latex}{\out{\hypertarget{method-CensoringModel-fit}{}}}
\subsection{Method \code{fit()}}{
Fit the censoring model using flipped event indicators.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CensoringModel$fit(data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{A data frame containing time, status, and covariates.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-CensoringModel-predict"></a>}}
\if{latex}{\out{\hypertarget{method-CensoringModel-predict}{}}}
\subsection{Method \code{predict()}}{
Predict the censoring survival function.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CensoringModel$predict(new_data, time.points = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{new_data}}{New covariate data for prediction.}

\item{\code{time.points}}{Optional time grid for survival curve evaluation.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list with predicted survival probabilities and time points.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-CensoringModel-clone"></a>}}
\if{latex}{\out{\hypertarget{method-CensoringModel-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CensoringModel$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
