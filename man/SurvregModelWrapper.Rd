% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_wrappers.R
\name{SurvregModelWrapper}
\alias{SurvregModelWrapper}
\title{Parametric Survival Model Wrapper}
\description{
Wrapper for parametric survival models via \code{survival::survreg()}.
}
\section{Super class}{
\code{\link[confsurv:SurvivalModelWrapper]{confsurv::SurvivalModelWrapper}} -> \code{SurvregModelWrapper}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{dist}}{A character string indicating the distribution used by the model (e.g., "weibull", "lognormal").}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-SurvregModelWrapper-new}{\code{SurvregModelWrapper$new()}}
\item \href{#method-SurvregModelWrapper-fit}{\code{SurvregModelWrapper$fit()}}
\item \href{#method-SurvregModelWrapper-predict_quantiles}{\code{SurvregModelWrapper$predict_quantiles()}}
\item \href{#method-SurvregModelWrapper-clone}{\code{SurvregModelWrapper$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="confsurv" data-topic="SurvivalModelWrapper" data-id="parse_formula"><a href='../../confsurv/html/SurvivalModelWrapper.html#method-SurvivalModelWrapper-parse_formula'><code>confsurv::SurvivalModelWrapper$parse_formula()</code></a></span></li>
<li><span class="pkg-link" data-pkg="confsurv" data-topic="SurvivalModelWrapper" data-id="predict"><a href='../../confsurv/html/SurvivalModelWrapper.html#method-SurvivalModelWrapper-predict'><code>confsurv::SurvivalModelWrapper$predict()</code></a></span></li>
<li><span class="pkg-link" data-pkg="confsurv" data-topic="SurvivalModelWrapper" data-id="predict_interpolate"><a href='../../confsurv/html/SurvivalModelWrapper.html#method-SurvivalModelWrapper-predict_interpolate'><code>confsurv::SurvivalModelWrapper$predict_interpolate()</code></a></span></li>
<li><span class="pkg-link" data-pkg="confsurv" data-topic="SurvivalModelWrapper" data-id="select_columns"><a href='../../confsurv/html/SurvivalModelWrapper.html#method-SurvivalModelWrapper-select_columns'><code>confsurv::SurvivalModelWrapper$select_columns()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SurvregModelWrapper-new"></a>}}
\if{latex}{\out{\hypertarget{method-SurvregModelWrapper-new}{}}}
\subsection{Method \code{new()}}{
Constructor for SurvregModelWrapper.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SurvregModelWrapper$new(use_covariates = NULL, dist = "weibull")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{use_covariates}}{Optional character vector of covariate names to use.}

\item{\code{dist}}{Distribution to use in the survreg model (e.g., \"weibull\", \"lognormal\").}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SurvregModelWrapper-fit"></a>}}
\if{latex}{\out{\hypertarget{method-SurvregModelWrapper-fit}{}}}
\subsection{Method \code{fit()}}{
Fit a parametric survival model using survreg.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SurvregModelWrapper$fit(formula, data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{formula}}{A formula of the form \code{Surv(time, status) ~ predictors}.}

\item{\code{data}}{A data.frame with survival outcome and covariates.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SurvregModelWrapper-predict_quantiles"></a>}}
\if{latex}{\out{\hypertarget{method-SurvregModelWrapper-predict_quantiles}{}}}
\subsection{Method \code{predict_quantiles()}}{
Predict survival time quantiles.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SurvregModelWrapper$predict_quantiles(new_data, probs = c(0.25, 0.5, 0.75))}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{new_data}}{A data.frame of covariates.}

\item{\code{probs}}{A numeric vector of quantile probabilities.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A data.frame of quantile predictions.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SurvregModelWrapper-clone"></a>}}
\if{latex}{\out{\hypertarget{method-SurvregModelWrapper-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SurvregModelWrapper$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
