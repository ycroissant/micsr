% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsq.R
\name{rsq}
\alias{rsq}
\alias{rsq.lm}
\alias{rsq.micsr}
\title{Coefficient of determination}
\usage{
rsq(x, type)

\method{rsq}{lm}(x, type = c("raw", "adj"))

\method{rsq}{micsr}(
  x,
  type = c("mcfadden", "cox_snell", "cragg_uhler", "aldrich_nelson", "veall_zimm",
    "estrella", "cor", "ess", "rss", "tjur", "mckel_zavo", "wald", "score", "lr")
)
}
\arguments{
\item{x}{fitted model}

\item{type}{the type of coefficient of determination}
}
\value{
a numeric scalar.
}
\description{
A generic function to compute different flavors of coefficients of determination
}
\examples{
pbt <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = 'probit')
rsq(pbt)
rsq(pbt, "estrella")
rsq(pbt, "veall_zimm")
}
