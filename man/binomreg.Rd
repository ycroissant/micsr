% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binomreg.R
\name{binomreg}
\alias{binomreg}
\alias{residuals.binomreg}
\alias{glance.binomreg}
\alias{predict.binomreg}
\title{Binomial regression}
\usage{
binomreg(
  formula,
  data,
  weights,
  subset,
  na.action,
  offset,
  contrasts = NULL,
  link = c("identity", "probit", "logit"),
  method = c("ml", "twosteps", "minchisq", "test"),
  robust = TRUE,
  start = NULL,
  ...
)

\method{residuals}{binomreg}(object, ..., type = c("deviance", "pearson", "response"))

\method{glance}{binomreg}(x, ...)

\method{predict}{binomreg}(object, ..., type = c("response", "link"), newdata = NULL)
}
\arguments{
\item{formula}{a symbolic description of the model}

\item{data}{a data frame,}

\item{subset, weights, na.action, offset, contrasts}{see \code{stats::lm},}

\item{link}{one of \code{"identity"}, \code{"probit"} and "\code{logit}" to fit
respectively the linear probability, the probit and the logit
model}

\item{method}{\code{"ml"} for maximum likelihood (the only relevant
method for a regression without instrumental variables),
\code{"twosteps"} for two-steps estimator, \code{"minchisq"} for minimum
chi-squared estimator and \code{"test"} to get the exogeneity test,}

\item{robust}{only when \code{method = "twosteps"}, should the robust
covariance matrix be computed?}

\item{start}{a vector of starting values}

\item{...}{further arguments}

\item{object, x, type}{a \code{binomreg} object and the type of residuals
for the \code{residuals} method}

\item{newdata}{a new data frame for the \code{predict} method}
}
\value{
an object of class \code{c("binomreg", "micsr")}, see
\code{micsr::micsr} for further details
}
\description{
A unified interface for binomial regression models, including
linear probability, probit and logit models
}
\examples{
pbt <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = 'probit')
lpm <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = 'identity')
summary(pbt, vcov = "opg")
}
\keyword{models}
