% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binomreg.R
\name{binomreg}
\alias{binomreg}
\alias{logLik.binomreg}
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
  link = c("lm", "probit", "logit"),
  start = NULL,
  ...
)

\method{logLik}{binomreg}(object, ..., type = c("model", "null"))

\method{residuals}{binomreg}(object, ..., type = c("deviance", "pearson", "response"))

\method{glance}{binomreg}(x, ...)

\method{predict}{binomreg}(object, ..., type = c("response", "link"), newdata = NULL)
}
\arguments{
\item{formula}{a symbolic description of the model, (for the count
component and for the selection equation).}

\item{data}{a data frame,}

\item{subset, weights, na.action, offset}{see \code{stats::lm},}

\item{link}{one of \code{"lm"}, \code{"probit"} and "\code{logit}" to fit
respectively the linear probability, the probit and the logit
model}

\item{start}{a vector of starting values, in this case, no
estimation}

\item{...}{further arguments}

\item{object, x, type}{a \code{binomreg} object and the type of
log-likelihood / residuals for the \code{logLik} / \code{residuals}
method}

\item{newdata}{a new data frame for the \code{predict} method}
}
\value{
an object of class \verb{c(}"binomreg", "micsr")\verb{, see }micsr::micsr` for further details
}
\description{
A unified interface perform binomial regression calling \code{lm} to fit
the linear-probability model and \code{glm} to fit the probit and the
logit model.
}
