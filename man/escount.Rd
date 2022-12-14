% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/escount.R
\name{escount}
\alias{escount}
\title{Endogenous switching and sample selection models for count data}
\usage{
escount(
  formula,
  data,
  subset,
  weights,
  na.action,
  offset,
  start = NULL,
  R = 16,
  hessian = FALSE,
  method = c("twosteps", "ml"),
  model = c("es", "ss")
)
}
\arguments{
\item{formula}{a Formula object which includes two responses (the
count and the binomial variables) and two sets of covariates
(for the count component and for the selection equation).}

\item{data}{a data frame,}

\item{subset, weights, na.action, offset}{see \code{stats::lm}}

\item{start}{an optional vector of starting values,}

\item{R}{the number of points for the Gauss-Hermite quadrature
approximation,}

\item{hessian}{if \code{TRUE}, the numerical hessian is computed,
otherwise the covariance matrix of the coefficients is computed
using the outer product of the gradient,}

\item{method}{one of \code{'ML'} for maximum likelihood estimation (the
default) or \code{'twosteps'} for the two-steps NLS method,}

\item{model}{one of \code{'es'} for endogenous switching (the default)
or \code{'ss'} for sample selection.}
}
\value{
an object of class \code{'escount'} and \code{'micsr'} which is a
list containing the following components: - coefficients: a
named vector of coefficients, - sigma, rho: the estimated
values of sigma and rho for the \code{twosteps} method, - residuals:
the vector of residuals, - fitted.values: the fitted values -
vcov: estimation of the covariance matrix of the estimators -
value: value of the objective function at convergence - model:
the model frame - call: the matched call - K: the number of
covariates - L: the number of instruments - df.residual: the
degrees of freedom of the regression - xlevels: a record of the
levels of the factors used in fitting - na.action: information
returned by \code{model.frame} on the sepcial handling of \code{NA}'s
}
\description{
Heckman's like estimator are extended to count data, using either
non-linear least squares, using a correction term computed from a
first-stage probit estimation or maximum likelihood
}
\references{
\insertRef{TERZ:98}{micsr}

\insertRef{GREE:01}{micsr}
}
\author{
Yves Croissant
}
