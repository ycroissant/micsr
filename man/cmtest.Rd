% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmtest.R
\name{cmtest}
\alias{cmtest}
\alias{cmtest.tobit}
\alias{cmtest.micsr}
\alias{cmtest.censReg}
\alias{cmtest.glm}
\alias{cmtest.weibreg}
\title{Conditional moments test}
\usage{
cmtest(
  x,
  test = c("normality", "reset", "heterosc", "skewness", "kurtosis"),
  powers = 2:3,
  heter_cov = NULL,
  opg = FALSE
)

\method{cmtest}{tobit}(
  x,
  test = c("normality", "reset", "heterosc", "skewness", "kurtosis"),
  powers = 2:3,
  heter_cov = NULL,
  opg = FALSE
)

\method{cmtest}{micsr}(
  x,
  test = c("normality", "reset", "heterosc", "skewness", "kurtosis"),
  powers = 2:3,
  heter_cov = NULL,
  opg = FALSE
)

\method{cmtest}{censReg}(
  x,
  test = c("normality", "reset", "heterosc", "skewness", "kurtosis"),
  powers = 2:3,
  heter_cov = NULL,
  opg = FALSE
)

\method{cmtest}{glm}(
  x,
  test = c("normality", "reset", "heterosc", "skewness", "kurtosis"),
  powers = 2:3,
  heter_cov = NULL,
  opg = FALSE
)

\method{cmtest}{weibreg}(
  x,
  test = c("normality", "reset", "heterosc", "skewness", "kurtosis"),
  powers = 2:3,
  heter_cov = NULL,
  opg = FALSE
)
}
\arguments{
\item{x}{a fitted model, currently a tobit model either fitted by
\code{AER::tobit}, \code{censReg::censReg} or \code{micsr::tobit1} or a probit model fitted by
\code{glm} with \code{family = binomial(link = "probit")} or by \code{micsr::binomreg} with \code{link = "probit"}}

\item{test}{the kind of test to be performed, either a normality
test (or separately a test that the skewness or kurtosis
are 0 and 3), a heteroscedasticity test or a reset
test,}

\item{powers}{the powers of the fitted values that should be used
in the reset test,}

\item{heter_cov}{a one side formula that indicates the covariates
that should be used for the heteroscedasticity test (by default
all the covariates used in the regression are used),}

\item{opg}{a boolean, if \code{FALSE} (the default), the analytic
derivatives are used, otherwise the outer product of the
gradient formula is used}
}
\value{
an object of class \code{"htest"} containing the following components:
\itemize{
\item data.mane: a character string describing the fitted model
\item statistic: the value of the test statistic
\item parameter: degrees of freedom
\item p.value: the p.value of the test
\item method: a character indicating what type of test is performed
}
}
\description{
Conditional moments tests for maximum likelihood estimators,
particularly convenient for the probit and the tobit model to test
relevance of functional form, omitted variables, heteroscedasticity and
normality.
}
\examples{
charitable$logdon <- with(charitable, log(donation) - log(25))
ml <- tobit1(logdon ~ log(donparents) + log(income) + education +
             religion + married + south, data = charitable)
cmtest(ml, test = "heterosc")
cmtest(ml, test = "normality", opg = TRUE)
}
\references{
\insertRef{NEWE:85}{micsr}

\insertRef{PAGA:VELL:89}{micsr}

\insertRef{TAUC:85}{micsr}

\insertRef{WELL:03}{micsr}
}
\author{
Yves Croissant
}
\keyword{htest}
