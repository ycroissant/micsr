# The `micsr` Package - Companion package to the forthcoming book "Microeconometrics With R", CRC/Chapman

## warning

- the example of ndvtest doesn't work anymor
- the ivldv vignette is highly commented

## About

`micsr` provides methods described in the book that are not available in R. This includes:

- testing functions:
  - `scoretest` for score (or Lagrange Multiplier) tests,
  - `cmtest` for conditional moment tests,
  - `sargantest` for Sargan-Hansen's test of overidentified moment
    conditions,
  - `haustest` for Hausman's test,
  - `ndvtest` for non-degenerate Vuong test of Shi,
  - `endogtest` tests for endogeneity (for probit and tobit model).
- function to estimate models:
  - `expreg` exponentional mean regression model,
  - `escount` endogenous switching model for count data,
  - `ivldv` instrumental variable estimators for probit and logit,
  - `tobit1` one-equation tobit model,
  - `binomreg` binomal variable regression,
  - `ordreg` ordered variable regression,
  - `bivprobit` bivariate probit model
- more than 30 data sets

# Installation

The `micsr` is not yet on CRAN. To install it, first install the `remotes` package and enter:

```
remotes::install_github("ycroissant/micsr", build_vignettes = TRUE)
```
