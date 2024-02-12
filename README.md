# **micsr**: Microeconometrics with R

## About

The micsr package is the companion package to the book
"Microeconometrics with R" (Chapman and Hall/CRC The R Series). It
includes function to estimate and to test models, miscellanous
tools and data sets:

`micsr` provides methods described in the book that are not available in R. This includes:

- testing functions:
  - `scoretest` for score (or Lagrange Multiplier) tests,
  - `cmtest` for conditional moment tests,
  - `sargan` for Sargan-Hansen's test of overidentified moment
    conditions,
  - `hausman` for Hausman's test,
  - `ndvuong` for non-degenerate Vuong test of Shi,
  - `ftest` for F test

- function to estimate models:
  - `binomreg` binomal variable regression,
  - `bivprobit` bivariate probit model
  - `clm` constrained linear regression,
  - `escount` endogenous switching model for count data,
  - `expreg` exponentional mean regression model,
  - `loglm` log-linear models
  - `tobit1` one-equation tobit model,
  - `ordreg` ordered variable regression,
  - `poisreg` poissong regression model,
  - `pscore`: matching models

- miscellanous tools

  - `gaze`: print a short summary of an object,
  - `dummy`: generate a set of dummy variables from a factor,
  - `newton`: Newton-Raphson optimization method, using the analytical gradient and hessian,
  - `mills`: compute the inverse mills ratio and its first two derivatives,
  - `stder`: extract the standard errors of a fitted model,
  - `npar`: extract the number of parameters in a fitted model.

- data sets:
  - `apples`: Apple production, Ivaldi and al. (1996), constrained
  linear model,
  - `birthwt`: Cigarette smoking and birth weigth, Mullahy
(1997), exponentional conditional mean regression model,
  - `charitable`: Intergenerational transmission of charitable
giving, Wilhem (2008), Tobit-1 model,
  - `cigmales`: Cigarettes consumption and smoking habits, 
Mullahy (1997), exponentional conditional mean regression mdodel,
  - `drinks`: Physician advice on alcohol consumption, Kenkel and
Terza (2001), endogenous switching model for count data,
  - `ferediv`: Foreign exchange derivatives use by large US bank
holding companies, Adkins (2012), instrumental variable probit
model,
  - `fin_reform`: Political economy of financial reforms, Abiad and
Mody (2005), ordered regression model,
  - `housprod`: Household production, Kerkhofs and Kooreman (2003),
bivariate probit model,
  - `mode_choice`: Choice between car and transit, Horowitz (1993),
probit model,
  - `trade_protection`: Lobying and trade protection, Atschke and
Sherlund (2006), instrumental variable Tobit-1 model,
  - `trips`: Determinants of household trip taking, Terza (1998),
endogenous switching model for count data,
  - `turnout`: Turnout in Texas liquor referenda, Coate and Conlin
(2004), non-degenerate Vuong test,
  - `twa`: Temporary help jobs and permanent employment, Ichino,
Mealli and Nannicini (2008), matching.

- vignettes:
  - charitable: Estimating the Tobit-1 model with the charitable
data set
  - escount: Endogenous switching or sample selection models for
count data
  - expreg: Exponentional conditional mean models with endogeneity
  - ndvvuong: Implementation of Shi's non-degeranate Vuong test



We tried to keep the sets of package on which **micsr** depends on as
small as possible. **micsr** depends on **Formula**, **generics**,
**Rdpack**, **knitr**, **sandwich** and on a subset of the
**tidyverse** metapackage (**ggplot2**, **dplyr**, **purrr**,
**tidyselect**, **magrittr**, **tibble**, **rlang**). We borrowed the
gaussian quadrature function from the **statmod** package (Smyth and
al., 2023), and the distribution function of quadratic forms in normal
variables from the **CompQuadForm** package (Duchesne and Lafaye,
2010).



## Installation

The `micsr` is not yet on CRAN. To install it, first install the `pak` package and enter:

```
pak::pkg_install("ycroissant/micsr")
```
