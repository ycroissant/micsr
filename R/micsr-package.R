#' **micsr** : Microeconometrics with R
#'
#' The micsr package is the companion package to the book
#' "Microeconometrics with R" (Chapman and Hall/CRC The R Series). It
#' includes function to estimate and to test models, miscellanous
#' tools and data sets:
#'
#' - functions to estimate models:
#'
#'   - `binomreg`: binomial regression models, Rivers and Vuong (1988),
#'   - `bivprobit`: bivariate probit model
#'   - `clm`: constrained linear models,
#'   - `escount`: endogenous switching and selection model for count data, Terza (1998),
#'   - `expreg`: exponential conditional mean models, Mullahy (1997),
#'   - `loglm`: log-linear models,
#'   - `ordreg`: ordered regression models,
#'   - `poisreg`: poisson models,
#'   - `pscore`: matching, Dehejia and Wahba (2002),
#'   - `tobit1`: tobit-1 model, Tobin (1958), Smith and Blundel (1986), Powel (1986).
#'
#' - functions for statistical tests and diagnostic:
#'
#'   - `cmtest`: conditional moment tests, Newey (1985), Tauchen (1985),
#'   - `ftest`: F statistic,
#'   - `hausman`: Hausman's test, Hausman (1978),
#'   - `ndvuong`: non-degenerate Vuong test, Vuong (1989), Shi (2015),
#'   - `rsq`: different flavors of R squared,
#'   - `sargan`: Sargan's test, Sargan (1958),
#'   - `scoretest`: score, or Lagrange multiplier test.
#'
#' - miscellanous tools
#'
#'   - `gaze`: print a short summary of an object,
## #'   - `dummy`: generate a set of dummy variables from a factor,
#'   - `newton`: Newton-Raphson optimization method, using the analytical gradient and hessian,
#'   - `mills`: compute the inverse mills ratio and its first two derivatives,
#'   - `stder`: extract the standard errors of a fitted model,
#'   - `npar`: extract the number of parameters in a fitted model.
#'
#' - data sets:
#' 
#'   - `apples`: Apple production, Ivaldi and al. (1996), constrained
#'   linear model,
#'   - `birthwt`: Cigarette smoking and birth weigth, Mullahy
#' (1997), exponentional conditional mean regression model,
#'   - `charitable`: Intergenerational transmission of charitable
#' giving, Wilhem (2008), Tobit-1 model,
#'   - `cigmales`: Cigarettes consumption and smoking habits, 
#' Mullahy (1997), exponentional conditional mean regression mdodel,
#'   - `drinks`: Physician advice on alcohol consumption, Kenkel and
#' Terza (2001), endogenous switching model for count data,
#'   - `ferediv`: Foreign exchange derivatives use by large US bank
#' holding companies, Adkins (2012), instrumental variable probit
#' model,
#'   - `fin_reform`: Political economy of financial reforms, Abiad and
#' Mody (2005), ordered regression model,
#'   - `housprod`: Household production, Kerkhofs and Kooreman (2003),
#' bivariate probit model,
#'   - `mode_choice`: Choice between car and transit, Horowitz (1993),
#' probit model,
#'   - `trade_protection`: Lobying and trade protection, Atschke and
#' Sherlund (2006), instrumental variable Tobit-1 model,
#'   - `trips`: Determinants of household trip taking, Terza (1998),
#' endogenous switching model for count data,
#'   - `turnout`: Turnout in Texas liquor referenda, Coate and Conlin
#' (2004), non-degenerate Vuong test,
#'   - `twa`: Temporary help jobs and permanent employment, Ichino,
#' Mealli and Nannicini (2008), matching.
#'
#' - vignettes:
#'
#'   - charitable: Estimating the Tobit-1 model with the charitable
#'   data set
#'   - escount: Endogenous switching or sample selection models for
#'   count data
#'   - expreg: Exponentional conditional mean models with endogeneity
#'   - ndvvuong: Implementation of Shi's non-degeranate Vuong test
#'
#' 
#' We tried to keep the sets of package on which **micsr** depends on
#' as small as possible. **micsr** depends on **Formula**,
#' **generics**, **Rdpack**, **knitr**, **sandwich** and on a subset
#' of the **tidyverse** metapackage (**ggplot2**, **dplyr**,
#' **purrr**, **tidyselect**, **magrittr**, **tibble**, **rlang**). We
#' borrowed the gaussian quadrature function from the **statmod**
#' package (Smyth and al., 2023), and the distribution function of
#' quadratic forms in normal variables from the **CompQuadForm**
#' package (Duchesne and Lafaye, 2010).
#'
#' @docType package
#' @keywords internal
#' @references
#' \insertRef{ABIA:MODY:05}{micsr}
#'
#' \insertRef{ADKI:12}{micsr}
#' 
#' \insertRef{COAT:CONL:04}{micsr}
#'
#' \insertRef{DEHE:WAHB:02}{micsr}
#' 
#' \insertRef{DUCH:LAFA:10}{micsr}
#'
#' \insertRef{HAUS:78}{micsr}
#'
#' \insertRef{ICHI:MEAL:NANN:08}{micsr}
#'
#' \insertRef{IVAL:LADO:OSSA:SIMI:96}{micsr}
#'
#' \insertRef{KENK:TERZ:01}{micsr}
#'
#' \insertRef{KERK:KOOR:03}{micsr}
#' 
#' \insertRef{MATS:SHER:06}{micsr}
#' 
#' \insertRef{MULL:97}{micsr}
#'
#' \insertRef{NEWE:85}{micsr}
#'
#' \insertRef{POWE:86}{micsr}
#'
#' \insertRef{RIVE:VUON:88}{micsr}
#'
#' \insertRef{SARG:58}{micsr}
#'
#' \insertRef{SHI:15}{micsr}
#'
#' \insertRef{SMIT:BLUN:86}{micsr}
#'
#' \insertRef{SMYT:CHEN:23}{micsr}
#' 
#' \insertRef{TAUC:85}{micsr}
#' 
#' \insertRef{TERZ:98}{micsr}
#'
#' \insertRef{TOBI:58}{micsr}
#'
#' \insertRef{VUON:89}{micsr}
#'
#' \insertRef{WILH:08}{micsr}
"_PACKAGE"
