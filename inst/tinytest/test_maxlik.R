library(dplyr)
library(survival)
.my_eps <- 1E-04

check_gh <- function(x) x$check_gradient$gradient < .my_eps &
                            x$check_gradient$hessian < .my_eps
set.seed(1)
# weibreg
z <- unemp_duration
z$w <- runif(nrow(z))
ph_no <- weibreg(Surv(duration / sd(duration), censored == "no") ~ gender + log(wage + 1),
                 z, mixing = FALSE, weights = w, model = "ph", check_gradient = TRUE, opt = "nr")
ph_yes <- update(ph_no, mixing = TRUE)
aft_no <- update(ph_no, model = "aft")
aft_yes <- update(aft_no, mixing = TRUE)
expect_true(check_gh(ph_no),  info = "ph")
expect_true(check_gh(ph_no),  info = "ph, mixed")
expect_true(check_gh(aft_no),  info = "aft")
expect_true(check_gh(aft_yes),  info = "aft, mixed")

# binomreg
z <- mode_choice
z$w <- runif(nrow(z))
#z$w <- rep(1, nrow(z))
pbt <- binomreg(mode ~ cost + ivtime + ovtime, data = z,
                link = 'probit', weights = w, check_gradient = TRUE)
lgt <- update(pbt, link = 'logit')
lin <- update(pbt, link = 'identity')
expect_true(check_gh(pbt), info = "probit")
expect_true(check_gh(lgt), info = "logit")
expect_true(check_gh(lin), info = "identity")


# ordreg
z <- fin_reform
z$w <- runif(nrow(z))
z$y <- as.numeric(as.factor(z$dindx))
z$y <- dplyr::case_when(z$y < 5 ~ 1, z$y == 5 ~ 2, z$y > 5 ~ 3)
z$e <- sample(c(1, 0), nrow(z), replace = TRUE, prob = c(.9, .1))
olgt <- ordreg(y ~ rhs1 + catchup, z, link = "logit", check_gradient = TRUE, weights = w)
opbt <- update(olgt, link = "probit")
expect_true(check_gh(olgt), info = "ordered logit")
expect_true(check_gh(opbt), info = "ordered probit")
survlgt <- ordreg(Surv(y, e) ~ rhs1 + catchup, z, link = "logit", check_gradient = TRUE)
survpbt <- update(survlgt, link = 'probit')
expect_true(check_gh(survlgt), info = "ordered logit with censored data")
expect_true(check_gh(survpbt), info = "ordered probit with censored data")

# Poisreg

z <- trips
z$w <- runif(nrow(z))
pois <- poisreg(trips ~ workschl + size + dist, z, check_gradient = TRUE, weights = w)
# problème avec la hessienne de nb1
nb1 <- update(pois, mixing = "gamma", vlink = "nb1")
lnr <- update(pois, mixing = "lognorm")
expect_true(check_gh(pois), info = "poisson")
expect_true(check_gh(nb1), info = "negbin 1")
expect_true(check_gh(lnr), info = "log-normal Poisson")


# tobit1
charitable$logdon <- with(charitable, log(donation) - log(25))
charitable$w <- runif(nrow(charitable))
charitable$w <- charitable$w / mean(charitable$w)
tbt1 <- tobit1(logdon ~ log(income) + religion, data = charitable,
               check_gradient = TRUE, weights = w)
tbt1two <- update(tbt1, right = 5)
tbt1twob <- update(tbt1, right = 5, left = 2)
c1 <- update(tbt1, sample = "truncated")
c1b <- update(tbt1, sample = "truncated", left = 1)

expect_true(check_gh(tbt1), info = "tobit censored")
expect_true(check_gh(tbt1two), info = "two-limits tobit censored")
expect_true(check_gh(tbt1twob), info = "two-limits user defined tobit censored")
expect_true(check_gh(c1), info = "tobit truncated")
expect_true(check_gh(c1), info = "user defined tobit truncated")

# ivldv


## inst <- ~ sic3 + k_serv + inv + engsci + whitecol + skill + semskill + cropland + 
##     pasture + forest + coal + petro + minerals + scrconc + bcrconc + scrcomp +
##     bcrcomp + meps + kstock + puni + geog2 + tenure + klratio + bunion
## trade_protection <- dplyr::mutate(micsr::trade_protection,
##                                  y = ntb / (1 + ntb),
##                                  x1 = vshipped / imports / elast,
##                                  x2 = cap * x1,
##                                  x3 = labvar)
## GH <- ivldv(Formula::as.Formula(y  ~  x1 + x2, inst), trade_protection,
##             method = "twosteps", model = "tobit") 
## Full <- ivldv(Formula::as.Formula(y ~ x1 + x2 + labvar, inst), trade_protection,
##               method = "twosteps", model = "tobit") 
## Short <- ivldv(Formula::as.Formula(y ~ x1 + I(x2 + labvar), inst),
##                  trade_protection, method = "twosteps", model = "tobit")
## bank_msq <- ivldv(federiv ~ eqrat + optval + bonus + ltass + linsown + linstown +
##                   roe + mktbk + perfor + dealdum + div + year | . - eqrat - bonus -
##                   optval + no_emp + no_subs + no_off + ceo_age + gap + cfa,
##                   data = federiv, method = "minchisq")
## bank_ml <- update(bank_msq, method = "ml")



## form_ord <- Surv(duration, end == "recall") ~ age + sex + educ + race +
##     nb + ui + marital + unemp + wifemp + homeowner + occupation + industry
## form_ord <- Surv(duration, end == "recall") ~ age + sex
## recall_ord <- ordreg(form_ord, recall, link = "cloglog", check_gradient = TRUE)

