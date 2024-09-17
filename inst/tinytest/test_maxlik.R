check_gh <- function(x) x$check_gradient$gradient < 1E-04 &
                      x$check_gradient$hessian < 1E-04
set.seed(1)
# weibreg
z <- unemp_duration
z$w <- runif(nrow(z))
ph_no <- weibreg(Surv(duration, censored == "no") ~ gender + age + log(wage + 1),
                 z, mixing = FALSE, weights = w, model = "ph", check_gradient = TRUE)
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
z$e <- sample(c(1, 0), nrow(z), replace = TRUE, prob = c(.9, .1))
olgt <- ordreg(factor(dindx) ~ rhs1 + catchup, z, link = "logit", check_gradient = TRUE, weights = w)
opbt <- update(olgt, link = "probit")
expect_true(check_gh(olgt), info = "ordered logit")
expect_true(check_gh(opbt), info = "ordered probit")
survlgt <- ordreg(Surv(dindx, e) ~ rhs1 + catchup, z, link = "logit", check_gradient = TRUE)
expect_true(check_gh(survlgt), info = "ordered logit with censored data")


# Poisreg

z <- trips
z$w <- runif(nrow(z))
pois <- poisreg(trips ~ workschl + size + dist + smsa + fulltime + distnod +
                   realinc + weekend + car, z, check_gradient = TRUE)
# problème avec la hessienne de nb1
nb1 <- update(pois, mixing = "gamma", vlink = "nb1")
lnr <- update(pois, mixing = "lognorm")
expect_true(check_gh(pois), info = "poisson")
expect_true(check_gh(nb1), info = "negbin 1")
expect_true(check_gh(lnr), info = "log-normal Poisson")

pois <- update(pois, weights = w)
nb1 <- update(nb1,  weights = w)
lnr <- update(lnr,  weights = w)
expect_true(check_gh(pois), info = "poisson with weights")
expect_true(check_gh(nb1), info = "negbin 1 with weights")
expect_true(check_gh(lnr), info = "log-normal Poisson with weights")


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

