#' Coefficient of determination
#'
#' A generic function to compute different flavors of coefficients of determination
#'
#' @name rsq
#' @param x fitted model
#' @param type the type of coefficient of determination
#' @return a numeric scalar
#' @importFrom stats model.response model.frame resid logLik
#' @examples
#' pbt <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = 'probit')
#' rsq(pbt)
#' rsq(pbt, "estrella")
#' rsq(pbt, "veall_zimm")
#' @export
rsq <- function(x, type){
    UseMethod("rsq")
}

#' @rdname rsq
#' @export
rsq.lm <- function(x, type = c("raw", "adj")){
    .type <- match.arg(type)
    if (.type == "raw") summary(x)$r.squared else summary(x)$adj.r.squared
}

#' @rdname rsq
#' @export
rsq.binomreg <- function(x, type = c("mcfadden", "cox_snell", "cragg_uhler",
                                     "aldrich_nelson", "veall_zimm", "estrella",
                                     "cor", "ess", "rss", "tjur",
                                     "mckel_zavo")){
    type <- match.arg(type)
    y <- model.response(model.frame(x))
    yb <- mean(y)
    logL <- as.numeric(logLik(x))
    if (type == "mcfadden"){
        r2 <- 1 - deviance(x) / deviance(x, type = "null")
    }
    if (type == "cox_snell"){
        # deviance = -2 * logLik
        logLik0 <- - deviance(x, type = "null") / 2
        r2 <- 1 - exp((logLik0 - logL) * 2 / nobs(x))
    }
    if (type == "cragg_uhler"){
        # deviance = -2 * logLik
        logLik0 <- - deviance(x, type = "null") / 2
        r2 <- (1 - exp((logLik0 - logL) * 2 / nobs(x))) /
            (1 - exp(logLik0 * 2 / nobs(x)))
    }
    if (type == "aldrich_nelson"){
        logLik0 <- - deviance(x, type = "null") / 2
        r2 <- 2 * (logL - logLik0) /
            (2 * (logL - logLik0) + nobs(x))
    }
    if (type == "veall_zimm"){
        logLik0 <- - deviance(x, type = "null") / 2
        r2 <- 2 * (logL - logLik0) /
            (2 * (logL - logLik0) + nobs(x))
        r2 <- r2 * (2 * logLik0 - nobs(x)) / (2 * logLik0)
    }
    if (type == "estrella"){
        logLik0 <- - deviance(x, type = "null") / 2
        logLikx <- as.numeric(logL)
        r2 <- 1 - exp(-2 / nobs(x) * logLik0 * log(logLikx - logLik0))
        r2 <- 1 - (logLikx / logLik0) ^ (- 2 / nobs(x) * logLik0)
    }
    if (type == "cor"){
        r2 <- sum( (fitted(x) - yb) * (y - yb)) ^ 2 /
            sum((fitted(x) - yb) ^ 2) / sum((y - yb) ^ 2)
    }
    if (type == "rss"){
        r2 <- 1 - sum(resid(x, "response") ^ 2) / sum( (y - yb) ^ 2)
    }
    if (type == "ess"){
        r2 <- sum( (fitted(x) - yb) ^ 2) / sum( (y - yb) ^ 2)
    }
    if (type == "tjur"){
        r2_model <- sum( (fitted(x) - yb) ^ 2) / sum( (y - yb) ^ 2)
        r2_resid <- 1 - sum(resid(x, "response") ^ 2) / sum( (y - yb) ^ 2)
        r2 <- (r2_model + r2_resid) / 2
        r2
    }
    if (type == "mckel_zavo"){
        lp <- x$linear.predictors
        hlp <- mean(lp)
        ESS <- sum( (lp- hlp) ^ 2)
        .model <- x$call$model
        RSS <- nobs(x) * ifelse(.model == "probit", 1, pi ^ 2 / 3)
        r2 <- ESS / (ESS + RSS)
    }
    r2
}
