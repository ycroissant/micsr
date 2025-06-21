#' `micsr` class
#'
#' The `micsr` class is intend to deal with a lot of different models
#' that are estimated in the `micsr` package. More specifically, some
#' models may be estimated using different estimation methods, like
#' maximum likelihood, GMM or two-steps estimators. Objects of class
#' `micsr` have an `est_method` item which is used by the different
#' methods in order to have a relevent behaviour for the different
#' methods.
#'
#' @name micsr
#' @param x,object an object which inherits the `micsr` class
#' @param formula a formula
#' @param subset,grep,fixed,invert,coef invert see `micsr::select_coef
#' @param vcov the method used to compute the covariance matrix of the
#'     estimators (only for the ML estimator), one of `hessian` (the
#'     opposite of the inverse of the hessian), `info` (the inverse of
#'     the opposite of the expected value of the hessian), `opg` (the
#'     outer product of the gradient)
#' @param digits,width see `print`
#' @param conf.int,conf.level see `broom:tidy.lm`
#' @param lhs,rhs see `Formula::model.frame.Formula`
#' @param type,omega,sandwich see `sandwich::sandwich`
#' @param covariates a set of covariates for the `effects` method,
#' @param se whether the standard errors sould be computed for
#'     predictions and slopes
#' @param shape the shape of the predictions for `mlogit` objects
#' @param newdata a new data frame to compute the predictions
#' #' @param se a boolean indicating whether the standard errors should
#'     be computed
#' @param k see `AIC`
#' @param sum return either the sum of the contributions or the vector
#'     of contribution
#' @param ... further arguments
#' @return
#' 
#' Objects of class `micsr` share a lot of common elements with `lm`:
#' `coefficients`, `residuals`, `fitted.values`, `model`, `terms`,
#' `df.residual`, `xlevels`, `na.action`, and `call`. `npar` is a
#' named vector containing the index of subset of coefficients, it is
#' used to print a subset of the results.  It also has a `est_method`
#' element and, depending of its value, contains further elements. In
#' particular, for model fitted by maximum likelihood, `value`
#' contains the individual contribution to the log-likelihood
#' function, `gradient` the individual contribution to the gradient,
#' `hessian` the hessian and `information` the information
#' matrix. `logLik` contains the log-likelihood values of the
#' proposed, null and saturated models. `tests` contains the values of
#' the test that all the coefficients of the covariates are 0, using
#' the three classical tests.
#'
#' The `llobs` function is provided as a generic to extract the
#' individual contributions to the log-likelihood
#' 
#' Specific methods have been writen for `micsr` objects: `nobs`,
#' `generics::tidy`, `generics::glance`, `sandwich::meat`,
#' `sandwich::estfun`, `predict`, `model.matrix`,
#' `Formula::model.part`.
#'
#' `logLik`, `BIC`, `AIC` and `deviance` methods have a `type`
#' argument to select theproposed, null or saturated model.
#'
#' `vcov` and `summary` methods have a `vcov` argument to select the
#' estimator of the covariance matrix, which can be either based on
#' the hessian, the gradient or the information.
#'
#' `vcov`, `summary` and `coef` have a subset argument to select only
#' a subset of the coefficients
#' @importFrom stats nobs deviance BIC AIC logLik predict deviance
#' @importFrom sandwich bread estfun vcovHC sandwich
#' @importFrom Formula model.part
#' @useDynLib micsr, .registration=TRUE
NULL

# @importFrom Rcpp evalCpp

#' @importFrom survival Surv
#' @export
survival::Surv

## #' @importFrom dplyr mutate
## #' @export
## dplyr::mutate

## #' @importFrom magrittr %>%
## #' @export
## magrittr::`%>%`

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom sandwich bread
#' @export
sandwich::bread

#' @importFrom sandwich estfun
#' @export
sandwich::estfun

#' @importFrom sandwich vcovHC
#' @export
sandwich::vcovHC


#' @rdname micsr
#' @export
llobs <- function(x, ...)
    UseMethod("llobs")

#' @importFrom Formula model.part
#' @export
Formula::model.part

pretty_nms <- function(x, subset = NA){
    .subset <- subset
    .cov_subsets <- c("covariates", "heterosc", "instruments")
    .multi <- length(intersect(.subset, .cov_subsets)) > 1
    if (! .multi){
        if ("covariates" %in% .subset) x[grep("covar_", x)] <- substr(x[grep("covar_", x)], 7, 100)
        if ("instruments" %in% .subset) x[grep("instr_", x)] <- substr(x[grep("instr_", x)], 7, 100)
    }
    x
}

#' @rdname micsr
#' @export
coef.micsr <- function(object, ..., subset = NA, fixed = FALSE,
                       grep  = NULL, invert = FALSE, coef = NULL){
    .sel <- select_coef(object, subset = subset, fixed = fixed, grep = grep,
                        coef = coef, invert = invert)
    .coef <- object$coefficients[.sel]
    names(.coef) <- pretty_nms(names(.coef), subset)
    .coef
}

#' @rdname micsr
#' @export
vcov.micsr <- function(object, ..., vcov = NULL, subset = NA, fixed = FALSE,
                       grep = NULL, invert = FALSE, coef = NULL){
    .vcov_method <- vcov
    .est_method <- object$est_method
    if (.est_method == "ml"){
        if (is.null(.vcov_method)){
            if (! is.null(object$info)){
                .vcov_method = "hessian"
            } else {
                if (! is.null(object$hessian)){
                    .vcov_method = "hessian"
                } else {
                    .vcov_method = "opg"
                }
            }
        } else {
            if (! .vcov_method %in% c("info", "hessian", "opg")){
                stop("irrelevant method")
            }
        }
        if (.vcov_method == "info"){
                if (is.null(object$info)){
                    stop("no information matrix available")
                } else {
                    .vcov <- object$info
                }
        }
        if (.vcov_method == "hessian"){
            if (is.null(object$hessian)){
                stop("no hessian matrix available")
            } else {
                .vcov <- - object$hessian
            }
        }
        if (.vcov_method == "opg") .vcov <- crossprod(object$gradient)
    } else {
        .vcov <- object$vcov
    }
    nms <- rownames(.vcov)
    .sel <- select_coef(object, subset = subset, fixed = fixed,
                        grep = grep, coef = coef, invert = invert)
    nms <- nms[.sel]
    .vcov <- .vcov[.sel, .sel, drop = FALSE]
    colnames(.vcov) <- rownames(.vcov) <- pretty_nms(nms, subset)
    solve(.vcov)
}

#' @rdname micsr
#' @export
summary.micsr <- function(object, ...,
                          vcov = c("hessian", "info", "opg"),
                          subset = NA, fixed = FALSE, grep = NULL, invert = FALSE, coef = NULL){

    .est_method <- object$est_method
    .vcov_method <- match.arg(vcov)
    .vcov <- vcov(object, subset = subset, vcov = .vcov_method,
                  fixed = fixed, grep = grep, coef = coef, invert = invert)
    std.err <- sqrt(diag(.vcov))
    b <- coef(object, subset = subset, fixed = fixed,
              grep = grep, coef = coef, invert = invert)
    z <- b / std.err
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    object$coefficients <- cbind(b, std.err, z, p)
    colnames(object$coefficients) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
    if (.est_method == "gmm"){
        object$sargan <- sargan(object)
    }
    structure(object, class = c("summary.micsr", "micrs"))
}

#' @rdname micsr
#' @export
coef.summary.micsr <- function(object, ...) object$coefficients

#' @rdname micsr
#' @export
print.micsr <- function (x, digits = max(3L, getOption("digits") - 3L), ...){
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L, 
                  quote = FALSE)
    cat("\n")
    invisible(x)
}

#' @rdname micsr
#' @export
print.summary.micsr <- function (x, digits = max(3, getOption("digits") - 2), width = getOption("width"), ...){
    .est_method <- x$est_method
    if (.est_method == "trimmed"){
        .est_method_lib <- "Trimmed estimation"
        gof_stat <- NULL
    }
    if (.est_method == "ml"){
        .est_method_lib <- "Maximum likelihood estimation"
        gof_stat <- "log-Likelihood"
    }
    if (.est_method == "lm"){
        .est_method_lib <- "Ordinary Least Squares"
        gof_stat <- NA
    }

    if (.est_method == "twosteps"){
        .est_method_lib <- "Two-steps estimation"
        gof_stat <- "deviance"
    }
    if (.est_method == "gmm"){
        .est_method_lib <- "General Method of Moments"
        gof_stat <- "quadratic form of the moments"
    }
    if (.est_method == "ols"){
        .est_method_lib <- "Ordinary Least Squares"
        gof_stat <- "Sum of Squared Residuals"
    }
    if (.est_method == "twostep"){
        .est_method_lib <- "Two-steps"
        gof_stat <- NA
    }
    if (.est_method == "minchisq"){
        .est_method_lib <- "Minimum chi-squared"
        gof_stat <- NA
    }
    
    cat(.est_method_lib, "\n", sep = "")
    printCoefmat(coef(x), digits = digits)
    if (! is.null(gof_stat)){
        cat("\n", gof_stat, ": ", format(x$logLik["model"], digits = digits), "\n\n", sep = "")
    }
    if (.est_method == "twosteps"){
        hrho <- x$coefficients
        cat("Estimated value of sigma: ", format(x$sigma, digits = digits), "\n", sep = "")
        cat("Implied value for rho   : ", format(x$rho, digits = digits), "\n", sep = "")
    }
    if (.est_method == "gmm"){
        print(x$sargan)
    }
}

#' @rdname micsr
#' @export
logLik.micsr <- function(object, ..., type = c("model", "null", "saturated"), sum = TRUE){
    .type <- match.arg(type)
    if (sum){
        .val <- object$logLik[.type]
        if (is.na(.val)) stop(paste("the ", .type, " log-likelihood is not available", sep = ""))
        .nobs <- nobs(object)
        .df <- switch(.type,
                      model = npar(object),
                      null = 1,
                      saturated = nobs(object))
        structure(.val, nobs = .nobs, df = .df, class = "logLik")
    }
    else object$value[, .type]
}



## logLik.micsr <- function(object, ..., type = c("model", "null", "saturated")){
##     .type <- match.arg(type)
##     .val <- object$logLik[.type]
##     if (is.na(.val)) stop(paste("the ", .type, " log-likelihood is not available", sep = ""))
##     .nobs <- nobs(object)
##     .df <- switch(.type,
##                   model = npar(object),
##                   null = 1,
##                   saturated = nobs(object))
##     structure(.val, nobs = .nobs, df = .df, class = "logLik")
## }



#' @rdname micsr
#' @export
BIC.micsr <- function(object, ..., type = c("model", "null")){
    if (object$est_method != "ml") NULL
    else{
        .type <- match.arg(type)
        ll <- logLik(object, type = .type)
        npar <- attr(ll, "df")
        N <- nobs(object)
        result <- - 2  *  as.numeric(ll) + log(N) * npar
    }
    result
}

#' @rdname micsr
#' @export
AIC.micsr <- function(object, ..., k = 2, type = c("model", "null")){
    if (object$est_method != "ml") NULL
    else{
        .type <- match.arg(type)
        ll <- logLik(object, type = .type)
        npar <- attr(ll, "df")
        N <- nobs(object)
        result <- - 2  *  as.numeric(ll) + k * npar
    }
    result
}

#' @rdname micsr
#' @export
deviance.micsr <- function(object, ..., type = c("model", "null")){
    .type <- match.arg(type)
    logLik_saturated <- unname(ifelse(is.na(object$logLik["saturated"]), 0,
                               object$logLik["saturated"]))
    - 2 * (unname(object$logLik[.type]) - logLik_saturated)
}


# !!!!! uniquement pour escount

## #' @rdname micsr
## #' @export
## predict.micsr <- function(object, ..., newdata = NULL){
##     .est_method <- object$est_method
##     if (is.null(newdata)) fitted(object)
##     else{
## #        .formula <- Formula(formula(paste(deparse(object$call$formula), collapse = "")))
##         .formula <- Formula(formula(object$terms))
##         mf <- model.frame(.formula, newdata, dot = "previous")
##         X <- model.matrix(.formula, mf, rhs = 1)
##         Z <- model.matrix(.formula, mf, rhs = 2)
##         d <- model.part(.formula, newdata, lhs = 2, drop = TRUE)
##         q <- 2 * d - 1
##         K <- ncol(X)
##         L <- ncol(Z)
##         if (.est_method == "ml"){
##             beta <- coef(object, subset = "covariates")
##             alpha <- coef(object, subset = "instruments")
##             theta <- unname(prod(coef(object, subset = "vcov")))
##         }
##         else{
##             beta <- object$coefficients[1:K]
##             theta <- object$coefficients[K + 1]
##             alpha <- coef(object$first)
##         }
##         aZ <- drop(Z %*% alpha)
##         bX <- drop(X %*% beta)
##         exp(bX + pnorm(q * (theta + aZ), log.p = TRUE) - pnorm(q * aZ, log.p = TRUE))
##     }
## }

## #' @rdname micsr
## #' @export
## predict.binomreg <- function(object, type = c("response", "link"), ..., newdata = NULL){
##     .type <- match.arg(type)
    
    

#' @rdname micsr
#' @method model.part micsr
#' @export
model.part.micsr <- function(object, ..., lhs = 1){
    .formula <- Formula(formula(paste(deparse(object$call$formula), collapse = "")))
    model.part(.formula, object$model, lhs = lhs, drop = TRUE, dot = "previous")
}


## model.matrix.micsr <- function(object, formula = NULL, ..., rhs = 1){
## #    .formula <- Formula(formula(paste(deparse(object$call$formula), collapse = "")))
##     if (is.null(formula)) .formula <- object$formula
##     else{ .formula <- Formula(formula)
##         .formula <- update(object$formula, formula)
##         print(.formula)
##     }
##     model.matrix(.formula, object$model, rhs = rhs, dot = "previous")
## }


#' @rdname micsr
#' @method model.matrix micsr
#' @export
model.matrix.micsr <- function(object, formula = NULL, ..., rhs = 1){
#    .formula <- Formula(formula(paste(deparse(object$call$formula), collapse = "")))
    ## if (is.null(formula)) .formula <- object$formula
    ## else{ .formula <- Formula(formula)
    ##     .formula <- update(object$formula, formula)
    ## }
    if (is.null(formula)) .formula <- object$terms
    else{.formula <- Formula(formula)
        .formula <- update(object$terms, formula)
    }
    cl <- object$call
    cl$formula <- .formula
    m <- match(c("formula", "data", "subset", "weights"),
               names(cl), 0L)
    cl <- cl[c(1L, m)]
    mf <- cl
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    X <- model.matrix(.formula, mf, rhs = rhs, dot = "previous")
    X
}

#' @rdname micsr
#' @method estfun micsr
#' @export
estfun.micsr <- function(x, ...){
    if (x$est_method == "twosteps") stop("no estfun method for two-steps models")
    x$gradient
}

#' @rdname micsr
#' @method vcovHC micsr
#' @export
vcovHC.micsr <- function(x, type, omega = NULL, sandwich = TRUE, ...){
    if (x$est_method == "twosteps") stop("no meat method for two-steps models")
    .meat <- crossprod(estfun(x)) / nobs(x)
    .bread <- bread(x)
    sandwich(x, bread. = .bread, meat. = .meat, ...)
}

#' @rdname micsr
#' @method bread micsr
#' @export
bread.micsr <- function(x, ...){
    if (x$est_method == "twosteps") stop("no meat method for two-steps models")
#    if (! is.null(x$info)) solve(x$info) * nobs(x)
#    else solve(- x$hessian) * nobs(x)
    solve(- x$hessian) * nobs(x)
}

#' @rdname micsr
#' @export
nobs.micsr <- function(object, ...){
    if (! is.null(object$model)) nrow(object$model)#length(object$residuals)
    else nrow(object$gradient)
}

#' @rdname micsr
#' @method llobs micsr
#' @export
llobs.micsr <- function(x, ...) x$value


#' @rdname micsr
#' @method llobs mlogit
#' @export
llobs.mlogit <- function(x, ...) log(fitted(x, outcome = TRUE))

#' @rdname micsr
#' @method tidy micsr
#' @export
tidy.micsr <- function(x, conf.int = FALSE, conf.level = 0.95, ...){
    result <- summary(x)$coefficients
    nms <- rownames(result)
    rownames(result) <- NULL
    result <- data.frame(term = nms, result)
    names(result) <- c("term", "estimate", "std.error", "statistic", "p.value")
    result
}

#' @rdname micsr
#' @method glance micsr
#' @export
glance.micsr <- function(x, ...){
    .est_method <- x$est_method
    N <- nobs(x)
    result <- data.frame(nobs = nobs(x))
    if (.est_method == "ml")
        result <- data.frame(nobs = nobs(x), logLik = logLik(x), AIC = AIC(x), BIC = BIC(x))
#    if (.est_method == "twosteps"){
#        result <- data.frame(nobs = nobs(x), impliedsigma = x$sigma)
#    }
    result
}

#' @rdname micsr
#' @method residuals micsr
#' @export
residuals.micsr <- function(object, ..., type = c("deviance", "pearson", "response")){
    type <- match.arg(type)
    y <- model.response(model.frame(object))
    hy <- fitted(object)
    wt <- model.weights(model.frame(object))
    if (is.null(wt)) wt <- rep(1, length(y))
    if (inherits(object, "binomreg")) sd_hy <- sqrt(hy * (1 - hy))
    if (inherits(object, "poisreg")) sd_hy <- sqrt(hy)
    if (type == "response") resid <- y - hy
    if (type == "pearson") resid <- (y - hy) / sd_hy * sqrt(wt)
    if (type == "deviance"){
        sgn <- ifelse(y > hy, 1, - 1)
        resid <- sgn * sqrt(2) * sqrt(wt) * 
            sqrt(logLik(object, type = "saturated", sum = FALSE) -
                 logLik(object, type = "model", sum = FALSE))
    }
    resid
}


residuals.micsr <- function(object, ..., type = c("deviance", "pearson", "response")){
    .type <- match.arg(type)
    y <- model.response(model.frame(object))
    mu <- fitted(object)
    wt <- model.weights(model.frame(object))
    if (is.null(wt)) wt <- rep(1, length(y))
    if (.type == "response") .resid <- y - mu
    if (.type == "pearson") .resid <- (y - mu) / sqrt(object$family$variance(mu)) * sqrt(wt)
    if (.type == "deviance") .resid <- sqrt(2) * ifelse(y > mu, 1, -1) * sqrt(wt) * 
                                 (object$value[, "saturated"] - object$value[, "model"])  ^ .5
    .resid
}

    
