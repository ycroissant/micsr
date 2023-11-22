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
#' @param x,object an object which inherits the `micsr` class,
#' @param formula a formula
#' @param subset a character which indicates which subset of
#'     coefficients should be extracted: one of `noinst` (all the
#'     coefficients except those corresponding to instrumental
#'     variables), `all`, `covar` (only the coefficients of the
#'     covariates), `inst` (only the coefficients of the instrumental
#'     variables) and `misc` (ony the "miscelanous" coefficients,
#'     typicaly a standard deviation or a coefficient of correlation),
#' @param vcov the method used to compute the covariance matrix of the
#'     estimators (only for the ML estimator), one of `hessian` (the
#'     opposite of the inverse of the hessian), `info` (the inverse of
#'     the opposite of the expected value of the hessian), `opg` (the
#'     outer product of the gradient)
#' @param digits,width see `base::print`
#' @param conf.int,conf.level see `broom:tidy.lm`
#' @param lhs,rhs see `Formula::model.frame.Formula`
#' @param type,omega,sandwich see `sandwich::sandwich`
#' @param newdata a new data frame to compute the predictions
#' @param k see `stats::AIC`
#' @param ... further arguments
#' @importFrom stats nobs deviance BIC AIC logLik predict deviance
#' @importFrom sandwich bread estfun vcovHC sandwich
#' @importFrom Formula model.part
#' @importFrom Rcpp evalCpp
#' @useDynLib micsr, .registration=TRUE
NULL

#' @importFrom dplyr mutate
#' @export
dplyr::mutate

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


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

#' @importFrom nonnest2 llcont
#' @export
nonnest2::llcont

#' @importFrom Formula model.part
#' @export
Formula::model.part

#select_coef <- function(object, which = c("noinst", "all", "covar", "inst", "misc")){ 
select_coef <- function(object, subset = NA){   
    # ancillary  : instruments, heteroscedasticity
    # covariates : 
    # misc       : standard deviations / coefficients of correlation / cholesky matrix
    # vcov       : - ancillary
    # all        : all coefficients
    .npar <- object$npar
    if (is.null(.npar) | is.null(attr(.npar, "default"))) idx <- 1:length(object$coefficients)
    else{
        if (length(subset) == 1){
            if (is.na(subset)) .subset = attr(.npar, "default")
            else{
                if (subset == "all") .subset <- names(.npar)
                else{
                    if (subset %in% names(.npar)) .subset <- subset
                    else stop("irrelevant subset argument")
                }
            }
        }
        else{
            .subset <- subset
            if (! all(.subset %in% names(.npar))) stop("irrelevant subset")
        }
        idx <- data.frame(subset = rep(names(.npar), times = .npar),
                          idx = 1:sum(.npar))
        .sel <- .subset
        idx <- subset(idx, subset %in% .sel)[["idx"]]
    }
    idx
}

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
coef.micsr <- function(object, ..., subset = NA){
    .subset <- subset
    is.na_subset <- (length(.subset) == 1) && is.na(.subset)
    if (is.na_subset) .subset <- attr(object$npar, "default")
    .est_method <- object$est_method
    .sel <- select_coef(object, .subset)
    .coef <- object$coefficients[.sel]
    names(.coef) <- pretty_nms(names(.coef), .subset)
    .coef
}

#' @rdname micsr
#' @export
vcov.micsr <- function(object, ..., vcov = c("info", "hessian", "opg"), subset = NA){
    .subset <- subset
    is.na_subset <- (length(.subset) == 1) && is.na(.subset)
    if (is.na_subset) .subset <- attr(object$npar, "default")
    .est_method <- object$est_method
    if (.est_method == "ml"){
        .vcov_method <- match.arg(vcov)
        if (.vcov_method == "info" & is.null(object$info)) .vcov_method = "hessian"
        if (.vcov_method == "hessian") .vcov <- solve(- object$hessian)
        if (.vcov_method == "opg") .vcov <- solve(crossprod(object$gradient))
        if (.vcov_method == "info") .vcov <- object$info
    }
    else .vcov <- object$vcov
    nms <- rownames(.vcov)
    .sel <- select_coef(object, subset = .subset)
    nms <- nms[.sel]
    .vcov <- .vcov[.sel, .sel, drop = FALSE]
    colnames(.vcov) <- rownames(.vcov) <- pretty_nms(nms, .subset)
    .vcov
}

#' @rdname micsr
#' @export
summary.micsr <- function(object, ...,
                          vcov = c("hessian", "info", "opg"),
                          subset = NA){
    .subset <- subset
    is.na_subset <- (length(.subset) == 1) && is.na(.subset)
    if (is.na_subset) .subset <- attr(object$npar, "default")
    .est_method <- object$est_method
    .vcov_method <- match.arg(vcov)
    .vcov <- vcov(object, subset = .subset, vcov = .vcov_method)
    std.err <- sqrt(diag(.vcov))
    b <- coef(object, subset = .subset)
    z <- b / std.err
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    object$coefficients <- cbind(b, std.err, z, p)
    colnames(object$coefficients) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
    if (.est_method == "gmm"){
        object$sargantest <- sargantest(object)
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
    if (.est_method == "ml"){
        .est_method_lib <- "Maximum likelihood estimation"
        gof_stat <- "log-Likelihood"
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
    cat("\n", gof_stat, ": ", format(x$value, digits = digits), "\n\n", sep = "")
    if (.est_method == "twosteps"){
        hrho <- x$coefficients
        cat("Estimated value of sigma: ", format(x$sigma, digits = digits), "\n", sep = "")
        cat("Implied value for rho   : ", format(x$rho, digits = digits), "\n", sep = "")
    }
    if (.est_method == "gmm"){
        print(x$sargantest)
    }
}

#' @rdname micsr
#' @export
logLik.micsr <- function(object, ...){
    if (object$est_method != "ml") NULL
    else object$value
}

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
    - 2 * as.numeric(logLik(object, type = .type))
}

#' @rdname micsr
#' @export
predict.micsr <- function(object, ..., newdata = NULL){
    .est_method <- object$est_method
    if (is.null(newdata)) fitted(object)
    else{
#        .formula <- Formula(formula(paste(deparse(object$call$formula), collapse = "")))
        .formula <- Formula(formula(object$terms))
        mf <- model.frame(.formula, newdata, dot = "previous")
        X <- model.matrix(.formula, mf, rhs = 1)
        old_options <- options(warn = -1)
        Z <- model.matrix(.formula, mf, rhs = 2)
        options(old_options)
        d <- model.part(.formula, newdata, lhs = 2, drop = TRUE)
        q <- 2 * d - 1
        K <- ncol(X)
        L <- ncol(Z)
        if (.est_method == "ml"){
            beta <- coef(object, subset = "covariates")
            alpha <- coef(object, subset = "instruments")
            theta <- unname(prod(coef(object, subset = "vcov")))
        }
        else{
            beta <- object$coefficients[1:K]
            theta <- object$coefficients[K + 1]
            alpha <- coef(object$first)
        }
        aZ <- drop(Z %*% alpha)
        bX <- drop(X %*% beta)
        exp(bX + pnorm(q * (theta + aZ), log.p = TRUE) - pnorm(q * aZ, log.p = TRUE))
    }
}

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
    if (is.null(formula)) .formula <- object$formula
    else{ .formula <- Formula(formula)
        .formula <- update(object$formula, formula)
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
    if (! is.null(x$info)) x$info * nobs(x)
    else solve(- x$hessian) * nobs(x)
}

#' @rdname micsr
#' @export
nobs.micsr <- function(object, ...) nrow(object$model)#length(object$residuals)

#' @rdname micsr
#' @method llcont micsr
#' @export
llcont.micsr <- function(x, ...) x$logLik


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

    
        
