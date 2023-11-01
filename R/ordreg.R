#' Ordered regression
#'
#' A unified interface perform binomial regression calling `lm` to fit
#' the linear-probability model and `glm` to fit the probit and the
#' logit model.
#'
#' @name ordreg
#' @param formula a symbolic description of the model, (for the count
#'     component and for the selection equation).
#' @param data a data frame,
#' @param subset,weights,na.action,offset see `stats::lm`,
#' @param model one of `"lm"`, `"probit"` and "`logit`" to fit
#'     respectively the linear probability, the probit and the logit
#'     model
#' @param link one of `probit` and `logit`
#' @param start a vector of starting values, in this case, no
#'     estimation
#' @param object,type a `ordreg` object and the type of log-likelihood
#'     for the `logLik` method
#' @param ... further arguments
#' @return an object of class `micsr`, see `micsr::micsr` for further
#'     details
#' @importFrom stats glm plogis
#' @importFrom Formula Formula
#' @examples
#' mod1 <- ordreg(factor(dindx) ~ rhs1 + catchup, fin_reform, link = "logit")
#' @export
ordreg <- function(formula, data, weights, subset, na.action, offset,
                   link = c("probit", "logit"), start = NULL, ...){
    .call <- match.call()
    .link <- match.arg(link)
    cl <- match.call(expand.dots = FALSE)
    .formula <- cl$formula <- Formula(formula)
    m <- match(c("formula", "data", "subset", "weights"),
               names(cl), 0L)
    # construct the model frame and components
    cl <- cl[c(1L, m)]
    mf <- cl
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)
    X <- X[, - 1, drop = FALSE]
    y <- model.response(mf)
    J <- length(unique(y))
    nms_y <- levels(y)
    nms_thr <- paste(nms_y[- J], nms_y[ - 1], sep = "|")
    y <- as.integer(y)
    K <- ncol(X)
    N <- length(y)
    .df.residual <- N - K
    if (.link == "probit"){
        F_link <- pnorm
        f_link <- dnorm
        e_link <- function(x) - x * f_link(x)
    }
    if (.link == "logit"){
        F_link <- function(x) exp(x) / (1 + exp(x))
        f_link <- function(x) exp(x) / (1 + exp(x)) ^ 2
        e_link <- function(x) (exp(x) - exp(3 * x)) / (1 + exp(x)) ^ 4
    }
    lnl <- function(param, gradient = FALSE, hessian = FALSE, sum = TRUE){
        beta <- param[1L:K]
        mu <- c(-100, param[(K + 1L):(K + J - 1L)], 100)
        bX <- as.numeric(crossprod(t(X), beta))
        z1 <- mu[y + 1] - bX               ; z2 <- mu[y] - bX
        F1 <- F_link(z1)                   ; F2 <- F_link(z2)
        Li <- F1 - F2
        lnL <- log(Li)
        if (sum) lnL <- sum(lnL)
        if (gradient){
            W <- matrix(0, nrow(X), J)
            # The first two cols are not estimated coefficients (-infty and 0)
            # so remove them
            W1 <- (col(W) == (y + 1))[, - 1, drop = FALSE]  ; W2 <- (col(W) == y)[, - 1, drop = FALSE]
            f1 <- f_link(mu[y + 1] - bX)            ; f2 <- f_link(mu[y] - bX)
            gmu <- (W1 * f1 - W2 * f2)
            gb <- - (f1 - f2) * X
            gradi <- cbind(gb, gmu) / Li
            .gradient <- gradi
            if (sum) .gradient <- apply(gradi, 2, sum)
            attr(lnL, "gradient") <- .gradient
        }
        if (hessian){
            e1 <- e_link(mu[y + 1] - bX)              ; e2 <- e_link(mu[y] - bX)                     
            M1 <- cbind(- X, W1)
            M2 <- cbind(- X, W2)
            H <- crossprod(M1 * e1 / Li, M1) - crossprod(M2 * e2 / Li, M2) - crossprod(gradi)
            attr(lnL, "hessian") <- H
        }
        lnL
    }
    fy <- prop.table(table(y))
    fyc <- qnorm(cumsum(fy))
    sup.coef <- fyc[1:(J-1)]
    .start <- c(rep(0, K), sup.coef)
    names(.start) <- c(colnames(X), nms_thr)
    .coefs <- newton(lnl, coefs = .start, direction = "max", ...)
    .lp <- drop(X %*% .coefs[1:K])
    .lnl_conv <- lnl(.coefs, gradient = TRUE, hessian = TRUE, sum = FALSE)
    .fitted.values <- cbind(0, F_link(outer(- .lp, .coefs[K + (1L:(J - 1))], "+")), 1)
    .fitted.values <- .fitted.values[, - 1, drop = FALSE] - .fitted.values[, - (J + 1L), drop = FALSE]
    .df.residual <- N - K - J - 1L
    colnames(.fitted.values) <- nms_y
    .npar <- c(covariates = K, threshold = J - 1L)
    attr(.npar, "default") <- "covariates"
    result <- list(coefficients = .coefs,
                   model = mf,
                   gradient = attr(.lnl_conv, "gradient"),
                   hessian = attr(.lnl_conv, "hessian"),
                   linear.predictors = .lp,
                   logLik = as.numeric(.lnl_conv),
                   fitted.values = .fitted.values,
                   df.residual = .df.residual,
                   est_method = "ml",
                   formula = formula,
                   npar = .npar,
                   value = sum(as.numeric(.lnl_conv)),
                   call = .call
                   )
    structure(result, class = c("ordreg", "binomreg", "micsr"))
}



ordreg2 <- function(formula, data, weights, subset, na.action, offset,
                   link = c("probit", "logit"), start = NULL, ...){
    .call <- match.call()
    .link <- match.arg(link)
    cl <- match.call(expand.dots = FALSE)
    .formula <- cl$formula <- Formula(formula)
    m <- match(c("formula", "data", "subset", "weights"),
               names(cl), 0L)
    # construct the model frame and components
    cl <- cl[c(1L, m)]
    mf <- cl
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)
    X <- X[, - 1, drop = FALSE]
    y <- model.response(mf)
    if (inherits(y, "Surv")){
        cens <- unname(y[, 2]) == 0
        y <- unname(y[, 1])
    } else cens <- rep(0, length(y))
    if (! inherits(y, "factor")) y <- factor(y)
    J <- length(unique(y))
    nms_y <- levels(y)
    nms_thr <- paste(nms_y[- J], nms_y[ - 1], sep = "|")
    y <- as.integer(y)
    K <- ncol(X)
    N <- length(y)
    .df.residual <- N - K
    if (.link == "probit"){
        F_link <- pnorm
        f_link <- dnorm
        e_link <- function(x) - x * f_link(x)
    }
    if (.link == "logit"){
        F_link <- function(x) exp(x) / (1 + exp(x))
        f_link <- function(x) exp(x) / (1 + exp(x)) ^ 2
        e_link <- function(x) (exp(x) - exp(3 * x)) / (1 + exp(x)) ^ 4
    }
    lnl <- function(param, gradient = FALSE, hessian = FALSE, sum = TRUE){
        beta <- param[1L:K]
        mu <- c(-100, param[(K + 1L):(K + J - 1L)], 100)
        bX <- as.numeric(crossprod(t(X), beta))
        z1 <- mu[y + 1] - bX               ; z2 <- mu[y] - bX
        F1 <- F_link(z1)                   ; F2 <- F_link(z2)
        q <- (1 - 2 * cens)
        u <- (1 - cens)
        Li <- q * F1 - u * F2 + (1 - u)
        Li <- pmax(1E-06, Li)
        lnL <- log(Li)
        if (sum) lnL <- sum(lnL)
        if (gradient | hessian){
            W <- matrix(0, nrow(X), J)
            # The first two cols are not estimated coefficients (-infty and 0)
            # so remove them
            W1 <- (col(W) == (y + 1))[, - 1, drop = FALSE]  ; W2 <- (col(W) == y)[, - 1, drop = FALSE]
            f1 <- f_link(mu[y + 1] - bX)                    ; f2 <- f_link(mu[y] - bX)
            gmu <- q * W1 * f1 - u * W2 * f2
            gb <- - (q * f1 - u * f2) * X
            gradi <- cbind(gb, gmu) / Li
            colnames(gradi) <- c(colnames(X), nms_thr)
            .gradient <- gradi
            if (sum) .gradient <- apply(gradi, 2, sum)
            attr(lnL, "gradient") <- .gradient
        }
        if (hessian){
            e1 <- e_link(mu[y + 1] - bX)              ; e2 <- e_link(mu[y] - bX)
            M1 <- cbind(- X, W1)
            M2 <- cbind(- X, W2)
            H <- crossprod(q * M1 * e1 / Li, M1) - crossprod(u * M2 * e2 / Li, M2) - crossprod(gradi)
            colnames(H) <- rownames(H) <- c(colnames(X), nms_thr)
            attr(lnL, "hessian") <- H
        }
        lnL
    }
    fy <- prop.table(table(y))
    fyc <- qnorm(cumsum(fy))
    sup.coef <- fyc[1:(J-1)]
    .start <- c(rep(0, K), sup.coef)
    names(.start) <- c(colnames(X), nms_thr)
    if (FALSE){
        # gradient
        numgrad <- numDeriv::grad(lnl, .start)
        angrad <- attr(lnl(.start, gradient = TRUE, sum = TRUE), "gradient")
        print(cbind(numgrad, angrad, za = numgrad - angrad))
        # hessian
        numhess <- numDeriv::hessian(lnl, .start)
        anhess <- attr(lnl(.start, hessian = TRUE, sum = TRUE), "hessian")
        colnames(numhess) <- rownames(numhess) <- colnames(anhess)
        print(numhess)
        print(anhess)
        print(summary(abs(as.numeric(numhess - anhess))))
    stop()
    }
    za_f <- function(x) -lnl(x, sum = TRUE)
    za_g <- attr(function(x) - lnl(x, gradient = TRUE, sum = TRUE), "gradient")
    za_h <- attr(function(x) lnl(x, hessian = TRUE), "hessian")
    za <- optim(.start, za_f, za_g, method = "BFGS", hessian = T);
    print(za)
    print(solve(-za$hessian))
    stop()
    .coefs <- newton(lnl, coefs = .start, direction = "max", ...)
    stop("lebon")
    .lp <- drop(X %*% .coefs[1:K])
    .lnl_conv <- lnl(.coefs, gradient = TRUE, hessian = TRUE, sum = FALSE)
    .fitted.values <- cbind(0, F_link(outer(- .lp, .coefs[K + (1L:(J - 1))], "+")), 1)
    .fitted.values <- .fitted.values[, - 1, drop = FALSE] - .fitted.values[, - (J + 1L), drop = FALSE]
    .df.residual <- N - K - J - 1L
    colnames(.fitted.values) <- nms_y
    .npar <- c(covariates = K, threshold = J - 1L)
    attr(.npar, "default") <- "covariates"
    result <- list(coefficients = .coefs,
                   model = mf,
                   gradient = attr(.lnl_conv, "gradient"),
                   hessian = attr(.lnl_conv, "hessian"),
                   linear.predictors = .lp,
                   logLik = as.numeric(.lnl_conv),
                   fitted.values = .fitted.values,
                   df.residual = .df.residual,
                   est_method = "ml",
                   formula = formula,
                   npar = .npar,
                   value = sum(as.numeric(.lnl_conv)),
                   call = .call
                   )
    structure(result, class = c("ordreg", "binomreg", "micsr"))
}



#' @rdname ordreg
#' @export
logLik.ordreg <- function(object, ..., type = c("model", "null")){
    .type <- match.arg(type)
    if (.type == "model")
        result <- structure(object$value, nobs = nobs(object), df = sum(object$npar), class = "logLik")
    if (.type == "null"){
        yb <- prop.table(table(model.response(model.frame(object))))
        J <- length(yb)
        result <- nobs(object) * sum(yb * log(yb))
        result <- structure(result, nobs = nobs(object), df = J - 1L, class = "logLik")
    }
    result
}

if (FALSE){
    library("tidyverse")
    load("../../katz.rda")
    katz <- katz %>% mutate(dur40 = pmin(duration, 40), out40 = outcome ,dur2 = ceiling(duration / 20) * 20)
    katz$out40[katz$duration > 40] <- "censored"
}
