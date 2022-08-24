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
    q <- 2 * y - 1
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
    .coefs <- newton(lnl, coefs = .start, direction = "max", trace = 0)
    .lp <- drop(X %*% .coefs[1:K])
    .lnl_conv <- lnl(.coefs, gradient = TRUE, hessian = TRUE, sum = FALSE)
    .fitted.values <- cbind(0, pnorm(outer(- .lp, .coefs[K + (1L:(J - 1))], "+")), 1)
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
