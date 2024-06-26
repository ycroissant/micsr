#' Ordered regression
#'
#' Maximum-likelihood estimation of a model for which the response is
#' ordinal
#' 
#' @name ordreg
#' @param formula a symbolic description of the model
#' @param data a data frame
#' @param subset,weights,na.action,offset see `lm`
#' @param link one of `probit` and `logit`
#' @param start a vector of starting values, in this case, no
#'     estimation
#' @param method optimization method
#' @param ... further arguments
#' @return an object of class `micsr`, see `micsr::micsr` for further
#'     details.
#' @importFrom stats glm plogis
#' @importFrom Formula Formula
#' @keywords models
#' @examples
#' mod1 <- ordreg(factor(dindx) ~ rhs1 + catchup, fin_reform, link = "logit")
#' library(survival)
#' ud <- unemp_duration %>%
#'       mutate(years = floor(duration / 365),
#'              years = ifelse(years == 6, 5, years))
#' mod2 <- ordreg(Surv(years, censored == "no") ~ gender + age + log(1 + wage), ud,
#'                link = "cloglog", method = "bfgs")
#' @export
ordreg <- function(formula, data, weights, subset, na.action, offset,
                   link = c("probit", "logit", "cloglog"),
                   method = c("bfgs", "nr"), start = NULL,  ...){
    .method <- match.arg(method)
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
    y <- model.response(mf)
    if (inherits(y, "Surv")){
        y <- as.matrix(model.response(mf))
        e <- as.logical(y[, 2])
        y <- as.factor(y[, 1])
    }
    else {
        e <- rep(1, length(y))
        if (! inherits(y, "factor")) y <- as.factor(y)
    }
    X <- model.matrix(mt, mf)
    X <- X[, - 1, drop = FALSE]
    if (compute_rank(X) < ncol(X)){
        .rank <- compute_rank(X)
        .ncol <- ncol(X)
        stop(paste("the rank of X = ", .rank, " < the number of columns = ", .ncol, sep = ""))
    }
    J <- length(unique(y))
    nms_y <- levels(y)
    nms_thr <- paste(nms_y[- J], nms_y[ - 1], sep = "|")
    y <- as.integer(y)
    Y <- model.matrix(~ factor(y) - 1)
    K <- ncol(X)
    N <- length(y)
    # the observations for the last class should be coded as
    # uncensored
    e[y == J] <- TRUE
    .df.residual <- N - K
    if (.link == "probit"){
        F_link <- pnorm
        f_link <- dnorm
        q_link <- qnorm
        e_link <- function(x) - x * f_link(x)
    }
    if (.link == "logit"){
        F_link <- function(x) exp(x) / (1 + exp(x))
        f_link <- function(x) exp(x) / (1 + exp(x)) ^ 2
        e_link <- function(x) (exp(x) - exp(3 * x)) / (1 + exp(x)) ^ 4
        q_link <- function(x) log(x / (1 - x))
    }
    if (.link == "cloglog"){
        F_link <- function(x) 1 - exp(- exp(x))
        f_link <- function(x) exp(x) * exp(- exp(x))
        e_link <- function(x) exp(-exp(x)) * (exp(x) - exp(2 * x))
        q_link <- function(x) log(- log(1 - x))
    }

    lnl <- function(param, gradient = FALSE, hessian = FALSE, sum = TRUE,
                    information = FALSE, X = X, y = y, e = e, opposite = FALSE){
        sgn <- ifelse(opposite, - 1, 1)
        beta <- param[1L:K]
        mu <- c(-100, param[(K + 1L):(K + J - 1L)], 100)
        bX <- as.numeric(crossprod(t(X), beta))
        z1 <- mu[y + 1] - bX               ; z2 <- mu[y] - bX
        F1 <- F_link(z1)                   ; F2 <- F_link(z2)
        Ei <- F1 - F2
        Ci <- 1 - F1
        Li <- ifelse(e == 1, Ei, Ci)
        lnL <- sgn * sum(log(Li))
        if (sum) lnL <- sum(lnL)
        if (gradient | hessian){
            W <- matrix(0, nrow(X), J)
            # The first two cols are not estimated coefficients (-infty and 0)
            # so remove them
            W1 <- (col(W) == (y + 1))[, - 1, drop = FALSE]  ; W2 <- (col(W) == y)[, - 1, drop = FALSE]
            f1 <- f_link(mu[y + 1] - bX)                    ; f2 <- f_link(mu[y] - bX)
            W1c <- matrix(FALSE, nrow(X), J - 1) ; W1c[, J - 1] <- TRUE ; W2c <- W1
            f1c <- rep(0, length(y))                        ; f2c <- f1
            W1 <- e * W1 + (1 - e) * W1c                    ; W2 <- e * W2 + (1 - e) * W2c
            f1 <- e * f1 + (1 - e) * f1c                    ; f2 <- e * f2 + (1 - e) * f2c
            gmu <- (W1 * f1 - W2 * f2)
            gb <- - (f1 - f2) * X
            gradi <- sgn * cbind(gb, gmu) / Li
            .gradient <- gradi
            if (sum) .gradient <- apply(gradi, 2, sum)
            attr(lnL, "gradient") <- .gradient
        }
        if (hessian){
            e1 <- e_link(mu[y + 1] - bX)              ; e2 <- e_link(mu[y] - bX)
            M1 <- cbind(- X, W1)
            M2 <- cbind(- X, W2)
            H <- crossprod(M1 * e1 / Li, M1) - crossprod(M2 * e2 / Li, M2) - crossprod(gradi)
            H <- sgn * H
            attr(lnL, "hessian") <- H
        }
        if (information){
            attr(lnL, "info") <- - H
        }
        lnL
    }
    sup.coef <- q_link(cumsum(prop.table(table(y))))[1: (J - 1)]
    .start <- c(rep(0.1, K), sup.coef)
    names(.start) <- c(colnames(X), nms_thr)
    .coefs <- maximize(lnl, .start, method = .method, trace = 0, X = X, y = y, e = e)
    .lp <- drop(X %*% .coefs[1:K])
    .lnl_conv <- lnl(.coefs, gradient = TRUE, hessian = TRUE, sum = FALSE,
                     X = X, y = y, e = e, information = TRUE)
    .null_intercept <- sup.coef
    fy <- prop.table(table(y))
    logLik_null <- N * sum(fy * log(fy))
    logLik_saturated <- 0
    logLik_model <- sum(.lnl_conv)
    .fitted.values <- cbind(0, F_link(outer(- .lp, .coefs[K + (1L:(J - 1))], "+")), 1)
    .fitted.values <- .fitted.values[, - 1, drop = FALSE] -
        .fitted.values[, - (J + 1L), drop = FALSE]
    .probabilities <- .fitted.values
    .fitted.values <- as.numeric(apply(Y * .fitted.values, 1, sum))
    .df.residual <- N - K - J - 1L
    colnames(.probabilities) <- nms_y
    .npar <- c(covariates = K, threshold = J - 1L)
    attr(.npar, "default") <- "covariates"
    # Null model
    .logLik <- c(model = logLik_model, null = logLik_null, saturated = logLik_saturated)
    null_coefs <- c(rep(0, ncol(X)), .null_intercept)
    lnl_null <- lnl(null_coefs, gradient = TRUE, hessian = TRUE,
                    information = TRUE, X = X, y = y, e = e)
    .null_gradient <- attr(lnl_null, "gradient")
    .null_info <- attr(lnl_null, "info")
    .lm <- drop(crossprod(.null_gradient, solve(.null_info, .null_gradient)))
    .model_info <- attr(.lnl_conv, "info")
#    .w <- drop(crossprod(.coefs[- 1], t(crossprod(.coefs[-1], .model_info[- 1, - 1]))))
    .vcov <- solve(.model_info)
    .w <- drop(crossprod(.coefs[1:K], t(crossprod(.coefs[1:K],
                                                  solve(.vcov[1:K, 1:K, drop = FALSE])))))
    .lr <- 2 * unname(.logLik["model"] - .logLik["null"])
    tests <- c(wald = .w, score = .lm, lr = .lr)

    result <- list(coefficients = .coefs,
                   model = mf,
                   gradient = attr(.lnl_conv, "gradient"),
                   hessian = attr(.lnl_conv, "hessian"),
                   linear.predictors = .lp,
                   logLik = .logLik,
                   probabilities = .probabilities,
                   fitted.values = .fitted.values,
                   df.residual = .df.residual,
                   est_method = "ml",
                   formula = formula,
                   npar = .npar,
                   tests = tests,
                   value = as.numeric(.lnl_conv),
                   call = .call
                   )
    structure(result, class = c("ordreg", "binomreg", "micsr"))
}


fitted.ordreg <- function(object, ..., type = c("outcome", "probabilities")){
    .type <- match.arg(type)
    ifelse(.type == "outcome",
           object$fitted.values,
           object$probabilities)
}

