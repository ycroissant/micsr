#' Ordered regression
#'
#' Maximum-likelihood estimation of a model for which the response is
#' ordinal
#' 
#' @name ordreg
#' @param formula a symbolic description of the model
#' @param data a data frame
#' @param subset,weights,na.action,offset,contrasts see `lm`
#' @param link one of `probit` and `logit`
#' @param start a vector of starting values,
#' @param opt optimization method
#' @param maxit maximum number of iterations
#' @param trace printing of intermediate result
#' @param check_gradient if `TRUE` the numeric gradient and hessian
#'     are computed and compared to the analytical gradient and
#'     hessian
#' @param object a `ordreg` object
#' @param type one of `"outcome"` or `"probabilities"` for the
#'     `fitted` method
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
#'                link = "cloglog", opt = "bfgs")
#' @export
ordreg <- function(formula, data, weights, subset, na.action, offset, contrasts = NULL, 
                   link = c("probit", "logit", "cloglog"),
                   start = NULL, 
                   opt = c("bfgs", "nr", "newton"),
                   maxit = 100, trace = 0, 
                   check_gradient = FALSE, ...){
    .opt <- match.arg(opt)
    .call <- match.call()
    .link <- match.arg(link)
    cl <- match.call(expand.dots = FALSE)
    .formula <- cl$formula <- Formula(formula)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
               names(cl), 0L)
    # construct the model frame and components
    cl <- cl[c(1L, m)]
    mf <- cl
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf)
    wt <- as.vector(model.weights(mf))
    if (is.null(wt)) wt <- 1 else wt <- wt #/ mean(wt)
    .offset <- model.offset(mf)
    if (inherits(y, "Surv")){
        y <- as.matrix(model.response(mf))
        e <- as.logical(y[, 2])
        y <- as.factor(y[, 1])
    }
    else {
        e <- rep(1, length(y))
        if (! inherits(y, "factor")) y <- as.factor(y)
    }
    X <- model.matrix(mt, mf, contrasts)
    if (colnames(X)[1] == "(Intercept)") X <- X[, - 1, drop = FALSE]
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
                    information = FALSE, X = X, y = y, e = e, weights,
                    opposite = FALSE){
        sgn <- ifelse(opposite, - 1, 1)
        J <- length(unique(y))
        K <- ncol(X)
        beta <- param[1L:K]
        mu <- c(-100, param[(K + 1L):(K + J - 1L)], 100)
        bX <- as.numeric(crossprod(t(X), beta))
        z1 <- mu[y + 1] - bX               ; z2 <- mu[y] - bX
        F1 <- F_link(z1)                   ; F2 <- F_link(z2)
        Li <- (2 * e - 1) * F1 - e * F2 + (1 - e)
        lnl <- sgn * log(Li)
        if (sum) lnl <- sum(lnl * weights)
        if (gradient | hessian){
            W <- matrix(0, nrow(X), J)
            W1 <- (col(W) == (y + 1))[, - 1, drop = FALSE]
            W2 <- (col(W) == y)[, - 1, drop = FALSE]
            M1 <- cbind(- X, W1)
            M2 <- cbind(- X, W2)
            # The first two cols are not estimated coefficients (-infty and 0)
            # so remove them
            f1 <- f_link(mu[y + 1] - bX)
            f2 <- f_link(mu[y] - bX)
            gradi <- (2 * e - 1) * f1 * M1 -  e * f2 * M2
            gradi <- sgn * gradi / Li
            .gradient <- gradi
            colnames(.gradient) <- names(param)
            if (sum) .gradient <- apply(gradi * weights, 2, sum)
            attr(lnl, "gradient") <- .gradient
        }
        if (hessian){
            e1 <- e_link(mu[y + 1] - bX)
            e2 <- e_link(mu[y] - bX)
            H <- crossprod((2 * e - 1) * weights * M1 * e1 / Li,  M1) -
                crossprod(e * weights * M2 * e2 / Li, M2) -
                crossprod(sqrt(weights) * gradi)
            colnames(H) <- rownames(H) <- names(param)
            attr(lnl, "hessian") <- sgn * H
        }
        if (information){
            attr(lnl, "info") <- - H
        }
        lnl
    }

    sup.coef <- q_link(cumsum(prop.table(table(y))))[1: (J - 1)]
    if (is.null(start)){
        .start <- c(rep(0.1, K), sup.coef)
    } else .start <- start
    nms <- c(colnames(X), nms_thr)
    names(.start) <- nms
    if (maxit > 0){
        .coefs <- maximize(lnl, start = .start, trace = trace, method = .opt, maxit = maxit,
                           X = X, y = y, e = e, weights = wt, ...)
    } else .coefs <- .start
    .linpred <- drop(X %*% .coefs[1:K])
    .lnl_conv <- lnl(.coefs, gradient = TRUE, hessian = TRUE, sum = FALSE,
                     X = X, y = y, e = e, weights = wt, information = TRUE)

    
    
    fun <- function(x) lnl(x, gradient = TRUE, hessian = TRUE, sum = TRUE,
                           X = X, y = y, e = e, weights = wt, information = TRUE)
    if(check_gradient) z <- check_gradient(fun, .coefs) else z <- NA

    .null_intercept <- sup.coef
    fy <- prop.table(table(y))

    lnl_model <- as.numeric(.lnl_conv)
    lnl_null <- sum(fy * log(fy))
    lnl_saturated <- 0
    values <- cbind(model = lnl_model, saturated = lnl_saturated, null = lnl_null)
    .logLik <- apply(values * wt, 2, sum)
    
    .fitted.values <- cbind(0, F_link(outer(- .linpred, .coefs[K + (1L:(J - 1))], "+")), 1)
    .fitted.values <- .fitted.values[, - 1, drop = FALSE] -
        .fitted.values[, - (J + 1L), drop = FALSE]
    colnames(.fitted.values) <- nms_y
    .df.residual <- N - K - J - 1L
    .npar <- c(covariates = K, threshold = J - 1L)
    attr(.npar, "default") <- "covariates"
    # Null model
    null_coefs <- c(rep(0, ncol(X)), .null_intercept)
    lnl_null <- lnl(null_coefs, gradient = TRUE, hessian = TRUE,
                    information = TRUE, X = X, y = y, e = e, weights = wt)
    .null_gradient <- attr(lnl_null, "gradient")
    .null_info <- attr(lnl_null, "info")
    .lm <- drop(crossprod(.null_gradient, solve(.null_info, .null_gradient)))
    .model_info <- attr(.lnl_conv, "info")
    .vcov <- solve(.model_info)
    .w <- drop(crossprod(.coefs[1:K], t(crossprod(.coefs[1:K],
                                                  solve(.vcov[1:K, 1:K, drop = FALSE])))))
    .lr <- 2 * unname(.logLik["model"] - .logLik["null"])
    tests <- c(wald = .w, score = .lm, lr = .lr)
    result <- list(coefficients = .coefs,
                   model = mf,
                   terms = mt,
                   value = values,
                   gradient = attr(.lnl_conv, "gradient"),
                   hessian = attr(.lnl_conv, "hessian"),
                   fitted.values = .fitted.values,
                   linear.predictors = .linpred,
                   logLik = .logLik,
                   tests = tests,
                   df.residual = .df.residual,
                   npar = .npar,
                   est_method = "ml",
                   call = .call,
                   na.action = attr(mf, "na.action"),
                   offset = .offset,
                   contrasts = attr(X, "contrasts"),
                   xlevels = .getXlevels(mt, mf),
                   check_gradient = z)
    structure(result, class = c("ordreg", "binomreg", "micsr"))
}


#' @rdname ordreg
#' @method fitted ordreg
#' @export
fitted.ordreg <- function(object, ..., type = c("outcome", "probabilities")){
    .type <- match.arg(type)
    .probs <- object$fitted.values
    if(.type == "probabilities") .probs
    else{
        y = model.response(model.frame(object))
        if (inherits(y, "Surv")) y <- y[, 1]
        Y <- model.matrix(~ factor(y) - 1, data.frame(y = y))
        as.numeric(apply(Y * .probs, 1, sum))
    }
}

