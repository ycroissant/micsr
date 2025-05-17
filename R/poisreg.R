#' Poisson regression
#'
#' A unified interface to perform Poisson, Negbin and log-normal Poisson models
#'
#' @name poisreg
#' @param formula a symbolic description of the model, (for the count
#'     component and for the selection equation)
#' @param data a data frame
#' @param subset,weights,na.action,offset,contrasts see `stats::lm`,
#' @param start a vector of starting values
#' @param mixing the mixing distribution, one of `"none"`, `"gamma"`
#'     and `"lognorm"`
#' @param vlink one of `"nb1"` and `"nb2"`
#' @param opt optimization method
#' @param maxit maximum number of iterations
#' @param trace printing of intermediate result
#' @param check_gradient if `TRUE` the numeric gradient and hessian
#'     are computed and compared to the analytical gradient and
#'     hessian
#' @param object a `poisreg` object
#' @param type the type of residuals for the `residuals` method
#' @param ... further arguments
#' @param vcov the covariance matrix estimator to use for the score
#'     test
#' @return an object of class `c("poisreg", "micsr")`, see
#'     `micsr::micsr` for further details.
#' @importFrom stats glm plogis qlogis
#' @importFrom Formula Formula
#' @keywords models
#' @examples
#' nb1 <- poisreg(trips ~ workschl + size + dist + smsa + fulltime + distnod +
#'                realinc + weekend + car, trips, mixing = "gamma", vlink = "nb1")
#' @export
poisreg <- function(formula, data, weights, subset, na.action, offset, contrasts = NULL, 
                    start = NULL, mixing = c("none", "gamma", "lognorm"),
                    vlink = c("nb1", "nb2"),
                    opt = c("bfgs", "nr", "newton"),
                    maxit = 100, trace = 0,
                    check_gradient = FALSE, ...){
    .call <- match.call()
    cl <- match.call(expand.dots = FALSE)
    .opt <- match.arg(opt)
    .formula <- cl$formula <- Formula(formula)
    .mixing <- match.arg(mixing)
    .vlink <- match.arg(vlink)
    # construct the model frame and components
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
               names(cl), 0L)
    cl <- cl[c(1L, m)]
    mf <- cl
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf, contrasts)
    wt <- as.vector(model.weights(mf))
    if (is.null(wt)) wt <- 1 else wt <- wt #/ mean(wt)
    .offset <- model.offset(mf)
    y <- model.response(mf)
    yb <- mean(y)
    K <- ncol(X)
    N <- length(y)
    .df.residual <- N - K
    .null_logLik <- N * yb * (log(yb) - 1) - sum(lfactorial(y))
    .sat_logLik <- sum((y * (log(y) - 1) - lfactorial(y))[y > 0])
    .null_intercept <- log(yb)
    
    # Poisson model
    if (.mixing == "none"){
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, info = FALSE,
                        opposite = FALSE, sum = TRUE, X, y, weights){
            sgn <- ifelse(opposite, - 1, + 1)
            linpred <- drop(X %*% coefs)
            mu <- exp(linpred)
            lnl <- - mu + y * log(mu) - lfactorial(y)
            lnl <- lnl * sgn
            if (gradient) grad <- (y - mu) * X
            if (hessian) hess <- -  crossprod(weights *  mu * X, X)
            if (info) information <- crossprod(weights * mu * X, X)
            if (sum){
                lnl <- sum(weights * lnl)
                if (gradient) grad <- apply(weights * grad, 2, sum)
            }
            if (gradient) attr(lnl, "gradient") <- grad * sgn
            if (hessian) attr(lnl, "hessian") <- hess * sgn
            if (info) attr(lnl, "info") <- information
            lnl
        }
        if (is.null(start)){
            start <- rep(0, K)
            names(start) <- colnames(X)
        }
        .npar <- c(covariates = K)
        attr(.npar, "default") <- "covariates"
        attr(.npar, "null") <- 1
    }
    # Negbin model
    if (.mixing == "gamma"){
        if (.vlink == 'nb1') k <- 1 else k <- 2
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, info = FALSE,
                        opposite = FALSE, sum = TRUE, X, y, weights){
            sgn <- ifelse(opposite, - 1, + 1)
            K <- ncol(X)
            beta <- coefs[1L:K]
            sigma <- coefs[K + 1L]
            bX <- as.numeric(crossprod(t(X), beta))
            l <- exp(bX)
            if (abs(sigma) >= 1E-10){
                v <- l ^ (2 - k) / sigma
                lnl <- lgamma(y + v) - lgamma(y + 1) - lgamma(v) + v * log(v) -
                    v * log(v + l) + y * log(l) - y * log(v + l)
                lnl <- sgn * lnl
                if (sum) lnl <- sum(weights * lnl)
                if (gradient){
                    lb <- l
                    vb <- (2 - k) * v
                    vs <- - v / sigma
                    Ll <- - (v + y) / (l + v) + y / l
                    Lv <-  log(v) + 1 - log(l + v) - (v + y) / (l + v) + digamma(y + v) - digamma(v)
                    gb <- (Ll + (2 - k) * v * Lv / l) * l
                    gs <- - v / sigma * Lv
                    gradi <-  cbind(gb * X, sigma = gs)
                    if (sum) gradi <- apply(weights * gradi, 2, sum)
                    attr(lnl, "gradient") <- gradi * sgn
                }
                if (hessian){
                    Lll <- (v + y) / (l + v) ^ 2 - y / l ^ 2
                    Llv <- (y - l) / (l + v) ^ 2
                    Lvv <- 1 / v - 1 / (l + v) + (y - l) / (l + v) ^ 2 + trigamma(y + v) - trigamma(v)
                    lbb <- l
                    vbb <- (2 - k) ^ 2 * l
                    vbs <- - (2 - k) * v / sigma
                    vss <- 2 * v / sigma ^ 2
                    ## Hbb <- Ll * lbb + Lv * vbb +
                    ##     vb * (Lvv * vb + Llv * lb) +
                    ##     lb * (Llv * vb + Lll * lb)                
                    vl <- vb / l
                    Hbb <- l * Lll + Ll + l * Llv * vl + (2 - k) * ( (Lv + v * Lvv) * vl + v * Llv)
                    Hbb <- Hbb * l
                    Hbv <- Lv * vbs + vb * (Lvv * vs ) + lb * (Llv * vs)
                    Hvv <- Lv * vss + vs * (Lvv * vs)
                    Hbb <- crossprod(weights * X * Hbb, X)
                    Hbv <- apply(weights * Hbv * X, 2, sum)
                    Hvv <- sum(weights * Hvv)
                    H <- rbind(cbind(Hbb, sigma = Hbv), sigma = c(Hbv, Hvv))
                    attr(lnl, "hessian") <- H * sgn
                }
            } else {
                linpred <- drop(X %*% beta)
                mu <- exp(linpred)
                lnl <- - mu + y * log(mu) - lfactorial(y)
                lnl <- sgn * lnl
                if (sum) lnl <- sum(weights * lnl)
                if (gradient){
                    gb <- (y - mu)
                    gs <- ( (y - mu) ^ 2- y) / (2 * mu ^ (2 - k))
                    gradi <-  cbind(gb * X, sigma = gs)
                    if (sum) gradi <- apply(weights * gradi, 2, sum)
                    attr(lnl, "gradient") <- gradi * sgn
                }
                if (hessian){
                    Hbb <- -  crossprod(weights *  mu * X, X)
                    Hss <- sum(weights * (y ^ 2 / 2 - mu * (y - mu) ^ 2)/ (mu ^ (4 - 2 * k)))
                    Hbs <- - apply(weights * ( y - mu + (2 - k) / 2 * ( (y - mu) ^ 2 - y) / mu) * X, 2, sum)
                    H <- rbind(cbind(Hbb, sigma = Hbs), sigma = c(Hbs, Hss))
                    attr(lnl, "hessian") <- H * sgn
                }
                if (info){
                    information <- rbind(cbind(crossprod(weights *  mu * X, X), sigma = rep(0, ncol(X))),
                                  sigma = c(rep(0, ncol(X)), sum(weights * (mu - 1) / (2 * mu ^ (3 - 2 * k)))))
                    attr(lnl, "info") <- information
                }
            }
            lnl
        }
        if (is.null(start)){
            start <- c(rep(0, K), 1E-08)
            names(start) <- c(colnames(X), "sigma")
        }
        .npar <- c(covariates = K, mixing = 1)
        attr(.npar, "default") <- c("covariates", "mixing")
        attr(.npar, "null") <- 1
    }
    if (.mixing == "lognorm"){
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, info = FALSE,
                        opposite = FALSE, sum = TRUE, X, y, weights){
            sgn <- ifelse(opposite, - 1, + 1)
            K <- ncol(X)
            beta <- coefs[1L:K]
            sigma <- coefs[K + 1L]
            bX <- X %*% beta
            gq <- gauss_hermite(40)
            lambda <- sapply(gq$nodes, function(r) bX + sqrt(2) * sigma * r)
            qr <- exp(- exp(lambda)) * exp(lambda * y)
            W <- matrix(gq$weights, nrow = nrow(qr), ncol = ncol(qr), byrow = TRUE)
            O <- matrix(gq$nodes, nrow = nrow(qr), ncol = ncol(qr), byrow = TRUE)
            q <- apply(W * qr, 1, sum)
            lnl <- - 0.5 * log(pi) - lfactorial(y) + log(q)
            lnl <- sgn * lnl
            if (sum) lnl <- sum(weights * lnl)
            if (gradient){
                q_nr <- exp(- exp(lambda)) * exp(lambda * y)
                dq_nr <- q_nr * (y - exp(lambda))
                q_n <-  apply(W * q_nr, 1, sum)
                dq_n <- apply(W * dq_nr, 1, sum)
                dsq_n <- apply(W * dq_nr * sqrt(2) * O, 1, sum)
                g <- cbind(dq_n / q_n * X, sigma = dsq_n / q_n)
                if (sum) g <- apply(weights * g, 2, sum)
                attr(lnl, "gradient") <- g * sgn
            }
            if (hessian){
                d2q_n <- apply(weights * (W * (dq_nr * (y - exp(lambda)) - q_nr * exp(lambda))), 1, sum)
                H <- crossprod(d2q_n * X / q_n, X) - crossprod(sqrt(weights) * dq_n / q_n * X)
                d2cq_n <- apply(weights * W * (dq_nr * (y - exp(lambda)) - q_nr * exp(lambda)) * sqrt(2) * O, 1, sum)
                d2sq_n <- apply(weights * W * (dq_nr * (y - exp(lambda)) - q_nr * exp(lambda)) * 2  * O ^ 2, 1, sum)
                H_gs <- apply(d2cq_n * X / q_n, 2, sum) - apply(weights * dsq_n * dq_n / q_n ^ 2 * X, 2, sum)
                H_ss <- sum(d2sq_n / q_n) - sum(weights * dsq_n ^ 2 / q_n ^  2)
                H <- rbind(cbind(H, sigma = H_gs), sigma = c(H_gs, H_ss))
                attr(lnl, "hessian") <- H * sgn
            }
            lnl
        }
        if (is.null(start)){
            start <- c(glm(y ~ X - 1, family = poisson)$coef, sigma  = 1E-1)
            names(start) <- c(colnames(X), "sigma")
        }
        .npar <- c(covariates = K, mixing = 1)
        attr(.npar, "default") <- c("covariates", "mixing")
        attr(.npar, "null") <- 1
    }

    if (maxit > 0){
        .coefs <- maximize(lnl, start = start, method = .opt, trace = trace, maxit = maxit,
                           X = X, y = y, weights = wt, ...)
    } else .coefs <- start
    ## if (.method == "bfgs"){
    ##     f <- function(x) - lnl(x, X = X, y = y, weights = wt)
    ##     g <- function(x) - attr(lnl(x, X = X, y = y, weights = wt, gradient = TRUE), "gradient")
    ##     .coefs <- optim(start, f, g, method = "BFGS")$par
    ## }
    ## else .coefs <- newton(lnl, X = X, y = y, weights = wt, trace = 1, coefs = start, direction = "max")

    ## f <- function(x) - lnl(x, X = X, y = y, weights = wt)
    ## g <- function(x) - attr(lnl(x, X = X, y = y, weights = wt, gradient = TRUE), "gradient")
#    .coefs <- optim(start, f, g, method = "BFGS")$par

    .linpred <- drop(X %*% .coefs[1:K])
    .fitted <- exp(.linpred)
    .lnl_conv <- lnl(.coefs, X = X, y = y, weights = wt, gradient = TRUE,
                     hessian = TRUE, info = TRUE, sum = FALSE, opposite = FALSE)

    fun <- function(x) lnl(x, X = X, y = y, weights = wt, gradient = TRUE,
                     hessian = TRUE, info = TRUE, sum = TRUE)
    if (check_gradient) z <- check_gradient(fun, .coefs) else z <- NA

    l_model <- as.numeric(.lnl_conv)
    l_saturated <- ifelse(y > 0, y * (log(y) - 1) - lfactorial(y), 0)
    l_null <- - yb + y * log(yb) - lfactorial(y)
    
    if (.mixing == "gamma"){
        if (.vlink == 'nb1') k <- 1 else k <- 2
        sigma <- as.numeric(.coefs["sigma"])
        l <- y
        v <- l ^ (2 - k) / sigma
        l_saturated <- lgamma(y + v) - lgamma(y + 1) - lgamma(v) + v * log(v) - v * log(v + l) +
            ifelse(y > 0, y * log(l), 0) - y * log(v + l)
    }
    
    values <- cbind(model = l_model, saturated = l_saturated, null = l_null)
    .logLik <- apply(values * wt, 2, sum)
    
    # Null model
    .coefs_0 <- .coefs
    .index_0 <- attr(.npar, "null")
    .coefs_0[.index_0] <- .null_intercept
    .coefs_0[- .index_0] <- 0
    if (.mixing == "gamma") .coefs_0["sigma"] <- 1E-7
    .lnl_0 <- lnl(.coefs_0, gradient = TRUE, hessian = TRUE, info = TRUE,
                  X = X, y = y, weights = wt)
    .g_0 <- attr(.lnl_0, "gradient")
    .info_0 <- attr(.lnl_0, "info")
    if (is.null(.info_0)) .info_0 <- - attr(.lnl_0, "hessian")
    .V_0 <- solve(.info_0)

    # The three statistics
    .lm <- drop(crossprod(.g_0, solve(.info_0, .g_0)))
    .info_model <- attr(.lnl_conv, "info")
    if (is.null(.info_model)) .info_model <- - attr(.lnl_conv, "hessian")
    .vcov <- solve(.info_model)[- .index_0, - .index_0, drop = FALSE]
    .w <- drop(crossprod(.coefs[- .index_0],
                         t(crossprod(.coefs[-  .index_0],
                                     solve(.vcov)))))
    .lr <- 2 * unname(.logLik["model"] - .logLik["null"])
    .rsq <- c(w = .w / (.w + N), lr = 1 - exp(-.lr / N), lm = .lm / N)
    tests <- c(wald = .w, score = .lm, lr = .lr)

    # families
    if (.mixing == "gamma"){
        env <- new.env(parent=.GlobalEnv)
        .sigma <- .coefs["sigma"]
        assign(".sigma", .sigma, envir=env)
        if (.vlink == "nb2") .variance <- function(mu) (1 + .sigma * mu) * mu
        if (.vlink == "nb1") .variance <- function(mu) (1 + .sigma)      * mu
        environment(.variance) <- env
    }
    if (.mixing == "lognorm"){
        env <- new.env(parent=.GlobalEnv)
        .sigma <- .coefs["sigma"]
        assign(".sigma", .sigma, envir=env)
        .variance <- function(mu) (1 + (exp(.sigma ^ 2) - 1) * exp(0.5 * .sigma ^ 2) * mu) *
                                      exp(0.5 * .sigma ^ 2) * mu
        environment(.variance) <- env
    }
    if (.mixing == "none"){
        .variance = function(mu) mu
    }
    .family <- structure(list(
        vlink  = .vlink,
        variance = .variance,
        family = "Poisson",
        link = "log",
        mixing = .mixing
    ),
    class = "family")
    
    result <- list(coefficients = .coefs,
                   model = mf,
                   terms = mt,
                   value = values,
                   gradient = attr(.lnl_conv, "gradient"),
                   hessian = attr(.lnl_conv, "hessian"),
                   info = attr(.lnl_conv, "info"),
                   fitted.values = .fitted,
                   linear.predictors = .linpred,
                   logLik = .logLik,
                   tests = tests,
                   df.residual = .df.residual,
                   npar = .npar,
                   est_method = "ml",
                   call = .call,
                   na.action = attr(mf, "na.action"),
                   weights = wt,
                   offset = .offset,
                   contrasts = attr(X, "contrasts"),
                   xlevels = .getXlevels(mt, mf),
                   family = .family,
                   check_gradient = z)
    structure(result, class = c("poisreg", "micsr"))
}

#' @rdname poisreg
#' @method scoretest poisreg
#' @export
scoretest.poisreg <- function(object, ..., vcov = NULL){
    .vcov <- ifelse(is.null(vcov), "opg", vcov)
    new <- list(...)[[1]]
    cls <- class(object)[1]
    if (inherits(new, "formula")){
        class(object) <- setdiff(class(object), "poisreg")
        scoretest(object, ...)
    } else {
        new <- update(object, ..., maxit = 0, start = c(coef(object), sigma = 0))
        gs <- sum(new$gradient[, "sigma"])
        if (.vcov == "opg") hm1 <- solve(crossprod(new$gradient))["sigma", "sigma"]
        if (.vcov == "hessian") hm1 <- solve(- new$hessian)["sigma", "sigma"]
        if (.vcov == "info") hm1 <- 1 / new$info["sigma", "sigma"]
#        hm1 <- - solve(new$hessian)["sigma", "sigma"]
        .statistic <- gs ^ 2 * hm1
        mu <- object$fitted
        k <- 2
        names(.statistic) <- "chisq"
        .parameter <- 1
        .method <- "score test"
        .pval <- pchisq(.statistic, df = 1, lower.tail = FALSE)
        .alternative <- "overdispersion"
        .data.name <- paste(deparse(formula(object)))
        structure(list(statistic = .statistic,
                       parameter = .parameter,
                       p.value = .pval,
                       method = .method,
                       data.name = .data.name),
                  class = "htest")
    }
}    
    

#' @rdname poisreg
#' @export
residuals.poisreg <- function(object, ..., type = c("deviance", "pearson", "response")){
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
