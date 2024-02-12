#' Poisson regression
#'
#' A unified interface to perform Poisson, Negbin and log-normal Poisson models
#'
#' @name poisreg
#' @param formula a symbolic description of the model, (for the count
#'     component and for the selection equation)
#' @param data a data frame
#' @param subset,weights,na.action,offset see `stats::lm`,
#' @param start a vector of starting values
#' @param mixing the mixing distribution, one of `"none"`, `"gamma"`
#'     and `"lognorm"`
#' @param vlink one of `"nb1"` and `"nb2"`
#' @param method the optimization method, one of `"newton"` and `"bfgs"`
#' @param ... further arguments
#' @return an object of class `c("poisreg", "micsr")`, see
#'     `micsr::micsr` for further details
#' @importFrom stats glm plogis qlogis
#' @importFrom Formula Formula
#' @keywords models
#' @examples
#' nb1 <- poisreg(trips ~ workschl + size + dist + smsa + fulltime + distnod +
#'                realinc + weekend + car, trips, mixing = "gamma", vlink = "nb1")
#' @export
poisreg <- function(formula, data, weights, subset, na.action, offset,
                    start = NULL, mixing = c("none", "gamma", "lognorm"),
                    method = c("bfgs", "newton"), vlink = c("nb1", "nb2"), ...){
    .call <- match.call()
    cl <- match.call(expand.dots = FALSE)
    .formula <- cl$formula <- Formula(formula)
    .method <- match.arg(method)
    .mixing <- match.arg(mixing)
    .vlink <- match.arg(vlink)
    # construct the model frame and components
    m <- match(c("formula", "data", "subset", "weights"),
               names(cl), 0L)
    cl <- cl[c(1L, m)]
    mf <- cl
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)
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
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, information = FALSE, sum = TRUE, X, y){
            linpred <- drop(X %*% coefs)
            mu <- exp(linpred)
            lnl <- - mu + y * log(mu) - lfactorial(y)
            if (gradient) grad <- (y - mu) * X
            if (hessian) hess <- -  crossprod( mu * X, X)
            if (information) info <- crossprod( mu * X, X)
            if (sum){
                lnl <- sum(lnl)
                if (gradient) grad <- apply(grad, 2, sum)
            }
            if (gradient) attr(lnl, "gradient") <- grad
            if (hessian) attr(lnl, "hessian") <- hess
            if (information) attr(lnl, "info") <- info
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
        if (vlink == 'nb1') k <- 1 else k <- 2
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, information = FALSE, sum = TRUE, X, y){
            K <- ncol(X)
            beta <- coefs[1L:K]
            sigma <- coefs[K + 1L]
            bX <- as.numeric(crossprod(t(X), beta))
            l <- exp(bX)
            v <- l ^ (2 - k) / sigma
            lnL <- lgamma(y + v) - lgamma(y + 1) - lgamma(v) + v * log(v) -
                v * log(v + l) + y * log(l) - y * log(v + l)
            if (sum) lnL <- sum(lnL)
            if (gradient){
                lb <- l
                vb <- (2 - k) * v
                vs <- - v / sigma
                Ll <- - (v + y) / (l + v) + y / l
                Lv <-  log(v) + 1 - log(l + v) - (v + y) / (l + v) + digamma(y + v) - digamma(v)
                gb <- (Ll * l + (2 - k) * v * Lv)
                gs <- - v / sigma * Lv
                gradi <-  cbind(gb * X, sigma = gs)
                if (sum) gradi <- apply(gradi, 2, sum)
                attr(lnL, "gradient") <- gradi
            }
            if (hessian){
                Lll <- (v + y) / (l + v) ^ 2 - y / l ^ 2
                Llv <- (y - l) / (l + v) ^ 2
                Lvv <- 1 / v - 1 / (l + v) + (y - l) / (l + v) ^ 2 + trigamma(y + v) - trigamma(v)
                lbb <- l
                vbb <- (2 - k) ^ 2 * l
                vbs <- - (2 - k) * v / sigma
                vss <- 2 * v / sigma ^ 2
                Hbb <- Ll * lbb + Lv * vbb +
                    vb * (Lvv * vb + Llv * lb) +
                    lb * (Llv * vb + Lll * lb)
                Hbv <- Lv * vbs +
                    vb * (Lvv * vs ) +
                    lb * (Llv * vs)
                Hvv <- Lv * vss +
                    vs * (Lvv * vs)
                Hbb <- crossprod(X * Hbb, X)
                Hbv <- apply(Hbv * X, 2, sum)
                Hvv <- sum(Hvv)
                H <- rbind(cbind(Hbb, sigma = Hbv), sigma = c(Hbv, Hvv))
                attr(lnL, "hessian") <- H
            }
            lnL
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
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, information = FALSE, sum = TRUE, X, y){
            K <- ncol(X)
            beta <- coefs[1L:K]
            sigma <- coefs[K + 1L]
            bX <- X %*% beta
#            gq <- statmod::gauss.quad(40, "hermite")
            gq <- gaussian_quad(40, "hermite")
            lambda <- sapply(gq$nodes, function(r) bX + sqrt(2) * sigma * r)
            qr <- exp(- exp(lambda)) * exp(lambda * y)
            W <- matrix(gq$weights, nrow = nrow(qr), ncol = ncol(qr), byrow = TRUE)
            O <- matrix(gq$nodes, nrow = nrow(qr), ncol = ncol(qr), byrow = TRUE)
            q <- apply(W * qr, 1, sum)
            lnl <- - 0.5 * log(pi) - lfactorial(y) + log(q)
            if (sum) lnl <- sum(lnl)
            if (gradient){
                q_nr <- exp(- exp(lambda)) * exp(lambda * y)
                dq_nr <- q_nr * (y - exp(lambda))
                q_n <-  apply(W * q_nr, 1, sum)
                dq_n <- apply(W * dq_nr, 1, sum)
                dsq_n <- apply(W * dq_nr * sqrt(2) * O, 1, sum)
                g <- cbind(dq_n / q_n * X, sigma = dsq_n / q_n)
                if (sum) g <- apply(g, 2, sum)
                attr(lnl, "gradient") <- g
            }
            if (hessian){
                d2q_n <- apply(W * (dq_nr * (y - exp(lambda)) - q_nr * exp(lambda)), 1, sum)
                H <- crossprod(d2q_n * X / q_n, X) - crossprod(dq_n / q_n * X)
                d2cq_n <- apply(W * (dq_nr * (y - exp(lambda)) - q_nr * exp(lambda)) * sqrt(2) * O, 1, sum)
                d2sq_n <- apply(W * (dq_nr * (y - exp(lambda)) - q_nr * exp(lambda)) * 2  * O ^ 2, 1, sum)
                H_gs <- apply(d2cq_n * X / q_n, 2, sum) - apply(dsq_n * dq_n / q_n ^ 2 * X, 2, sum)
                H_ss <- sum(d2sq_n / q_n) - sum(dsq_n ^ 2 / q_n ^  2)
                H <- rbind(cbind(H, sigma = H_gs), sigma = c(H_gs, H_ss))
                attr(lnl, "hessian") <- H
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

    ## cat("_________________\n")
    ## cat("numerical gradient\n")
    ## print(numDeriv::grad(lnl, start, X = X, y = y))
    ## cat("_________________\n")
    ## cat("_________________\n")
    ## cat("numerical hessian\n")
    ## print(numDeriv::hessian(lnl, start, X = X, y = y))
    ## cat("_________________\n")
    ## va <- lnl(start, X = X, y = y, gradient = TRUE, hessian = TRUE)
    ## cat("_________________\n")
    ## cat("analytical gradient\n")
    ## print(attr(va, "gradient"))
    ## cat("_________________\n")
    ## cat("_________________\n")
    ## cat("analytical hessian\n")
    ## print(attr(va, "hessian"))
    ## cat("_________________\n")

    if (.method == "bfgs"){
        f <- function(x) - lnl(x, X = X, y = y)
        g <- function(x) - attr(lnl(x, X = X, y = y, gradient = TRUE), "gradient")
        .coefs <- optim(start, f, g, method = "BFGS")$par
    }
    else .coefs <- newton(lnl, X = X, y = y, trace = 1, coefs = start, direction = "max")

    .linpred <- drop(X %*% .coefs[1:K])
    .mu <- exp(.linpred)
    .lnl_conv <- lnl(.coefs, X = X, y = y, gradient = TRUE, hessian = TRUE, info = TRUE, sum = FALSE)
    .fitted <- dpois(y, .mu)
    .logLik <- c(model = sum(as.numeric(.lnl_conv)),
                 saturated = .sat_logLik,
                 null = .null_logLik)
    
    # Null model
    .coefs_0 <- .coefs
    .index_0 <- attr(.npar, "null")
    .coefs_0[.index_0] <- .null_intercept
    .coefs_0[- .index_0] <- 0
    if (.mixing == "gamma") .coefs_0["sigma"] <- 1E-7
    .lnl_0 <- lnl(.coefs_0, gradient = TRUE, hessian = TRUE, information = TRUE, X = X, y = y)
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
    tests <- c(w = .w, lm = .lm, lr = .lr)
    result <- list(coefficients = .coefs,
                   model = mf,
                   gradient = attr(.lnl_conv, "gradient"),
                   hessian = attr(.lnl_conv, "hessian"),
                   info = attr(.lnl_conv, "info"),
                   linear.predictors = .linpred,
                   logLik = .logLik,
                   fitted.values = .fitted,
                   df.residual = .df.residual,
                   est_method = "ml",
                   formula = formula,
                   npar = .npar,
                   value = as.numeric(.lnl_conv),
                   tests = tests,
                   call = .call
                   )
    structure(result, class = c("poisreg", "micsr"))
}


