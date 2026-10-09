#' Weibull regression model for duration data
#'
#' The Weibull model is the most popular model for duration data. This
#' function enables the estimation of this model with two alternative
#' (but equivalent) parametrization: the Accelerate Failure Time and
#' the Proportional Hazard. Moreover heterogeneity can be introduced,
#' which leads to the Gamma-Weibull model
#' @name weibreg
#' @param formula a symbolic description of the model
#' @param data a data frame
#' @param subset,weights,na.action,offset,contrasts see `stats::lm`,
#' @param start a vector of starting values
#' @param model one of `"aft"` or `"ph"`
#' @param opt the optimization method
#' @param robust a boolean if `TRUE`, the log of the shape and the
#'     variance parameters are estimated
#' @param maxit maximum number of iterations
#' @param trace an integer
#' @param mixing if `TRUE`, the Gamma-Weibull model is estimated
#' @param check_gradient if `TRUE` the numeric gradient and hessian
#'     are computed and compared to the analytical gradient and
#'     hessian
#' @param ... further arguments
#' @param x,object a `weibreg` object
#' @param vcov the covariance matrix estimator to use for the score
#'     test
#' @return an object of class `c("weibreg", "micsr")`, see
#'     `micsr::micsr` for further details.
#' @importFrom stats glm plogis nlm
#' @importFrom Formula Formula
#' @importFrom numDeriv hessian grad
#' @importFrom survival Surv
#' @keywords models
#' @examples
#' library(survival)
#' wz <- weibreg(Surv(duration, censored == "no") ~ gender + age + log(wage + 1),
#'          unemp_duration, mixing = TRUE, model = "ph")
#' @export
weibreg <- function(formula, data, weights, subset, na.action, offset, contrasts = NULL,
                    model = c("aft", "ph"), opt = c("bfgs", "newton", "nr"), start = NULL, maxit = 100,
                    robust = TRUE, trace = 0, mixing = FALSE, check_gradient = FALSE, ...){
    .start <- start
    .robust <- robust
    .mixing <- mixing
    .model <- match.arg(model)
    .opt <- match.arg(opt)
    .call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    .formula <- mf$formula# <- Formula(formula)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    # construct the model frame and components
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- as.matrix(model.response(mf))
    wt <- as.vector(model.weights(mf))
    if (is.null(wt)) wt <- 1 else wt <- wt / mean(wt)
    .offset <- model.offset(mf)
    e <- as.logical(y[, 2])
    y <- y[, 1]
    X <- model.matrix(mt, mf, contrasts)
    start_null <- c(ifelse(.model == "aft", 1, -1) * mean(log(y)), ifelse(robust, 0, 1))
    lnl_fun <- ifelse(mixing, lnl_mweibull, lnl_weibull)

    if (is.null(.start)){
        start_lm <- coef(lm(log(y) ~ X - 1, subset = e))
        names(start_lm) <- colnames(X)
        if (.model == "ph") start_lm <- - start_lm
        .start <- c(start_lm, shape = ifelse(robust, 0, 1))
        if (mixing) .start <- c(.start, var =  ifelse(robust, 0, 1))
    }
    if (maxit > 0){
        coefs <- maximize(lnl_fun, start = .start, model = .model, method = .opt,
                          X = X, y = y, e = e, weights = wt, robust = .robust, trace = trace)
    } else {
        coefs <- .start
    }
    names(coefs) <- names(.start)
    if (.robust){
        coefs[ncol(X) + 1] <- exp(coefs[ncol(X) + 1])
        if (.mixing) coefs[ncol(X) + 2] <- exp(coefs[ncol(X) + 2])
    }
    
    lnl_conv <- lnl_fun(coefs, sum = FALSE, opposite = FALSE, model = .model,
                        X = X, y = y, e = e, weights = wt,
                        robust = FALSE, trace = trace)
    
    if (check_gradient){
        fun <- function(x) lnl_fun(x, sum = TRUE, opposite = FALSE, model = .model,
                                   X = X, y = y, e = e, weights = wt,
                                   robust = FALSE, trace = trace)    
        z <- check_gradient(fun, coefs)
    } else z <- NA

    .linpred <- drop(X %*% coefs[1:ncol(X)])
    .fitted <- exp(.linpred + lgamma(1 + 1 / coefs["shape"]))
    result <- list(coefficients = coefs,
                   model = mf,
                   terms = mt,
                   value = as.numeric(lnl_conv),
                   gradient = attr(lnl_conv, "gradient"),
                   hessian = attr(lnl_conv, "hessian"),
                   info = attr(lnl_conv, "info"),
                   linear.predictor = .linpred,
                   fitted.values = .fitted,
                   logLik = c(model = sum(as.numeric(lnl_conv))),#, null = lnl_null),
                   npar = structure(c(covariates = ncol(X), vcov = ifelse(mixing, 2, 1)),
                                    default = c("covariates", "vcov")),
                   est_method = "ml",
                   call = .call,
                   na.action = attr(mf, "na.action"),
                   weights = wt,
                   offset = .offset,
                   contrasts = attr(X, "contrasts"),
                   xlevels = .getXlevels(mt, mf),
                   check_gradient = z)
    structure(result, class = c("weibreg", "micsr", "lm"))
}
    

lnl_weibull <- function(param, sum = TRUE, gradient = TRUE, hessian  = TRUE, info = TRUE, opposite = TRUE,
                        model = c("aft", "ph"), X, y, e, weights, robust = FALSE, ...){
    sgn <- ifelse(opposite, -1, + 1)
    d <- e
    beta <- param[1:ncol(X)]
    
    delta <- param[ncol(X) + 1]
    bX <- as.numeric(X %*% beta)
    .model <- match.arg(model)
    if (robust) delta <- exp(delta)
    if (.model == "aft"){
        A <- exp(delta * (log(y) - bX))
        lnl <- d * (log(delta) - delta * bX + (delta - 1) *  log(y)) - A
    }
    if (.model == "ph"){
        A <- y ^ delta * exp(bX)
        lnl <- d * (log(delta) + bX + (delta - 1) *  log(y)) - A
    }
    lnl <- sgn * lnl
    if (sum) lnl <- sum(lnl * weights)
    
    if (gradient){
        if (.model == "aft"){
            g_b <- (A - d) *  delta
            g_d <- d * (1 / delta + (log(y) - bX)) - (log(y) - bX) * A
        }
        if (.model == "ph"){
            g_b <- d - A
            g_d <- d / delta + (d - A) * log(y)
        }
        if (robust) g_d <- g_d * delta
        grad <- cbind(X * g_b, shape = g_d)
        grad <- sgn * grad
        if (sum) grad <- apply(grad * weights, 2, sum)
        attr(lnl, "gradient") <- grad
    }
    
    if (hessian){
        if (.model == "aft"){
            h_bb <- - A * delta ^ 2
            h_bd <- - d  + A * (1 + delta * (log(y) - bX))
            h_dd <- - d / delta ^ 2 - (log(y) - bX) ^ 2 * A
        }
        if (.model == "ph"){
            h_bb <- - A
            h_bd <- - A * log(y)
            h_dd <- - d / delta ^ 2 - log(y) ^ 2 * A
        }
        if (robust) h_dd <- h_dd * delta ^ 2 + g_d
        if (robust) h_bd <- h_bd * delta
        h_bb <- crossprod(weights * h_bb * X, X)
        h_bd <- apply(h_bd * X * weights, 2, sum)
        h_dd <- sum(weights * h_dd)
        hessian <- rbind(cbind(h_bb, shape = h_bd),
                         shape = c(h_bd, h_dd))
        attr(lnl, "hessian") <- hessian * sgn
    }
    if (info){
        if (.model == "aft"){
            emc <- 0.5772156649015328606065120
            I1 <- (1 - emc) / delta
            I2 <- (pi ^ 2 / 6 - 2 * emc + emc ^ 2) / delta ^ 2
            i_bb <- delta ^ 2 * crossprod(weights * X, X)
            i_bd <- - apply( (- d + 1 + delta * I1) * X * weights, 2, sum)
            i_dd <- sum( (d / delta ^ 2 + I2) * weights)
            info <- rbind(cbind(i_bb, shape = i_bd),
                             shape = c(i_bd, i_dd))
            attr(lnl, "info") <- info
        }
    }
    lnl
}

lnl_mweibull <- function(param, sum = TRUE, gradient = TRUE, hessian  = TRUE, info = TRUE, opposite = TRUE,
                         X, y, e, weights, robust = TRUE, model = c("aft", "ph"), ...){
    .model <- match.arg(model)
    sgn <- ifelse(opposite, -1, + 1)
    beta <- param[1:ncol(X)]
    delta <- param[ncol(X) + 1]
    if (robust) delta <- exp(delta)
    rho <- param[ncol(X) + 2]
    if (robust) rho <- exp(rho)
    bX <- as.numeric(X %*% beta)
    d <- e
    loglambda <- log(y) - bX
    if (abs(rho) > 1E-8){
        if (.model == "ph"){
            A <- exp(bX) * y ^ delta
            lnl <- d * (log(delta) + (delta - 1) * log(y) +         bX) - (1 / rho + d) * log(1 + rho * A)
        }
        if (.model == "aft"){
            A <- (y / exp(bX)) ^ delta
            lnl <- d * (log(delta) + (delta - 1) * log(y) - delta * bX) - (1 / rho + d) * log(1 + rho * A)
        }
        if (opposite) lnl <- - lnl
        if (sum) lnl <- sum(lnl * weights)
        if (gradient){
            if (.model == "ph"){
                g_b <- d - (1 / rho + e) * rho * A / (1 + rho * A)
                g_d <- d * (1 / delta + log(y)) - (1 / rho + d) * log(y) * rho * A / (1 + rho * A)
                g_r <- 1 / rho ^ 2 * log(1 + rho * A) - (1 / rho + d) * A / (1 + rho * A)
            }
            if (.model == "aft"){
                g_b <- - d * delta + (1 / rho + d) * (rho * delta * A) / (1 + rho *  A)
                g_d <- d * (1 / delta + loglambda) - (1 / rho + d) * (rho * loglambda * A) / (1 + rho * A)
                g_r <- log(1 + rho * A) / rho ^ 2 - (1 / rho + d) * A / (1 + rho * A)
            }
            if (robust){
                g_d <- g_d * delta
                g_r <- g_r * rho
            }
            grad <- cbind(g_b * X, shape = g_d, var =  g_r)
            if (opposite) grad <- - grad
            if (sum) grad <- apply(grad * weights, 2, sum)
            attr(lnl, "gradient") <- grad
        }
        if (hessian){
            if (.model == "ph"){
                h_bb <- - (1 / rho + d)                 * rho * A / (1 + rho * A) ^ 2
                h_bd <- - (1 / rho + d)                 * rho * A / (1 + rho * A) ^ 2 * log(y)
                h_br <-  (A - d)                              * A / (1 + rho * A) ^ 2
                h_dd <- - d / delta ^ 2 - (1 / rho + d) * rho * A / (1 + rho * A) ^ 2 * log(y) ^ 2
                h_dr <-  (A - d)                              * A / (1 + rho * A) ^ 2 * log(y)
                h_rr <- - 2 / rho ^ 3 * log(1 + rho * A) + 2 / rho ^ 2 * A / (1 + rho * A) +
                    (1 / rho + d) * (A / (1 + rho * A)) ^ 2
            }
            if (.model == "aft"){
                h_bb <- - (1 / rho + d) * delta ^  2 * (rho * A) / (1 + rho * A) ^ 2
                h_bd <- - d + (1 / rho + d) * (rho * A) / ( 1 + rho * A) * (1 +  delta * loglambda / (1 + rho * A))
                h_br <-  delta * A * (d - A) / (1 + rho * A) ^ 2
                h_dd <- - d / delta ^ 2 - (1 / rho + d) * loglambda ^ 2 * rho * A / (1 + rho * A) ^ 2
                h_dr <- loglambda * A * (A - d) / (1 + rho * A) ^ 2
                h_rr <- - 2 * log(1 + rho * A) / rho ^ 3 + (2 * A + 3 * rho * A ^ 2 + d * rho ^ 2 * A ^ 2) / (rho ^ 2 * (1 + rho * A) ^ 2)
            }
            if (robust){
                h_dd <- h_dd * delta ^ 2 + g_d
                h_rr <- h_rr * rho  ^ 2 + g_r
                h_bd <- h_bd * delta
                h_br <- h_br * rho
                h_dr <- h_dr * rho * delta
            }
            h_bb <- crossprod(weights * h_bb * X, X)
            h_bd <- apply(weights * h_bd * X, 2, sum)
            h_br <- apply(weights * h_br * X, 2, sum)
            hessian <- rbind(cbind(h_bb, shape = h_bd, var = h_br),
                             shape = c(h_bd, sum(weights * h_dd), sum(weights * h_dr)),
                             var = c(h_br, sum(weights * h_dr), sum(weights * h_rr)))
            if (opposite) hessian <- - hessian
            attr(lnl, "hessian") <- hessian
        }
    } else {
        if (.model == "ph"){
            A <- exp(bX) * y ^ delta * rho
            lnl <- e * (log(delta) + (delta - 1) * log(y) + bX) - (1 / rho + e) * log(1 + A)
        }
        if (.model == "aft"){
            lnl <- d * (log(delta) - delta * loglambda + (delta - 1) * log(y)) - (1 / rho + d) * log(1 + delta * A)
        }
        if (opposite) lnl <- - lnl
        if (sum) lnl <- sum(lnl * weights)

        if (gradient){
            g_b <- - d * delta  + delta * A
            g_s <- d * (1 / delta + loglambda) - loglambda * A
            g_r <- 1 / 2 * A ^ 2 - 2 / 3 * rho * A ^ 3 - d * (A - rho * A ^ 2 + rho ^ 2 * A ^ 3)
            grad <- cbind(g_b * X, shape = g_s, var =  g_r)
            if (opposite) grad <- - grad
            if (sum) grad <- apply(grad * weights, 2, sum)
            attr(lnl, "gradient") <- grad
        }
        if (hessian){
            h_bb <- - crossprod(weights * delta ^ 2 *  A * X, X)
            h_bd <- apply(weights * (- d + A * (1 + delta * loglambda)) * X, 2, sum)
            h_br <- apply(weights * delta * A * (d - A) * X, 2, sum)
            h_dd <- sum(weights * (- d / delta ^ 2 - loglambda ^ 2 * A))
            h_dr <- sum(weights * loglambda * A * (A - d))
            h_rr <- sum(weights * (d * A ^ 2 - 2 / 3 * A ^ 3))
            hessian <- rbind(cbind(h_bb, shape = h_bd, var = h_br),
                             shape = c(h_bd, h_dd, h_dr),
                             var = c(h_br, h_dr, h_rr))
            if (opposite) hessian <- - hessian
            attr(lnl, "hessian") <- hessian
        }
        if (info){
            emc <- 0.5772156649015328606065120
            I1 <- (1 - emc) / delta
            I2 <- (pi ^ 2 / 6 - 2 * emc + emc ^ 2) / delta ^ 2
            i_bb <- - delta ^ 2 * crossprod(weights * X, X)
            i_bs <- apply(weights * (- d + (1 + delta * I1)) * X, 2, sum)
            i_br <- apply(weights * delta * (d - 2) * X, 2, sum)
            i_ss <- sum(weights * (- d / delta ^ 2 - I1))
            i_sr <- sum(weights * (- I2 * d + (3 - 2 * emc) / delta))
            i_rr <- 2 * sum(weights * (d - 2))
            .info <- - rbind(cbind(i_bb, shape = i_bs, var = i_br),
                             shape = c(i_bs, i_ss, i_sr),
                             var = c(i_br, i_sr, i_rr))
            
            attr(lnl, "info") <- .info
        }
    }
    lnl
}

#' @rdname weibreg
#' @export
gres <- function(x){
    .mixing <- x$call$mixing
    if (is.null(.mixing)) .mixing <- FALSE
    .type <- ifelse(is.null(x$call$model), "aft", x$call$model)
    delta <- coef(x)["shape"]
    y <- model.response(model.frame(x))
    linpred <- x$linear.predictor
    e <- y[, 2]
    is.censored <- 1 - e
    y <- y[, 1]
    if (! .mixing){
        if (.type == "aft") gr <- exp( delta * (log(y) - linpred)) + is.censored
        if (.type == "ph") gr <-  exp(linpred) * y ^ delta         + is.censored
    }
    else{
        varia <- coef(x)["var"]
        if (.type == "aft") gr <- log(1 + varia * exp( delta * (log(y) - linpred))) / varia + is.censored
        if (.type == "ph")  gr <- log(1 + varia * exp(linpred) * y ^ delta)         / varia + is.censored
    }
    gr
}
    

#' @rdname weibreg
#' @method scoretest weibreg
#' @export
scoretest.weibreg <- function(object, ..., vcov = NULL){
    .vcov <- ifelse(is.null(vcov), "opg", vcov)
    new <- list(...)[[1]]
    cls <- class(object)[1]
    if (inherits(new, "formula")){
        class(object) <- setdiff(class(object), "weibreg")
        scoretest(object, ...)
    } else {
        new <- update(object, maxit = 0, start = c(coef(object), var = 0), mixing = TRUE, robust = FALSE)
        gs <- sum(new$gradient[, "var"])
        if (.vcov == "opg") hm1 <- solve(crossprod(new$gradient))["var", "var"]
        if (.vcov == "hessian") hm1 <- solve(- new$hessian)["var", "var"]
        if (.vcov == "info") hm1 <- 1 / new$info["var", "var"]
#        hm1 <- - solve(new$hessian)["var", "var"]
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


## lnl_cloglog <- function(param, gradient = TRUE, sum = TRUE, opposite = TRUE, hessian = TRUE){
##     sgn <- ifelse(opposite, -1, + 1)
##     beta <- param[1:ncol(X)]
##     sigma <- param[ncol(X) + 1]
##     bX <- as.numeric(X %*% beta)
##     resid <- log(y) - bX
##     prob <- exp(sigma * resid) / (1 + exp(sigma * resid))
##     ln_hazard <- log(sigma) + sigma * resid - log(y) - log(1 + exp(sigma * resid))
##     ln_surv <- - log(1 + exp(sigma * resid))
##     lnl <- (e * ln_hazard + ln_surv) * sgn
##     if (sum) lnl <- sum(lnl)
##     if (gradient){
##         g_sigma <- e * (1 / sigma + resid) - (1 + e) * resid * prob
##         g_beta <- - e * sigma + (1 + e) * sigma * prob
##         grad <- cbind(g_beta * X, g_sigma) * sgn
##         if (sum) grad <- apply(grad, 2, sum)
##         attr(lnl, "gradient") <- grad
##     }
##     if (hessian){
##         P_s <-   resid * prob * (1 - prob)
##         P_b <- - sigma * prob * (1 - prob)
##         h_ss <- - e / sigma ^ 2 - (1 + e) * resid * P_s
##         h_bs <- - e + (1 + e) * prob + (1 + e) * sigma * P_s
##         h_bb <-  (1 + e) * sigma * P_b
##         h_bb <- crossprod(h_bb * X, X)
##         h_bs <- apply(h_bs * X, 2, sum)
##         h_ss <- sum(h_ss)
##         attr(lnl, "hessian") <- rbind(cbind(h_bb, h_bs),
##                                       c(h_bs, h_ss)) * sgn
##     }
##     lnl
## }
    


## scoretest.weibreg <- function(object, ..., vcov = NULL){
##     x <- object
##     if (is.null(vcov)){
##         if (is.null(x$info)) .vcov <- "hessian" else .vcov <- "info"
##     }
##     else .vcov <- vcov
##     objects <- list(object, ...)
##     sup_objects <- list(...)
##     sup_objects_1 <- sup_objects[[1]]
##     if (inherits(sup_objects_1, "weibreg") | inherits(sup_objects_1, "formula")){
##         update_covariates <- TRUE
##     } else update_covariates <- FALSE
##     if (update_covariates & length(sup_objects) != 1)
##         stop("Two models should be provided")
##     if (update_covariates){
##         if(! inherits(objects[[2]], "formula")) objects[[2]] <- formula(objects[[2]])
##         newform <- objects[[2]]
##         nX <- model.matrix(x, newform)
##         oX <- model.matrix(x)
##         new_nms <- colnames(nX)
##         old_nms <- colnames(oX)
##         all_old_nms <- names(coef(x))
##         sup_coef <- setdiff(all_old_nms, old_nms)
##         L <- length(new_nms) - length(old_nms)
##         if(! all(old_nms %in% new_nms))
##             stop("the old model is not nested in the new one")
##         start <- coef(x)[c(new_nms, sup_coef)]
##         start[is.na(start)] <- 0
##         names(start) <- c(new_nms, sup_coef)
##         y <- update(x, formula = newform, start = start, maxit = 0)
##         .data.name <- paste(deparse(newform))
##     } else {
##         if (length(sup_objects) != 1)
##             stop("only a mixing argument should be provided")
##         if (names(sup_objects)[1] != "mixing")
##             stop("a mixing argument should be provided")
##         start <- c(coef(x), var = 0)
##         y <- update(x, start = start, mixing = TRUE, maxit = 0, robust = FALSE)
##         L <- 1
##         print(y$coefficients)
##         print(apply(y$gradient, 2, sum));stop()
##         .data.name <- "mixing = TRUE"
##     }
##     .gradient <- apply(y$gradient, 2, sum)
##     if (.vcov == "hessian") information <- - y$hessian
##     if (.vcov == "opg") information <- crossprod(y$gradient)
##     if (.vcov == "info"){
##         if (is.null(y$info)) stop("no information matrix estimate available")
##         else information <- y$info
##     }
##     .stat <- drop(crossprod(.gradient, solve(information, .gradient)))
##     structure(list(statistic = c(chisq = .stat),
##                    p.value = pchisq(.stat, lower.tail = FALSE, df = L),
##                    parameter = c(df = L),
##                    method = "score test",
##                    data.name = .data.name,
##                    alternative = "the constrained model is rejected"),
##               class = "htest")    
## }
