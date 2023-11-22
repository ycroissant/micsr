#' Binomial regression
#'
#' A unified interface perform binomial regression calling `lm` to fit
#' the linear-probability model and `glm` to fit the probit and the
#' logit model.
#'
#' @name binomreg
#' @param formula a symbolic description of the model, (for the count
#'     component and for the selection equation).
#' @param data a data frame,
#' @param subset,weights,na.action,offset see `stats::lm`,
#' @param link one of `"lm"`, `"probit"` and "`logit`" to fit
#'     respectively the linear probability, the probit and the logit
#'     model
#' @param start a vector of starting values, in this case, no
#'     estimation
#' @param object,x,type a `binomreg` object and the type of
#'     log-likelihood / residuals for the `logLik` / `residuals`
#'     method
#' @param ... further arguments
#' @param newdata a new data frame for the `predict` method
#' @return an object of class `c(`"binomreg", "micsr")`, see
#'     `micsr::micsr` for further details
#' @importFrom stats glm plogis
#' @importFrom Formula Formula
#' @examples
#' pbt <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = 'probit')
#' lpm <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = 'lm')
#' summary(pbt, vcov = "opg")
#' @export
binomreg <- function(formula, data, weights, subset, na.action, offset,
                     link = c("lm", "probit", "logit"), start = NULL, ...){
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
    y <- model.response(mf)
    q <- 2 * y - 1
    K <- ncol(X)
    N <- length(y)
    .df.residual <- N - K
    .null_deviance <- - 2 * N * (mean(y) * log(mean(y)) + (1 - mean(y)) * log(1 - mean(y)))
    if (.link == "lm"){
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, information = FALSE, sum = TRUE, X, y){
            K <- ncol(X)
            N <- length(y)
            beta <- coefs[1:K]
            sig <- coefs[K + 1]
            linpred <- drop(X %*% beta)
            lnl <- dnorm(y, mean = linpred, sd = sig, log = TRUE)
            if (gradient){
                grad <- cbind((y - linpred) / sig ^ 2 * X, - 1 / sig + (y - linpred) ^ 2 / sig ^ 3)
            }
            if (hessian){
                hess_bb <- - crossprod(X) / sig ^ 2
                hess_ss <- - 2 / sig ^ 2
                hess <- rbind(cbind(hess_bb, rep(0, K)),
                              c(rep(0, K), hess_ss))
            }
            if (information){
                info_bb <- crossprod(X) / sig ^ 2
                info_ss <- 2 / sig ^ 2
                info <- rbind(cbind(info_bb, rep(0, K)),
                              c(rep(0, K), info_ss))
            }
            if (sum){
                lnl <- sum(lnl)
                if (gradient) grad <- apply(grad, 2, sum)
            }
            if (gradient) attr(lnl, "gradient") <- grad
            if (hessian) attr(lnl, "hessian") <- hess
            if (information) attr(lnl, "ifno") <- info
            lnl
        }
    }
                
    if (.link == "probit"){
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, information = FALSE, sum = TRUE, X, y){
            linpred <- drop(X %*% coefs)
            q <- 2 * y - 1
            lnl <- pnorm(q * linpred, log.p = TRUE)
            if (gradient) grad <- q * mills(q * linpred) * X
            if (hessian) hess <-    -      crossprod(sqrt(- dmills(q * linpred)) * X)
            if (information) info <- solve(crossprod(sqrt(mills(linpred) * mills(- linpred)) * X))
            if (sum){
                lnl <- sum(lnl)
                if (gradient) grad <- apply(grad, 2, sum)
            }
            if (gradient) attr(lnl, "gradient") <- grad
            if (hessian) attr(lnl, "hessian") <- hess
            if (information) attr(lnl, "info") <- info
            lnl
        }
    }

    if (.link == "logit"){
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, information = FALSE, sum = TRUE, X, y){
            linpred <- drop(X %*% coefs)
            q <- 2 * y - 1
            lnl <- plogis(q * linpred, log.p = TRUE)
            if (gradient) grad <- (y - exp(linpred) / (1 + exp(linpred))) * X
            if (hessian) hess <-     -     crossprod(sqrt( exp(linpred) / (1 + exp(linpred)) ^ 2) * X)
            if (information) info <- solve(crossprod(sqrt( exp(linpred) / (1 + exp(linpred)) ^ 2) * X))
            if (sum){
                lnl <- sum(lnl)
                if (gradient) grad <- apply(grad, 2, sum)
            }
            if (gradient) attr(lnl, "gradient") <- grad
            if (hessian) attr(lnl, "hessian") <- hess
            if (information) attr(lnl, "info") <- info
            lnl
        }
    }
    if (is.null(start)){
        start <- rep(0, K)
        names(start) <- colnames(X)
        if (.link != "lm") .coefs <- newton(lnl, X = X, y = y, trace = 0, coefs = start, direction = "max")
        else{
            .coefs <- drop(solve(crossprod(X), crossprod(X, y)))
            .sigma <- sqrt(mean((y - drop(X %*% .coefs) ^ 2)))
            .coefs <- c(.coefs, .sigma)
        }
    }
    else .coefs <- start
    if (.link != "lm") .linpred <- drop(X %*% .coefs)
    else .linpred <- drop(X %*% .coefs[- (ncol(X) + 1)])
    .lnl_conv <- lnl(.coefs, X = X, y = y, gradient = TRUE, hessian = TRUE, info = TRUE, sum = FALSE)
    if (.link == "logit") .fitted <- plogis(.linpred)
    if (.link == "probit") .fitted <- pnorm(.linpred)
    if (.link == "lm") .fitted <- .linpred
    .npar <- c(covariates = K)
    if (.link == "lm") .npar <- c(covariates = K, vcov = 1)
    attr(.npar, "default") <- "covariates"
    .logLik <- structure(sum(as.numeric(.lnl_conv)), nobs = length(y), df = length(.coefs), class = "logLik")
    result <- list(coefficients = .coefs,
                   model = mf,
                   gradient = attr(.lnl_conv, "gradient"),
                   hessian = attr(.lnl_conv, "hessian"),
                   info = attr(.lnl_conv, "info"),
                   linear.predictors = .linpred,
                   logLik = as.numeric(.lnl_conv),
                   fitted.values = .fitted,
                   df.residual = .df.residual,
                   est_method = "ml",
                   formula = formula,
                   npar = .npar,
                   value = .logLik,
                   call = .call
                   )
    structure(result, class = c("binomreg", "micsr"))
}

#' @rdname binomreg
#' @export
logLik.binomreg <- function(object, ..., type = c("model", "null")){
    .type <- match.arg(type)
    if (.type == "model")
        result <- structure(object$value, nobs = nobs(object), df = sum(object$npar), class = "logLik")
    if (.type == "null"){
        yb <- mean(model.response(model.frame(object)))
        result <- nobs(object) * (yb * log(yb) + (1 - yb) * log(1 - yb))
        result <- structure(result, nobs = nobs(object), df = 1, class = "logLik")
    }
    result
}

#' @rdname binomreg
#' @export
residuals.binomreg <- function(object, ..., type = c("deviance", "pearson", "response")){
    type <- match.arg(type)
    y <- model.response(model.frame(object))
    hy <- fitted(object)
    sd_hy <- sqrt(hy * (1 - hy))
    if (type == "response") resid <- y - hy
    if (type == "pearson") resid <- (y - hy) / sd_hy
    if (type == "deviance") resid <- (2 * y - 1) * sqrt(- 2  * object$logLik)
    resid
}

#' @rdname binomreg
#' @method glance binomreg
#' @export
glance.binomreg <- function(x, ...){
    .est_method <- x$est_method
    N <- nobs(x)
    result <- data.frame(nobs = nobs(x))
    if (.est_method == "ml"){
        result <- data.frame(null.deviance = deviance(x, type = "null"),
                             df.null = nobs(x) - 1,
                             logLik = logLik(x),
                             AIC = AIC(x),
                             BIC = BIC(x),
                             deviance = deviance(x, type = "model"),
                             df.residual = df.residual(x),
                             nobs = nobs(x))
    }
#    if (.est_method == "twosteps"){
#        result <- data.frame(nobs = nobs(x), impliedsigma = x$sigma)
#    }
    result
}

#' @rdname binomreg
#' @export
predict.binomreg <- function(object, ..., type = c("response", "link"), newdata = NULL){
    .type <- match.arg(type)
    .link <- object$call$link
    if (is.null(newdata)){
        if (.type == "response") result <- object$fitted
        if (.type == "link") result <- object$linear.predictors
    }
    else{
        .formula <- Formula(object$formula)
        mf <- model.frame(.formula, newdata, dot = "previous")
        X <- model.matrix(.formula, mf, rhs = 1)
        y <- model.part(.formula, newdata, lhs = 1, drop = TRUE)
        result <- drop(X %*% coef(object))
        cum_fun <- switch(.link,
                          "logit" = plogis,
                          "probit" = pnorm,
                          "lm" = function(x) x)
        if (.type == "response") result <- cum_fun(result)
    }
    result
}

## #' @rdname binomreg
## #' @export
## get_predict.binomreg <- function(model,
##                                  newdata = insight::get_data(model),
##                                  vcov = NULL,
##                                  conf_level = 0.95,
##                                  type = "response",
##                                  ...) {

##     out <- stats::predict(model, type = type, newdata = newdata)
##     out <- data.frame(rowid = seq_len(length(out)), predicted = out)
##     return(out)
## }

# ajout dans le fichier type_dictionary.R
# ajout du fichier methods_binomreg
# sanity_model.R ajouter dans la liste des modèles supportés


## binomreg <- function(formula, data, weights, subset, na.action, offset, model = c("lm", "probit", "logit"), ...){
##     .formula <- Formula(formula)
##     .call <- match.call()
##     .model <- match.arg(model)
##     cl <- match.call(expand.dots = FALSE)
##     m <- match("model", names(cl), 0L)
##     cl <- cl[- m]
##     if (.model == "lm") cl[[1L]] <- as.name("lm")
##     else cl[[1L]] <- as.name("glm")
##     if (.model == "probit") dist <- binomial(link = 'probit')
##     if (.model == "logit") dist <- binomial(link = 'logit')
##     if (.model != "lm") cl$family <- dist
##     fitted_model <- eval(cl, parent.frame())
##     mf <- model.frame(fitted_model)
##     tm <- terms(fitted_model)
##     .coefs <- coef(fitted_model)
##     X <- model.matrix(tm, mf)
##     y <- model.response(mf)
##     q <- 2 * y - 1
##     .linpred <- drop(X %*% .coefs)
##     K <- ncol(X)
##     N <- length(y)
    
##     if (model == "probit"){
##         .gradient <- q * mills(q * .linpred) * X
##         .hessian <- - crossprod(sqrt(- dmills(q * .linpred)) * X)
##         .info <- crossprod(sqrt(mills(.linpred) * mills(- .linpred)) * X)
##         .logLik <- pnorm(q * .linpred, log.p = TRUE)
##     }
##     if (model == "logit"){
##         .gradient <- (y - exp(.linpred) / (1 + exp(.linpred))) * X
##         .hessian <- - crossprod(sqrt( exp(.linpred) / (1 + exp(.linpred)) ^ 2) * X)
##         .info <- - .hessian
##         .logLik <- y * .linpred - log(1 + exp(.linpred))
##     }
##     if (model == "lm"){
##         .sigma <- sqrt(deviance(fitted_model) / length(y))
##         .gradient <- cbind( (y - .linpred) * X / .sigma ^ 2,
##                            - 1 / .sigma + (y - .linpred) ^ 2 / .sigma ^ 3)
##         .hessian_bb <- - crossprod(X) / .sigma ^ 2
##         .hessian_bs <- rep(0, length(.coefs))
##         .hessian_ss <- - 2 * length(y) / .sigma ^ 2
##         .hessian <- rbind(cbind(.hessian_bb, sigma = .hessian_bs),
##                           sigma = c(t(.hessian_bs), .hessian_ss))
##         .info <- - .hessian
##         .logLik <- dnorm(.linpred, log = TRUE)
##         .coefs <- c(.coefs, sigma = .sigma)
##     }
##     .value <- structure(sum(.logLik), nobs = length(y), df = length(.coefs), class = "logLik")
##     .npar <- c(covariates = K)
##     if (model == "lm") .npar <- c(.npar, vcov = 1)
##     result <- list(coefficients = .coefs,
##                    fitted.values = .linpred,
##                    gradient = .gradient,
##                    logLik = .logLik,
##                    hessian = .hessian,
##                    value = .value,
##                    est_method = "ml",
##                    call = .call,
##                    info = .info,
##                    model = fitted_model$model,
##                    npar = .npar,
##                    formula = .formula)
##     structure(result, class  = "micsr")
## }

    
