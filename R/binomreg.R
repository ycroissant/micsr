#' Binomial regression
#'
#' A unified interface for binomial regression models, including
#' linear probability, probit and logit models
#'
#' @name binomreg
#' @param formula a symbolic description of the model
#' @param data a data frame,
#' @param subset,weights,na.action,offset,contrasts see `stats::lm`,
#' @param link one of `"identity"`, `"probit"` and "`logit`" to fit
#'     respectively the linear probability, the probit and the logit
#'     model
#' @param method `"ml"` for maximum likelihood (the only relevant
#'     method for a regression without instrumental variables),
#'     `"twosteps"` for two-steps estimator, `"minchisq"` for minimum
#'     chi-squared estimator and `"test"` to get the exogeneity test,
#' @param start a vector of starting values
#' @param robust only when `method = "twosteps"`, should the robust
#'     covariance matrix be computed?
#' @param opt optimization method
#' @param maxit maximum number of iterations
#' @param trace printing of intermediate result
#' @param object,x,type a `binomreg` object and the type of residuals
#'     for the `residuals` method
#' @param check_gradient if `TRUE` the numeric gradient and hessian
#'     are computed and compared to the analytical gradient and
#'     hessian
#' @param ... further arguments
#' @param newdata a new data frame for the `predict` method
#' @return an object of class `c("binomreg", "micsr")`, see
#'     `micsr::micsr` for further details
#' @importFrom stats glm plogis qlogis dlogis
#' @importFrom Formula Formula model.part
#' @importFrom dplyr case_when
#' @keywords models
#' @examples
#' pbt <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = 'probit')
#' lpm <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = 'identity')
#' summary(pbt, vcov = "opg")
#' @export
binomreg <- function(formula, data, weights, subset, na.action, offset, contrasts = NULL,
                     link = c("identity", "probit", "logit"),
                     method = c("ml", "twosteps", "minchisq", "test"),
                     start = NULL,
                     robust = TRUE,
                     opt = c("newton", "nr", "bfgs"),
                     maxit = 100, trace = 0,
                     check_gradient = FALSE, ...){
    .opt <- match.arg(opt)
    .method <- match.arg(method)
    .call <- match.call()
    .link <- match.arg(link)
    .linkfun <- switch(.link,
                       probit = function(mu) qnorm(mu),
                       logit = function(mu) qlogis(mu))
    .linkinv <- switch(.link,
                       probit = function(mu) pnorm(mu),
                       logit = function(mu) plogis(mu))
    .mu.eta <- switch(.link,
                      probit = function(mu) dnorm(mu),
                      logit = function(mu) dlogis(mu))

    mf <- match.call(expand.dots = FALSE)
    .formula <- Formula(formula)
    if (length(.formula)[2] == 2){
        mf$model <- "probit"
        mf$method <- .method
        mf[[1L]] <- as.name("ivldv")#quote(micsr::ivldv())
        result <- eval(mf, parent.frame())
        result$call <- .call
        return(result)
    } else {
        if (.method != "ml")
            stop("with a one-part formula, the only relevant method is ml")
    }
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    # construct the model frame and components
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf)
    wt <- as.vector(model.weights(mf))
    if (is.null(wt)) wt <- 1 else wt <- wt# / mean(wt)
    .offset <- model.offset(mf)
    X <- model.matrix(mt, mf, contrasts)
    if (compute_rank(X) < ncol(X)){
        .rank <- compute_rank(X)
        .ncol <- ncol(X)
        stop(paste("the rank of X = ", .rank, " < the number of columns = ", .ncol, sep = ""))
    }
    yb <- mean(y)
    q <- 2 * y - 1
    K <- ncol(X)
    N <- length(y)
    .df.residual <- N - K
    .null_logLik <- N * (yb * log(yb) + (1 - yb) * log(1 - yb))
    .sat_logLik <- 0
    .null_deviance <- - 2 * .null_logLik
        
    if (.link == "identity"){
        .null_intercept <- yb
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, information = FALSE,
                        opposite = FALSE, sum = TRUE, X, y, weights){
            sgn <- ifelse(opposite, - 1, + 1)
            K <- ncol(X)
            N <- length(y)
            beta <- coefs[1:K]
            sig <- coefs[K + 1]
            linpred <- drop(X %*% beta)
            lnl <- dnorm(y, mean = linpred, sd = sig, log = TRUE)
            if (gradient){
                grad <- cbind((y - linpred) / sig ^ 2 * X,
                              - 1 / sig + (y - linpred) ^ 2 / sig ^ 3)
            }
            if (hessian){
                hess_bb <- - crossprod(sqrt(weights) * X) / sig ^ 2
                hess_ss <- - 2 / sig ^ 2
                hess_ss <- 1 / sig ^ 2 * sum(weights) - 3 * sum(weights * (y - linpred) ^  2) / sig ^ 4
                hess <- rbind(cbind(hess_bb, sigma = rep(0, K)),
                              sigma = c(rep(0, K), hess_ss))
            }
            if (information){
                info_bb <- crossprod(sqrt(weights) * X) / sig ^ 2
                info_ss <- 2 / sig ^ 2
                info <- rbind(cbind(info_bb, rep(0, K)),
                              c(rep(0, K), info_ss))
            }
            if (sum){
                lnl <- sum(weights * lnl)
                if (gradient) grad <- apply(weights * grad, 2, sum)
            }
            if (gradient) attr(lnl, "gradient") <- grad * sgn
            if (hessian) attr(lnl, "hessian") <- hess * sgn
            if (information) attr(lnl, "info") <- info
            lnl * sgn
        }
    }
                
    if (.link == "probit"){
        .null_intercept <- qnorm(yb)
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, information = FALSE,
                        opposite = FALSE, sum = TRUE, X, y, weights){
            sgn <- ifelse(opposite, +1, -1)
            linpred <- drop(X %*% coefs)
            q <- 2 * y - 1
            lnl <- pnorm(q * linpred, log.p = TRUE)
            if (gradient) grad <-  q * mills(q * linpred) * X
            if (hessian) hess <- -   crossprod(sqrt(- weights *  mills(q * linpred, 1)) * X)
            if (information) info <- crossprod(sqrt(weights * mills(linpred) * mills(- linpred)) * X)
            if (sum){
                lnl <- sum(lnl * weights)
                if (gradient) grad <- apply(grad * weights, 2, sum)
            }
            if (gradient) attr(lnl, "gradient") <- grad * sgn
            if (hessian) attr(lnl, "hessian") <- hess * sgn
            if (information) attr(lnl, "info") <- info
            lnl * sgn
        }

        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, information = FALSE,
                        opposite = FALSE, sum = TRUE, X, y, weights){
            sgn <- ifelse(opposite, - 1, + 1)
            linpred <- drop(X %*% coefs)
            q <- 2 * y - 1
            lnl <- pnorm(q * linpred, log.p = TRUE)
            if (gradient) grad <-  q * mills(q * linpred) * X
            if (hessian) hess <- -   crossprod(sqrt(- weights *  mills(q * linpred, 1)) * X)
            if (information) info <- crossprod(sqrt(weights * mills(linpred) * mills(- linpred)) * X)
            if (sum){
                lnl <- sum(lnl * weights)
                if (gradient) grad <- apply(grad * weights, 2, sum)
            }
            if (gradient) attr(lnl, "gradient") <- grad * sgn
            if (hessian) attr(lnl, "hessian") <- hess * sgn
            if (information) attr(lnl, "info") <- info
            lnl * sgn
        }

    }

    if (.link == "logit"){
        .null_intercept <- qlogis(yb)
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, information = FALSE,
                        opposite = FALSE, sum = TRUE, X, y, weights){
            sgn <- ifelse(opposite, - 1, + 1)
            linpred <- drop(X %*% coefs)
            q <- 2 * y - 1
            lnl <- plogis(q * linpred, log.p = TRUE)
            if (gradient) grad <- (y - exp(linpred) / (1 + exp(linpred))) * X
            if (hessian) hess <-  -   crossprod(sqrt(weights * exp(linpred) / (1 + exp(linpred)) ^ 2) * X)
            if (information) info <- crossprod(sqrt(weights * exp(linpred) /  (1 + exp(linpred)) ^ 2) * X)
            if (sum){
                lnl <- sum(lnl * weights)
                if (gradient) grad <- apply(grad * weights, 2, sum)
            }
            if (gradient) attr(lnl, "gradient") <- grad * sgn
            if (hessian) attr(lnl, "hessian") <- hess * sgn
            if (information) attr(lnl, "info") <- info
            lnl * sgn
        }
    }
    ## if (is.null(start)){
    ##     start <- rep(0, K)
    ##     names(start) <- colnames(X)
    ##     if (.link != "identity"){
    ##         .coefs <- newton(lnl, X = X, y = y, weights = wt,
    ##                          trace = 0, coefs = start, direction = "max")
    ##     }
    ##     else{
    ##         .coefs <- drop(solve(crossprod(sqrt(wt) * X), crossprod(wt * X, y)))
    ##         .sigma <- sqrt(mean(wt * (y - drop(X %*% .coefs) ^ 2)))
    ##         .coefs <- c(.coefs, sigma = .sigma)
    ##     }
    ## }    
    ## else .coefs <- start

    if (is.null(start)){
        start <- rep(0, K)
    }
    names(start) <- colnames(X)
    if (maxit > 0){
        if (.link != "identity"){
##             .coefs <- newton(lnl, X = X, y = y, weights = wt,
##                              trace = 1, coefs = start, direction = "max")
## stop()
            .coefs <- maximize(lnl, start = start, method = .opt, trace = trace, maxit = maxit,
                               X = X, y = y, weights = wt, ...)
        }
        else{
            .coefs <- drop(solve(crossprod(sqrt(wt) * X), crossprod(wt * X, y)))
            .sigma <- sqrt(mean(wt * (y - drop(X %*% .coefs) ^ 2)))
            .coefs <- c(.coefs, sigma = .sigma)
        }
    }   
    else .coefs <- start
    
    if (.link != "identity") .linpred <- drop(X %*% .coefs)
    else .linpred <- drop(X %*% .coefs[- (ncol(X) + 1)])


    .lnl_conv <- lnl(.coefs, X = X, y = y, weights = wt, gradient = TRUE,
                     hessian = TRUE, info = TRUE, sum = FALSE)
    fun <- function(x) lnl(x, X = X, y = y, weights = wt, gradient = TRUE,
                           hessian = TRUE, info = TRUE, sum = TRUE)
    if (check_gradient) z <- check_gradient(fun, .coefs) else z <- NA

    l_model <- as.numeric(.lnl_conv)
    l_saturated <- 0
    l_null <- yb * log(yb) + (1 - yb) * log(1 - yb)
    values <- cbind(model = l_model, saturated = l_saturated, null = l_null)
    .logLik <- apply(values * wt, 2, sum)
    
    if (.link == "logit") .fitted <- plogis(.linpred)
    if (.link == "probit") .fitted <- pnorm(.linpred)
    if (.link == "identity") .fitted <- .linpred
    .npar <- c(covariates = K)
    if (.link == "identity") .npar <- c(covariates = K, vcov = 1)
    attr(.npar, "default") <- "covariates"
    
    # Null model
    null_coefs <- rep(0, ncol(X))
    names(null_coefs) <- colnames(X)
    null_coefs["(Intercept)"] <- .null_intercept
    if (.link == "identity") null_coefs <- c(null_coefs, sigma = sqrt(mean((y - mean(y)) ^ 2)))
    lnl_null <- lnl(null_coefs, gradient = TRUE, hessian = TRUE, information = TRUE,
                    X = X, y = y, weights = wt)
    .null_gradient <- attr(lnl_null, "gradient")
    .null_info <- attr(lnl_null, "info")
    .lm <- drop(crossprod(.null_gradient, solve(.null_info, .null_gradient)))
    .model_info <- attr(.lnl_conv, "info")
    .vcov <- solve(.model_info)
    .wald <- drop(crossprod(.coefs[- 1],
                            t(crossprod(.coefs[-1],
                                        solve(.vcov[- 1, - 1, drop = FALSE])))))
    .lr <- 2 * unname(.logLik["model"] - .logLik["null"])
    tests <- c(wald = .wald, score = .lm, lr = .lr)

    gres <- (y - pnorm(.linpred)) * dnorm(.linpred) / (pnorm(.linpred) * (1 - pnorm(.linpred)))
    .family <- structure(list(
        link = .link,
        linkfun = .linkfun,
        linkinv = .linkinv,
        mu.eta = .mu.eta,
        variance = function(mu) mu * (1 - mu),
        family = "binomial"),
        class = "family")
    
    result <- list(coefficients = .coefs,
                   model = mf,
                   terms = mt,
                   value = values,
                   gradient = attr(.lnl_conv, "gradient"),
                   hessian = attr(.lnl_conv, "hessian"),
                   info = attr(.lnl_conv, "info"),
                   fitted.values = .fitted,
                   residuals = gres,
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
                   check_gradient = z,
                   family = .family)
    structure(result, class = c("binomreg", "micsr"))
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



# ajout dans le fichier type_dictionary.R
# ajout du fichier methods_binomreg
# sanity_model.R ajouter dans la liste des modèles supportés

