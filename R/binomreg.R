#' Binomial regression
#'
#' A unified interface for binomial regression models, including
#' linear probability, probit and logit models
#'
#' @name binomreg
#' @param formula a symbolic description of the model
#' @param data a data frame,
#' @param subset,weights,na.action,offset see `stats::lm`,
#' @param link one of `"identity"`, `"probit"` and "`logit`" to fit
#'     respectively the linear probability, the probit and the logit
#'     model
#' @param method `"ml"` for maximum likelihood (the only relevant
#'     method for a regression without instrumental variables),
#'     `"twosteps"` for two-steps estimator, `"minchisq"` for minimum
#'     chi-squared estimator and `"test"` to get the exogeneity test
#' @param start a vector of starting values
#' @param object,x,type a `binomreg` object and the type of residuals
#'     for the `residuals` method
#' @param ... further arguments
#' @param newdata a new data frame for the `predict` method
#' @return an object of class `c("binomreg", "micsr")`, see
#'     `micsr::micsr` for further details
#' @importFrom stats glm plogis qlogis
#' @importFrom Formula Formula model.part
#' @keywords models
#' @examples
#' pbt <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = 'probit')
#' lpm <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = 'identity')
#' summary(pbt, vcov = "opg")
#' @export
binomreg <- function(formula, data, weights, subset, na.action, offset,
                     link = c("identity", "probit", "logit"),
                     method = c("ml", "twosteps", "minchisq", "test"),
                     start = NULL, ...){
    .method <- match.arg(method)
    .call <- match.call()
    .link <- match.arg(link)
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
    m <- match(c("formula", "data", "subset", "weights"),
               names(mf), 0L)
    # construct the model frame and components
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)
    y <- model.response(mf)
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
            if (information) attr(lnl, "info") <- info
            lnl
        }
    }
                
    if (.link == "probit"){
        .null_intercept <- qnorm(yb)
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, information = FALSE, sum = TRUE, X, y){
            linpred <- drop(X %*% coefs)
            q <- 2 * y - 1
            lnl <- pnorm(q * linpred, log.p = TRUE)
            if (gradient) grad <- q * mills(q * linpred) * X
            if (hessian) hess <- -   crossprod(sqrt(- dmills(q * linpred)) * X)
            if (information) info <- crossprod(sqrt(mills(linpred) * mills(- linpred)) * X)
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
        .null_intercept <- qlogis(yb)
        lnl <- function(coefs, gradient = FALSE, hessian = FALSE, information = FALSE, sum = TRUE, X, y){
            linpred <- drop(X %*% coefs)
            q <- 2 * y - 1
            lnl <- plogis(q * linpred, log.p = TRUE)
            if (gradient) grad <- (y - exp(linpred) / (1 + exp(linpred))) * X
            if (hessian) hess <-  -   crossprod(sqrt( exp(linpred) / (1 + exp(linpred)) ^ 2) * X)
            if (information) info <- crossprod(sqrt( exp(linpred) / (1 + exp(linpred)) ^ 2) * X)
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
        if (.link != "identity") .coefs <- newton(lnl, X = X, y = y, trace = 0, coefs = start, direction = "max")
        else{
            .coefs <- drop(solve(crossprod(X), crossprod(X, y)))
            .sigma <- sqrt(mean((y - drop(X %*% .coefs) ^ 2)))
            .coefs <- c(.coefs, sigma = .sigma)
        }
    }
    else .coefs <- start
    if (.link != "identity") .linpred <- drop(X %*% .coefs)
    else .linpred <- drop(X %*% .coefs[- (ncol(X) + 1)])
    .lnl_conv <- lnl(.coefs, X = X, y = y, gradient = TRUE, hessian = TRUE, info = TRUE, sum = FALSE)
    if (.link == "logit") .fitted <- plogis(.linpred)
    if (.link == "probit") .fitted <- pnorm(.linpred)
    if (.link == "identity") .fitted <- .linpred
    .npar <- c(covariates = K)
    if (.link == "identity") .npar <- c(covariates = K, vcov = 1)
    attr(.npar, "default") <- "covariates"
    .logLik <- structure(sum(as.numeric(.lnl_conv)), nobs = length(y), df = length(.coefs), class = "logLik")
    .logLik <- c(model = sum(as.numeric(.lnl_conv)),
                 saturated = .sat_logLik,
                 null = .null_logLik)

    # Null model
    null_coefs <- rep(0, ncol(X))
    names(null_coefs) <- colnames(X)
    null_coefs["(Intercept)"] <- .null_intercept
    if (.link == "identity") null_coefs <- c(null_coefs, sigma = sqrt(mean((y - mean(y)) ^ 2)))
    lnl_null <- lnl(null_coefs, gradient = TRUE, hessian = TRUE, information = TRUE, X = X, y = y)
    .null_gradient <- attr(lnl_null, "gradient")
    .null_info <- attr(lnl_null, "info")
    .lm <- drop(crossprod(.null_gradient, solve(.null_info, .null_gradient)))
    .model_info <- attr(.lnl_conv, "info")
#    .w <- drop(crossprod(.coefs[- 1], t(crossprod(.coefs[-1], .model_info[- 1, - 1]))))
    .vcov <- solve(.model_info)
    .w <- drop(crossprod(.coefs[- 1], t(crossprod(.coefs[-1], solve(.vcov[- 1, - 1, drop = FALSE])))))
    .lr <- 2 * unname(.logLik["model"] - .logLik["null"])
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
    structure(result, class = c("binomreg", "micsr"))
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
    if (type == "deviance") resid <- (2 * y - 1) * sqrt(- 2  * object$value)
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

# ajout dans le fichier type_dictionary.R
# ajout du fichier methods_binomreg
# sanity_model.R ajouter dans la liste des modèles supportés

