#' Instrumental variable estimation for exponential conditional mean
#' models
#'
#' Exponential conditional mean models are particularly useful for
#' non-negative responses (including count data). Least squares and one or two steps
#' IV estimators are available
#' 
#' @name expreg
#' @param formula a two-part right hand side formula, the first part
#'     describing the covariates and the second part the instruments
#' @param data a data frame,
#' @param subset,weights,na.action,offset see `stats::lm`
#' @param method one of `"gmm"` (the default), `"iv"` or `ls`.
#' @param error one of `"mult"` (the default) or `"add"` in order to
#'     get a model with respectively a multiplicative or an additive
#'     error
#' @param ... further arguments
#' @return An object of class `"micsr"` which is a list
#'     containing the following components:
#' - coefficients: a named vector of coefficients,
#' - residuals: the vector of residuals,
#' - fitted.values: the fitted values
#' - vcov: estimation of the covariance matrix of the estimators
#' - value: value of the objective function at convergence
#' - model: the model frame
#' - call: the matched call
#' - K: the number of covariates
#' - L: the number of instruments
#' - df.residual: the degrees of freedom of the regression
#' - xlevels: a record of the levels of the factors used in fitting
#' - na.action: information returned by `model.frame` on the sepcial handling of `NA`'s
#' @importFrom stats .getXlevels coef glm model.matrix model.response
#'     nobs optim pchisq pnorm poisson printCoefmat vcov model.frame
#' @importFrom Formula Formula
#' @keywords models
#' @references \insertRef{MULL:97}{micsr}
#' @author Yves Croissant
#' @examples
#' cigmales <- dplyr::mutate(cigmales,
#'                           age2 = age ^ 2, educ2 = educ ^ 2, educage = educ * age,
#'                           age3 = age ^ 3, educ3 = educ ^ 3)
#' expreg(cigarettes ~ habit + price + restaurant + income + age + age2 + educ + educ2 +
#'                      famsize + race | . - habit + reslgth + lagprice + age3 + educ3 + educage,
#'                      data = cigmales)
#' expreg(birthwt ~ cigarettes + parity + race + sex | parity + race + sex +
#'                   edmother + edfather + faminc + cigtax, data = birthwt)
#' @export
expreg <- function(formula,
                   data,
                   subset,
                   weights,
                   na.action,
                   offset,
                   method = c("iv", "gmm", "ls"),
                   error = c("mult", "add"),
                   ...){
    .est_method <- match.arg(method)
    .call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    .formula <- mf$formula <- Formula(formula)
    n_rhs <- length(.formula)[2]
    if (.est_method %in% c("iv", "gmm") & n_rhs == 1)
        stop("iv or gmm estimation requires a two-part formula")
    if (.est_method == "ls" & n_rhs == 2)
        warning("the second rhs of the Formula is irrelevant for ls estimation")
    .error <- match.arg(error)
    m <- match(c("formula", "data", "subset", "weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf$dot <- "previous"
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    # Extract the element of the model
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    offset <- model.offset(mf)
    X <- model.matrix(.formula, mf, rhs = 1)
    K <- ncol(X)
    names_X <- colnames(X)
    if (n_rhs > 1) Z <- model.matrix(.formula, mf, rhs = 2)
    else Z <- X
    y <- model.response(mf)
    N <- length(y)
    start <- coef(glm(y ~ X - 1, family = poisson))
#    old_options <- options( warn = -1 )

    objfun <- function(param, error, W){
        mu <- exp(drop(X %*% param))
        N <- length(y)
        if (error == "add") e <- y - mu
        else e <- y / mu - 1
        em <- drop(crossprod(Z, e)) / N
        drop(crossprod(crossprod(W, em), em))
    }

    objgr <- function(param, error, W){
        mu <- exp(drop(X %*% param))
        N <- length(y)
        if (error == "add"){e <- y - mu; g <- mu}
        else{e <- y / mu - 1; g <- (e + 1)}
        em <- drop(crossprod(Z, e)) / N
        - 2 * drop(crossprod(g * X, Z) %*% W %*% em) / N
    }

    W <- solve(crossprod(Z)) * N
    results <- optim(par = start,
                     fn = objfun,
                     gr = objgr,
                     method = "BFGS",
                     error = .error,
                     W = W)
    .coef <- results$par
    .value <- results$value

    mu <- exp(drop(X %*% .coef))
    if (.error == "add") e <- y - mu else e <- y / mu - 1

    if (.est_method == "ls") .value <- sum(e ^ 2)
    
    if (.est_method == "gmm"){
        EM <- Z * e
        W <- solve(crossprod(EM)) * N
        results <- optim(par = .coef,
                         fn = objfun,
                         gr = objgr,
                         method = "BFGS",
                         error = .error,
                         W = W)
        .coef <- results$par
        .value = results$value
        .vcov <- solve(crossprod(X, Z) %*% W %*% crossprod(Z, X)) * N
    }
    else{
        sig2 <- mean(e ^ 2)
        .vcov <- sig2 * solve(crossprod(X, Z) %*% W %*% crossprod(Z, X)) * N
    }
    names(.coef) <- rownames(.vcov) <- colnames(.vcov) <- names_X
    mu <- exp(drop(X %*% .coef))    
    result <- list(coefficients = .coef,
                   residuals = y - mu,
                   fitted.values = mu,
                   vcov = .vcov,
                   value = .value,
                   model = mf,
                   call = .call,
                   K = ncol(X),
                   formula = .formula,
                   terms = mt,
                   npar = c(covariates = K),
                   df.residual = N - ncol(X),
                   xlevels = .getXlevels(mt, mf),
                   na.action = attr(mf, "na.action"),
                   est_method = .est_method
                   )
    if (.est_method %in% c("iv", "gmm"))  result$L = ncol(Z)
#    options(old_options)
    structure(result, class = c("micsr"))
}

