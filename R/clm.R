#' Constrained least squares
#'
#' Compute the least squares estimator using linear constrains on the
#' coefficients.
#'
#' @name clm
#' @param x a linear model fitted by `lm`,
#' @param R a matrix of constrains (one line for each constrain, one
#'     column for each coefficient),
#' @param q an optional vector of rhs values (by default a vector of
#'     0)
#' @param object a `clm` object for the `summary` and the `vcov`
#'     methods
#' @param \dots further arguments
#' @return an object of class `clm` which inherits from class `lm`
#' @export
clm <- function(x, R, q = NULL){
    .call <- x$call
    if (is.null(q)) q <- rep(0, nrow(R))
    X <- model.matrix(x)
    y <- model.response(model.frame(x))
    beta_nc <- coef(x)
    excess <- R %*% beta_nc - q
    XpXm1 <- solve(crossprod(X))
    beta_c <- beta_nc - drop(XpXm1 %*% t(R) %*% solve(R %*% XpXm1 %*% t(R)) %*% excess)
    names(beta_c) <- names(beta_nc)
    .fitted <- drop(X %*% beta_c)
    .resid <- y - .fitted
    .df.residual <- length(y) - length(beta_c)
    .s2 <- sum(.resid ^ 2) / .df.residual
    structure(list(coefficients = beta_c, model = model.frame(x),
                   residuals = .resid, fitted.values = .fitted,
                   df.residual = .df.residual,
                   call = .call, terms = terms(x), R = R, q = q),
              class = c("clm", "lm"))
}

#' @rdname clm
#' @export
vcov.clm <- function(object, ...) vcov(summary(object))

#' @rdname clm
#' @export
summary.clm <- function(object, ...){
    .sigma <- sqrt(sum(object$residuals ^ 2) / object$df.residual)
    XpXm1 <- solve(crossprod(model.matrix(object$terms, object$model)))
    y <- model.response(model.frame(object))
    N <- length(y)
    Kp1 <- length(coef(object))
    ybar <- mean(y)
    ESS <- sum( (object$fitted.values - ybar) ^ 2)
    RSS <- sum(object$residuals ^ 2)
    TSS <- sum( (y - ybar) ^ 2)
    .r.squared <- ESS / TSS
    .adj.r.squared <- 1 - (1 - .r.squared) * (N - 1) / object$df.residual
    .fstatistic <- c(value = ESS / .sigma ^ 2 / object$df.residual,
                     numdf = N - object$df.residual + 1,
                     dendf = object$df.residual)
    .cov.unscaled <- (XpXm1 - XpXm1 %*% t(object$R) %*%
                      solve(object$R %*% XpXm1 %*%
                            t(object$R)) %*% object$R %*% XpXm1)
    .aliased <- rep(FALSE, Kp1)
    names(.aliased) <- names(coef(object))
    std.err <- .sigma * sqrt(diag(.cov.unscaled))
    b <- coef(object)
    z <- b / std.err
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    .coef <- cbind(b, std.err, z, p)
    colnames(.coef) <- c("Estimate", "Std. Error", 
                         "z-value", "Pr(>|z|)")
    structure(list(call = object$call,
                   terms = object$terms,
                   residuals = object$residuals,
                   coefficients = .coef,
                   aliased = .aliased,
                   sigma = .sigma,
                   df = c(Kp1, object$df.residual, Kp1),
                   r.squared = .r.squared,
                   adj.r.squared = .adj.r.squared,
                   fstatistic = .fstatistic,
                   cov.unscaled = .cov.unscaled),
              class = c("summary.lm", "lm"))
}
