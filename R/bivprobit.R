#' Bivariate probit
#'
#' Estimation of bivariate probit models by maximum likelihood
#'
#' @name bivprobit
#' @param formula a symbolic description of the model, (for the count
#'     component and for the selection equation).
#' @param data a data frame,
#' @param subset,weights,na.action,offset see `stats::lm`,
#' @param type for the `logLik` method
#' @param ... further arguments
#' @param object a `bivprobit` object
#' @param type for the `logLik` method
#' @return an object of class `micsr`, see `micsr::micsr` for further details
#' @examples
#' pins <- bivprobit(doctor | privateins ~ privateins + size + smsa + age + sex + educ + log(wage) |
#'                   . - privateins + pluriloc + nbemp, private_ins)
#' @export
bivprobit <- function (formula, data, weights, subset, na.action, offset, ...) {
    .call <- match.call()
    cl <- match.call(expand.dots = FALSE)
    .formula <- cl$formula <- Formula(formula)
    m <- match(c("formula", "data", "subset", "weights"), names(cl), 
        0L)
    cl <- cl[c(1L, m)]
    mf <- cl
    mf[[1L]] <- as.name("model.frame")
    mf$dot <- "previous"
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    old_options <- options(warn = -1)
    X1 <- model.matrix(.formula, mf, rhs = 1)
    X2 <- model.matrix(.formula, mf, rhs = 2)
    options(old_options)
    y1 <- model.part(.formula, mf, lhs = 1, drop = TRUE)
    y2 <- model.part(.formula, mf, lhs = 2, drop = TRUE)
    q1 <- 2 * y1 - 1
    q2 <- 2 * y2 - 1
    K1 <- ncol(X1)
    K2 <- ncol(X2)
    .npar <- c(eq1 = K1, eq2 = K2, vcov = 1)
    attr(.npar, "default") <- c("eq1", "eq2", "vcov")
    N <- length(y1)
    pbt1 <- glm(y1 ~ X1 - 1, family = binomial(link = "probit"))
    pbt2 <- glm(y2 ~ X2 - 1, family = binomial(link = "probit"))

    lnl <- function(param, gradient = FALSE, hessian = FALSE, 
                    sum = TRUE, robust = TRUE, X1, X2, q1, q2) {
        beta1 <- param[1:ncol(X1)]
        beta2 <- param[(ncol(X1) + 1):(ncol(X1) + ncol(X2))]
        v <- param[ncol(X1) + ncol(X2) + 1]
        if (robust) {
            fun_rho <- function(x) atan(x) * 2/pi
            grad_rho <- function(x) 2/pi/(1 + x^2)
            hess_rho <- function(x) -4/pi * x/(1 + x^2)^2
        }
        else {
            fun_rho <- function(x) x
            grad_rho <- function(x) x
            hess_rho <- function(x) x
        }
        rho <- fun_rho(v)
        w1 <- q1 * drop(X1 %*% beta1)
        w2 <- q2 * drop(X2 %*% beta2)
        r <- rho * q1 * q2
        PHI2 <- pbnorm(w1, w2, r)
        phi2 <- exp(-0.5 * (w1 ^ 2 + w2 ^ 2 - 2 * r * w1 * w2)/(1 - 
            r ^ 2))/(2 * pi * sqrt(1 - r ^ 2))
        result <- log(PHI2)
        if (sum) result <- sum(result)
        if (gradient) {
            g1 <- dnorm(w1) * pnorm((w2 - r * w1)/sqrt(1 - r ^ 2))
            g2 <- dnorm(w2) * pnorm((w1 - r * w2)/sqrt(1 - r ^ 2))
            grad1 <- q1 * g1 / PHI2
            grad2 <- q2 * g2 / PHI2
            gradr <- q1 * q2 * phi2 / PHI2
            gradv <- gradr * grad_rho(v)
            .gradient <- cbind(grad1 * X1, grad2 * X2, gradv)
            if (sum) .gradient <- apply(.gradient, 2, sum)
            attr(result, "gradient") <- .gradient
        }
        if (hessian) {
            delta <- 1/sqrt(1 - r^2)
            g1 <- dnorm(w1) * pnorm((w2 - r * w1)/sqrt(1 - r ^ 2))
            g2 <- dnorm(w2) * pnorm((w1 - r * w2)/sqrt(1 - r ^ 2))
            v1 <- delta * (w2 - r * w1)
            v2 <- delta * (w1 - r * w2)
            h_11 <- -w1 * g1/PHI2 - r * phi2/PHI2 - g1 ^ 2 / PHI2 ^ 2
            h_22 <- -w2 * g2/PHI2 - r * phi2/PHI2 - g2 ^ 2 / PHI2 ^ 2
            h_12 <- q1 * q2 * (phi2/PHI2 - g1 * g2/PHI2^2)
            h_1r <- q2 * phi2/PHI2 * (r * delta * v1 - w1 - g1/PHI2) * 
                hess_rho(v) * (-1)
            h_2r <- q1 * phi2/PHI2 * (r * delta * v2 - w2 - g2/PHI2) * 
                hess_rho(v) * (-1)
            gradr <- q1 * q2 * phi2/PHI2
            gradv <- gradr * grad_rho(v)
            h_rr <- phi2/PHI2 * (delta^2 * r * (1 - delta^2 * 
                (w1^2 + w2^2 - 2 * r * w1 * w2)) + delta^2 * 
                w1 * w2 - phi2/PHI2) * grad_rho(v)^2 + gradr * 
                hess_rho(v)
            H1 <- cbind(crossprod(h_11 * X1, X1), crossprod(h_12 * 
                X1, X2), apply(h_1r * X1, 2, sum))
            H2 <- cbind(crossprod(h_12 * X2, X1), crossprod(h_22 * 
                X2, X2), apply(h_2r * X2, 2, sum))
            Hr <- c(apply(h_1r * X1, 2, sum), apply(h_2r * X2, 
                2, sum), sum(h_rr))
            H <- rbind(H1, H2, Hr)
            attr(result, "hessian") <- H
        }
        result
    }
    start_1 <- coef(pbt1)
    start_2 <- coef(pbt2)
    names(start_1) <- paste("eq1", colnames(X1), sep = "_")
    names(start_2) <- paste("eq2", colnames(X2), sep = "_")
    .start <- c(start_1, start_2, rho = 0)
    # use of optim, newton not reliable
    f <- function(x) - lnl(x, sum = TRUE, robust = TRUE, gradient = FALSE, hessian = FALSE, X1 = X1, X2 = X2, q1 = q1, q2 = q2)
    g <- function(x) - attr(lnl(x, sum = TRUE, robust = TRUE, gradient = TRUE, hessian = FALSE, X1 = X1, X2 = X2, q1 = q1, q2 = q2), "gradient")
    w <- optim(.start, f, g, method = "BFGS")
    .coefs <- w$par
    .coefs_null <- newton(lnl, c(qnorm(mean((q1 + 1) / 2)), qnorm(mean((q2 + 1) / 2)), 0),
                          trace = FALSE, direction = "max",
                          X1 = matrix(rep(1, N), ncol = 1), X2 = matrix(rep(1, N), ncol = 1), q1 = q1, q2 = q2)
    lnl_null <- lnl(.coefs_null, gradient = FALSE, sum = TRUE, 
        robust = TRUE, X1 = matrix(rep(1, N), ncol = 1), X2 = matrix(rep(1, 
            N), ncol = 1), q1 = q1, q2 = q2)
print(.coefs_null)
    .coefs[K1 + K2 + 1] <- atan(.coefs[K1 + K2 + 1]) * 2 / pi
    .lnl_conv <- lnl(.coefs, sum = FALSE, gradient = TRUE, hessian = TRUE, 
        robust = FALSE, X1 = X1, X2 = X2, q1 = q1, q2 = q2)
    result <- list(coefficients = .coefs,
                   model = mf,
                   value = sum(as.numeric(.lnl_conv)), 
                   logLik = as.numeric(.lnl_conv),
                   logLik_null = lnl_null, 
                   gradient = attr(.lnl_conv, "gradient"),
                   hessian = attr(.lnl_conv, "hessian"),
                   est_method = "ml",
                   npar = .npar,
                   formula = formula, 
                   call = .call)
    structure(result, class = c("bivprobit", "micsr"))
}


#' @rdname bivprobit
#' @export
logLik.bivprobit <- function(object, ..., type = c("model", "null")){
    .type <- match.arg(type)
    if (.type == "model")
        result <- structure(object$value, nobs = nobs(object), df = sum(object$npar), class = "logLik")
    if (.type == "null"){
        result <- structure(object$logLik_null, nobs = nobs(object), df = 3, class = "logLik")
    }
    result
}
