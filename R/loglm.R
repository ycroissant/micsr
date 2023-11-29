#' Log-linear model
#'
#' Estimation of log-linear model; the estimation is done by `lm`, but
#' the correct log-likelihood related quantities are returned
#'
#' @name loglm
#' @param formula,data see `lm`
#' @return
#' An object of class `micsr`
#' @author Yves Croissant
#' @examples
#' lm(dist ~ log(speed), cars)
#' loglm(dist ~ log(speed), cars)
#' @export
loglm <- function(formula, data){
    cl <- match.call()
    mf <- model.frame(formula, data)
    mt <- attr(mf, "terms")
    X <- model.matrix(formula, mf)
    K <- ncol(X)
    y <- model.response(mf)
    ofs <- model.offset(mf)
    if (is.null(ofs)) ofs <- 0
    N <- length(y)
    .form <- update(formula, log(.) ~ .)
    lm_reg <- lm(.form, data)
    sig <- sqrt(mean(resid(lm_reg) ^ 2))
    mu <- drop(X %*% coef(lm_reg)) + ofs
    z <- (log(y) - mu) / sig
    .lnl <-  dnorm(z, log = TRUE) - log(sig) - log(y)
    .grad_beta <- z / sig
    .grad_sig <- - 1 / sig + z ^ 2 / sig
    .gradient <- cbind(.grad_beta * X, sigma = .grad_sig)
    .hessian_beta_beta <- - crossprod(X / sig)
    .hessian_beta_sig <- - 1 / sig ^ 2 * apply(z * X, 2, sum)
    .hessian_sig_sig <- N / sig ^ 2 - 3 / sig ^ 2 * sum(z ^ 2)
    .hessian <- rbind(cbind(.hessian_beta_beta, sigma = rep(0, K)),
                      sigma = c(rep(0, K), .hessian_sig_sig))
    .npar <- structure(c(covariates = K, sigma = 1),
                       default = c("covariates", "sigma"))
    structure(list(coefficients = c(coef(lm_reg), sigma = sig),
                   model = mf,
                   gradient = .gradient,
                   hessian = .hessian,
                   linear.predictor = mu,
                   logLik = .lnl,
                   fitted.values = fitted(lm_reg),
                   residuals = resid(lm_reg),
                   est_method = "ml",
                   formula = formula,
                   npar = .npar,
                   value = sum(.lnl),
                   call = cl,
                   terms = mt),
              class = c("micsr", "lm"))
                   
}
