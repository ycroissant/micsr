#' Generalized production function
#'
#' Log-likelihood function for the generalized production function of
#' Zellner and Revankar (1969)
#'
#' @name zellner_revankar
#' @param theta the vector of parameters
#' @param y the vector of response
#' @param Z the matrix of covariates
#' @param sum if `FALSE`, a vector of individual contributions to the
#'     likelihood and the matrix of individual contributions to the
#'     gradient are returned, if `TRUE` a log-likelihood scalar and a
#'     gradient vector are returned
#' @param gradient if `TRUE`, the gradient is returned as an attribute
#' @param hessian if `TRUE`, the hessian is returned as an attrubute
#' @param repar if `TRUE`, the likelihood is parametrized such that
#'     the constant return to scale hypothesis implies that two
#'     coefficients are 0
#' @author Yves Croissant
#' @references \insertRef{ZELL:REVA:69}{micsr}
#' @return a function.
#' @export
zellner_revankar <- function(theta, y, Z, sum = FALSE, gradient = TRUE,
                         hessian = TRUE, repar = TRUE){
    K <- ncol(Z) - 1
    N <- length(y)
    .gamma <- theta[1:(K + 1)]
    .lambda <- theta[K + 2]
    .sigma <- theta[K + 3]
    if (repar){
        Z[, 2:K] <- Z[, 2:K] - Z[, K + 1]
    }
    mu <- drop(Z %*% .gamma)
    if (repar) mu <- mu + Z[, K + 1]
    eps <- log(y) + .lambda * y - mu
    l <- - log(.sigma) + log(1 + .lambda * y) - log(y) +
        dnorm(eps / .sigma, log = TRUE)
    if (sum) l <- sum(l)
    if (gradient){
        g_beta <- eps / .sigma ^ 2 * Z
        g_lam <- y / (1 + .lambda * y) - eps / .sigma ^ 2 * y
        g_sigma <- - 1 / .sigma + eps ^ 2 / .sigma ^ 3
        .gradient <- cbind(g_beta, lambda = g_lam, sigma = g_sigma)
        if (sum) .gradient <- apply(.gradient, 2, sum)
        attr(l, "gradient") <- .gradient
    }
    if (hessian){
        h_bb <- - 1 / .sigma ^ 2 * crossprod(Z)
        h_bg <- 1 / .sigma ^ 2 * apply(Z * y, 2, sum)
        h_bs <- - 2 / .sigma ^ 3 * apply(Z * eps, 2, sum)
        h_gg <- - sum(y ^ 2 / (1 + .lambda * y) ^ 2 +
                      y ^ 2 / .sigma ^ 2)
        h_gs <- 2 / .sigma ^ 3 * sum(eps * y)
        h_ss <- N / .sigma ^ 2 - 3 / .sigma ^ 4 * sum(eps ^ 2)
        H <- rbind(cbind(h_bb, h_bg, h_bs),
                   c(h_bg, h_gg, h_gs),
                   c(h_bs, h_gs, h_ss))
        dimnames(H) <- list(names(theta), names(theta))
        attr(l, "hessian") <- H
    }
    l
}
