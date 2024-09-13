lnliv_ldv <- function(param, X1, X2, W, y, sum = TRUE, gradient = FALSE, right, model){
    L <- ncol(X1) + ncol(W)
    G <- ncol(W)
    K <- ncol(X1) + ncol(X2)
    N <- length(y)
    X <- cbind(X1, X2)
    Z <- cbind(X1, W)
    beta <- param[1:L]
    rho <- param[L + (1:G)]
    pi2 <- param[L + G + 1:(G * K)]
    PI2 <- matrix(pi2, nrow = G, byrow = TRUE)
    chol <-  param[L + G + G * K + 1:(G * (G + 1) / 2)]
    Chol <- matrix(0, G, G)
    Chol[lower.tri(Chol, diag = TRUE)] <- chol
    if (model == "tobit") sigma <- param[L + G + G * K + G * (G + 1) / 2 + 1]
    if (model == "probit") sigma <- 1
    d <- sqrt(sigma ^ 2 - sum(rho ^ 2))
    Wres <- (W %*% Chol - X %*% t(PI2))
    lp <- as.numeric(Z %*% beta + Wres %*% rho)
    q <- 2 * y - 1
    pos <- as.numeric(y > 0)
    # get the position of the diagonal elements of the Cholesky decomposition matrix
    if (G > 1) idx <- cumsum(c(1, G:2)) else idx <- 1
    # the log of the marginal density
    lnl_g <- - 1 / 2 * G * log(2 * pi) + sum(log(chol[idx])) - 1 / 2 * apply(Wres ^ 2, 1, sum)
    # the log of the conditional density
    if (model == "probit") lnl_f <- pnorm(q * lp / d, log.p = TRUE)
    if (model == "tobit") lnl_f <- (1 - pos) * pnorm(- lp / d, log.p = TRUE) +
                              pos * (- log(d) + dnorm( (y - lp) / d, log = TRUE))
    lnl <- lnl_f + lnl_g
    if (sum) lnl <- sum(lnl)
    if (gradient){
        # get the vector of the columns of Wres (1, 2 2, 3 3 3, etc.)
        V <- matrix(1:G, nrow = G, ncol = G) ; idx <- V[lower.tri(t(V), diag = TRUE)]
        lnl_g_pi <- Reduce("cbind", lapply(1:G, function(i) Wres[, i] * X))
        lnl_g_c2 <- diag(G)[lower.tri(diag(G), diag = TRUE)] / chol
        lnl_g_c2 <- matrix(lnl_g_c2, byrow = TRUE, ncol = length(lnl_g_c2), nrow = N)
        lnl_g_c3 <- - Wres[, rep(1:G, times = G:1)] * W[, idx]
        lnl_g_c <- lnl_g_c2 + lnl_g_c3
        if (model == "probit"){
            mu <- mills(q * lp / d)
            lnl_f_beta <- mu * q / d * Z
            lnl_f_rho <- mu * q / d * (Wres + t(rho) %x% (lp / d ^ 2))
            lnl_f_pi <- - mu * q / d  * t(rho) %x% X
            lnl_f_c <- mu * q / d * t(t(W[, idx]) * rep(rho, times = G:1))
        }
        if (model == "tobit"){
            mu <- mills(- lp / d)
            lnl_f_beta <- (- (1 - pos) * mu / d +
                       pos * (y - lp) / d ^ 2) * Z
            lnl_f_rho <- - (1 - pos) * mu / d * (Wres + t(rho) %x% (lp / d ^ 2)) +
                pos / d ^ 2 * (Wres * (y - lp) + t(rho) %x% (1 - (y - lp) ^ 2 / d ^ 2))
            lnl_f_pi <-  (1 - pos) * mu / d  * t(rho) %x% X +
                - pos / d ^ 2 * t(rho) %x% ((y - lp) * X)
            lnl_f_c1 <- - (1 - pos) * mu / d * t(t(W[, idx]) * rep(rho, times = G:1))
            lnl_f_c2 <- pos / d ^ 2 * (y - lp) * t(t(W[, idx]) * rep(rho, times = G:1))
            lnl_f_c <- lnl_f_c1 + lnl_f_c2
            lnl_f_s <- (1 - pos) * mu * lp / d ^ 3 * sigma -
                pos * sigma / d ^ 2 * (1 - (y - lp) ^ 2 / d ^ 2)
        }
        lnl_beta <- lnl_f_beta
        lnl_rho <- lnl_f_rho
        lnl_pi <- lnl_f_pi + lnl_g_pi
        lnl_c <- lnl_f_c + lnl_g_c
        if (model == "tobit") lnl_s <- lnl_f_s
        gradObs <- cbind(lnl_beta, lnl_rho, lnl_pi, lnl_c)
        colnames(gradObs) <- names(param)
        if (model == "tobit") gradObs <- cbind(gradObs, lnl_s)
        if (sum) gradObs <- apply(gradObs, 2, sum)
        attr(lnl, "gradient") <- gradObs
    }
    lnl
}

