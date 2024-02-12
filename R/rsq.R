#' Coefficient of determination
#'
#' A generic function to compute different flavors of coefficients of determination
#'
#' @name rsq
#' @param x fitted model
#' @param type the type of coefficient of determination
#' @return a numeric scalar
#' @importFrom stats model.response model.frame resid logLik
#' @examples
#' pbt <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = 'probit')
#' rsq(pbt)
#' rsq(pbt, "estrella")
#' rsq(pbt, "veall_zimm")
#' @export
rsq <- function(x, type){
    UseMethod("rsq")
}

#' @rdname rsq
#' @export
rsq.lm <- function(x, type = c("raw", "adj")){
    .type <- match.arg(type)
    if (.type == "raw") summary(x)$r.squared else summary(x)$adj.r.squared
}

#' @rdname rsq
#' @export
rsq.micsr <- function(x, type = c("mcfadden", "cox_snell", "cragg_uhler",
                                     "aldrich_nelson", "veall_zimm", "estrella",
                                     "cor", "ess", "rss", "tjur",
                                  "mckel_zavo", "w", "lm", "lr")){
    type <- match.arg(type)
    y <- model.response(model.frame(x))
    if (! is.numeric(y)) y <- as.numeric(y)
    yb <- mean(y)
    logL <- as.numeric(logLik(x))
    logL0 <- - deviance(x, type = "null") / 2
    df.model <- df.residual(x)
    K <- as.numeric(x$npar) - 1
    N <- nobs(x)
    if (! inherits(x, "ordreg")){
        ESS <- sum((fitted(x) - yb) ^ 2)
        TSS <- sum( (y - yb) ^ 2)
        RSS <- sum(resid(x, type = "response") ^ 2)
    }
    LR <- 2 * (logL - logL0)
    LR_star <- - 2 * logL0
    if (type %in% c("w", "lm", "lr")){
        if (! inherits(x, "micsr")) stop("only implemented for micsr objects")
        if (type == "w") r2 <- x$tests["w"] / (N + x$tests["w"])
        if (type == "lm") r2 <- x$tests["lm"] / N
        if (type == "lr") r2 <- 1 - exp(- x$tests["lr"] / N)
    }
    if (type == "f") r2 <- K * F / (K * F + df.model) 
    if (type == "mcfadden") r2 <- 1 - logL / logL0
    if (type == "cox_snell") r2 <- 1 - exp(- LR / N)
    if (type == "cragg_uhler") r2 <- (1 - exp(- LR / N)) / (1 - exp(- LR_star / N))
    if (type == "aldrich_nelson")  r2 <- LR / (LR + N)
    if (type == "veall_zimm") r2 <- (LR / (LR + N)) / (LR_star / (LR_star + N))
    if (type == "estrella") r2 <- 1 - (logL / logL0) ^ (- 2 / N * logL0)
    if (type == "cor") r2 <- sum( (fitted(x) - yb) * (y - yb)) ^ 2 /
                           sum((fitted(x) - yb) ^ 2) / sum((y - yb) ^ 2)
    if (type == "rss")  r2 <- 1 - RSS / TSS
    if (type == "ess") r2 <- ESS / TSS
    if (type == "tjur"){
        r2_model <- ESS / TSS
        r2_resid <- 1 - RSS / TSS
        r2 <- (r2_model + r2_resid) / 2
        r2
    }
    if (type == "mckel_zavo"){
        lp <- x$linear.predictors
        hlp <- mean(lp)
        ESS <- sum( (lp- hlp) ^ 2)
        .link <- x$call$link
        RSS <- nobs(x) * ifelse(.link == "probit", 1, pi ^ 2 / 3)
        r2 <- ESS / (ESS + RSS)
    }
    r2
}
