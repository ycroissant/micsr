#' Conditional moments test
#'
#' Conditional moments tests for maximum likelihood estimators,
#' particularly convenient for the probit and the tobit model to test
#' relevance of functional form, omitted variables, heteroscedasticity and
#' normality.
#'
#' @name cmtest
#' @param x a fitted model, currently a tobit model either fitted by
#'     `AER::tobit`, `censReg::censReg` or `micsr::tobit1` or a probit model fitted by
#'     `glm` with `family = binomial(link = "probit")` or by `micsr::binomreg` with `link = "probit"`
#' @param test the kind of test to be performed, either a normality
#'     test (or separately a test that the skewness or kurtosis
#'     are 0 and 3), a heteroscedasticity test or a reset
#'     test,
#' @param powers the powers of the fitted values that should be used
#'     in the reset test,
#' @param heter_cov a one side formula that indicates the covariates
#'     that should be used for the heteroscedasticity test (by default
#'     all the covariates used in the regression are used),
#' @param opg a boolean, if `FALSE` (the default), the analytic
#'     derivatives are used, otherwise the outer product of the
#'     gradient formula is used
#' @return an object of class `"htest"` containing the following components:
#' - data.mane: a character string describing the fitted model
#' - statistic: the value of the test statistic
#' - parameter: degrees of freedom
#' - p.value: the p.value of the test
#' - method: a character indicating what type of test is performed
#' @importFrom stats binomial dnorm pnorm model.frame model.response model.matrix coef glm formula family as.formula pchisq
#' @importFrom Rdpack reprompt
#' @author Yves Croissant
#' @keywords htest
#' @references
#'
#' \insertRef{NEWE:85}{micsr}
#'
#' \insertRef{PAGA:VELL:89}{micsr}
#' 
#' \insertRef{TAUC:85}{micsr}
#'
#' \insertRef{WELL:03}{micsr}
#' 
#' @examples
#' charitable$logdon <- with(charitable, log(donation) - log(25))
#' ml <- tobit1(logdon ~ log(donparents) + log(income) + education +
#'              religion + married + south, data = charitable)
#' cmtest(ml, test = "heterosc")
#' cmtest(ml, test = "normality", opg = TRUE)
#' @export
cmtest <- function(x, test = c("normality", "reset", "heterosc",
                               "skewness", "kurtosis"),
                   powers = 2:3, heter_cov = NULL, opg = FALSE){    
    
    UseMethod("cmtest")
}

#' @rdname cmtest
#' @export
cmtest.tobit <- function(x, test = c("normality", "reset", "heterosc",
                                     "skewness", "kurtosis"),
                         powers = 2:3, heter_cov = NULL, opg = FALSE){
    test <- match.arg(test)
    param <- c(coef(x), x$scale)
    X <- model.matrix(x)
    mf <- model.frame(x)
    y <- x$y[, 1]
    result <- cmtest_tobit(param, X, y, mf, test = test,
                           powers = powers, heter_cov = heter_cov, opg = opg)
    .data.name <- paste(deparse(formula(x)))
    if (length(.data.name) > 1) .data.name <- paste(.data.name[1], "...")
    result$data.name <- .data.name
    result
}

#' @rdname cmtest
#' @export
cmtest.micsr <- function(x, test = c("normality", "reset", "heterosc",
                                     "skewness", "kurtosis"),
                         powers = 2:3, heter_cov = NULL, opg = FALSE){
    test <- match.arg(test)
    param <- coef(x)
    X <- model.matrix(x)
    mf <- model.frame(x)
    y <- model.response(mf)
    if (! inherits(x, "tobit1") & ! inherits(x, "binomreg"))
        stop("a tobit1 or a binomreg object is required")
    else{
        if (inherits(x, "tobit1"))
            result <- cmtest_tobit(param, X, y, mf, test = test,
                                   powers = powers, heter_cov = heter_cov, opg = opg)
        else result <- cmtest_probit(param, X, y, mf, test = test,
                                     powers = powers, heter_cov = heter_cov, opg = opg)
    }
    .data.name <- paste(deparse(formula(x)))
    if (length(.data.name) > 1) .data.name <- paste(.data.name[1], "...")
    result$data.name <- .data.name
    result
}


#' @rdname cmtest
#' @export
cmtest.censReg <- function(x, test = c("normality", "reset", "heterosc",
                               "skewness", "kurtosis"),
                         powers = 2:3, heter_cov = NULL, opg = FALSE){
    test <- match.arg(test)
    param <- coef(x)
    mf <- model.frame(x)
    param[length(coef(x))] <- exp(param[length(coef(x))])
    y <- unname(model.response(mf))
    X <- model.matrix(x)
    result <- cmtest_tobit(param, X, y, mf, test = test,
                           powers = powers, heter_cov = heter_cov, opg = opg)
    .data.name <- paste(deparse(formula(x)))
    if (length(.data.name) > 1) .data.name <- paste(.data.name[1], "...")
    result$data.name <- .data.name
    result
}

lnl_cens <- function(param, X, y, sum = TRUE, gradient = FALSE, hessian = FALSE){
    mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
    K <- length(param) - 1
    beta <- param[1:K]
    sig <- param[K + 1]
    yb <- ifelse(y > 0, 1, 0)
    bX <- as.numeric(X %*% beta)
    bxs <- bX / sig
    mls <- mills(- bxs)
    e <- y - bX
    lnl <-  pnorm(- bxs, log.p = TRUE) * (1 - yb) +
        (- log(sig) - 1 / 2 * log(2 * pi) - 1 / 2 * (e / sig) ^ 2) * yb
    if (sum) lnl <- sum(lnl)
    if (gradient){
        grad_beta <-  - sig * mls * (1 - yb) + e * yb
        grad_sigma <- (sig * mls * bX) * (1 - yb) + (- sig ^ 2 + e ^ 2) * yb
        grad <- cbind(grad_beta * X / sig ^ 2, grad_sigma / sig ^ 3)
        if (sum) grad <- apply(grad, 2, sum)
        attr(lnl, "gradient") <- grad
    }
    if (hessian){
        X1 <- X[y > 0, ]
        X0 <- X[y == 0, ]
        H_bb <- (- crossprod(X0 * (mls * (mls - bxs))[y == 0], X0) - crossprod(X1) ) / sig ^ 2
        H_bs <- - 2 / sig ^ 3 * apply(e[y > 0] * X1, 2, sum) +
            1 / sig ^ 2 * apply( (mls * (1 + bxs * (mls - bxs)))[y == 0]  * X0, 2, sum)
        H_ss <- sum(y > 0) / sig ^ 2 - 3  / sig ^ 4 * sum(e[y > 0] ^ 2) -
            2 / sig ^ 3 * sum( (mls * bX)[y == 0]) -
            1 / sig ^ 4 * sum( (mls * (mls - bxs)  * bX ^ 2)[y == 0])
        H <- rbind(cbind(H_bb, H_bs),
                   c(H_bs, H_ss))
        attr(lnl, "hessian") <- H
    }
    lnl
}

grad_cens <- function(param, X, y)
    attr(lnl_cens(param, X, y, gradient = TRUE), "gradient")

hess_cens <- function(param, X, y)
    attr(lnl_cens(param, X, y, hessian = TRUE), "hessian")


cmtest_tobit <- function(param, X, y, mf, test = c("normality", "reset", "heterosc",
                                                      "skewness", "kurtosis"),
                         powers = 2:3, heter_cov = NULL, opg = FALSE){    
    S <- attr(lnl_cens(param, X, y, sum = FALSE, gradient = TRUE, hessian = FALSE), "gradient")
    N <- length(y)
    has.intercept <- any(colnames(X) == "(Intercept)")    
    K <- length(param) - 1 - has.intercept
    bX <- as.numeric(X %*% param[1 : (K + 1)])
    sigma <- unname(param[K + 2])
    z <- bX / sigma
    e <- y - bX
    ones <- rep(1, length(y))
    Ipos <- as.numeric(y > 0)
    mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
    d_mills <- function(x) - mills(x) * (mills(x) + x)
    mls <- mills(- z)
    d_mls <- d_mills(- z)
    if (test %in% c("normality", "skewness", "kurtosis")){
        if (test %in% c("normality", "skewness")){
            psi3 <- Ipos * e ^ 3 - (1 - Ipos) * (z ^ 2 + 2) * sigma ^ 3 * mls
            if (! opg){
                m3_z <- - 3 * Ipos * sigma * e ^ 2 - (1 - Ipos) * 2 * z * sigma ^ 3 * mls +
                    (1 - Ipos) * (z ^ 2 + 2) * sigma ^ 3 * d_mls
                m3_s <- - 3 * Ipos * e ^ 2 * z - 3 * (1 - Ipos) * (z ^ 2 + 2) * sigma ^ 2 * mls
                m3_beta  <- m3_z * (1   / sigma)
                m3_sigma <- m3_z * (- z / sigma) + m3_s
                d_m3 <- cbind(m3_beta * X, m3_sigma)
            }
        }
        if (test %in% c("normality", "kurtosis")){
            psi4 <- Ipos * (e ^ 4 - 3 * sigma ^ 4) + (1 - Ipos) * (z ^ 2 + 3) * sigma ^ 4 * mls * z
            if (! opg){
                m4_z <- - 4 * sigma * Ipos * e ^ 3 + (1 - Ipos) * (3 * z ^ 2 + 3) * sigma ^ 4 * mls -
                    (1 - Ipos) * (z ^ 3 + 3 * z) * sigma ^  4 * d_mls
                m4_s <- Ipos * (- 4 * z * e ^ 3 - 12 * sigma ^ 3) +
                    4 * (1 - Ipos) * sigma ^ 3 * (z ^ 3 + 3 * z) * mls
                m4_beta  <- m4_z * (1   / sigma)
                m4_sigma <- m4_z * (- z / sigma) + m4_s
                d_m4 <- cbind(m4_beta * X, m4_sigma)
            }
        }
        if (test == "normality"){
            M <- cbind(psi3, psi4)
            if (! opg)  W <- cbind(apply(d_m3, 2, sum), apply(d_m4, 2, sum))
            test.name <- "Conditional Expectation Test for Normality"
        }
        if (test == "skewness"){
            M <- matrix(psi3, ncol = 1)
            if (! opg) W <- apply(d_m3, 2, sum)
            test.name <- "Conditional Expectation Test for Skewness"
        }
        if (test == "kurtosis"){
            M <- matrix(psi4, ncol = 1)
            if (! opg) W <- apply(d_m4, 2, sum)
            test.name <- "Conditional Expectation Test for Kurtosis"
        }
    }
    if (test == "reset"){
        gres <- - sigma * (1 - Ipos) * mls + Ipos * e
        form <- as.formula(paste(" ~ 0",
                                 paste("I( bX ^ ", powers, ")", sep = "", collapse = " + "), sep = " + "))
        gres <- Ipos * e - sigma * (1 - Ipos) * mls
        Ys <- model.matrix(form) 
        M <-Ys  *  gres
        if (! opg){
            gres_z <- - Ipos * sigma + sigma * (1 - Ipos) * d_mls
            gres_s <- - Ipos * z + (1 - Ipos) * mls
            gres_beta <- gres_z / sigma
            gres_sigma <- gres_s - gres_z * (z / sigma)
            W <- crossprod(gres_beta * X, Ys)
            W <- rbind(W, sigma = apply(gres_sigma * Ys, 2, sum))
        }
        test.name <- "Reset test"
    }
    if (test == "heterosc"){
        if (is.null(heter_cov)) XH <- X
        else XH <- model.matrix(heter_cov, mf)
        if (colnames(XH)[1] == "(Intercept)") XH <- XH[, -1, drop = FALSE]
        rho <- Ipos * ( e ^ 2 - sigma ^ 2) + (1 - Ipos) * sigma ^ 2 * mls * z
        M <- rho * XH
        if (! opg){
            rho_z <- Ipos * (- 2 * e * sigma) + (1 - Ipos) * sigma ^ 2 * (- z * d_mls + mls)
            rho_s <- Ipos * (- 2 * e * z - 2 * sigma) + (1 - Ipos) * (2 * sigma * mls * z)
            rho_beta <- rho_z / sigma
            rho_sigma <- rho_z * (- z / sigma) + rho_s
            W <- cbind(crossprod(X, rho_beta * XH))
            W <- rbind(W, sigma = apply(rho_sigma * XH, 2, sum))
        }
        test.name <- "Heteroscedasticity Test"
    }
    if (opg){
        W <- t(S) %*% M
        F <- crossprod(S)
    }
    else{
        F <- - attr(lnl_cens(param, X, y, sum = TRUE,
                             gradient = FALSE, hessian = TRUE), "hessian")
        W <- - W
    }
    Q <- crossprod(M - S %*% solve(F, W))
    m <- apply(M, 2, sum)
    stat <- as.numeric(m %*% solve(Q, m))
    df <- ncol(M)
    if (df > 1){
        stat <- c("chisq" = stat)
        param <- c("df" = df)
        pval <- pchisq(stat, df, lower.tail = FALSE)
    }
    else{
        stat <- c("z" = sqrt(stat))
        param <- NULL
        pval <- pnorm(stat, lower.tail = FALSE) * 2
    }
    structure(list(statistic = stat,
                   parameter = param,
                   p.value = pval,
                   method = test.name),
              class = "htest")
}

#' @rdname cmtest
#' @export
cmtest.glm <- function(x, test = c("normality", "reset", "heterosc",
                                   "skewness", "kurtosis"),
                       powers = 2:3, heter_cov = NULL, opg = FALSE){
    .family <- family(x)
    if (! (.family$family == "binomial" & .family$link == "probit"))
        stop("cmtest only implemented for probit models")
    test <- match.arg(test)
    param <- coef(x)
    mf <- model.frame(x)
    y <- unname(model.response(mf))
    X <- model.matrix(x)
    result <- cmtest_probit(param, X, y, mf, test = test,
                            powers = powers, heter_cov = heter_cov, opg = opg)
    .data.name <- paste(deparse(formula(x)))
    if (length(.data.name) > 1) .data.name <- paste(.data.name[1], "...")
    result$data.name <- .data.name
    result
}
    
cmtest_probit <- function(param, X, y, mf, test = c("normality", "reset", "heterosc",
                                                      "skewness", "kurtosis"),
                          powers = 2:3, heter_cov = NULL, opg = FALSE){
    test <- match.arg(test)
    bX <- z <- drop(X %*% param)
    N <- length(y)
    has.intercept <- any(colnames(X) == "(Intercept)")    
    K <- length(param) - 1 - has.intercept
    gres <- dnorm(z) / pnorm(z) / (1 - pnorm(z)) * (y - pnorm(z))
    gres_z <- - (z + dnorm(z) * (1 - 2 * pnorm(z)) / (pnorm(z) * (1 - pnorm(z)))) * gres -
        dnorm(z) ^ 2 / (pnorm(z) * (1 - pnorm(z)))
    S <- gres * X
    if (test %in% c("normality", "skewness", "kurtosis")){
        if (test %in% c("normality", "skewness")){
            psi3 <- z ^ 2 * gres
            if (! opg){
                psi3_z <- 2 * z * gres + z ^ 2 * gres_z
                d_m3 <- psi3_z * X
            }
        }
        if (test %in% c("normality", "kurtosis")){
            psi4 <- z ^ 3 * gres
            if (! opg){
                psi4_z <- 3 * z ^ 2 * gres + z ^ 3 * gres_z
                d_m4 <- psi4_z * X
            }
        }
        if (test == "normality"){
            M <- cbind(psi3, psi4)
            if (! opg)  W <- cbind(apply(d_m3, 2, sum), apply(d_m4, 2, sum))
            test.name <- "Conditional Expectation Test for Normality"
        }
        if (test == "skewness"){
            M <- matrix(psi3, ncol = 1)
            if (! opg) W <- apply(d_m3, 2, sum)
            test.name <- "Conditional Expectation Test for Skewness"
        }
        if (test == "kurtosis"){
            M <- matrix(psi4, ncol = 1)
            if (! opg) W <- apply(d_m4, 2, sum)
            test.name <- "Conditional Expectation Test for Kurtosis"
        }
    }
    if (test == "reset"){
        form <- as.formula(paste(" ~ 0",
                                 paste("I( bX ^ ", powers, ")", sep = "", collapse = " + "), sep = " + "))
        Ys <- model.matrix(form) 
        M <- Ys  *  gres
        if (! opg){
            W <- crossprod(gres_z * X, Ys)
        }
        test.name <- "Reset test"
    }
    if (test == "heterosc"){
        if (is.null(heter_cov)) XH <- X
        else XH <- model.matrix(heter_cov, mf)
        if (colnames(XH)[1] == "(Intercept)") XH <- XH[, -1, drop = FALSE]
        rho <- z * gres
        M <- rho * XH
        if (! opg){
            rho_z <- gres + z * gres_z
            W <- crossprod(X, rho_z * XH)
        }
        test.name <- "Heteroscedasticity Test"
    }
    if (opg){
        W <- t(S) %*% M
        F <- crossprod(S)
    }
    else{
        F <- - crossprod(gres_z * X, X)
        W <- - W
    }
    Q <- crossprod(M - S %*% solve(F, W))
    m <- apply(M, 2, sum)
    stat <- as.numeric(m %*% solve(Q, m))
    df <- ncol(M)
    if (df > 1){
        stat <- c("chisq" = stat)
        param <- c("df" = df)
        pval <- pchisq(stat, df, lower.tail = FALSE)
    }
    else{
        stat <- c("z" = sqrt(stat))
        param <- NULL
        pval <- pnorm(stat, lower.tail = FALSE) * 2
    }
    structure(list(statistic = stat,
                   parameter = param,
                   p.value = pval,
                   method = test.name),
              class = "htest")
}
