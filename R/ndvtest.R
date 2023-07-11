#' Shi test for non-nested models
#' 
#' The Shi test correct the bias of the Vuong test
#' 
#' 
#' @aliases ndvtest
#' @param x a first fitted model,
#' @param y a second fitted model,
#' @param size the size of the test,
#' @param pval should the p-value be computed ?
#' @param nested a boolean, `TRUE` for nested models,
#' @param vartest a boolean, if `TRUE`, the variance test is computed,
#' @param ndraws the number of draws for the simulations,
#' @param diffnorm a creuser,
#' @param seed the seed,
#' @param numbers a user provided matrix of random numbers
#' @param nd a boolean, if `TRUE` (the default) the non-degenarate Vuong test is computed,
#' @param print.level the level of details to be printed,
#' @param \dots further arguments,
#' @param object an object of class `maxLik2`.
#' @return an object of class \code{"htest"}
#' @importFrom Rdpack reprompt
#' @importFrom stats pnorm
#' @importFrom CompQuadForm davies
#' @seealso the classical Vuong test is implemented in `pscl::vuong` and `nonnest2::vuongtest`.
#' @references
#'
#' \insertRef{VUON:89}{micsr}
#'
#' \insertRef{SHI:15}{micsr}
#'
#' @importFrom stats ecdf logLik optimize qnorm quantile rnorm uniroot
#' @keywords htest
#' @examples
#' 
#' # A poisson model example from the nonnest2 man page
#' data("housing", package = "MASS")
#' house1 <- glm(Freq ~ Infl + Type + Cont, family = poisson, data = housing)
#' house2 <- glm(Freq ~ Infl + Sat,         family = poisson, data = housing)
#' nonnest2::vuongtest(house1, house2)
#' ndvtest(house1, house2)
#' data("bioChemists", package = "pscl")
#' bio1 <- glm(art ~ fem + mar + phd + ment, family=poisson, data=bioChemists)
#' bio2 <- pscl::hurdle(art ~ fem + mar + phd + ment, data=bioChemists)
#' bio3 <- pscl::zeroinfl(art ~ fem + mar + phd + ment, data=bioChemists)
#' nonnest2::vuongtest(bio3, bio2)
#' ndvtest(bio3, bio2)
#' @export
ndvtest <- function(x, y, size = 0.05, pval = TRUE,
                    nested = FALSE, vartest = FALSE,
                    ndraws = 1E04, diffnorm = 0.1, seed = 1,
                    numbers = NULL, nd = TRUE,
                    print.level = 0){
    
    # if pval is TRUE, size is adjusted, otherwise it is fixed
    data.name <- c(
        paste(deparse(substitute(x))),
        paste(deparse(substitute(y)))
    )
    data.name <- paste(data.name, collapse = "-")
    set.seed(seed)

    # the two models are renamed f and g like in Vuong paper, the
    # generics llcont, estfun and bread are used to extract the
    # relevant features of the log-likelihood
    f <- x                    ;  g  <- y    
    N <- length(llcont(f))
    K <- ncol(estfun(f)) + ncol(estfun(g))
    
    # Compute the LR statistic as an average and its variance
    LR <- as.numeric(logLik(f) - logLik(g)) / N
    w2 <- sum((llcont(f)- llcont(g)) ^ 2) / N - LR ^ 2#(LR / N) ^ 2
    
    # solveA is a block-diagonal matrix containing (H_f / N) ^ -1 and
    # - (H_g / N) ^ - 1. The bread function returns - (H / N) ^ - 1 =
    # - N x H ^ -1 = - N x vcov
    solveA <- bdiag(- bread(f), bread(g))
    # B binds the column of G_f and - G_g
    B <- cbind(estfun(f), - estfun(g))
    # substract the mean of the gradient (usefull if the gradient is
    # not null at the optimum)
    B <- t(t(B) - apply(B, 2, mean))
    # B is the estimation of the variance of the gradient
    B <- crossprod(B) / N
    eigen_B <- eigen(B)
    sqrtB <- eigen_B$vectors %*% diag(sqrt(abs(eigen_B$values))) %*% t(eigen_B$vectors)
    V <- sqrtB %*% solveA %*% sqrtB
    V <- eigen(V, symmetric = TRUE)$values
    # Compute the 4 matrices that are used to compute empirically the
    # distribution of the Vuong statistic
    if (! nested){
        # rho is a K vector containing 0 except at the position of the
        # higher eigen value of V (in absolute value), where a 1 is
        # inserted
        rho <- as.numeric((abs(V) - max(abs(V))) == 0)
        rho <- rho / sqrt(sum(rho))
    }
    else{
        # rho is irrelevant for nested models (a vector of K 0)
        rho <- rep(0, length(V))
    }
    # Sigma is the covariance matrix of the K + 1 normal random normal
    # draws
    Sigma <- rbind(c(1, rho),
                   cbind(rho, diag(K)))
    # To get the square root of Sigma, just insert a line of 0 at the
    # position of the highest eigen value
    Sigma[c(FALSE, as.logical(rho)), ] <- 0

    # if numbers is not null, it is a user provided matrix of random
    # draws. Otherwise rnorm is used.
    if (is.null(numbers)) Z <- matrix(rnorm(ndraws * (K + 1)), ndraws)
    else Z <- numbers

    # then transform the random numbers according to the relevant
    # covariance matrix
    Z <- Z %*% Sigma
    ZL <- Z[,   1]
    ZP <- Z[, - 1]

    # the scaled sum of the eigen values
    trVsq <- sum(V ^ 2)
    Vnmlzd <- V / sqrt(trVsq)

    # the 4 scalars used to compute random draws from the same
    # distribution as the one of the modified Vuong statistic
    A1 <- ZL
    A2 <- as.numeric(apply(ZP, 1, function(x) sum(Vnmlzd * x ^ 2) / 2) -
                     sum(Vnmlzd) / 2)
    A3 <- as.numeric(apply(ZP, 1, function(x) sum(Vnmlzd * rho * x)))
    A4 <- as.numeric(apply(ZP, 1, function(x) sum(x ^ 2 * Vnmlzd ^ 2)))

    Tmod <- function(sigma, cst){
        # R draws in the distribution of the Vuong statistic
        num <- sigma * A1 - A2
        denom <- sigma ^ 2 - 2 * sigma * A3 + A4 + cst
        num / sqrt(denom)
    }

    quant <- function(sigma, cst, size)
        # compute the empirical quantile of level 1 - size in the
        # distribution of the Vuong statistic
        as.numeric(quantile(abs(Tmod(sigma, cst)), 1 - size))
        
    sigstar <- function(cst, size)
        # compute the value of sigma which maximize the empirical
        # quantile of level 1 - size of the distribution of the Vuong
        # statistic
        optimize(function(x) quant(x, cst, size), c(0, 5), maximum = TRUE)$maximum
    
    seekpval <- function(size){
        # for a given size, compute the constant so that the critical
        # value equals the target
        c.value <- 0
        cv.value <- quant(sigstar(0, size), 0, size)
        # what is the interest of that equation
        cv.normal <- max(qnorm(1 - size / 2), quantile(abs(ZL), 1 - size))
        cv.normal <- qnorm(size / 2, lower.tail = FALSE)
        cv.target <- cv.normal + diffnorm
        if (cv.value < cv.target){
            cv.value <- max(cv.value, cv.normal)
        }
        else {
            froot <- function(cst) quant(sigstar(cst, size), cst, size) - cv.target
            zo <- uniroot(froot, c(0, 10))
            c.value <- zo$root
            cv.value <- quant(sigstar(c.value, size), c.value, size)
        }
        LRmod <- LR + sum(V) / (2 * N)
        w2mod <- w2 + c.value * sum(V ^ 2) / N
        Tnd <-  sqrt(N) * LRmod / sqrt(w2mod)
        pvalue <- 1 - ecdf(abs(Tmod(sigstar(c.value, size), c.value)))(abs(Tnd))
        list(pvalue = pvalue, stat = Tnd,
             constant = c.value, cv = cv.value)
    }

    if (vartest){
        # compute the variance test and quit"
        W <- eigen(crossprod(B, - solveA), only.values = TRUE)$values
        stat <- N * w2
        pvalue <- CompQuadForm::davies(stat, W ^ 2)$Qq
        results <- list(statistic = c(w2 = w2),
                        p.value = pvalue,
                        data.name = data.name,
                        alternative = "positive variance",
                        method = "variance test")
    }
    else{
        if (nested){
            # for the nested model, the statistic is the classical
            # log-likelihood ratio and the distribution is a weighted
            # chi^2
            Tvuong <- 2 * LR * N
            W <- eigen(crossprod(B, - solveA), only.values = TRUE)$values
            pval_vuong <- CompQuadForm::davies(Tvuong, W)$Qq
            if (nd){
                # for the non-degenarate version, just modify the
                # numerator, and compute the one sided p-value
                LRmod <- LR + sum(V) / (2 * N)
                Tnd <- sqrt(N) * (LRmod) / sqrt(w2)
            if (Tnd < 0) pvalue <- ecdf(Tmod(0, 0))(Tnd)
                else pvalue <- 1 - ecdf(Tmod(0, 0))(Tnd)
                results <- list(statistic = c(z = Tnd),
                                method = "Non-degenerate Vuong test for nested models",
                                p.value = pvalue,
                                data.name = data.name,
                                alternative = "different models",
                                parameters = c(
                                    vuong_stat = Tvuong,
                                    vuong_p.value = pval_vuong))
            }
            else{
                results <- list(statistic = c(wchisq = Tvuong),
                                method = "Vuong test for nested models",
                                p.value = pval_vuong,
                                data.name = data.name,
                                alternative = "different models")
            }
        }
        else{
            Tvuong <-sqrt(N) * LR / sqrt(w2)
            pval_vuong <- pnorm(abs(Tvuong), lower.tail = FALSE) * 2
            if (nd){
                if (! pval){
                    results <- seekpval(size)
                    results <- list(statistic = c(z = as.numeric(results$stat)),
                                    method = "Non-degenerate Vuong test for non-nested models",
                                    data.name = data.name,
                                    alternative = "different models",
                                    parameters = c(
                                        size = size,
                                        vuong_stat = Tvuong,
                                        constant = results$constant,
                                        'crit-value' = results$cv,
                                        'sum e.v.' = sum(V),                                        
                                        'vuong_p.value' = pval_vuong))
                }
                else{
                    froot <- function(alpha) seekpval(alpha)$pvalue - alpha
                    if (froot(1E-100) < 0) results <- list(pvalue = 0, stat = Inf, constant = 0, cv = 0)
                    else{
                        pvalue <- uniroot(froot, c(1E-100, 1-1E-100))
                        results <- seekpval(pvalue$root)
                    }
                    results <- list(statistic = c(z = as.numeric(results$stat)),
                                    method = "Non-degenerate Vuong test for non-nested models",
                                    p.value = results$pvalue,
                                    data.name = data.name,
                                    alternative = "different models",
                                    parameters = c(
                                        vuong_stat = Tvuong,
                                        constant = results$constant,
                                        'sum e.v.' = sum(V),                                        
                                        'vuong_p.value' = pval_vuong))
                }
            }
            else{
                results <- list(statistic = c(z = Tvuong),
                                method = "Vuong test for non-nested models",
                                p.value = pval_vuong,
                                data.name = data.name,
                                alternative = "different models")
            }
        }
    }
    structure(results, class = "htest")
}

# a function to construct block-diagonal matrix (can't remember from
# which package it is borrowed)
bdiag <- function(...){
  if (nargs() == 1)
    x <- as.list(...)
  else
    x <- list(...)
  n <- length(x)
  if(n == 0) return(NULL)
  x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
              stop("Zero-length component in x"))
  d <- array(unlist(lapply(x, dim)), c(2, n))
  rr <- d[1,]
  cc <- d[2,]
  rsum <- sum(rr)
  csum <- sum(cc)
  out <- array(0, c(rsum, csum))
  ind <- array(0, c(4, n))
  rcum <- cumsum(rr)
  ccum <- cumsum(cc)
  ind[1,-1] <- rcum[-n]
  ind[2,] <- rcum
  ind[3,-1] <- ccum[-n]
  ind[4,] <- ccum
  imat <- array(1:(rsum * csum), c(rsum, csum))
  iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
                                               (y[3]+1):y[4]], imat=imat)
  iuse <- as.vector(unlist(iuse))
  out[iuse] <- unlist(x)
  return(out)
} 

#' @rdname ndvtest
#' @method llcont maxLik2
#' @export
llcont.maxLik2 <- function(x, ...) x$lnl


#' @rdname ndvtest
#' @method bread maxLik2
#' @export
bread.maxLik2 <- function(x, ...) - solve(x$hessian) * length(x$lnl)

#' @rdname ndvtest
#' @method estfun maxLik2
#' @export
estfun.maxLik2 <- function(x, ...) x$gradient

#' @rdname ndvtest
#' @method logLik maxLik2
#' @export
logLik.maxLik2 <- function(object, ...) sum(object$lnl)

#' @importFrom sandwich estfun
#' @export
sandwich::estfun

#' @importFrom sandwich bread
#' @export
sandwich::bread

#' @importFrom sandwich meat
#' @export
sandwich::meat

#' @importFrom nonnest2 llcont
#' @export
nonnest2::llcont

#' Simulated pdfs for the Vuong statistics using linear models
#'
#' This function can be used to reproduce the examples given of Shi (2015) which illustrate the fact that the distribution of the Vuong statistic may be very different from a standard normal
#' 
#' @aliases sim_lm
#' @param N sample size
#' @param R the number of replications
#' @param Kf the number of covariates for the first model
#' @param Kg the number of covariates for the second model
#' @param a the share of the variance of `y` explained by the two
#'     competing models
#' @return a numeric of length `N` containing the values of the Vuong statistic
#'
#' @importFrom stats lm.fit
#' 
#' @references
#'
#' \insertRef{SHI:15}{micsr}
#'
#' @examples
#' sim_lm(N = 100, R = 10, Kf = 10, Kg = 2, a = 0.5)
#' @export
sim_lm <- function(N = 1E03, R = 1E03, Kf = 15, Kg = 1, a = 0.125){
    Zf <- array(rnorm(R * Kf * N, sd = a / sqrt(Kf)), dim = c(N, Kf, R))
    Zg <- array(rnorm(R * Kg * N, sd = a / sqrt(Kg)), dim = c(N, Kg, R))
    eps <- matrix(rnorm(R * N, sd = sqrt(1 - a ^ 2)), N, R)
    y <- apply(Zf, c(1, 3), sum) + apply(Zg, c(1, 3), sum) + eps
    get_lnl <- function(x){
        res <- x$residuals
        N <- length(res)
        s2 <- sum(res ^ 2) / N
        - 1 / 2 * log(2 * pi) - 1 / 2 * log(s2) - 1 / 2 * res ^ 2 / s2
    }
    res_f <- sapply(1:R, function(i) get_lnl(lm.fit(Zf[,, i], y[, i])))
    res_g <- sapply(1:R, function(i) get_lnl(lm.fit(Zg[,, i, drop = FALSE], y[, i])))
    LR <- apply(res_f - res_g, 2, mean)
    w2 <- apply( (res_f - res_g) ^ 2, 2, mean) - LR ^ 2
    Vuong <- sqrt(N) * LR / sqrt(w2)
    Vuong
}



#' Turnout
#'
#' Turnout in Texas liquor referenda
#'
#' @name turnout
#' @docType data
#' @keywords dataset
#'
#' @format a list of three fitted models:
#' - group: the group-rule-utilitarian model,
#' - intens: the intensity model,
#' - sur: the reduced form SUR model.
#'
#' @description these three models are replication in R of stata's
#'     code available on the web site of the American Economic
#'     Association. The estimation is complicated by the fact that
#'     some linear constraints are imposed. The estimation was
#'     performed using the `maxLik` package. As the Hessian is near
#'     singular, the `bread` method for `maxLik` which use the `vcov`
#'     method returns an error. Therefore, we use a new `maxLik2`
#'     class and write specific `llcont`, `estfun` and `bread` methods
#'     for this class.
#'
#' @source
#' [American Economic Association data archive](https://www.aeaweb.org/aer/).
#' @references
#' \insertRef{COAT:CONL:04}{micsr}
#'
#' @examples
#' \dontrun{
#' ndvtest(turnout$group, turnout$intens)
#' ndvtest(turnout$group, turnout$sur)
#' ndvtest(turnout$intens, turnout$sur)
#' }
NULL
