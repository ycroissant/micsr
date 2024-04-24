#' Endogenous switching and sample selection models for count data
#'
#' Heckman's like estimator for count data, using either
#' maximum likelihood or a two-steps estimator
#' 
#' @name escount
#' @param formula a `Formula` object which includes two responses (the
#'     count and the binomial variables) and two sets of covariates
#'     (for the count component and for the selection equation)
#' @param data a data frame,
#' @param subset,weights,na.action,offset see `stats::lm`
#' @param start an optional vector of starting values,
#' @param R the number of points for the Gauss-Hermite quadrature
#' @param hessian if `TRUE`, the numerical hessian is computed,
#'     otherwise the covariance matrix of the coefficients is computed
#'     using the outer product of the gradient
#' @param method one of `'ML'` for maximum likelihood estimation (the
#'     default) or `'twosteps'` for the two-steps NLS method
#' @param model one of `'es'` for endogenous switching (the default)
#'     or `'ss'` for sample selection
#' @return an object of class `c("escount,micsr)"`, see `micsr::micsr` for further details.
#' @importFrom stats .getXlevels coef glm model.matrix model.response
#'     model.frame model.offset model.weights formula fitted nobs
#'     optim pchisq pnorm poisson printCoefmat vcov binomial dnorm
#'     dpois
#' @keywords models
#' @importFrom Formula Formula model.part
#' @author Yves Croissant
#' @references \insertRef{TERZ:98}{micsr}
#' 
#' \insertRef{GREE:01}{micsr}
#' @importFrom Rdpack reprompt
#' @importFrom Formula Formula model.part
#' @examples
#' trips_2s <- escount(trips | car ~ workschl + size + dist + smsa + fulltime + distnod +
#' realinc + weekend + car | . - car - weekend + adults, data = trips, method = "twosteps")
#' trips_ml <- update(trips_2s, method = "ml")
#' @export
escount <- function(formula,
                    data,
                    subset,
                    weights,
                    na.action,
                    offset,
                    start = NULL,
                    R = 16,
                    hessian = FALSE,
                    method = c("twosteps", "ml"),
                    model = c("es", "ss")){
    start_provided <- ! is.null(start)
    model <- match.arg(model)
    ss <- ifelse(model == "ss", 1, 0)
    .est_method <- match.arg(method)
    # An efficient function for computing the inverse Mills ratio
    mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
    lmills <- function(x) dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE)
    
    # Compute the model frame
    .call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    .formula <- mf$formula <- Formula::Formula(formula)
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
    names_X <- colnames(X)
    Z <- model.matrix(.formula, mf, rhs = 2)
    names_Z <- colnames(Z)
    y <- model.part(.formula, mf, lhs = 1, drop = TRUE)
    d <- model.part(.formula, mf, lhs = 2, drop = TRUE)
    q <- 2 * d - 1
    N <- length(y)
    K <- ncol(X)
    L <- ncol(Z)

    # Compute starting values (beta and theta, using the mills ratio
    # approximation
    binom <- glm(d ~ Z - 1, family = binomial(link = 'probit'))
    bZ <- drop(Z %*% coef(binom))

    if (model == "ss"){
        mills_corr <- lmills(bZ)
        count <- glm(y ~ X + mills_corr - 1, family = poisson, subset = d == 1)
        .start <- unname(coef(count))
    }
    else{
        mills_corr <- lmills(q * bZ)
        count <- glm(y ~ X + mills_corr - 1, family = poisson)
        count <- glm(y ~ X - 1, family = poisson)        
        .start <- c(unname(coef(count)), 0)
    }
    
        # Compute the expected values
    ## Exp_y <- function(param, X, Z){
    ##     alpha <- param[1:L]
    ##     beta <- param[L + c(1:K)]
    ##     theta <- param[K + L + 1]
    ##     bX <- drop(X %*% beta)
    ##     aZ <- drop(Z %*% alpha)
    ##     log_Ey <- bX +
    ##         pnorm(q * (theta + aZ), log.p = TRUE) -
    ##         pnorm(q * aZ, log.p = TRUE)
    ##     exp(log_Ey)
    ## }

    obs_count <- (d * ss + (1 - ss))
    # pos is 1 for d = 1 if ss = 1 and 1 for all the observations
    # otherwise
    if (.est_method == "ml"){
        # The likelihood function for the Terza model, approximated using
        # gauss-hermite quadratures
        lnl <- function(param, gradient = FALSE, sum = FALSE, R = 16){
#            rn <- statmod::gauss.quad(R, kind = "hermite")
            rn <- gaussian_quad(R, "hermite")
            alpha <- param[1:L]
            beta <- param[L + c(1:K)]
            beta <- param[1:K]
            alpha <- param[K + c(1:L)]
            sig <- param[K + L + 1]
            rho <- param[K + L + 2]
            bX <- drop(X %*% beta)
            aZ <- drop(Z %*% alpha)
            # l is the contrubution of an individual to the likelihood for
            # a given value of the random deviate, g_## are the
            # corresponding derivatives
            l <- function(eps){
                mu <- exp(bX + sig * sqrt(2) * eps)
                ln_mu <- bX + sig * sqrt(2) * eps
                v <- (aZ + rho * sqrt(2) * eps) / sqrt(1 - rho ^ 2)
                lnl_count <- obs_count * (- exp(ln_mu) + y * ln_mu - lfactorial(y))
                lnl_binom <- pnorm(q * v, log.p = TRUE)
                exp(lnl_count + lnl_binom)
            }
            g_beta <- function(eps){
                mu <- exp(bX + sig * sqrt(2) * eps)
                (y - mu) * obs_count
            }
            g_alpha <- function(eps){
                v <- (aZ + rho * sqrt(2) * eps) / sqrt(1 - rho ^ 2)
                q * mills(q * v) / sqrt(1 - rho ^ 2)
            }
            g_sig <- function(eps){
                mu <- exp(bX + sig * sqrt(2) * eps)
                (sqrt(2) * eps * (y - mu)) * obs_count
            }
            g_rho <- function(eps){
                v <- (aZ + rho * sqrt(2) * eps) / sqrt(1 - rho ^ 2)
                (sqrt(2) * eps * sqrt(1 - rho ^ 2) +
                 rho * (aZ + rho * sqrt(2) * eps) /
                 sqrt(1 - rho ^ 2)) / (1 - rho ^ 2) *
                    (2 * d - 1) * mills(q * v)
            }
            # the contributions of an individual to the likelihood are
            # computed for the R values of eps, giving a matrix of
            # dimension N x R
            l_nr <- sapply(rn$nodes, l)
            # the contributions of an individual to the likelihood are
            # computed as a linear combination (with relevant weights) of
            # the R values
            l_n <- apply(t(l_nr) * rn$weights, 2, sum) / sqrt(pi)
            # the contribution of the log-likelihood function
            lnl <- log(l_n)
            if (sum) lnl <- sum(lnl)
            if (gradient){
                # Compute the N x R matrix for all the components of the
                # gradient
                g_beta_r <- sapply(rn$nodes, g_beta) * l_nr
                g_alpha_r <- sapply(rn$nodes, g_alpha) * l_nr
                g_sig_r <- sapply(rn$nodes, g_sig) * l_nr
                g_rho_r <- sapply(rn$nodes, g_rho) * l_nr
                # Get a vector as linear combinations of the columns of
                # the matrices
                g_beta_r <- apply(t(g_beta_r) * rn$weights, 2, sum) / sqrt(pi)
                g_alpha_r <- apply(t(g_alpha_r) * rn$weights, 2, sum) / sqrt(pi)
                g_sig_r <- apply(t(g_sig_r) * rn$weights, 2, sum) / sqrt(pi)
                g_rho_r <- apply(t(g_rho_r) * rn$weights, 2, sum) / sqrt(pi)
                # Compute the matrix (N x (K + L + 2)) of the individual
                # contributions to the gradient
                grad <- cbind(g_beta_r * X, g_alpha_r * Z, g_sig_r, g_rho_r) / l_n
                if (sum) grad <- apply(grad, 2, sum)
                attr(lnl, "gradient") <- grad
            }
            lnl
        }
        # functions for the log-likelihood and its gradient relevant for
        # optim

        fn <- function(param) - lnl(param, gradient = FALSE, sum = TRUE)
        gr <- function(param) - attr(lnl(param, gradient = TRUE, sum = TRUE), "gradient")
        gr_obs <- function(param) attr(lnl(param, gradient = TRUE, sum = FALSE), "gradient")
        fn_obs <- function(param) lnl(param, gradient = FALSE, sum = FALSE)

        # if starting values are not provided, use the 2-steps
        # estimator
        if (! start_provided){
            .ncall <- .call
            .ncall$method <- "twosteps"
            tsp <- eval(.ncall, parent.frame())
            .beta <- tsp$coefficients
            .beta <- .beta[- length(.beta)]
            .sigma <- tsp$sigma
            .rho <- tsp$rho
            .start <- c(.beta, coef(binom), sigma = .sigma, rho = .rho)
        }
        else .start <- start
        names(.start) <- c(paste("covar", names_X, sep = "_"), paste("instr", names_Z, sep = "_"), "sigma", "rho")
        z <- optim(.start, fn, gr, method = "BFGS", hessian = TRUE, control = list(trace = 0))
        .coef <- z$par
        .logLik <- z$value
        .gradient <- gr_obs(.coef)
        .hessian <- - z$hessian
        .lnl <- fn_obs(.coef)
        logLik_model <- c(model = sum(.lnl))
        
        bX <- drop(X %*% .coef[1:K])
        aZ <- drop(Z %*% .coef[K + c(1:L)])
        theta <- unname(prod(.coef[K + L + (1:2)]))
        .fitted <- exp(bX + pnorm(q * (theta + aZ), log.p = TRUE) -
                       pnorm(q * aZ, log.p = TRUE))
        .resid <- y - .fitted
        .npar <- c(covariates = K, instruments = L, vcov = 2)
        attr(.npar, "default") <- c("covariates", "vcov")
        result <- list(coefficients = .coef,
                       residuals = .resid,
                       fitted.values = .fitted,
                       model = mf,
                       terms = mt,
                       logLik = logLik_model,
                       value = .lnl,
                       npar = .npar,
                       gradient = .gradient,
                       hessian = .hessian,
                       df.residual = N - K - L,
                       xlevels = .getXlevels(mt, mf),
                       na.action = attr(mf, "na.action"),
                       call = .call,
                       est_method = "ml"
                       )
    }
    if (.est_method == "twosteps"){
        if (model == "ss"){
            oy <- y; oX <- X; oZ <- Z; oq <- q
            pos <- d > 0
            y <- y[pos]
            X <- X[pos, ]  
            Z <- Z[pos, ]
            q <- q[pos]
        }
        aZ <- drop(Z %*% coef(binom)) # ssr returns the sum of squares residuals as a function of
                                        # beta and theta (Z %*% alpha is taken from the first-step
                                        # probit regression). It optionnaly returns the gradient
        ssr <- function(param, gradient = FALSE){
            beta <- param[1:K]
            theta <- param[K + 1]
            bX <- drop(X %*% beta)
            log_Ey <- bX + pnorm(q * (theta + aZ), log.p = TRUE) -
                pnorm(q * aZ, log.p = TRUE)
            Ey <- exp(log_Ey)
            e <- (y - Ey)
                                    result <- sum(e ^ 2)
            if (gradient){
                grad_beta <- Ey * X
                grad_theta <- Ey * mills(q * (theta + aZ)) * q
                gradient <- cbind(grad_beta, theta = grad_theta)
                gradient <- - 2 * e * gradient
                                        attr(result, "gradient") <- apply(gradient, 2, sum)
            }
            result
        }
        # Two function two separately compute the objective function and the gradient"
        obj_ssr <- function(param) ssr(param, gradient = FALSE)
        grad_ssr <- function(param) attr(ssr(param, gradient = TRUE), "gradient")
        # Compute the minimum
        names(.start) <- c(colnames(X), "theta")
        conv_nls <- optim(.start, obj_ssr, grad_ssr, method = "BFGS", control = list(maxit = 1000))
        .value <- conv_nls$value
        # Compute the linear predictors for the fitted coefficients
        hbeta <- conv_nls$par[1:K]
        htheta <- unname(conv_nls$par[K + 1])
        halpha <- coef(binom)
        if (model == "ss"){
            y <- oy; X <- oX; Z <- oZ; q <- oq
        }
        bX <- drop(X %*% hbeta)
        aZ <- drop(Z %*% halpha)
        # Estimate sigma using the previous fitted coefficient using
        # gaussian-quadrature in a one-parameter likelihood function

        # l is returns the individual contribution to the likelihood
        # for given eps and a value of theta
        l <- function(eps, sig){
            v <- (aZ + (htheta / sig) * sqrt(2) * eps) / sqrt(1 - htheta ^ 2 / sig ^ 2)
            mu <- exp(bX - 0.5 * sig ^ 2 + sqrt(2) * sig * eps)
            ln_mu <- bX - 0.5 * sig ^ 2 + sqrt(2) * sig * eps
            lnl_count <- - exp(ln_mu) + y * ln_mu - lfactorial(y)
            lnl_binom <- pnorm(q * v, log.p = TRUE)
            exp(lnl_count * d + lnl_binom)
        }

        # objFun compute the log-likelihood function as a function of
        # sig only
#        rn <- statmod::gauss.quad(R, kind = "hermite")
        rn <- gaussian_quad(R, kind = "hermite")

        sig <- abs(2 * htheta)
        objFun <- function(sig){
            le <- function(eps) l(eps, sig)
            l_nr <- sapply(rn$nodes, le)
            l_n <- apply(t(l_nr) * rn$weights, 2, sum) / sqrt(pi)
            lnl_n <- log(l_n)
            - sum(lnl_n)
        }
        # Compute the ML estimator of sigma
        conv_sig <- optim(sig, objFun, method = "BFGS")
        hsigma <- conv_sig$par
        hrho <- htheta / hsigma
        if (hsigma < 0){
            hsigma <- - hsigma
            hrho <- - hrho
        }
        h_param <- c(coef(binom), conv_nls$par)
        # covariance matrix of the estimator: robust to
        # heteroscedasticity and taking into account the fact that
        # alpha is estimated in the first step
        if (model == "ss"){
            pos <- d > 0
            y <- y[pos]
            X <- X[pos, ]  
            Z <- Z[pos, ]
            q <- q[pos]
            aZ <- aZ[pos]
        }
        Ey <- exp(bX + pnorm(q * (htheta + aZ), log.p = TRUE) -
            pnorm(q * aZ, log.p = TRUE))
        e <- (y - Ey)
        G1 <- (Ey * cbind(X, mills(aZ + htheta) * q))
        G2 <- (Ey * (mills(aZ + htheta) - mills(aZ)) * q * Z)
        V_binom <- summary(binom)$cov.scaled
        G1G1I <- solve(crossprod(G1))
        .vcov <- G1G1I %*% (crossprod(G1 * e) +
                            crossprod(G1, G2) %*% V_binom %*%
                            crossprod(G2, G1)) %*% G1G1I
        result <- list(coefficients = conv_nls$par,
                       sigma = hsigma,
                       rho = hrho,
                       vcov = .vcov,
                       residuals = e,
                       fitted.values = Ey,
                       model = mf,
                       terms = mt,
                       value = .value,
                       npar = c(covariates = K, vcov = 1),
                       df.residual = N - ncol(X),
                       xlevels = .getXlevels(mt, mf),
                       na.action = attr(mf, "na.action"),
                       call = .call,
                       first = binom,
                       est_method = "twosteps"
                       )
    }
    result <- structure(result, class = c("escount", "micsr"))
}

