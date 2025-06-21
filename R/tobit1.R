## Marche pas si on fournit des valeurs initiales


#' Truncated response model
#'
#' Estimation of models for which the response is truncated, either on
#' censored or truncated samples using OLS, NLS, maximum
#' likelihood, two-steps estimators or trimmed estimators
#'
#' @name tobit1
#' @param formula a symbolic description of the model; if two right
#'     hand sides are provided, the second one described the set of
#'     instruments if `scedas` is `NULL`, which is the
#'     default. Otherwise, the second part indicates the set of
#'     covariates for the variance function
#' @param data,subset,weights,na.action,offset,contrasts see `lm`
#' @param start an optional vector of starting values
#' @param left,right left and right truncation points for the response
#'     The default is respectively 0 and +Inf which corresponds to the
#'     most classic (left-zero truncated) tobit model
#' @param scedas the functional form used to specify the conditional
#'     variance, either `"exp"` or `"pnorm"`
#' @param sample either `"censored"` (the default) to estimate the
#'     censored (tobit) regression model or `"truncated"` to estimated
#'     the truncated regression model
#' @param method one of `"ml"` for maximum likelihood, `"lm"` for
#'     (biased) least squares estimators, `"twostep"` for two-steps
#'     consistent estimators, `"trimmed"` for symetrically censored
#'     estimator, `"minchisq"` and `"test"`. The last two are only
#'     relevant for instrumental variable estimation (when the formula
#'     is a two-parts formula and `scedas` is `NULL`)
#' @param opt optimization method
#' @param maxit maximum number of iterations
#' @param trace printing of intermediate result
#' @param check_gradient if `TRUE` the numeric gradient and hessian
#'     are computed and compared to the analytical gradient and
#'     hessian
#' @param object a `tobit1` object
#' @param ... further arguments
#' @importFrom stats binomial coef dnorm glm lm model.matrix
#'     model.response pnorm sigma df.residual fitted logLik
#'     model.frame printCoefmat residuals terms vcov nobs
#'     model.weights .getXlevels predict delete.response predict
#'     update
#' @keywords models
#' @return An object of class `c("tobit1", "micsr")`, see
#'     `micsr::micsr` for further details.
#' @author Yves Croissant
#' @references \insertRef{POWE:86}{micsr}
#' @examples
#' charitable$logdon <- with(charitable, log(donation) - log(25))
#' ml <- tobit1(logdon ~ log(donparents) + log(income) + education +
#'              religion + married + south, data = charitable)
#' scls <- update(ml, method = "trimmed")
#' tr <- update(ml, sample = "truncated")
#' nls <- update(tr, method = "nls")
#' @export
tobit1 <- function(formula, data, subset, weights, na.action, offset, contrasts = NULL,
                   start = NULL, left = 0, right = Inf,
                   scedas = NULL,
                   sample = c("censored", "truncated"),
                   method = c("ml", "lm", "twostep", "trimmed",
                              "nls", "minchisq", "test"),
                   opt = c("bfgs", "nr", "newton"),                  
                   maxit = 100, trace = 0,
                   check_gradient = FALSE,
                   ...){
    .call <- match.call()
    .method <- match.arg(method)
    .sample <- match.arg(sample)
    zerotrunc <- ifelse(left == 0 & is.infinite(right) & (right > 0), TRUE, FALSE)
    .scedas <- scedas
    if (! is.null(scedas)){
        if (! scedas %in% c("exp", "pnorm"))
            stop("scedas should be either equal to exp or to pnorm")
        .scedas <- .scedas
    }
    mf <- cl <- match.call(expand.dots = FALSE)
    .formula <- Formula(formula)
    if (length(.formula)[2] == 2 & is.null(.scedas)){
        mf$model <- "tobit"
        if (! .method %in% c("twostep", "minchisq", "ml", "test"))
            stop("method should be one of twostep, minchisq, ml and test")
        mf$method <- .method
        mf[[1L]] <- as.name("ivldv")#quote(micsr::ivldv())
        result <- eval(mf, parent.frame())
        result$call <- .call
        return(result)
    } else {
        if (! .method %in% c("twostep", "ml", "nls", "lm", "trimmed"))
            stop("method should be one of twostep, ml, nls and lm")
    }

    .formula <- mf$formula <- Formula(formula)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    # construct the model frame and components
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    if (length(.formula)[2] > 1){
        X <- model.matrix(.formula, mf, rhs = 1, contrasts.arg = contrasts)
        formh <- formula(.formula, rhs = 2)
        if (attr(terms(formh), "intercept") == 1L) formh <- update(formh, . ~ . - 1)
        Z <- model.matrix(formh, mf)
    }
    else{
        X <- model.matrix(.formula, mf, contrasts)
        Z <- NULL
    }
    K <- ncol(X)
    y <- model.response(mf)
    N <- length(y)
    wt <- model.weights(mf)
    if (is.null(wt)) wt <- rep(1, N) else wt <- wt #/ mean(wt)
    .offset <- model.offset(mf)
    # identify the untruncated observations
    P <- as.numeric(y > left & y < right)
    Plog <- as.logical(P)
    share_untr <- mean(Plog)
    
    # check whether the sample is censored or truncated (only when the
    # model is estimated)
    is_cens_smpl <- any(P == 0)
    if (maxit > 0){
        if (.sample == "censored" & ! is_cens_smpl)
            stop("the tobit model requires a censored sample")
        if (.method != "twostep" & is_cens_smpl & .sample == "truncated"){
            X <- X[Plog, ]
            y <- y[Plog]
            wt <- wt[Plog]
        }
    }

    # compute the starting values if they are not provided
    if (is.null(start)){
        init_lm <- cl
        init_lm <- cl[c(1L, m)]
        init_lm[[1L]] <- as.name("lm")
        if (length(.formula)[2] > 1)
            init_lm$formula <- formula(.formula, rhs = 1)
        init_lm <- eval(init_lm, parent.frame())
        coefs_init <- c(coef(init_lm), sigma = sigma(init_lm))
    } else {
        coefs_init <- start
    }
    # lm estimator, biased
    if (.method == "lm"){
        if (.sample == "truncated") result <- lm(y ~ X - 1, subset = Plog)
        else result <- lm(y ~ X - 1)
        coefs <- coef(result)
        names(coefs) <- colnames(X)
        result <- list(coefficients = coefs,
                       linear.predictor = as.numeric(X %*% coefs),
                       fitted.values = fitted(result),
                       residuals = residuals(result),
                       df.residual = df.residual(result),
                       vcov = vcov(result),
                       npar = structure(c(covariates = length(coefs)), default = "covariates"),
                       logLik = NA,
                       fomula = .formula,
                       model = model.frame(result),
                       terms = terms(model.frame(result)),
                       call = .call,
                       est_method = "lm")
    }

    # Two-steps estimator a la Heckman
    if (.method == "twostep"){
        if (! is_cens_smpl)
            stop("2 steps estimator requires a censored sample")
        pbt <- glm(P ~ X - 1, family = binomial(link = 'probit'))
        Xbs <- pbt$linear.predictor
        mls <- mills(Xbs)
        if (.sample == "truncated"){
            Z <- cbind(X, sigma = mls)
            result <- lm(y ~ Z - 1, subset = Plog)
            coefs <- coef(result)
            names(coefs) <- colnames(Z)
        }
        else{
            Z <- cbind(X * pnorm(Xbs), sigma = dnorm(Xbs))
            result <- lm(y ~ Z - 1)
            coefs <- coef(result)
            names(coefs) <- colnames(Z)
        }
        # consistent estimate of the covariance matrix
        sigma <- coefs["sigma"]
        ZZM1 <- solve(crossprod(Z))
        MRP <- mills(Xbs) * Xbs + mills(Xbs) ^ 2 
        SIGMA <- (1 - MRP)
        A <-  crossprod(Z * sqrt(SIGMA))
        DX <- X * MRP
        DXZ <- crossprod(DX, Z)
        Q <- t(DXZ) %*% vcov(pbt) %*% DXZ
        V <- sigma ^ 2 * ZZM1 %*% (A + Q) %*% ZZM1
        result <- list(coefficients = coefs,
                       linear.predictor = as.numeric(Z %*% coefs),
                       fitted.values = fitted(result),
                       residuals = residuals(result),
                       df.residual = df.residual(result),
                       formula = .formula,
                       vcov = V,
                       logLik = NA,
                       npar = structure(c(covariates = length(coefs)), default = "covariates"),
                       model = model.frame(result),
                       terms = terms(model.frame(result)),
                       call = .call,
                       est_method = "twostep")
    }

    # trimmed estimator
    if (.method == "trimmed"){
        if (! zerotrunc) stop("trimmed estimator only implemented for simple tobit")
        trim_trunc <- function(param, X, y, sum = TRUE, gradient = FALSE, hessian = TRUE){
            bX <- as.numeric(X %*% param)
            f <- (y - pmax(1 / 2 * y, bX)) ^ 2
            if (sum) f <- sum(f)
            if (gradient | hessian) ymin <- pmin(y, 2 * bX)
            if (gradient){
                grad <- - 2 * (y < (2 * bX)) * (y - bX) * X
                if (sum) grad <- apply(grad, 2, sum)
                attr(f, "gradient") <- grad
            }
            if (hessian) attr(f, "hessian")   <- 2 * crossprod( (y < (2 * bX)) * X, X)
            f
        }
        trim_cens <- function(param, X, y, sum = TRUE, gradient = FALSE, hessian = FALSE){
            sgn <- function(x) ifelse(x > 0, 1, -1)
            bX <- as.numeric(X %*% param)
            f <- (bX < 0) * (y ^ 2 / 2) +
                (bX > 0 & y < (2 * bX)) * ((y - bX) ^ 2)+
                (bX > 0 & y > (2 * bX)) * (y ^ 2 / 2 - bX ^ 2)
            if (sum) f <- sum(f)
            if (gradient | hessian) ymin <- pmin(y, 2 * bX)
            if (gradient){
                grad <- - 2 * (bX > 0)* (ymin - bX) * X
                if (sum) grad <- apply(grad, 2, sum)
                attr(f, "gradient") <- grad
            }
            if (hessian) attr(f, "hessian")   <- 2 * crossprod( (bX > 0) * sgn(2 * bX - y) * X, X)
            f
        }
        coefs <- coefs_init[1:(length(coefs_init) - 1)]
        if (.sample == "truncated") coefs <- newton(trim_trunc, coefs, trace = trace, X = X, y = y)
        else coefs <- newton(trim_cens, coefs, trace = trace, X = X, y = y)
        vcov_trim <- function(coefs){
            # only relevant for censored samples
            bX <- as.numeric(X %*% coefs)
            N <- sum(y > 0 & y < (2 * bX))
            N <- length(y)
            C_mat <- crossprod( (y > 0 & y < (2 * bX)) * X, X) / N
            D_mat <- crossprod( (bX > 0) * (pmin(y, 2 * bX) - bX) ^ 2 * X, X) / N
            solve(C_mat) %*% D_mat %*% solve(C_mat) / N
        }
        bX <- as.numeric(X %*% coefs)
        status <- rep(NA, nrow(X))
        status[y == 0] <- "left-cens"
        status[y > 2 * bX & bX > 0] <- "right-trimmed"
        status[bX < 0 & y > 0] <- "neg-linpred"
        .df.residual <- nrow(X) - length(coefs)
        result <- list(coefficients = coefs,
                       linear.predictor = bX,
                       fitted.values = bX,
                       residuals = y - bX,
                       df.residual = .df.residual,
                       formula = .formula,
                       vcov = vcov_trim(coefs),
                       npar = structure(c(covariates = length(coefs)), default = "covariates"),
                       logLik = NA,
                       model = mf,
                       terms = NA,
                       status = status,
                       call = .call,
                       est_method = "trimmed")
    }
    
    # non-linear least-squares
    if (.method == "nls"){
        if (! zerotrunc | .sample != "truncated")
            stop("nls estimator only implemented for simple truncated model")
        if (.sample == "truncated"){
            fun_nls <-  function(x, gradient = FALSE, hessian = FALSE){
                beta <- x[1:ncol(X)]
                sigma <- x[ncol(X) + 1]
                e <- as.numeric(y - X %*% beta - sigma * mills(X %*% beta / sigma))
                bXs <- as.numeric(X %*% beta) / sigma
                if (gradient | hessian){
                    e_beta <- - (1 + mills(bXs, 1))
                    e_sigma <- mills(bXs, 1) * bXs - mills(bXs)
                    e_beta_beta <- - mills(bXs, 2) / sigma
                    e_beta_sigma <- mills(bXs, 2) * bXs / sigma
                    e_sigma_sigma <- - mills(bXs, 2) * bXs ^ 2 / sigma
                    grad_beta <-  (- 2 * e * e_beta) * X
                    grad_sigma <- - 2 * e * e_sigma
                }
                if (hessian){
                    hess_beta_beta <- crossprod(- 2 * (e * e_beta_beta + e_beta ^ 2) * X, X)
                    hess_beta_sigma <- apply(   - 2 * (e * e_beta_sigma + e_beta * e_sigma) * X, 2, sum)
                    hess_sigma_sigma <- sum(    - 2 * (e * e_sigma_sigma + e_sigma ^ 2))
                }
                f <- sum(e ^ 2)
                if (gradient){
                    g <- - cbind(grad_beta, grad_sigma)
                    attr(f, "gradient") <- g
                }
                if (hessian){
                    h <- - rbind(cbind(hess_beta_beta, hess_beta_sigma),
                                 c(hess_beta_sigma, hess_sigma_sigma))
                    attr(f, "hessian") <- h
                }
                f
            }
            

            coefs <- coefs_init
            coefs <- newton(fun_nls, coefs, trace = trace)
            f <- fun_nls(coefs, gradient = TRUE, hessian = TRUE)
            h <- attr(f, "hessian")
            g_n <- attr(f, "gradient")
            B <- crossprod(g_n)
            A <- h
            .vcov <- solve(A) %*% B %*% solve(A)# / length(y)
        }
        .terms <-  terms(mf)
        attr(.terms, ".Environment") <- NULL
        bX <- as.numeric(X %*% coefs[1:ncol(X)])
        result <- list(coefficients = coefs,
                       linear.predictor = bX,
                       fitted.values = bX,
                       residuals = y - bX,
                       formula = .formula,
                       df.residual = length(y) - length(coefs),
                       vcov = .vcov,
                       npar = structure(c(covariates = length(coefs)), default = "covariates"),
                       logLik = NA,
                       model = mf,
                       terms = .terms,
                       call = .call,
                       est_method = "nls")
    }
    
    # Maximum likelihood estimator
    if (.method == "ml"){
        coefs_init[1:K] <- coefs_init[1:K] / coefs_init[K + 1]
        coefs_init[K + 1] <- 1 / coefs_init[K + 1]
        coefs <- newton(lnl_tp_olsen, coefs_init, trace = trace, X = X, y = y, wt = wt,
                        sum = FALSE, left = left, right = right,
                        direction = "max", sample = .sample, maxit = maxit)
        coefs[1:K] <- coefs[1:K] / coefs[K + 1]
        coefs[K + 1] <- 1 / coefs[K + 1]
        if (is.null(Z)){
            lnl_conv <- lnl_tp(coefs, X = X, y = y, wt = wt,
                               sum = FALSE, gradient = TRUE, hessian = TRUE,
                               left = left, right = right, sample = .sample)
            npar <- structure(c(covariates = ncol(X), vcov = 1), default = c("covariates", "vcov"))
        }
        else{
            sup_coef <- rep(0, ncol(Z))
            names(sup_coef) <- paste("sig_", colnames(Z), sep = "")
            coefs <- c(coefs, sup_coef)
            coefs <- newton(lnl_tp, coefs, trace = TRUE, X = X, y = y, wt = wt,
                            Z = Z, scedas = .scedas,
                            left = left, right = right, direction = "max", sample = .sample)
            lnl_conv <- lnl_tp(coefs, X = X, y = y, wt = wt,
                               scedas = .scedas, Z = Z,
                               sum = FALSE, gradient = TRUE, hessian = TRUE,
                               left = left, right = right, sample = .sample)
            npar <- structure(c(covariates = ncol(X), scedas = ncol(Z) + 1), default = c("covariates"))
        }
        fun <- function(x) lnl_tp(x, X = X, y = y, wt = wt,
                                  scedas = .scedas, Z = Z,
                                  sum = TRUE, gradient = TRUE, hessian = TRUE,
                                  left = left, right = right, sample = .sample)
        if(check_gradient) z <- check_gradient(fun, coefs) else z <- NA

        .hessian <- attr(lnl_conv, "hessian")
        .info <- attr(lnl_conv, "info")
        .gradient <- attr(lnl_conv, "gradient")
        .logLik <- sum(as.numeric(lnl_conv))
        beta <- coefs[1:K]
        sigma <- coefs[K + 1]
        .linpred <- as.numeric(X %*% beta)
        h <- .linpred / sigma
        Ppos <- pnorm(h)
        Epos <- .linpred + sigma * mills(h)
        .fitted <- data.frame(y = y, Ppos = Ppos, Epos = Epos, lp = .linpred)
        class(.fitted) <- c("tbl_df", "tbl", "data.frame")
#        .vcov <- solve(- .hessian)
#        dimnames(.vcov) <- list(names(coefs), names(coefs))
        .terms <-  terms(mf)
        attr(.terms, ".Environment") <- NULL

        # tobit without covariates (mu and sigma), only for the standard tobit
        if (left == 0 & is.infinite(right) & .sample == "censored" & maxit > 0){
            coef0 <- function(x){
                mu2 <- mean(x[x > 0] ^ 2)
                yb <- mean(x[x > 0])
                pr <- mean(x > 0)
                alpha <- qnorm(pr)
                h <- (alpha + dnorm(alpha) / pnorm(alpha)) / yb
                za <- function(h){
                    alpha <-  (h * mu2 - 1 / h) / yb
                    pr * alpha - pr * h * yb + (1 - pr) * mills(- alpha)
                }
                .sig <- 1 / h
                h_min <- 1 / (.sig * 10)
                h_max <- 1 / (.sig / 2)
                h_conv <- uniroot(za, c(0, h_max), tol = 1E-10)$root
                alpha_conv <- (h_conv * mu2 - 1 / h_conv) / yb
                .sig <- 1 / h_conv
                .mu <- .sig * alpha_conv
                .lnL <- (1 - pr) * pnorm(- alpha_conv, log.p = TRUE) - pr * log(2 * pi) / 2 +
                    pr * log(h_conv) - pr / 2 * h_conv ^ 2 * mu2 - pr / 2 * alpha_conv ^ 2 +
                    pr * alpha_conv * h_conv * yb
                c(mu = .mu, sigma = .sig, lnl = .lnL)
            }
            values <- cbind(model     = as.numeric(lnl_conv),
                            saturated = (- 1 / 2 * log(2 * pi)) * (y > 0) + 0 * (y == 0),
                            null      = coef0(y)["lnl"])
            .logLik <- apply(values * wt, 2, sum)
        } else {
            values <- as.numeric(lnl_conv)
            .logLik <- c(model = sum(as.numeric(lnl_conv)))
        }

        d <- ifelse(y > 0, 1, 0)

        gres <- - (1 - d) * sigma * dnorm(h) / (1 - pnorm(h)) + d * (y - .linpred)
        gres <- - (y == 0) * sigma * dnorm(h) / (1 - pnorm(h)) + (y > 0) * (y - .linpred)

        env <- new.env(parent = .GlobalEnv)
        .sigma <- coefs["sigma"]
        assign(".sigma", .sigma, envir = env)
        .variance <- function(mu) .sigma ^ 2 * (1 + mills(mu / .sigma, deriv = 1)) * pnorm(mu / .sigma)
        .gres <- function(mu) - (y == 0) * .sigma * dnorm(mu / .sigma) / (1 - pnorm(mu / .sigma)) + (y > 0) * (y - mu)
        environment(.variance) <- env
        
        .variance
        .family <- structure(list(
            link    = "identity",
            linkfun = function(mu) mu,
            linkinv = function(eta) eta,
            mu.eta  = function(eta) 1,
            variance = .variance),
            class = "family")
        
        result <- list(coefficients = coefs,
                       model = mf,
                       terms = mt,
                       value = values,
                       gradient = .gradient,
                       hessian = .hessian,
                       info = .info,
                       fitted.values = .fitted,
                       linear.predictor = .linpred,
                       logLik = .logLik,
#                       residuals = y - .fitted$Epos,
                       df.residual = length(y) - length(coefs),
                       npar = npar,
                       est_method = "ml",
                       call = .call,
                       na.action = attr(mf, "na.action"),
                       weights = wt,
                       offset = .offset,
                       contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(mt, mf),
                       check_gradient = z,
                       family = .family,
                       residuals = gres)
#                       formula = .formula,
#                       vcov = .vcov,
    }
    structure(result, class = c("tobit1", "micsr"))
}


#' @rdname tobit1
#' @export
fitted.tobit1 <- function(object, ...) object$fitted.values$Ppos * object$fitted.values$Epos
#fitted.tobit1 <- function(object, ...) object$fitted.values$lp
