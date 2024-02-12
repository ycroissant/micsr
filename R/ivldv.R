#' Instrumental variable estimators for limited dependent variable
#'
#' Estimation of simultaneous-equation models when the response is
#' binomial or censored
#'
#' @name ivldv
#' @param x on object returned by `ivldv`
#' @param formula a symbolic description of the model,
#' @param data a data frame,
#' @param subset,weights,na.action,offset see `lm`,
#' @param method one of `"ml"` for maximum likelihood, "twosteps"` and
#'     `"minchisq"`
#' @param model one of `"probit"` or `"tobit"`,
#' @param robust a boolean, if `TRUE`, a consistent estimation of the
#'     covariance of the coefficients is used for the 2-steps method,
#' @param left,right left and right limits of the dependent
#'     variable. The default is respectively 0 and +Inf which
#'     corresponds to the most classic (left-zero truncated) tobit
#'     model,
#' @param trace a boolean (the default if `FALSE`) if `TRUE` some
#'     information about the optimization process is printed,
#' @param ... further arguments
#' @importFrom Formula Formula
#' @importFrom stats coef dnorm lm model.matrix model.response pnorm
#'     fitted model.frame residuals optim update glm binomial pchisq
#'     vcov terms .getXlevels delete.response logLik printCoefmat
#' @return An object of class `c('ivldv', 'lm')`
#' @author Yves Croissant
#' @keywords models
#' @references
#' \insertRef{SMIT:BLUN:86}{micsr}
#'
#' \insertRef{RIVE:VUON:88}{micsr}
#' @importFrom Rdpack reprompt
#' @export
#' @examples
#' inst <- ~ sic3 + k_serv + inv + engsci + whitecol + skill + semskill + cropland + 
#'     pasture + forest + coal + petro + minerals + scrconc + bcrconc + scrcomp +
#'     bcrcomp + meps + kstock + puni + geog2 + tenure + klratio + bunion
#' trade_protection <- dplyr::mutate(micsr::trade_protection,
#'                                  y = ntb / (1 + ntb),
#'                                  x1 = vshipped / imports / elast,
#'                                  x2 = cap * x1,
#'                                  x3 = labvar)
#' GH <- ivldv(Formula::as.Formula(y  ~  x1 + x2, inst), trade_protection,
#'             method = "twosteps", model = "tobit") 
#' Full <- ivldv(Formula::as.Formula(y ~ x1 + x2 + labvar, inst), trade_protection,
#'               method = "twosteps", model = "tobit") 
#' Short <- ivldv(Formula::as.Formula(y ~ x1 + I(x2 + labvar), inst),
#'                  trade_protection, method = "twosteps", model = "tobit")
#' bank_msq <- ivldv(federiv ~ eqrat + optval + bonus + ltass + linsown + linstown +
#'                   roe + mktbk + perfor + dealdum + div + year | . - eqrat - bonus -
#'                   optval + no_emp + no_subs + no_off + ceo_age + gap + cfa,
#'                   data = federiv, method = "minchisq")
#' bank_ml <- update(bank_msq, method = "ml")
#' bank_2st <- update(bank_msq, method = "twosteps")
ivldv <- function(formula,
                  data,
                  subset = NULL,
                  weights = NULL,
                  na.action,
                  offset,
                  method = c("twosteps", "minchisq", "ml", "test"),
                  model = c("probit", "tobit"),
                  robust = TRUE,
                  left = 0,
                  right = Inf,
                  trace = 0,
                  ...){
    rm_intercept <- function(x){
        pos <- match("(Intercept)", colnames(x))
        if (! is.na(pos)) x[, - pos, drop = FALSE]
        else x
    }
    .est_method <- match.arg(method)
    compute_test <- (.est_method == "test")
    model <- match.arg(model)
    if (compute_test){
        .est_method <- "twosteps"
        robust <- FALSE
    }
    .call <- match.call()
    mf <- match.call(expand.dots = TRUE)
    has.int <- TRUE # a voir plus tard
    .formula <- mf$formula <- Formula(formula)
    m <- match(c("formula", "data", "subset", "weights"),
               names(mf), 0L)

    # construct the model frame and components
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf$dot <- "previous"
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    # Extract the element of the model
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    offset <- model.offset(mf)
    # response and its transformation
    y <- model.response(mf)
    q <- 2 * y - 1
    # model matrices
    ZZ <- rm_intercept(model.matrix(Formula(formula), mf, rhs = 1))  # L covariates (X1 + Y)
    X <- rm_intercept(model.matrix(Formula(formula), mf, rhs = 2))  # K exogenous variables (X1 + X2)
    names_X <- colnames(X)
    W <- ZZ[, setdiff(colnames(ZZ), colnames(X)), drop = FALSE]       # G endogenous variables (W)
    X1 <- ZZ[, intersect(colnames(ZZ), colnames(X)), drop = FALSE]    # K1 exogenous covariates (X1)
    X2 <- X[, setdiff(colnames(X), colnames(X1)), drop = FALSE]    # K1 exogenous covariates (X1)
    ZZ <- cbind(X1, W)
    names_ZZ <- colnames(ZZ)
    if (has.int){
        names_X <- c("(Intercept)", names_X)
        names_ZZ <- c("(Intercept)", names_ZZ)
    }
    N <- length(y)
    K1 <- ncol(X1)
    K <- ncol(X)
    G <- ncol(W)
    K2 <- K - K1
    L <- K1 + G
    if (has.int){
        X1 <- cbind("(Intercept)" = 1, X1)
        L <- L + 1
    }

    # estimation of the reduced forms for the endogenous variables
    step_1 <- lm(W ~ X)
    Wres <- residuals(step_1)
    Yhat <- fitted(step_1)
    Yres <- residuals(step_1)
    if (model == "probit") step_2 <- glm(y ~ ZZ + Wres, family = binomial(link = 'probit'))
    if (model == "tobit")  step_2 <- tobit1(y ~ ZZ + resid(step_1), left = left, right = right)

    lp <- step_2$linear.predictor
    if (.est_method == "twosteps"){
        .coef <- coef(step_2)
        nms_coef <- c("(Intercept)", colnames(ZZ), paste("rho", colnames(W), sep = "_"))
        if (model == "tobit") nms_coef <- c(nms_coef, "sigma")
        S_epsilon <- crossprod(Wres) / N
        # compute the consistent covariance matrix if required
        U <- cbind("(Intercept)" = 1, ZZ, Wres)
        X <- cbind("(Intercept)" = 1, X)
        if (model == "tobit"){
            .sigma <- .coef["sigma"]
            z <- lp / .sigma
            phi <- dnorm(z)
            Phi <- pnorm(z)
            a11 <- - (phi * z - phi ^ 2 / (1 - Phi) - Phi) / .sigma ^ 2
            a12 <-  (phi * z ^ 2 + phi - z * phi ^ 2 / (1 - Phi)) /  .sigma ^ 2
            a22 <- - (phi * z ^ 3 + phi * z - phi ^ 2 * z ^ 2 / (1 - Phi) - 2 * Phi) / .sigma ^ 2
            UAU <- rbind(cbind(  crossprod(U * a11, U), crossprod(U, a12)),
                         cbind(t(crossprod(U, a12)),    sum(a22)))
            UAX <- rbind(crossprod(U * a11, X), apply(X * a12, 2, sum))
        }
        if (model == "probit"){
            dmls <- dmills(q * lp)
            UAU <- crossprod(U * sqrt(- dmls))
            UAX <- crossprod(U * (- dmls), X)
        }
        if (robust){
            .rho <- coef(step_2)[K1 + 1 + (1:G)]
            delta <- as.numeric(t(.rho) %*% S_epsilon %*% .rho)
            Q <- solve(UAU) %*% UAX %*% solve(crossprod(X)) %*% t(UAX) %*% solve(UAU)
            .vcov <- solve(UAU) + delta * Q
        }
        else .vcov <- solve(UAU)
        .vcov_raw <- solve(UAU)
        names(.coef) <- colnames(.vcov) <- rownames(.vcov) <- nms_coef
    }
    
    if (compute_test){
        .data.name <- paste(deparse(substitute(formula)))[1]
        .data.name <- paste(.data.name, "...")
        coef_res <- grep("rho", names(.coef))
        B <- .coef[coef_res]
        V <- .vcov[coef_res, coef_res, drop = FALSE]
        stat <- as.numeric(crossprod(B,solve(V, B)))
        pval <- pchisq(stat, df = length(coef_res), lower.tail = FALSE)
        res <- list(statistic    = c(chisq = stat),
                    p.value      = pval,
                    parameter    = c(df = length(coef_res)),
                    method       = "Smith-Blundell / Rivers-Vuong test",
                    data.name    = .data.name,
                    alternative  = "endogeneity")
        class(res) <- "htest"
        return(res)
    }
    
    if (.est_method == "minchisq"){
        # Newey_2. probit /tobit with all the exogenous and the residuals of the first step
        if (model == "probit"){
            newey_2 <- glm(y ~ X + Wres, family = binomial(link = 'probit'))
            # this is not the default vcov computed by glm
            vcov_step_2 <- solve(crossprod(sqrt(- dmills(q * drop(cbind(1, X, Yres) %*% coef(newey_2)))) * cbind(1, X, Yres)))
        }
        if (model == "tobit"){
            newey_2 <- tobit1(y ~ X + Wres, left = left, right = right)
#            vcov_step_2 <- vcov(newey_2)
            vcov_step_2 <- newey_2$vcov
        }
        coef_newey_2 <- newey_2$coefficients
        alpha_hat <- coef_newey_2[1:(K + 1)]
        coef_res_w <- coef_newey_2[K + 1 + (1:G)]
        alpha_vcov <- vcov_step_2[1:(K + 1), 1:(K + 1)]
        
        # Newey_3. probit/tobit with the covariates and the residuals of the first step
        if (model == "probit") newey_3 <- glm(y ~ fitted(step_1) + X1 + resid(step_1) - 1, family = binomial(link = 'probit'))
        if (model == "tobit") newey_3 <- tobit1(y ~ fitted(step_1) + X1 + resid(step_1) - 1, left = left, right = right)
        coef_newey_3 <- newey_3$coefficients
#        coef_fit_w <- coef_newey_3[1 + (1:G)] bug 2024-01-03: mauvais indexage
        coef_fit_w <- coef_newey_3[1:G]
        rho_hat <- coef_res_w - coef_fit_w

        # Newey_4. computation of the augmented matrix of covariance
        W_rho <- drop(W %*% rho_hat)
        newey_4 <- lm(W_rho ~ X)
        SIGMA <- vcov(newey_4)
        OMEGA <- SIGMA + alpha_vcov
        
        # Newey_5. Computation of the estimator and its variance
        d <- solve(crossprod(cbind(1, X)), crossprod(cbind(1, X), cbind(1, ZZ)), tol = 1E-20)
        .vcov <- solve(t(d) %*% solve(OMEGA, tol = 1E-20) %*% d, tol = 1E-20)
        .coef <- drop(.vcov %*% t(d) %*% solve(OMEGA, alpha_hat, tol = 1E-20))
        names(.coef)[1] <- colnames(.vcov)[1] <- rownames(.vcov)[1] <- "(Intercept)"
    }
    if (.est_method == "ml"){
        nms_param <- list(covar = names_X, ancil = names_ZZ, vcov = c("sigma", "rho"))
        # computation of the starting values
        nu <- step_2$linear.predictor
        # Compute the covariance matrix of the SUR estimator and takes its Cholesky decomposition
        SIGMA <- crossprod(Wres) / N
        CHOL <- t(chol(solve(SIGMA)))
        chol <- CHOL[lower.tri(CHOL, diag = TRUE)]
        za <- outer(colnames(W), colnames(W), paste, sep= "|")
        nms_chol <- za[lower.tri(za, diag = TRUE)]
#        nms_chol <- paste("chol", 1:length(chol), sep = "_")
        # Other coefficients
        PI2 <- t(coef(step_1))
        beta <- coef(step_2)[1:L]
        rho <- coef(step_2)[L + (1:G)]
        nms_beta <- paste("covar", c(colnames(X1), colnames(W)), sep = "_")
        nms_beta <- c(colnames(X1), colnames(W))
        nms_rho <- paste("rho", colnames(W), sep = "_")
        # vector of estimated covariances
        covs <- drop(crossprod(Wres, nu)) / N
        if (model == "tobit") .sigma <- coef(step_2)["sigma"] else .sigma <- 1
        d <- sqrt(.sigma ^ 2 - drop((t(covs) %*% solve(SIGMA) %*% covs)))
        # transform the starting coefficients of the probit / tobit estimator
        beta <- beta * d
        rho <- drop(t(CHOL) %*% covs)
        PI2 <- drop(t(CHOL) %*% PI2)
        pi2 <- as.numeric(t(PI2))
        if (ncol(W) > 1){
            colnames(PI2) <- names_X
            nms_pi2 <- as.character(outer(colnames(PI2), rownames(PI2), paste, sep = "_"))
            nms_pi3 <- as.character(t(outer(rownames(PI2), colnames(PI2), paste, sep = "_")))
            nms_pi2 <- paste("ancil", nms_pi3, sep = "_")
        }
        else nms_pi2 <- paste(colnames(W), names(PI2), sep = "_")
        .start <- c(beta, rho, pi2, chol)
        nms_coef <- c(nms_beta, nms_rho, nms_pi2, nms_chol)
        if (model == "tobit"){
            .start <- c(.start, .sigma)
            nms_coef <- c(nms_coef, "sigma")
        }
        names(.start) <- nms_coef

        ## check_gradient <- function(f, g, param){
        ##     comp_grad <- cbind(analytic = gr(param), numeric = numDeriv::grad(func, param), diff = gr(param) - numDeriv::grad(func, param))
        ##     rownames(comp_grad) <- names(param)
        ##     comp_grad
        ## }
            
        func_obs <- function(param)    lnliv_ldv(param, X1 = X1, X2 = X2, W = W, y = y, sum = FALSE, gradient = FALSE, right = right, model = model)
        func <- function(param)    -      lnliv_ldv(param, X1 = X1, X2 = X2, W = W, y = y, sum = TRUE, gradient = FALSE, right = right, model = model)
        gr <- function(param)      - attr(lnliv_ldv(param, X1 = X1, X2 = X2, W = W, y = y, sum = TRUE, gradient = TRUE, right = right, model = model), "gradient")
        grObs <- function(param)    attr(lnliv_ldv(param, X1 = X1, X2 = X2, W = W, y = y, sum = FALSE, gradient = TRUE, right = right, model = model), "gradient")

        trace <- 0
        .optim <- optim(.start, func, gr, method = "BFGS",
                        hessian = TRUE, control = list(trace = trace, maxit = 1E04))
        .logLik <- structure(- .optim$value, nobs = length(y), df = length(.optim$par),
                             class = "logLik")
        .llcount <- func_obs(.optim$par)
        .coef <- .optim$par
        .gradient <- grObs(.coef)
    }
    .terms <-  terms(mf)
    attr(.terms, ".Environment") <- NULL
    result <- list(coefficients = .coef,
#                   linear.predictor = linear.predictor,
#                   fitted.values = .fitted,
#                   residuals = y - .fitted,
                   df.residual = N - length(.coef),
                   model = mf,
                   terms = .terms,
                   call = .call,
                   xlevels = .getXlevels(mt, mf),
                   na.action = attr(mf, "na.action"),
                   est_method = method)

    if (.est_method == "ml"){
        result$hessian <- - .optim$hessian
        result$gradient <- .gradient
        result$value <- .logLik
        .npar <- c(covariates = K1 + G + has.int,
                   resid = G,
                   instruments = (K1 + K2 + has.int) * G,
                   chol = G * (G + 1) / 2)
        attr(.npar, "default") <- c("covariates", "resid")
        if (model == "tobit"){
            .npar <- c(.npar, vcov = 1)
            attr(.npar, "default") <- c("covariates", "resid", "vcov")
        }
    }                    
    if (.est_method != "ml"){
        result$vcov <- .vcov
    }
    if (.est_method == "twosteps"){
        .npar <- structure(c(covariates = K1 + G + has.int, resid = G),
                           default =  c("covariates", "resid"))
        if (model == "tobit"){
            .npar <- c(.npar, vcov = 1)
            attr(.npar, "default") <- c("covariates", "resid", "vcov")
        }
    }
        
    if (.est_method == "minchisq"){
        .npar <- c(covariates = K1 + G + has.int)
    }
    result$npar <- .npar
    structure(result, class = c("ivldv", "micsr"))
}
#' @rdname ivldv
#' @export
endogtest <- function(x, ...) UseMethod("endogtest")

#' @rdname ivldv
#' @export
endogtest.formula <- function(x, ..., data, model = c("probit", "tobit")){
    model <- match.arg(model)
    cl <- match.call()
    cl[[1]] <- as.name("ivldv")
    names(cl)[2] <- "formula"
    cl$compute_test <- TRUE
    eval(cl, parent.frame())
}

#' @rdname ivldv
#' @export
endogtest.ivldv <- function(x, ...){
    update(x, compute_test = TRUE)
}


