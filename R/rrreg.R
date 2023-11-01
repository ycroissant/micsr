rr <- function(x, model = c("probit", "tobit")){
    model <- match.arg(model)
    delta <- qnorm(mean(x))
    if (model == "probit"){
        psi_1 <- dnorm(delta)
        psi_2 <- mean(x) - dnorm(delta) * delta
    }
    if (model == "tobit"){
        hat_Phi <- mean(x > 0)
        psi_1 <- hat_Phi
        hat_phi <- dnorm(qnorm(hat_Phi))
        s2 <- mean((x - mean(x)) ^ 2) /
            (hat_Phi - (hat_phi - qnorm(hat_Phi) * (1 - hat_Phi)) * (hat_phi + qnorm(hat_Phi) * hat_Phi))
        psi_2 <- sqrt(s2) * hat_phi
    }
    (x - psi_2) / psi_1
}

#' Recenter and rescale regression
#'
#' A method of estimation for probit and tobit models which consist on
#' recentering and rescaling the response so that linear estimators can then be applied
#'
#' @name rrreg
#' @param formula a symbolic description of the model,
#' @param data a data frame,
#' @param subset,weights,na.action,offset see `lm`,
#' @param method the estimation method, one of `"ols"`, `"iv"` and `"gmm"`
#' @param model the model to estimate, one of `"probit"` and `"tobit"`
#' @param ... further arguments
#' @return an object of class `micsr`
#' @importFrom stats qnorm
#' @references \insertRef{IWAT:01}{micsr}
#' @examples
#' library("dplyr")
#' data("PSID1976", package = "AER")
#' Mroz <- PSID1976 %>% as_tibble %>%
#'       mutate(participation = as.numeric(participation == "yes"),
#'              nwincome = (fincome - hours * wage) / 1000, exper2 = I(experience ^ 2))
#' rrols <- rrreg(participation ~ age + education + youngkids + oldkids | .
#'                - education + meducation + feducation, data = Mroz)
#' rriv <- update(rrols, method = "iv")
#' rrgmm <- update(rrols, method = "gmm")
#' # Adkins (2012) : Bank
#' rrols <- rrreg(federiv ~ eqrat + optval + bonus + ltass + linsown + linstown +
#'                    roe + mktbk + perfor + dealdum + div + year |
#'                    . - eqrat - bonus - optval +
#'                    no_emp + no_subs + no_off + ceo_age + gap + cfa,
#'                data = federiv, method = "ols")
#' rriv <- update(rrols, method = "iv")
#' rrgmm <- update(rriv, method = "gmm")
#' @export

rrreg <- function(formula,
                  data,
                  subset = NULL,
                  weights = NULL,
                  na.action,
                  offset,
                  method = c("ols", "iv", "gmm"),
                  model = c("probit", "tobit"),
                  ...){
    .model <- match.arg(model)
    .est_method <- match.arg(method)
    .call <- match.call()
    cl <- match.call(expand.dots = TRUE)
    cl$formula <- .formula <- Formula(formula)
    .method <- match.arg(method)
    m <- match(c("formula", "data", "subset", "weights"),
               names(cl), 0L)
    # construct the model frame and components
    cl <- cl[c(1L, m)]
    mf <- cl
    mf[[1L]] <- as.name("model.frame")
    mf$dot <- "previous"
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    # response and its transformation
    y <- rr(model.response(mf), model = .model)
    N <- length(y)
    X <- model.matrix(.formula, mf, rhs = 1)
    if (.est_method == "ols"){
        ols <- lm(y ~ X - 1)
        .vcov <- vcov(ols)
        .coef <- coef(ols)
    }
    else{
        Z <- model.matrix(.formula, mf, rhs = 2, dot = "previous")
        Xh <- fitted(lm(X ~ Z - 1))
        second <- lm(y ~ Xh - 1)
        if (.est_method == "iv"){
            .coef <- coef(second)
            .vcov <- vcov(second)
        }
        if (.est_method == "gmm"){
            e <- y - drop(X %*% coef(second))
            ZeZ <- crossprod(sqrt(e ^ 2) * Z)
            ZZ <- crossprod(Z)
            XZ <- crossprod(X, Z)
            Zy <- crossprod(Z, y)
            .vcov <- solve(XZ %*% solve(ZeZ) %*% t(XZ))
            .coef <- drop(.vcov %*% XZ %*% solve(ZeZ) %*% Zy)
            e <- y - drop(X %*% .coef)
            Ze <- drop(crossprod(Z, e))
            .value <- drop(t(Ze) %*% solve(ZeZ) %*% Ze)
        }
    }
    colnames(.vcov) <- rownames(.vcov) <- names(.coef) <- colnames(X)
    result <- list(coefficients = .coef,
                   vcov = .vcov,
                   call = .call,
                   est_method = .est_method,
                   model = mf,
                   K = ncol(X))
    if (.est_method == "gmm") result$value <- .value
    if (.est_method != "ols") result$L <- ncol(Z)
    structure(result, class = "micsr")
}
    
