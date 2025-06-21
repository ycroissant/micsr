#' Hausman test
#'
#' Hausman test; under the null both models are consistent but one of
#' them is more efficient, under the alternative, only one model is
#' consistent
#'
#' @name hausman
#' @param x the first model,
#' @param y the second model
#' @param omit a character containing the effects that are removed from the test
#' @param ... further arguments
#' @return an object of class `"htest"`.
#' @keywords htest
#' @author Yves Croissant
#' @importFrom stats pchisq
#' @references
#' \insertRef{HAUS:78}{micsr}
hausman <- function(x, y, omit = FALSE, ...)
    UseMethod("hausman")

#' @rdname hausman
#' @export
hausman.ivreg <- function(x, y, omit = FALSE, ...){
    .stat <- summary(x)$diagnostics["Wu-Hausman", ]
    .parameter <- .stat[1]
    .statistic <- .stat[3]
    names(.statistic) <- "chisq"
    .pval <- .stat[4]
    .method <- "Hausman Test"
    .data.name <- paste(deparse(formula(x)))
    rval <- list(statistic = .statistic,
                 parameter = .parameter,
                 p.value = .pval,
                 method = .method,
                 data.name = .data.name)
    structure(rval, class = "htest")
}

#' @rdname hausman
#' @export
hausman.micsr <- function(x, y, omit = NULL, ...){
    .data.name <- paste(paste(deparse(substitute(x))), "vs",
                       paste(deparse(substitute(y))))
    nms_x <- names(coef(x))
    nms_y <- names(coef(y))
    nms <- intersect(nms_x, nms_y)
    unkn <- setdiff(omit, nms)
    if (length(unkn)){
        effects_string <- ifelse(length(unkn) == 1, "effect", "effects")
        stop(cat(paste(paste(unkn, collapse = ", "), effects_string, "unknown\n"), sep = ""))
    }
    nms <- setdiff(nms, omit)
    delta <- coef(x)[nms] - coef(y)[nms]
    V <- vcov(x)[nms, nms] - vcov(y)[nms, nms]
    stat <- as.numeric(crossprod(solve(V, delta), delta))
    pval <- pchisq(stat, df = length(delta), lower.tail = FALSE)
    res <- list(statistic    = c(chisq = stat),
                p.value      = pval,
                parameter    = c(df = length(delta)),
                method       = "Hausman Test",
                data.name    = .data.name,
 #             null.value  = null.value,
                alternative  = "the second model is inconsistent")
    class(res) <- "htest"
    return(res)
}

#' F statistic
#'
#' Extract the F statistic that all the parameters except the
#' intercept are zero. Currently implemented only for models fitted by `lm` or `ivreg::ivreg`.
#' @name ftest
#' @param x a fitted object
#' @param covariate the covariate for which the test should be performed for the `ivreg` method
#' @param ... further arguments
#' @return an object of class `"htest"`.
#' @importFrom stats pf
#' @keywords htest
#' @export
ftest <- function(x, ...){
    UseMethod("ftest")
}

#' @rdname ftest
#' @export
ftest.lm <- function(x, ...){
    .fstat <- summary(x)$fstatistic
    .statistic <- .fstat["value"]
    names(.statistic) <- "F"
    .parameter <- .fstat[2:3]
    .method <- "F test"
    .data.name <- paste(deparse(formula(x)))
    .pval <- pf(.statistic, .parameter[1], .parameter[2], lower.tail = FALSE)
    rval <- list(statistic = .statistic,
                 parameter = .parameter,
                 p.value = .pval,
                 method = .method,
                 data.name = .data.name)
    structure(rval, class = "htest")
}
                 

#' @rdname ftest
#' @export
ftest.ivreg <- function(x, ..., covariate = NULL){
    .fstat <- summary(x)$diagnostics
    wi_rows <- grep("Weak instruments", rownames(.fstat))
    .fstat <- .fstat[wi_rows, , drop = FALSE]
    if (! is.null(covariate)) .covariate <- grep(covariate, rownames(.fstat))
    else .covariate <- 1
    if (length(.covariate) > 1) stop("several covariates")
    if (nrow(.fstat) > 1) .fstat <- .fstat[.covariate, ]
    .statistic <- .fstat["statistic"]
    .parameter <- .fstat[1:2]
    .method <- "F test"
    names(.statistic) <- "F"
    .data.name <- paste(deparse(formula(x)))
    .pval <- pf(.statistic, .parameter[1], .parameter[2], lower.tail = FALSE)
    rval <- list(statistic = .statistic,
                 parameter = .parameter,
                 p.value = .pval,
                 method = .method,
                 data.name = .data.name)
    structure(rval, class = "htest")
}    
    
#' Score test
#'
#' Score test, also knowned as Lagrange multiplier tests
#'
#' @name scoretest
#' @param object the first model,
#' @param ... for the `micsr` method, it should be the formula for the
#'     "large" model or an object from which a formula can be
#'     extracted
#' @param vcov an optional covariance matrix
#' @return an object of class `"htest"`.
#' @keywords htest
#' @author Yves Croissant
#' @importFrom stats pchisq
#' @examples
#' mode_choice <- transform(mode_choice, cost = cost * 8.42)
#' mode_choice <- transform(mode_choice, gcost = (ivtime + ovtime) * 8 + cost)
#' pbt_unconst <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = "probit")
#' pbt_const <- binomreg(mode ~ gcost, data = mode_choice, link = "logit")
#' scoretest(pbt_const , . ~ . + ivtime + ovtime)
#' @export
scoretest <- function(object, ...){
    UseMethod("scoretest")
}

#' @rdname scoretest
#' @export
scoretest.default <- function(object, ...){
    new <- list(...)[[1]]
    cls <- class(object)[1]
    nmodels <- length(new)
    if (! inherits(new, 'formula') | ! inherits(new, cls))
        stop("the updating argument doesn't have a correct class")
    if (inherits(new, cls)){
        ncoefs <- names(coef(new))
        new <- formula(formula(new))
    }
    else ncoefs <- names(coef(update(object, new, iterlim = 0)))
    start <- numeric(length = length(ncoefs))
    names(start) <- ncoefs
    supcoef <- ! ncoefs %in% names(coef(object))
    start[names(coef(object))] <- coef(object)
    newmodel <- update(object, new, start= start, iterlim = 0)
    data.name <- paste(deparse(formula(newmodel)))
    alt.hyp <- "unconstrained model"
    if (is.matrix(newmodel$gradient))
        gradvect <- apply(newmodel$gradient, 2, sum) else gradvect <- newmodel$gradient
    stat <- - sum(gradvect * solve(newmodel$hessian, gradvect))
    names(stat) <- "chisq"
    df <- c(df = length(coef(newmodel)) - length(coef(object)))
    pval <- pchisq(stat, df = df, lower.tail = FALSE)
    result <- list(statistic = stat,
                   parameter = df,
                   p.value = pval,
                   data.name = data.name,
                   method = "score test",
                   alternative = alt.hyp
                   )
    class(result) <- 'htest'
    result
}

#' @rdname scoretest
#' @export
scoretest.micsr <- function(object, ..., vcov = NULL){
    x <- object
    if (is.null(vcov)){
        if (is.null(x$info)) .vcov <- "hessian" else .vcov <- "info"
    }
    else .vcov <- vcov
    objects <- list(object, ...)
    if (length(objects) != 2) stop("Two models should be provided")
    if(! inherits(objects[[2]], "formula")) objects[[2]] <- formula(objects[[2]])
    newform <- objects[[2]]
    nX <- model.matrix(x, newform)
    new_nms <- colnames(nX)
    old_nms <- names(coef(x))
    L <- length(new_nms) - length(old_nms)
    if(! all(old_nms %in% new_nms))
        stop("the old model is not nested in the new one")
    start <- coef(x)[new_nms]
    start[is.na(start)] <- 0
    names(start) <- new_nms
    y <- update(x, formula = newform, start = start, maxit = 0)
    .gradient <- apply(y$gradient, 2, sum)
    if (.vcov == "hessian") information <- - y$hessian
    if (.vcov == "opg") information <- crossprod(y$gradient)
    if (.vcov == "info"){
        if (is.null(y$info)) stop("no information matrix estimate available")
        else information <- y$info
    }
    .stat <- drop(crossprod(.gradient, solve(information, .gradient)))
    structure(list(statistic = c(chisq = .stat),
                   p.value = pchisq(.stat, lower.tail = FALSE, df = L),
                   parameter = c(df = L),
                   method = "score test",
                   data.name = paste(deparse(newform)),
                   alternative = "the constrained model is rejected"),
              class = "htest")    
}

#' Sargan test for GMM models
#'
#' When a IV model is over-identified, the set of all the empirical
#' moment conditions can't be exactly 0. The test of the validity of
#' the instruments is based on a quadratic form of the vector of the
#' empirical moments
#' 
#' @name sargan
#' @param object a model fitted by GMM
#' @param ... further arguments
#' @return an object of class `"htest"`.
#' @keywords htest
#' @examples
#' cigmales <- cigmales |>
#'        transform(age2 = age ^ 2, educ2 = educ ^ 2,
#'                  age3 = age ^ 3, educ3 = educ ^ 3,
#'                  educage = educ * age)
#' gmm_cig <- expreg(cigarettes ~ habit + price + restaurant + income + age + age2 +
#'                  educ + educ2 + famsize + race | . - habit + age3 + educ3 +
#'                  educage + lagprice + reslgth, data = cigmales,
#'                  twosteps = FALSE)
#' sargan(gmm_cig)
#' @export
sargan <- function(object, ...){
    UseMethod("sargan")
}

#' @rdname sargan
#' @export
sargan.ivreg <- function(object, ...){
    .fstat <- summary(object)$diagnostics["Sargan", ]
    .statistic <- .fstat["statistic"]
    .parameter <- .fstat[1]
    .method <- "Sargan test"
    names(.statistic) <- "chisq"
    .data.name <- paste(deparse(formula(object)))
    .pval <- pchisq(.statistic, .parameter, lower.tail = FALSE)
    rval <- list(statistic = .statistic,
                 parameter = .parameter,
                 p.value = .pval,
                 method = .method,
                 data.name = .data.name)
    structure(rval, class = "htest")
}

#' @rdname sargan
#' @export
sargan.micsr <- function(object, ...){
    if (! object$est_method %in% c("iv", "gmm")) stop("Sargan test only relevant for GMM models")
    .stat <- nobs(object) * object$value
    .df <- object$L - object$K
    .pval <- pchisq(.stat, .df, lower.tail = FALSE)
    .data <- deparse(object$call$data)
    structure(list(statistic = c(chisq = .stat),
                   parameter = c(df = .df),
                   method = "Sargan Test",
                   p.value = .pval,
                   data.name = .data,
                   alternative = "the moment conditions are not valid"),
              class = "htest")
}

