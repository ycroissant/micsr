#' @importFrom stats sd

check_gradient <- function(f, coefs){
    object <- f(coefs)
    anal_grad <- attr(object, "gradient")
    anal_hess <- attr(object, "hessian")
    num_grad <- numDeriv::grad(f, coefs)
    num_hess <- numDeriv::hessian(f, coefs)
    dimnames(num_hess) <- dimnames(anal_hess)
    max_grad <- max(abs(anal_grad - num_grad))
    if (any(is.na(num_hess))) message("some NA in the numerical hessian")
    max_hess <- max(abs(anal_hess[! is.na(num_hess)] - num_hess[! is.na(num_hess)]) / sd(anal_hess))
    structure(list(gradient = max_grad, hessian = max_hess), class = "check_gradient",
              gradient = cbind(num_grad, anal_grad), anal_hess = anal_hess, num_hess = num_hess)
}

#' Compute the inverse Mills ratio and its first two derivatives
#'
#' The inverse Mills ratio is used in several econometric models,
#' especially different flavours of tobit model.
#' @name mills
#' @param x a numeric
#' @param deriv one of 0 (the default, returns the inverse Mills
#'     ratio), 1 (the first derivative) and 2 (the second derivative)
#' @keywords misc
#' @return a numeric.
#' @export
mills <- function(x, deriv = 0){
    if (! deriv %in% 0:2) stop("irrelevant value of deriv")
    .mills <-  exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
    if (deriv == 1) .mills <- - .mills * (x + .mills)
    if (deriv == 2) .mills <- .mills * ( (x + .mills) * (x + 2 * .mills) - 1)
    .mills
}

dmills <- function(x) - mills(x) * (x + mills(x))
d2mills <- function(x) mills(x) * ( (x + mills(x)) * (x + 2 * mills(x)) - 1)

#' Newton-Raphson  method for numerical optimization
#'
#' The Newton-Raphson method use the gradient and the hessian of a
#' function. For well behaved functions, it is extremely accurate.
#' @name newton
#' @param fun the function to optimize
#' @param coefs a vector of starting values
#' @param trace if positive or true, some information about the
#'     computation is printed
#' @param direction either `"min"` or `"max"`
#' @param tol the tolerance
#' @param maxit maximum number of iterations
#' @param ... further arguments, passed to fun
#' @keywords misc
#' @return a numeric vector, the parameters at the optimum of the
#'     function.
#' @export
newton <- function(fun, coefs, trace = 0, direction = c("min", "max"), tol = sqrt(.Machine$double.eps), maxit = 500, ...){
    if (maxit == 0) return(coefs)
    if (trace){
        cat("Initial values of the coefficients:\n")
    }
    direction <- match.arg(direction)
    i <- 1
    eps <- 10
    while (abs(eps) > tol){
        f <- fun(coefs, gradient = TRUE, hessian = TRUE, ...)
        g <- attr(f, "gradient")
        if (is.matrix(g)) g <- apply(g, 2, sum)
        h <- attr(f, "hessian")
        if (direction == "max"){
            f <- - f
            g <- - g
            h <- - h
        }
        lambda <- 1
        newcoefs <- coefs - as.numeric(solve(h, g))
        as_scalar <- function(x) sum(as.numeric(x))
        while (as_scalar(- fun(newcoefs, ...)) > as_scalar(f)){
            lambda <- lambda / 2
            if(trace) cat(paste("function is increasing, lambda set to:", lambda, "\n"))
            newcoefs <- coefs - lambda * as.numeric(solve(h, g))
        }
        eps <- as.numeric(crossprod(solve(h, g), g))
                if (trace) cat(paste("iteration:", i, "criteria:", round(eps, 5), "\n"))
        i <- i + 1
        if (i > maxit) stop("max iter reached")
        coefs <- newcoefs
    }
    coefs
}

#' Extract the standard errors of estimated coefficients
#'
#' The standard errors are a key element while presenting the results
#' of a model. They are the second column of the table of coefficient
#' and are used to compute the t/z-value. `stder` enables to retrieve
#' easily the vector of standard errors, either from a fitted model or
#' from a matrix of covariance
#' 
#' @name stder
#' @param x a fitted model or a matrix of covariance
#' @param vcov a function that computes a covariance matrix, or a character
#' @param subset,grep,fixed,invert invert see `micsr::select_coef
#' @param ... further arguments
#' @return a numeric vector
#' @keywords misc
#' @export
stder <- function(x, vcov, subset = NA, fixed = FALSE, grep = NULL, invert = FALSE, ...) UseMethod("stder")

#' @rdname stder
#' @export
stder.default <- function(x, vcov = NULL, subset = NA, fixed = FALSE, grep = NULL, invert = FALSE, ...){
    .vcov <- vcov
    if (is.matrix(x)) std <- sqrt(diag(x))
    else{
        if (! is.null(.vcov)){
            if (is.character(.vcov)){
                if (! inherits(x, "micsr"))
                    stop("object should be of class micsr")
                std <- sqrt(diag(vcov(x, vcov = .vcov, subset = subset, fixed = fixed, grep = grep, invert = invert)))
            }
            if (is.function(.vcov)){
                std <- sqrt(diag(.vcov(x, ...)))
            }
        }
        else  std <- sqrt(diag(vcov(x)))
    }
    std
}

## #' Transform a factor in a set of dummy variables
## #'
## #' The normal way to store cathegorical variables in R is to use
## #' factors, each modality being a level of this factor. Sometimes
## #' however, is is more convenient to use a set of dummy variables.
## #'
## #' @name dummy
## #' @param x a data frame
## #' @param ...  series of the data frame, should be factors
## #' @param keep a boolean, if `TRUE`, the original series is kept in
## #'     the data frame,
## #' @param prefix an optional prefix for the names of the computed
## #'     dummies,
## #' @param ref a boolean, if `TRUE`, a dummy is created for all the
## #'     levels, including the reference level
## #' @return a data frame
## #' @importFrom dplyr pull
## #' @importFrom tidyselect eval_select
## #' @importFrom rlang expr
## #' @keywords misc
## #' @examples
## #' charitable %>% dummy(religion, education)
## #' @export
## dummy <- function (x, ..., keep = FALSE, prefix = NULL, ref = FALSE) {
## #    error_call <- dplyr:::dplyr_error_call()
##     loc <- eval_select(expr(c(...)), data = x)
##     for (y in names(loc)){
##         if (inherits(x[[y]], "factor")){
##             levs <- x %>% pull({{ y }}) %>% levels
##             if (! ref) levs <- levs[- 1]
##             for (i in rev(levs)){
##                 if (! is.null(prefix)) nmi <- paste(prefix, i, sep = "") else nmi <- i
##                 x <- x %>% mutate({{ nmi }} := ifelse({{ y }} == i, 1, 0), .after = {{ y }})
##             }
##             if (! keep) x <- x %>% dplyr::select(- {{ y }})
##             x
##         }
##     }
##     x
## }

#' Number of parameters of a fitted model
#'
#' The number of observation of a fitted model is typically obtained
#' using the `nobs` method. There is no such generics to extract the
#' same information about the number of parameters. `npar` is such a
#' generic and has a special method for `micsr` objects with a
#' `subset` argument that enables to compute the number of parameters
#' for a subset of coefficients. The default method returns the length
#' of the vector of coefficients extracted using the `coef` function.
#'
#' @name npar
#' @param x a fitted model
#' @param subset a character indicating the subset of coefficients
#'     (only relevant for `micsr` models).
#' @return an integer.
#' @keywords misc
#' @author Yves Croissant
#' @export
npar <- function(x, subset = NULL)
    UseMethod("npar")

#' @rdname npar
#' @export
npar.default <- function(x, subset = NULL){
    length(coef(x))
}

#' @rdname npar
#' @export
npar.micsr <- function(x, subset = NULL){
    result <- x$npar
    if (! is.null(subset)){
        if (! is.character(subset)) stop("subset should be a character")
        if (any(! subset %in% names(result))) stop("unknown subset")
        result <- result[subset]
    }
    sum(as.numeric(result))
}

compute_rank <- function(x){
    abs_eigen_value <- abs(eigen(crossprod(x), only.values = TRUE)$values)
    sum(abs_eigen_value > sqrt(.Machine$double.eps))
}


maximize <- function(x, start, method = c("bfgs", "nr", "newton"), trace = 0, maxit = 100, ...){
    .method <- match.arg(method)
    if (method == "bfgs"){
        fun <- function(param) x(param, gradient = FALSE, hessian = FALSE, opposite = TRUE, ...)
        grad <- function(param) attr(x(param, gradient = TRUE, hessian = FALSE, opposite = TRUE, ...), "gradient")
        result <- optim(start, fun, grad, method = "BFGS", control = list(trace = trace, maxit = maxit))$par
    }
    if (method == "nr"){
        result <- nlm(x, start, print.level = trace, iterlim = maxit,
                      gradient = TRUE, hessian = TRUE, sum = TRUE, opposite = TRUE, ...)$estimate
    }
    if (method == "newton"){
        result <- newton(x, coefs = start, direction = "max", maxit = maxit, trace = trace, ...)
    }
    names(result) <- names(start)
    result
}

#' Compute quadratic form
#'
#' Compute quadratic form of a vector with a matrix, which can be the
#' vector of coefficients and the covariance matrix extracted from a
#' fitted model
#'
#' @name quad_form
#' @param x a numeric vector or a fitted model
#' @param m a square numeric matrix
#' @param inv a boolean, if `TRUE` (the default), the quadratic form
#'     is computed using the inverse of the matrix
#' @param subset a subset of the vector and the corresponding subset
#'     of the matrix
#' @param vcov if `NULL` the `vcov` method is used, otherwise it can
#'     be a function or, for `micsr` objects, a character
#' @param \dots arguments passed to `vcov` if it is a function
#' @export
quad_form <- function(x, m = NULL, inv = TRUE, subset = NULL, vcov = NULL, ...){
    .sub <- subset
    if (is.list(x)){
        if (is.null(coef(x))) stop("x should be a fitted model")
        .vcov <- vcov
        .coef <- coef(x)
        if (! is.null(.vcov)){
            if (is.character(.vcov)){
                if (! inherits(x, "micsr"))
                    stop("object should be of class micsr")
                .vcov <- vcov(x, vcov = .vcov)
            }
            if (is.function(.vcov)) .vcov <- .vcov(x, ...)
        }
        else .vcov <- vcov(x)
        if (is.null(.sub)) .sub <- 1:length(.coef)
        x <- .coef[.sub]
        m <- .vcov[.sub, .sub]
    }
    if (is.null(m)) stop("a matrix should be provided")
    if (length(m) == 1) m <- matrix(m, 1, 1)
    if (! is.matrix(m)  | ! is.numeric(m)) stop("the m argument should be a numeric matrix")
    if (nrow(m) != ncol(m)) stop("the matrix should be square")
    if (! is.numeric(x)) stop("the first argument should be numeric")
    if (is.matrix(x)){
        if (! (nrow(x) == 1 | ncol(x) == 1))
            stop("the first argument shouldn't be a matrix")
        else x <- drop(x)
    }
    if (length(x) != nrow(m)) stop("the length of the vector should be equal to the dimensions of the matrix")
    .sub <- subset
    if (is.null(.sub)) .sub <- 1:length(x)
    x <- x[.sub]
    m <- m[.sub, .sub]
    if (inv) qf <- as.numeric(crossprod(x, solve(m, x)))
    else qf <- as.numeric(crossprod(x, crossprod(m, x)))
    qf
}

#' select a subset of coefficients
#'
#' `micsr` objects have a `rpar` element which is vector of integers
#' with names that indicates the kind of the coefficients. For
#' example, if the 6 first coefficients are covariates parameters and
#' the next 3 parameters that define the distribution of the errors,
#' `npar` will be c(covariates = 6, vcov = 3). It has an attribute
#' which indicates the subset of coefficients that should be selected
#' by default. `select_coef` has a `subset` argument (a character
#' vector) and returns a vector of integers which is the position of
#' the coefficients to extract.
#'
#' @name select_coef
#' @param object a fitted model
#' @param subset a character vector, the type of parameters to extract
#' @param fixed if `TRUE`, the fixed parameters are selected
#' @param grep a regular expression
#' @param invert should the coefficients that **don't** match the
#'     pattern should be selected ?
#' @param coef a vector of coefficients
#' @return a numeric vector
#' @export
select_coef <- function(object, subset = NA, fixed = FALSE,
                        grep = NULL, invert = FALSE, coef = NULL){
    # ajouter subset = NULL de manière à ce que tous les coefficients soient sélectionnés
    .grep <- grep
    .coef <- coef
    .npar <- object$npar
    .names <- names(object$coefficients)
    .fixed <- attr(object$coefficients, "fixed")
    if (is.null(.fixed)) .fixed <- rep(FALSE, sum(.npar))
    if (is.null(.npar)){
        .npar <- structure(c(covariates = length(object$coefficients)),
                           default = "covariates")
    }

    # subset
    if (is.null(subset)) .subset <- names(.npar)
    if (length(subset) == 1){
        if (is.na(subset)) .subset = attr(.npar, "default")
        else{
            if (subset == "all") .subset <- names(.npar)
            else{
                if (subset %in% names(.npar)) .subset <- subset
                else stop("irrelevant subset argument")
            }
        }
    } else {
        .subset <- subset
        if (! all(.subset %in% names(.npar))) stop("irrelevant subset")
    }
    idx <- data.frame(subset = rep(names(.npar), times = .npar),
                      idx = 1:sum(.npar),
                      fixed = .fixed)
    .sel <- .subset
    idx <- subset(idx, subset %in% .sel)
    if (! fixed) idx <- subset(idx, ! fixed)
    idx <- idx$idx
    names(idx) <- .names[idx]

    # coef
    if (! is.null(.coef) | ! is.null(.grep)){
        z1 <- z2 <- numeric(0)
        if (! is.null(.coef)){
            if (! any(.coef %in% names(idx))) stop("unknown coefficient")
            z1 <- match(.coef, names(idx))
        }
        if (! is.null(.grep)){
            z2 <- grep(.grep, names(idx), invert = invert)
        }
        z <- unique(c(z1, z2))
        idx <- idx[z]
    }
    idx
}


# coef:          un vecteur de coefficients
# subset:        un vecteur de noms de sous-ensembles de coefficients
# grep / invert: une expression rationelle et un booleun pour inverser

