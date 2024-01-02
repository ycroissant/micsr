#' Compute the inverse Mills ratio and its first derivatives
#'
#' The inverse Mills ratio is used in several econometric models,
#' especially different flavours of tobit model.
#' @name mills
#' @param x a numeric
#' @param deriv one of 0 (the default, returns the inverse Mills
#'     ratio), 1 (the first derivative) and 2 (the second derivative)
#' @return a numeric
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

#' Newton-Rawlphson  method for numerical optimization
#'
#' The Newton-Rawlphson method use the gradient and the hessian of a
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
#' @return a numeric vector, the parameters at the optimum of the
#'     function
#' @export
newton <- function(fun, coefs, trace = 0, direction = c("min", "max"), tol = sqrt(.Machine$double.eps), maxit = 500, ...){
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
#' and are used to compute the t/z-value. `stderr` enables to retrieve
#' easily the vector of standard errors, either from a fitted model or
#' from a matrix of covariance
#' 
#' @name stder
#' @param x a fitted model or a matrix of covariance
#' @param .vcov a function that computes a covariance matrix, or a character
#' @param ... further arguments
#' @return a numeric vector
#' @export
stder <- function(x, .vcov, ...) UseMethod("stder")

#' @rdname stder
#' @export
stder.default <- function(x, .vcov = NULL, ...){
    if (is.matrix(x)) std <- sqrt(diag(x))
    else{
        if (! is.null(.vcov)){
            if (is.character(.vcov)){
                if (! inherits(x, "micsr"))
                    stop("object should be of class micsr")
                std <- sqrt(diag(vcov(x, vcov = .vcov)))
            }
            if (is.function(.vcov)){
                std <- sqrt(diag(.vcov(x, ...)))
            }
        }
        else  std <- sqrt(diag(vcov(x)))
    }
    std
}

#' Transform a factor in a set of dummy variables
#'
#' The normal way to store cathegorical variables in R is to use
#' factors, each modality being a level of this factor. Sometimes
#' however, is is more convenient to use a set of dummy variables.
#'
#' @name dummy
#' @param x a data frame
#' @param ...  series of the data frame, should be factors
#' @param keep a boolean, if `TRUE`, the original series is kept in
#'     the data frame,
#' @param prefix an optional prefix for the names of the computed
#'     dummies,
#' @param ref a boolean, if `TRUE`, a dummy is created for all the
#'     levels, including the reference level
#' @return a data frame
#' @importFrom tidyselect eval_select
#' @importFrom rlang expr
#' @examples
#' charitable %>% dummy(religion, education)
#' @export
dummy <- function (x, ..., keep = FALSE, prefix = NULL, ref = FALSE) {
#    error_call <- dplyr:::dplyr_error_call()
    loc <- eval_select(expr(c(...)), data = x)
    for (y in names(loc)){
        if (inherits(x[[y]], "factor")){
            levs <- x %>% pull({{ y }}) %>% levels
            if (! ref) levs <- levs[- 1]
            for (i in rev(levs)){
                if (! is.null(prefix)) nmi <- paste(prefix, i, sep = "") else nmi <- i
                x <- x %>% mutate({{ nmi }} := ifelse({{ y }} == i, 1, 0), .after = {{ y }})
            }
            if (! keep) x <- x %>% dplyr::select(- {{ y }})
            x
        }
    }
    x
}
## z <- scan("~/YvesPro2/treatment_effect/spacial/ERTU:KOCH:07/sp_solow/data-ek/data91.txt") %>% matrix(nrow = 91, byrow = TRUE) %>% .[, 1:2]
## pts <- st_multipoint(z) %>% st_sfc %>% st_cast(to = "POINT") %>% st_set_crs(st_crs(sps))

## b <- sps %>% na.omit %>% arrange(iso3) %>% add_column(pts = pts)


## dist <- st_distance(b$pts, b$point) %>% diag %>% units::set_units(km)
## b <- b %>% add_column(dist)


## dd <- scan("~/YvesPro2/Ecdat2/in/databrut/sp_solow/data-ek/data91.txt") %>%
##     matrix(ncol = 7, byrow = TRUE) %>%
##     as_tibble %>%
##     set_names(c("long", "lat", "lny60", "lny95", "gy", "lns", "lnngd")) %>%
##     bind_cols(read_csv("~/YvesPro2/Ecdat2/in/databrut/sp_solow/countries.txt")) %>% 
##     mutate(gdp60 = exp(lny60), gdp95 = exp(lny95), saving = exp(lns), labgwth = exp(lnngd) - 0.05) %>%
##     select(- (lny60:lnngd)) %>%
##     relocate(gdp60:labgwth, .before = 4) %>%
##     as_tibble %>%
##     mutate(code = ifelse(code == "ZAR", "COD", code)) %>%
##     relocate(code, name, .before = 1)
## long <- dd$long
## lat <- dd$lat
## longlat <- matrix(c(dd$long, dd$lat), nrow = 91) %>% st_multipoint %>% st_sfc(crs = 4326) %>% st_cast("POINT")
## dd <- st_sf(dd, longlat)
## dd <- dd %>%
##     mutate(lns = log(saving), lnngd = log(labgwth + 0.05),
##                     growth = (log(gdp95) - log(gdp60)) / 35)

## cts <- countries(include = "Hong Kong") %>% st_set_geometry("point") %>% select(iso3)
## dd <- dd %>% left_join(as_tibble(cts), join_by(code == iso3))# %>% st_set_geometry("point")
## sl <- lm(log(gdp95) ~ lns + lnngd, dd)
## d <- dnearneigh(dd, 0, Inf)
## W1 <- nb2listwdist(d, dd, type = "idw", 
##                    style = "W", zero.policy = TRUE, alpha = 2)
## W2 <- nb2listwdist(d, dd, type = "exp", 
##                    style = "W", zero.policy = TRUE, alpha = 2)
## moran.test(resid(sl), W1)
## lm.morantest(sl, W1)
## lm.morantest(sl, W2)


## sp_solow3 <- sp_solow2 %>% mutate(lns = lns + lnngd)
## un <- loglm(gdp95 ~ lns + lnngd, sp_solow2)
## deux <- loglm(gdp95 ~ lns, sp_solow3)

#' Number of parameters of a fitted model
#'
#' The number of observation of a fitted model is typically obtained
#' using the `nobs` method. There is no such generics to extract the
#' same information about the number of parameters. `npar` is such a
#' generic and has a special method for `micsr` objects with a
#' `subset` argument that enables to compute the number of parameters
#' for a subset of coefficients. The default method returns the length
#' of the vector of coefficients extracted using the `coef` function
#'
#' @name npar
#' @param x a fitted model
#' @param subset a character indicating the subset of coefficients
#'     (only relevant for `micsr` models).
#' @return an integer
#' @author Yves Croissant
#' @export
npar <- function(x, subset = NULL)
    UseMethod("npar")

npar.default <- function(x){
    length(coef(x))
}
npar.micsr <- function(x, subset = NULL){
    result <- x$npar
    if (! is.null(subset)){
        if (! is.character(subset)) stop("subset should be a character")
        if (any(! subset %in% names(result))) stop("unknown subset")
        result <- result[subset]
    }
    sum(as.numeric(result))
}


