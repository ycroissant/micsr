mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
dmills <- function(x) - mills(x) * (x + mills(x))
d2mills <- function(x) mills(x) * ( (x + mills(x)) * (x + 2 * mills(x)) - 1)

newton <- function(fun, coefs, trace, direction = c("min", "max"), tol = sqrt(.Machine$double.eps), ...){
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
        if (i > 500) stop("max iter reached")
        coefs <- newcoefs
    }
    coefs
}
