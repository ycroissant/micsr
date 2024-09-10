#' gaussian quadrature
#'
#' Calculate nodes and weights for Gaussian quadrature. This function
#' is a copy of `statmod::gauss.quad`
#'
#' @name gaussian_quad
#' @param n the number of nodes and weights
#' @param kind kind of Gaussian quadrature, one of `"legendre"`,
#'     `"chebyshev1"`, `"chebyshev2"`, `"hermite"`, `"jacobi"` or
#'     `"laguerre"`
#' @param alpha,beta see `statmod::gauss.quad`
#' @return a list containing the nodes and the weights
#' @export
gaussian_quad <- function(n, kind = "legendre", alpha = 0, beta = 0)
#	Calculate nodes and weights for Gaussian quadrature.
#	Adapted from Netlib routine gaussq.f
#	Gordon Smyth, Walter and Eliza Hall Institute
#	Suggestion from Stephane Laurent 6 Aug 2012
#	Created 4 Sept 2002. Last modified 28 Aug 2016.
{
	n <- as.integer(n)
	if(n<0L) stop("need non-negative number of nodes")
	if(n==0L) return(list(nodes=numeric(0L), weights=numeric(0L)))
	kind <- match.arg(kind,c("legendre","chebyshev1","chebyshev2","hermite","jacobi","laguerre"))
	i <- 1L:n
	i1 <- i[-n]
	switch(kind, legendre={
		lnmuzero <- log(2)
		a <- rep_len(0,n)
		b <- i1/sqrt(4*i1^2-1)
	}, chebyshev1={
		lnmuzero <- log(pi)
		a <- rep_len(0,n)
		b <- rep_len(0.5,n-1L)
		b[1] <- sqrt(0.5)
	}, chebyshev2={
		lnmuzero <- log(pi/2)
		a <- rep_len(0,n)
		b <- rep_len(0.5,n-1L)
	}, hermite={
		lnmuzero <- log(pi)/2
		a <- rep_len(0,n)
		b <- sqrt(i1/2)
	}, jacobi={
		ab <- alpha+beta
#		muzero <- 2^(ab+1) * gamma(alpha+1) * gamma(beta+1) / gamma(ab+2)
		lnmuzero <- (ab+1)*log(2) + lgamma(alpha+1) + lgamma(beta+1) - lgamma(ab+2)
		a <- i
		a[1] <- (beta-alpha)/(ab+2)
		i2 <- i[-1]
		abi <- ab+2*i2
		a[i2] <- (beta^2-alpha^2)/(abi-2)/abi
		b <- i1
		b[1] <- sqrt(4*(alpha+1)*(beta+1)/(ab+2)^2/(ab+3))
		i2 <- i1[-1]
		abi <- ab+2*i2
		b[i2] <- sqrt(4*i2*(i2+alpha)*(i2+beta)*(i2+ab)/(abi^2-1)/abi^2)
	}, laguerre={
		a <- 2*i-1+alpha
		b <- sqrt(i1*(i1+alpha))
		lnmuzero <- lgamma(alpha+1)
	})
	b <- c(b,0)
	z <- rep_len(0,n)
	z[1] <- 1
	ierr <- 0L
	out <- .Fortran("gausq2",n,as.double(a),as.double(b),as.double(z),ierr,PACKAGE="micsr")
	x <- out[[2]]
	w <- out[[4]]
	w <- exp(lnmuzero + 2*log(abs(w)))
	list(nodes=x,weights=w)
}

mydavies <- function(q, lambda, h = rep(1, length(lambda)), delta = rep(0, length(lambda)),
                     sigma = 0, lim = 10000, acc = 0.0001) {
  r <- length(lambda)
  if (any(delta < 0)) stop("All non centrality parameters in 'delta' should be positive!") 
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  
    out <- .C("qfc", lambdas = as.double(lambda), noncentral = as.double(delta),
              df = as.integer(h), r = as.integer(r), sigma = as.double(sigma),
              q = as.double(q), lim = as.integer(lim), acc = as.double(acc),
              trace = as.double(rep(0, 7)), ifault = as.integer(0), res = as.double(0), PACKAGE = "micsr")
  out$res <- 1 - out$res
  if (out$res > 1) warning("Consider playing with 'lim' or 'acc'.")
  return(list(trace = out$trace, ifault = out$ifault, Qq = out$res))
}

punorm <- function(z){
    z <- as.double(z)
    lg <- as.integer(length(z))
    prob <- double(lg)
    ans <- .Fortran("MYUNORM", prob, z, lg)[[1]]
    ans
}
    
pbnorm <- function(z1, z2, rho){
    z1 <- as.double(z1)
    z2 <- as.double(z2)
    al <- as.integer(length(z1))
    rho <- as.double(rho)
    prob <- double(al)
    ans <- .Fortran("MYBNORM", prob, z1, z2, rho, al)[[1]]
    ans
}

ptnorm <- function(z1, z2, z3, rho){
    # rho23 must be the largest corr coefficient in absolute value
    z <- cbind(z1, z2, z3)
    mxr <- which.max(abs(rho))
    if (mxr == 1){
        posh <- 1:3
        posr <- 1:3
    }
    if (mxr == 2){
        posh <- c(2, 1, 3)
        posr <- c(1, 3, 2)
    }
    if (mxr == 3){
        posh <- c(3, 1, 2)
        posr <- c(2, 3, 1)
    }
    z <- z[, posh, drop = FALSE]
    rho <- rho[posr]
    z1 <- as.double(z[, 1])
    z2 <- as.double(z[, 2])
    z3 <- as.double(z[, 3])
    al <- as.integer(nrow(z))
    rho <- as.double(rho)
    prob <- double(al)
    ans <- .Fortran("MYTNORM", prob, z1, z2, z3, rho, al)[[1]]
    ans
    return(ans)
}
