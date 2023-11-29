mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
dmills <- function(x) - mills(x) * (x + mills(x))
d2mills <- function(x) mills(x) * ( (x + mills(x)) * (x + 2 * mills(x)) - 1)

newton <- function(fun, coefs, trace = 0, direction = c("min", "max"), tol = sqrt(.Machine$double.eps), ...){
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

#' Short print of the summary of an object
#'
#' `print` and `print.summary` methods often returns long input, which
#' is suitable for the console, but to verbal for a printed output
#' like a book or an article written using quarto. `sight` is a generic
#' function which prints a short output
#' @name sight
#' @param x an object,
#' @param ... further arguments for the different methods,
#' @param first_stage a boolean for the `rdrobust::rdrobust` method,
#'     if `TRUE` the results of the first stage estimation,
#' @param coef the coefficients to be printed
#' @param digits the number of digits for the `lm` and the `ivreg`
#'     methods
#' @param signif.stars a boolean indicating whether the stars should
#'     be printed
#' @return invisibly its first argument
#' @export
sight <- function(x, ...) UseMethod("sight")


#' @rdname sight
#' @export
sight.lm <- function(x, ..., coef = NULL,
                       digits = max(3L, getOption("digits") - 3L), 
                       signif.stars = FALSE){
  .coef <- coef
  if (is.null(coef)){
    if (names(coef(x)[1]) == "(Intercept)") .coef <- 2:length(coef(x))
    else .coef <- 1:length(coef(x))
  }
  .x <- coef(summary(x))[.coef, , drop = FALSE]
  printCoefmat(.x, digits = digits, signif.stars = signif.stars, 
               na.print = "NA", ...)
  invisible(x)
}

#' @rdname sight
#' @export
sight.ivreg <- function(x, ..., coef,
                          digits = max(3L, getOption("digits") - 3L), 
                          signif.stars = getOption("show.signif.stars")){
  sight.lm(x, ..., coef = coef, digits = digits, signif.stars = signif.stars)
}

#' @rdname sight
#' @export
sight.rdrobust <- function(x, ..., first_stage = FALSE){
    cat(paste("Bandwidth (N used)     : ", sprintf("%3.3f", x$bws[1, 1]), " (", sum(x$N_h), ")\n", sep = ""))

    if (! is.null(x$tau_T) & first_stage){
        coef_first <- x$tau_T[1]
        se_first <- x$se_T[1]
        tstat_first <- x$t_T[c(1, 3)]
        pvals_first <- 2 * pnorm(tstat_first, lower.tail = FALSE)
        cat("First-stage estimate\n")
        cat("====================\n")
        cat(paste("coefficient (se)       : ",
                  sprintf("%3.3f", coef_first), " (", sprintf("%3.3f", se_first), ")\n", sep = ""))
        cat(paste("Conv. stat (p-value)   : ",
                  sprintf("%3.3f", tstat_first[1]), " (", sprintf("%3.3f", pvals_first[1]), ")\n", sep = ""))
        cat(paste("Robust. stat (p-value) : ",
                  sprintf("%3.3f", tstat_first[2]), " (", sprintf("%3.3f", pvals_first[2]), ")\n" , sep = ""))
    }
    coef <- x$Estimate[, 1]
    se <- x$Estimate[, 3]
    tstat <- x$Estimate[1:2] / x$Estimate[3:4]
    pvals <- 2 * pnorm(abs(tstat), lower.tail = FALSE)
    if (first_stage){
        cat("Treatment effect estimate\n")
        cat("=========================\n")
    }
    cat(paste("coefficient (se)       : ",
              sprintf("%3.3f", coef), " (", sprintf("%3.3f", se), ")\n", sep = ""))
    cat(paste("Conv. stat (p-value)   : ",
              sprintf("%3.3f", tstat[1]), " (", sprintf("%3.3f", pvals[1]), ")\n", sep = ""))
    cat(paste("Robust. stat (p-value) : ",
              sprintf("%3.3f", tstat[2]), " (", sprintf("%3.3f", pvals[2]), ")\n" , sep = ""))
    invisible(x)
}

#' @rdname sight
#' @export
sight.CJMrddensity <- function(x, ...){
    nobs <- c(x$N$eff_left, x$N$eff_right)
    stat <- x$test$t_jk
    pval <- x$test$p_jk
    bw <- c(x$h$left, x$h$right)
    cat(paste("Bandwith (left-right)     : ", sprintf("%3.3f", bw[1]), "-", sprintf("%3.3f", bw[2]), "\n", sep = ""))
    cat(paste("Observations (left-right) : ", nobs[1], "-", nobs[2], "\n", sep = ""))
    cat(paste("Statistic (p-value)       : ", sprintf("%3.3f", stat), " (", sprintf("%3.3f", pval), ")\n", sep = ""))
    invisible(x)
}

#' @rdname sight
#' @export
sight.htest <- function(x, ..., digits = 3){
    # name_raw = raw (mean-sd)
    # name = stat
    # df = param
    # pvalue = pval
    paste_element <- function(x, op = " = "){
        paste(names(x), x, sep = op)
    }
    paste_raw <- function(x, std, mu = NULL){
        if (is.null(mu)) paste(sprintf(.sprint, x), " (",
                               sprintf(.sprint, std), ")",
                               sep = "")
        else paste(sprintf(.sprint, x), " (",
                   sprintf(.sprint, mu), ", ",
                   sprintf(.sprint, std), ")",
                   sep = "")
    }
    .method <- x$method
    .nms <- names(x$statistic)
    .stat <- x$statistic
    .param <- x$parameter
    .pval <- x$p.value
    .sprint <- "%3.4f"
    .sprint <- paste("%3.", digits, "f", sep = "")
    .raw <- NULL
    if (.method %in% c("Welch Two Sample t-test", " Two Sample t-test")){
        .name_raw <- "diff"
        .stderr <- x$stderr
        .raw <- .stat * .stderr
        .raw <- paste_raw(.raw, .stderr)
        names(.raw) <- .name_raw
        .raw <- paste_element(.raw)
    }
    if (.method == "Moran I test under randomisation"){
        .nms <- "z"
        .name_raw <- "Moran I"
        .stderr <- sqrt(x$estimate[3])
        .mu <- x$estimate[2]
        .raw <- x$estimate[1]
        .raw <- paste_raw(.raw, .stderr, .mu)
        names(.raw) <- .name_raw
        .raw <- paste_element(.raw)
    }
        if (.method == "Geary C test under randomisation"){
        .nms <- "z"
        .name_raw <- "Geary C"
        .stderr <- sqrt(x$estimate[3])
        .mu <- x$estimate[2]
        .raw <- x$estimate[1]
        .raw <- paste_raw(.raw, .stderr, .mu)
        names(.raw) <- .name_raw
        .raw <- paste_element(.raw)
    }
    if (.method == "Global Moran I for regression residuals"){
        .nms <- "z"
        .name_raw <- "Moran I"
        .stderr <- sqrt(x$estimate[3])
        .mu <- x$estimate[2]
        .raw <- x$estimate[1]
        .raw <- paste_raw(.raw, .stderr, .mu)
        names(.raw) <- .name_raw
        .raw <- paste_element(.raw)
    }
    # degrees of freedom
    if (! is.null(.param)){
        .param <- sprintf("%0.0f", .param)
        names(.param) <- "df"
        .param <- paste_element(.param, ": ")
    }
    # p-value
    .pval <- sprintf(.sprint, .pval)
    names(.pval) <- "pval"
    .pval <- paste_element(.pval)
    # statistic
    .stat <- sprintf(.sprint, .stat)
    names(.stat) <- .nms
    .stat <- paste_element(.stat)
    cat(paste(paste(c(.raw, .stat, .param, .pval), collapse = ", "),
              "\n", sep = ""))
    invisible(x)
}

#' @rdname sight
#' @export
sight.LMtestlist <- function(x, ..., digits = 3) sight(x[[1]])



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
