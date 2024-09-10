#' Short print of the summary of an object
#'
#' `print` and `print.summary` methods often returns long input, which
#' is suitable for the console, but too verbal for a printed output
#' like a book or an article written using quarto. `gaze` is a generic
#' function which prints a short output
#' @name gaze
#' @param x an object,
#' @param ... further arguments for the different methods,
#' @param first_stage a boolean for the `rdrobust::rdrobust` method,
#'     if `TRUE` the results of the first stage estimation are printed
#' @param coef the coefficients to be printed
#' @param digits the number of digits for the `lm` and the `ivreg`
#'     methods
#' @param signif.stars a boolean indicating whether the stars should
#'     be printed
#' @return returns invisibly its first argument
#' @keywords misc
#' @examples
#' t.test(extra ~ group, sleep) %>% gaze
#' lm(dist ~ poly(speed, 2), cars) %>% gaze
#' lm(dist ~ poly(speed, 2), cars) %>% gaze(coef = "poly(speed, 2)2")
#' @export
gaze <- function(x, ...) UseMethod("gaze")

#' @rdname gaze
#' @export
gaze.lm <- function(x, ..., coef = NULL,
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

#' @rdname gaze
#' @export
gaze.micsr <- function(x, ..., coef = NULL,
                        digits = max(3L, getOption("digits") - 3L), 
                        signif.stars = FALSE){
  gaze.lm(x, ..., coef = coef, digits = digits, signif.stars = signif.stars)
}

#' @rdname gaze
#' @export
gaze.ivreg <- function(x, ..., coef = NULL,
                       digits = max(3L, getOption("digits") - 3L),
                       signif.stars = FALSE){
#                          signif.stars = getOption("show.signif.stars")){
  gaze.lm(x, ..., coef = coef, digits = digits, signif.stars = signif.stars)
}

#' @rdname gaze
#' @export
gaze.mlogit <- function(x, ..., coef = NULL,
                       digits = max(3L, getOption("digits") - 3L),
                       signif.stars = FALSE){
  gaze.lm(x, ..., coef = coef, digits = digits, signif.stars = signif.stars)
}


#' @rdname gaze
#' @export
gaze.rdrobust <- function(x, ..., first_stage = FALSE){
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

#' @rdname gaze
#' @export
gaze.CJMrddensity <- function(x, ...){
    nobs <- c(x$N$eff_left, x$N$eff_right)
    stat <- x$test$t_jk
    pval <- x$test$p_jk
    bw <- c(x$h$left, x$h$right)
    cat(paste("Bandwith (left-right)     : ", sprintf("%3.3f", bw[1]), "-", sprintf("%3.3f", bw[2]), "\n", sep = ""))
    cat(paste("Observations (left-right) : ", nobs[1], "-", nobs[2], "\n", sep = ""))
    cat(paste("Statistic (p-value)       : ", sprintf("%3.3f", stat), " (", sprintf("%3.3f", pval), ")\n", sep = ""))
    invisible(x)
}

#' @rdname gaze
#' @export
gaze.htest <- function(x, ..., digits = 3){
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
        if (length(.param) == 2) .param <- paste(.param, collapse = "-")#.param <- paste("(", paste(.param, collapse = "-"), ")", sep = "")
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

#' @rdname gaze
#' @export
gaze.anova <- function(x, ..., digits = 3){
    if ("Chisq" %in% names(x)){
        .stat <- c("Chisq" = x$Chisq[2])
        .pval <- x[["Pr(>Chisq)"]][2]
        .param <- abs(x$Df[2])
    }
    if ("F" %in% names(x)){
        .stat <- c("F" = x$F[2])
        .pval <- x[["Pr(>F)"]][2]
        if ("RSS" %in% names(x)) .param <- c(x$Df[2], min(x$RSS))
        if ("Res.Df" %in% names(x)) .param <- c(x$Df[2], min(x$Res.Df))
        
    }
    .result <- structure(list(method = "anova", statistic = .stat,
                              parameter = .param, p.value = .pval), class = "htest")
    gaze(.result, ..., digits = 3)
}
    

#' @rdname gaze
#' @export
gaze.LMtestlist <- function(x, ..., digits = 3){
    .sprint <- paste("%6.", digits, "f", sep = "")
    if (length(x) == 1) gaze(x[[1]], ..., digits= digits)
    else{
        .name <- unname(sapply(x, function(x) names(x$statistic)))
        .nchar <- nchar(.name)
        .nchar <- max(.nchar) - .nchar
        .blanks <- sapply(.nchar, function(x) paste(rep(" ", x), collapse = ""))
        .name <- paste(.name, .blanks, sep = "")
        for (i in 1:length(x)){            
            cat(.name[i], ": ", sprintf(.sprint, unname(x[[i]]$statistic)), ", df = ",
                      unname(x[[i]]$parameter), ", p-value = ", sprintf(.sprint, unname(x[[i]]$p.value)), "\n", sep = "")
        }
    }
    invisible(x)
}
        
#' @rdname gaze
#' @export
gaze.RStestlist <- function(x, ..., digits = 3){
    .sprint <- paste("%6.", digits, "f", sep = "")
    if (length(x) == 1) gaze(x[[1]], ..., digits= digits)
    else{
        .name <- unname(sapply(x, function(x) names(x$statistic)))
        .nchar <- nchar(.name)
        .nchar <- max(.nchar) - .nchar
        .blanks <- sapply(.nchar, function(x) paste(rep(" ", x), collapse = ""))
        .name <- paste(.name, .blanks, sep = "")
        for (i in 1:length(x)){            
            cat(.name[i], ": ", sprintf(.sprint, unname(x[[i]]$statistic)), ", df = ",
                      unname(x[[i]]$parameter), ", p-value = ", sprintf(.sprint, unname(x[[i]]$p.value)), "\n", sep = "")
        }
    }
    invisible(x)
}
        
