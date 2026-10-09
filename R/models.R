#' Print a list of models
#'
#' Provide a class for a list of models and print a model summary
#' table, with a choice of uncertainty and goodness of fit measures
#'
#' The function just construct a named list of models, the `print`
#' method prints a model summary
#'
#' @name models
#' @param ... each argument is a model, the name being the name that
#'     will appear on the table
#' @param x a `models` object
#' @param digits the number of digits
#' @param statistic the uncertainty measure, one of `"se"` (standard
#'     error), `"t"` the t-statistic, `"pval"` the probability values
#'     and `"stars"`
#' @param gof a vector of goodness of fit measures that may contain
#'     `"nobs"`, `"logLik"`, `"AIC"` and `"BIC"`
#' @param brackets a character indicating how the statistic is
#'     presented and separated from the coefficient, may be `"(" (the
#'     default for `"se"`, `"t"` and `"pval"`), `" "` (the default for
#'     `"stars"`, `"["` or `"-"
#' @return an object of class `models`, which is simply a named list
#'     of fitted models
#' @author Yves Croissant
#' @examples
#' ols <- lm(trips ~ car + dist + realinc, trips)
#' pois <- glm(trips ~ car + dist + realinc, family = poisson, trips)
#' nb2 <- poisreg(trips ~ car + realinc,  mixing = "gamma", vlink = "nb2", trips)
#' mdls <- models(OLS = ols, Poisson = pois, NB2 = nb2)
#' mdls
#' mdls |> print(digits = 2, statistic = "stars", brackets = "[", gof = c("nobs", "BIC"))
#' @export 
models <- function(...){
    structure(list(...), class = "models")
}

#' @rdname models
#' @export
print.models <- function(x, ..., 
                         digits = 3,
                         statistic = c("se", "t", "pval", "stars"),
                         gof = NULL,
                         brackets = NULL){
    statistic <- match.arg(statistic)
    if (is.null(brackets)){
        if (statistic %in% c("se", "t", "pval")) brackets <- c(" (", ")")
        if (statistic == "stars") brackets = c(" ", "")
    }
    if (length(brackets) == 1){
        .brackets <- brackets
        if (.brackets == "(") brackets <- c(" (", ")")
        if (.brackets == "[") brackets <- c(" [", "]")
        if (.brackets == "-") brackets <- c(" - ", "")
    }
    .models <- x
    av_gof <- c("nobs", "logLik", "AIC", "BIC")
    if (is.null(gof)){
        gof <- av_gof
    } else {
        if (! all(gof %in% av_gof)) stop("undefined gof")
    }
    
    get_st <- function(x, statistic){
        if (statistic == "se") .pos <- 2
        if (statistic == "t") .pos <- 3
        if (statistic %in% c("pval", "stars")) .pos <- 4
        x <- coef(summary(x))
        nms <- rownames(x)
        x <- x[, .pos]
        if (statistic == "stars"){
            s <- rep("", length(x))
            s[x < 0.1] <- "."
            s[x < 0.05] <- "*"
            s[x < 0.01] <- "**"
            s[x < 0.001] <- "***"
            x <- s
        }
        names(x) <- nms
        x
    }
    trms <- unique(Reduce("c", lapply(x, function(x) names(coef(x)))))
    z <- lapply(x, function(x) coef(x)[trms])
    st <- lapply(x, function(x) get_st(x, statistic)[trms])
    f_st <- function(.coef, .st){
        nas <- is.na(.coef)
        x <- paste(format(.coef, digits = digits),
                   brackets[1],
                   format(.st, digits = digits),
                   brackets[2],
                   sep = "")
        x[nas] <- ""
        x
    }
    x <- mapply(f_st, z, st)
    rownames(x) <- trms
    for (i in gof){
        agof <- sapply(.models, as.name(i)) |> format(digits = digits)
        x <- rbind(x, agof)
        rownames(x)[nrow(x)] <- i
    }
    print(as.data.frame(x))
    invisible(x)
}

