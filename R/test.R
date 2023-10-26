#' Hausman test
#'
#' Hausman test 
#'
#' @name haustest
#' @param x the first model,
#' @param y the second model
#' @param omit a character containing the effects that are removed from the test
#' @return a list with class `'htest'` containing the following components:
#' - data.mane: a character string describing the fitted model
#' - statistic: the value of the test statistic
#' - parameter: degrees of freedom
#' - p.value: the p.value of the test
#' - method: a character indicating what type of test is performed
#' - alternative: a character indicating the alternative hypothesis
#' @keywords htest
#' @author Yves Croissant
#' @importFrom stats pchisq
#' @references
#' \insertRef{HAUS:78}{micsr}
#' @examples
#' charitable$logdon <- with(charitable, log(donation) - log(25))
#' char_form <- logdon ~ log(donparents) + log(income) +
#'     education + religion + married + south
#' ml <- tobit1(char_form, data = charitable)
#' scls <- update(ml, method = "trimmed")
#' haustest(scls, ml, omit = "(Intercept)")
#' @export
haustest <- function(x, y, omit = NULL){
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


#' Score test
#'
#' Score test, also knowned as Lagrange multiplier tests
#'
#' @name scoretest
#' @param x the first model,
#' @param y the second model
#' @param vcov omit a character containing the effects that are
#'     removed from the test
#' @param ... further arguments
#' @return a list with class `'htest'` containing the following
#'     components: - data.mane: a character string describing the
#'     fitted model - statistic: the value of the test statistic -
#'     parameter: degrees of freedom - p.value: the p.value of the
#'     test - method: a character indicating what type of test is
#'     performed - alternative: a character indicating the alternative
#'     hypothesis
#' @keywords htest
#' @author Yves Croissant
#' @importFrom stats pchisq
#' @examples
#' mode_choice <- mode_choice %>%
#'    mutate(cost = cost * 8.42,
#'           gcost = (ivtime + ovtime) * 8 + cost)
#' pbt_unconst <- binomreg(mode ~ cost + ivtime + ovtime, data = mode_choice, link = "probit")
#' pbt_const <- binomreg(mode ~ gcost, data = mode_choice, link = "logit")
#' scoretest(pbt_const , . ~ . + ivtime + ovtime)

#' @export
scoretest <- function(x, y, ...){
    UseMethod("scoretest")
}



#' @rdname scoretest
#' @export
scoretest.micsr <- function(x, y, ..., vcov = NULL){
    if (is.null(vcov)){
        if (is.null(x$info)) .vcov <- "hessian" else .vcov <- "info"
    }
    else .vcov <- vcov
    newform <- y
    nX <- model.matrix(x, y)
    new_nms <- colnames(nX)
    old_nms <- names(coef(x))
    L <- length(new_nms) - length(old_nms)
    if(! all(old_nms %in% new_nms))
        stop("the old model is not nested in the new one")
    start <- coef(x)[new_nms]
    start[is.na(start)] <- 0
    names(start) <- new_nms
    y <- update(x, formula = y, start = start)
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
