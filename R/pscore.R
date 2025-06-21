#' Propensity scores
#'
#' Propensity scores estimation, using an algorithm that checks the
#' balancing hypothesis using strata and enable the estimation of the
#' treatment effect using stratification methods
#'
#' @name pscore
#' @param formula a Formula object; the left-hand side should contain
#'     two variables (`x1 + x2`), where x1 is the group variable and
#'     x2 the outcome. The group variable can be either a dummy for
#'     treated individuals or a factor with levels `"treated"` and
#'     `"control"`
#' @param data a data frame
#' @param maxiter the maximum number of iterations
#' @param tol stratas are cut in halves as long as the hypothesis of
#'     equal means is rejected at the `tol` level,
#' @param link the link for the binomial glm estimation, either
#'     `"logit"` or `"probit"`
#' @param x,object a `"pscore"` or a `"summary.pscore"` object
#' @param var_equal to compute the variance of the ATET, variances can
#'     be computed at the class/group level (`var_equal = "none"`), at
#'     the class level (`var_equal = "group"`), at the group level
#'     (`var_equal = "strata"`) or globally (`var_equal = "both"`)
#' @param step for the `print.summary` method, the step of the test to
#'     be printed: one of `"all"` (the default), `strata`,
#'     `covariates` and `atet`
#' @param smpl the sample to use, either the whole sample (`smpl =
#'     "total"`) or the sample with common support (`smpl = "cs"`)
#' @param digits number of digits for the `print` methods
#' @param ... further arguments
#' @keywords models
#' @importFrom stats na.omit quantile var aggregate reshape t.test
#' @references \insertRef{DEHE:WAHB:02}{micsr}
#'
#' \insertRef{BECK:ICHI:02}{micsr}
#' @return an object of class `"pscore"`, with the following elements:
#'
#' - `strata`: a tibble containing the stratas, the frequencies, the
#' means and the variances of the propensity scores for treated and
#' controled observations
#' - `cov_balance`: a tibble containing the results of the balancing
#' tests for every covariate; the results for the class with the
#' lowest p-value is reported
#' - `unchecked_cov`: a character vector containing the names of the
#' covariates for which the balancing test could be computed
#' - `model`: a tibble containing the original data, with
#' supplementary columns: `.gp` for the groups, `.resp` for the
#' outcome and `.cls` for the stratas
#' - `pscore`: the glm model fitted to compute the propensity scores
#' @examples
#' data_tuscany <- twa |>
#'                 subset(region == "Tuscany") |>
#'                 transform(dist2 = dist ^ 2,
#'                 livselfemp = I((city == "livorno") * (occup == "selfemp")),
#'                 perm = ifelse(outcome == "perm", 1, 0))
#' formula_tuscany <- perm + group ~ city + sex + marital + age +
#'    loc + children + educ + pvoto + training +
#'    empstat + occup + sector + wage + hour + feduc + femp + fbluecol +
#'    dist + dist2 + livselfemp
#' pscore(formula_tuscany, data_tuscany)
#' @export
pscore <- function(formula, data, maxiter = 4, tol = 0.005, link = c("logit", "probit")){
    .link <- match.arg(link)
    .maxiter <- maxiter
    .tol <- tol
    .model <- data
    # QAD fix, now y1 + y2 ~ ... instead of y1 | y2 ~ ...
    .formula <- paste(deparse(formula), collapse = "")
    .formula <- sub("\\+", "\\|", .formula)
    .formula <- Formula(formula(.formula))
    group_name <- paste(.formula[[2]][[3]])
    response_name <- paste(.formula[[2]][[2]])
    .model[[".gp"]] <- .model[[group_name]]
    .model[[".resp"]] <- .model[[response_name]]
    if (is.numeric(data[[group_name]]))
        .model[[".gp"]] <- factor(data[[group_name]], levels = 0:1,
                                 labels = c("control", "treated"))
    outcome_name <- paste(.formula[[2]][[3]])
    pscore <- glm(formula(.formula, lhs = 2), data = .model, family = binomial(link = .link))
    na_rows <- attr(model.frame(pscore), "na.action")
    if (! is.null(na_rows)) .model <- .model[- na_rows, ]
    .model$pscore <- fitted(pscore)
    min_sup <- min(subset(.model, .model$.gp == "treated")$pscore)
    max_sup <- max(subset(.model, .model$.gp == "treated")$pscore)

    mateub <- transform(.model, pscore = pscore + 2)
    .model <- .model |> transform(.cs = (pscore <= max_sup & pscore >= min_sup))
    scores <- subset(.model, .model$.cs)$pscore
    bks <- c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
    qtles <- quantile(scores, bks)
    .go_on <- TRUE
    .iter <- 0
    bks <- c(seq(0, 1, 0.2))
    while(.go_on){
        .iter <- .iter + 1
        da <- subset(.model, .model$.cs) |>
            transform(cls = cut(scores, bks, right = FALSE))
        z <- aggregate(.resp ~ .gp + cls, data = da,
                       FUN = function(x) c(n = length(x), mean = mean(x), var = var(x)),
                       simplify = TRUE)
        z <- cbind(z[, 1:2], n = z[[3]][, 1], mean = z[[3]][, 2], var = z[[3]][, 3])
        z <- reshape(z, direction = "wide", timevar = ".gp",
                     v.names = c("n", "mean", "var"), idvar = "cls",
                     sep = "_")
        z$cls <- as.character(z$cls)
        freq <- z[order(z$cls), c("cls", "n_control", "n_treated", "mean_control",
                               "mean_treated", "var_control", "var_treated")]
        no_cs <- is.na(freq$n_control) | is.na(freq$n_treated)
        r <- lapply(freq$cls[! no_cs],
                    function(x) t.test(pscore ~ .gp, var.equal = TRUE,
                                       data = subset(da, da$cls == x)))
        r <- lapply(r, function(x) c(x$estimate, x$p.value))
        r <- Reduce("rbind", r)
        dimnames(r) <- list(freq$cls[! no_cs], c("ps_control", "ps_treated", "p.val"))
        r <- cbind(cls = rownames(r), as.data.frame(r))
        rownames(r) <- NULL
        r <- merge(r, freq, by = "cls")
        r <- r[order(r$cls), ]
        if (any(na.omit(r$p.val) < .tol) & .iter <= .maxiter){
            pb <- which(r$p.val < .tol)
            sup_bks <- (bks[pb + 1] + bks[pb]) / 2
            bks <- sort(c(bks, sup_bks))
        }
        else .go_on <- FALSE
    }
   .model <- transform(.model, .cls = cut(pscore, bks, right = FALSE))
    strata <- r
    X <- model.matrix(.formula, data = model.frame(.formula, subset(.model, .model$.cs))) |>
        as.data.frame()
    X <- cbind(X, da[, c("pscore", "cls", ".gp")])

    xnms <- names(X)[- na.omit(match(c("(Intercept)", "cls", "pscore", ".gp"), names(X)))]
    cov_balance <- data.frame(name = character(0), classe = character(0), pvalue = numeric(0))
    unchecked_cov <- c()
    for (aname in xnms){
        X$covariate <- X[[aname]]
        r2 <- try(lapply(strata$cls[! no_cs],
                        function(x) t.test(covariate ~ .gp, var.equal = TRUE,
                                           data = subset(X, X$cls == x))), silent = TRUE)
        if (length(r2[[1]]) > 1){
            r2 <- lapply(r2, function(x) c(x$estimate, x$p.value))
            r2 <- Reduce("rbind", r2)
            dimnames(r2) <- list(strata$cls[! no_cs], c("control", "treated", "p.val"))
            r2 <- cbind(cls = strata$cls[! no_cs], as.data.frame(r2))
            rownames(r2) <- NULL
            r2 <- r2[order(r2$p.val), ]
            r2 <- r2[1, ]
            if (! is.na(r2$p.val)){
                cov_balance <- rbind(cov_balance,
                                     list(name = aname, classe = r2$cls, pvalue = r2$p.val))
            }
        } else {
            unchecked_cov <- c(unchecked_cov, aname)
        }


    }
    class(strata) <- class(cov_balance) <- class(.model) <- c("tbl_df", "tbl", "data.frame")
    structure(list(strata = strata, cov_balance = cov_balance,
                   unchecked_cov = unchecked_cov, model = .model,
                   pscore = pscore),
              class = "pscore")
}

#' @rdname pscore
#' @export
summary.pscore <- function(object, ...){
    # computation of the average treatment effect to the treated
    strata <- object$strata
    .model <- object$model
    no_cs <- is.na(strata$n_control) | is.na(strata$n_treated)
    rg_cs <- rg(object, smpl = "cs")#.model |> filter(.cs) |> pull(pscore) |> range
    rg_tot <- rg(object, smpl = "total")#.model |> pull(pscore) |> range
    nobs_tot <- nobs(object, smpl = "total")#.model |> pull(.gp) |> table |> as.numeric
    nobs_cs <- nobs(object, smpl = "cs")#.model |> filter(.cs) |> pull(.gp) |> table |> as.numeric    
    if (any(no_cs)){
        strata_save <- strata
        .model_save <- .model
        cls_no_cs <- strata$cls[no_cs]
        strata <- subset(strata, ! strata$.cls %in% cls_no_cs)
        .model <- subset(.model, ! strata$.cls %in% cls_no_cs)
    }

    ATET <- with(strata, sum(n_treated / sum(n_treated) *
                             (mean_treated - mean_control)))
    Nt <- sum(strata$n_treated)
    g <- strata$n_treated / strata$n_control
    f <- strata$n_treated / sum(strata$n_treated)
    Vg <- as.numeric(tapply(.model$.resp, .model$.gp, var))
    Vc <- Vg[1]
    Vt <- Vg[2]
    Vb <- tapply(.model$.resp, .model$.cls, var)
    Vgp <- aggregate(.resp ~ .gp + .cls, data = .model, FUN = var)
    Vgp <- reshape(Vgp, direction = "wide", timevar = ".gp", idvar = ".cls", v.names = ".resp")
    names(Vgp)[c(2:3)] <- c("control", "treated")
    Vtb <- Vgp$treated
    Vcb <- Vgp$control
    V <- var(.model$.resp)
    sd_both   <- sqrt(V * (1 + sum(f * g)) / Nt)
    sd_group  <- sqrt( (Vt + Vc * sum(f * g)) / Nt)
    sd_strata <- sqrt( sum(f * Vb * (1 + g)) / Nt)
    sd_none   <- sqrt( sum(f * (Vtb + g * Vcb) / Nt))
    atet <- c(atet = ATET, sd_both = sd_both, sd_group = sd_group,
              sd_strata = sd_strata, sd_none = sd_none)
    if (any(no_cs)){
        strata <- strata_save
        .model <- .model_save
    }
    strata <- strata[, c("cls", "n_treated", "n_control",
               "ps_treated", "ps_control", "p.val",
               "mean_treated", "mean_control")]
    names(strata)[c(1, 6)] <- c("strata", "p.value")
    structure(list(strata = strata, pscore = object$pscore,
                   cov_balance = object$cov_balance,
                   cov_unchecked = object$cov_unchecked,
                   range = list(total = rg_tot, cs = rg_cs),
                   nobs = list(total = nobs_tot, cs = nobs_cs),
                   atet = atet),
              class = "summary.pscore")
}

#' @rdname pscore
#' @export
print.pscore <- function(x, ..., digits = getOption("digits"),
                         var_equal = c("none", "strata", "group", "both")){
    var_equal <- match.arg(var_equal)
    x <- summary(x)
    cat("Number of strata", nrow(x$strata), "\n")
    cat("pscore range  :", paste(signif(rg(x, smpl = "total"), digits), collapse = " - "), "\n")
    cat("common support:", paste(signif(rg(x, smpl = "cs"), digits), collapse = " - "), "\n")
    cat("untreated obs used:", nobs(x, smpl = "cs")[1], "out of", nobs(x, smpl = "total")[1], "\n")
    cat("ATET: ", signif(mean(x), digits), " (", signif(stdev(x, var_equal = var_equal), digits), ")\n", sep = "")
    invisible(x)
}
    
#' @rdname pscore
#' @export
print.summary.pscore <- function(x, ...,
                                 digits = getOption("digits"),
                                 step = c("all", "strata", "covariates", "atet")){
    .step <- match.arg(step)
    .format <- match.arg(format)
    if (.step == "all"){
        cat("=========================================================\n")
        cat("Propensity score computation using Deheja-Wahba algorithm\n")
        cat("=========================================================\n\n")
        cat("----------------------------------------------------\n")
        cat("first-step: estimation of the participation equation\n")
        cat("----------------------------------------------------\n")
        cat("Binomial model with a", family(x$pscore)$link, "link\n")
        cat("Number of strata  :", nrow(x$strata), "\n")
        cat("pscore range      :", signif(x$range$total[1], digits), "-", signif(x$range$total[2], digits), "\n")
        cat("Common support    :", signif(x$range$cs[1], digits), "-", signif(x$range$cs[2], digits), "\n")
        cat("untreated obs used:", x$nobs$cs[1], "out of", x$nobs$total[1], "\n")
        cat("(p-value for the t-test of equal means)\n")
    }
    if (.step %in% c("all", "strata"))
        ## x$strata[, 1:6] |>
        ##     knitr::kable(digits = c(0, 0, 0, 3, 3, 3, 1, 1), format = .format) |>
        ##     print()
        print(x$strata[, 1:6], n = Inf, digits = digits)
        if (.step == "all"){
        cat("\n")
        cat("----------------------------------------------------\n")
        cat("second-step: checking the balance for the covariates\n")
        cat("----------------------------------------------------\n")
    }
    if (.step %in% c("all", "covariates")){
        cat("(for every covariate, strata with the lowest p-value)\n")
#        x$cov_balance |> knitr::kable(digits = digits, format = .format) |> print()
        print(x$cov_balance, n = Inf, digits = digits)
        if (! is.null(x$unchecked_cov)){
            cat("\nBalance checking impossible for:", paste(x$unchecked_cov, collapse = ", "), "\n")
        }
        cat("\n\n")
    }
    if (.step == "all"){
        cat("-----------------------------------\n")
        cat("third-step: computation of the ATET\n")
        cat("-----------------------------------\n")
    }
    if (.step %in% c("all", "atet")){
        cat("ATET                        :", signif(x$atet["atet"], digits), "\n")
        cat("sd: equal variance\n")
        cat("  - within groups           :", signif(x$atet["sd_group"], digits), "\n")
        cat("  - within strata           :", signif(x$atet["sd_strata"], digits), "\n")
        cat("  - within groups and strata:", signif(x$atet["sd_both"], digits), "\n")
        cat("  - no                      :", signif(x$atet["sd_none"], digits), "\n")
    }
}
    
#' @rdname pscore
#' @export
nobs.pscore <- function(object, ..., smpl = c("total", "cs")){
    .smpl <- match.arg(smpl)
    if (.smpl == "total") result <- object$model$.gp |> table() |> as.numeric()
    if (.smpl == "cs") result <- subset(object$model, object$model$.cs)$.gp |>
                                                         table() |> as.numeric()
    result
}

#' @rdname pscore
#' @export
nobs.summary.pscore <- function(object, ..., smpl = c("total", "cs")){
    .smpl <- match.arg(smpl)
    if (.smpl == "total") result <- object$nobs$total
    if (.smpl == "cs") result <- object$nobs$cs
    result
}


#' @rdname pscore
#' @export
rg <- function(object, ...) UseMethod("rg")

#' @rdname pscore
#' @export
rg.pscore <- function(object, ..., smpl = c("total", "cs")){
    .smpl <- match.arg(smpl)
    if (.smpl == "total") result <- range(object$model$pscore)
    if (.smpl == "cs") result <- range(subset(object$model, object$model$.cs)$pscore)
    result
}

#' @rdname pscore
#' @export
rg.summary.pscore <- function(object, ..., smpl = c("total", "cs")){
    .smpl <- match.arg(smpl)
    if (.smpl == "total") result <- object$range$total
    if (.smpl == "cs") result <- object$range$cs
    result
}

#' @rdname pscore
#' @export
stdev <- function(object, ...) UseMethod("stdev")

#' @rdname pscore
#' @export
mean.pscore <- function(x, ..., var_equal = c("none", "strat", "group", "both")){
    summary(x)$atet["atet"] |> unname()
}

#' @rdname pscore
#' @export
mean.summary.pscore <- function(x, ...){
    x$atet["atet"] |> unname()
}

#' @rdname pscore
#' @export
stdev.pscore <- function(object, ..., var_equal = c("none", "strata", "group", "both")){
    var_equal <- match.arg(var_equal)
    object <- summary(object)$atet
    switch(var_equal,
           none   = unname(object["sd_none"]),
           both   = unname(object["sd_both"]),
           strata = unname(object["sd_strata"]),
           group  = unname(object["sd_group"]))
}

#' @rdname pscore
#' @export
stdev.summary.pscore <- function(object, ..., var_equal = c("none", "strata", "group", "both")){
    var_equal <- match.arg(var_equal)
    object <- object$atet
    switch(var_equal,
           none   = unname(object["sd_none"]),
           both   = unname(object["sd_both"]),
           strata = unname(object["sd_strata"]),
           group  = unname(object["sd_group"]))
}
