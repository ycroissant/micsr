is_count <- function(x) is.integer(x) &
                            (length(unique(x)) > 2 |
                             (length(unique(x)) == 2 ) &
                             ! all(sort(unique(x)) %in% c(0L, 1L)))
is_dummy <- function(x) is.integer(x) &
                            (length(unique(x)) == 2  &
                             all(sort(unique(x)) %in% c(0L, 1L)))
is_factor <- function(x) inherits(x, "factor")
is_dbl <- function(x) is.numeric(x) & ! is_count(x) & ! is_dummy(x)

covar_type <- function(x){
        type <- NA
        if (is_count(x))  type <- "count"
        if (is_dummy(x))  type <- "dummy"
        if (is_factor(x)) type <- "fct"
        if (is_dbl(x))    type <- "dbl"
        type
    }

is_definite_positive <- function(x){
    isdefpos <- TRUE
    if (any(is.na(x))){
        isdefpos <- FALSE
    } else {
        .eigen <- eigen(x, only.values = TRUE)$values
        if (any(is.complex(.eigen))){
            isdefpos <- FALSE
        } else {
            isdefpos <- all(.eigen > 1E-08)
        }
    }
    isdefpos
}

get_data <- function(x) eval(x$call$data, parent.frame())

get_response <- function(x){
    .form <- formula(x)
    if (! inherits(.form, "Formula")) .form <- Formula(.form)
    z <- formula(.form, lhs = 1, rhs = 0)
    deparse(z[[2]])
}
    

get_covariates <- function(x, rhs = 1){
    .form <- formula(x)
    if (! inherits(.form, "Formula")) .form <- Formula(.form)
    z <- formula(.form, lhs = 0, rhs = rhs)
    names(get_all_vars(z, get_data(x)))
}

add_a_slope <- function(coefs, object, covar, level = NULL, type, newdata, alt = NULL, se = TRUE){
    .effect = switch(type,
                     dbl   = "dy/dx",
                     dummy = "1 - 0",
                     count = "+ 1",
                     fct   = "contrast")
    .estimate <- estimate(coefs = coefs, model = object, covar = covar, level = level,
                          type = type, newdata = newdata, alt = alt)
    if (se){
        .std.er <- std.er(model = object, covar = covar, level = level,
                          type = type, newdata = newdata, alt = alt)
    }
    .term <- covar
    if (type == "fct") .term <- paste(covar, level, sep = "_")
    if (inherits(object, "mlogit")){
        ids <- newdata[[idx_name(object, 1, 1)]]
        alts <- newdata[[idx_name(object, 2, 1)]]
        n_id <- length(unique(ids))
        n_alt <- length(unique(alts))
        ids <- rep(unique(ids), each = n_alt)
        alts <- rep(levels(alts), n_id)
        .term <- ifelse(is.null(alt), covar, paste(covar, alt, sep = "_"))
        aslp <- data.frame(term = .term,
                           id = ids,
                           alt = alts,
                           effect = .effect,
                           estimate = .estimate)
        
    } else {
        aslp <- data.frame(term = .term,
                           effect = .effect,
                           estimate = .estimate)
    }
    if (se) aslp[["std.error"]] <- .std.er
    aslp
}


# get the fitted values for a vector of parameters, a model and a data
# frame
predict_est <- function(coefs, model, newdata = NULL, shape = c("long", "wide")){
    .shape <- match.arg(shape)
    if (is.null(newdata)){
        preds <- update(model, start = coefs, maxit = 0, iterlim = 0)
    } else {
        if (inherits(data, "dfidx_mlogit")){
            class(data) <- setdiff(class(data), "dfidx_mlogit")
        }
        preds <- update(model, start = coefs, maxit = 0, data = newdata, iterlim = 0)
    }
    if (inherits(model, "mlogit")){
        if (.shape == "wide"){
            result <- preds$probabilities
        } else {
            result <- as.numeric(t(preds$probabilities))
        }
    } else {
        result <- fitted(preds)
    }
}

# get the standard deviations of the predictions using the delta
# method
predict_se <- function(model, newdata = NULL){
    jacob <- numDeriv::jacobian(predict_est, model$coefficients, model = model, newdata = newdata)
    sqrt(apply(jacob, 1, function(x) x %*% vcov(model, subset = "all") %*% x))
}

# predict method for micsr objects, returns either a data frame with
# the estimate and its standard error or a numeric vector of
# predictions. The result is of class predict

#' @rdname micsr
#' @export
predict.micsr <- function(object, ..., se = TRUE,
                          newdata = NULL, shape = c("long", "wide")){
    .shape <- match.arg(shape)
    if (.shape == "wide") se <- FALSE
    .preds <- predict_est(object$coefficients, object,
                          newdata = newdata, shape = .shape)
    if (se){
        .se <- predict_se(object, newdata = newdata)
        result <- structure(data.frame(estimate = .preds, std.error = .se),
                            class = c("predict", "tbl_df", "tbl", "data.frame"))
    } else {
        result <- structure(.preds, class = c("predict", "numeric"))
    }
    if (inherits(object, "mlogit") & .shape == "long"){
        result <- dfidx(cbind(model.frame(object)$idx, result), position = 1)
        class(result) <- c("tbl_dfidx", "dfidx", "tbl_df", "tbl", "data.frame")
        
    }
    result
}
    
# compute the slope for one covariate using numerical derivative of
# the covariate
estimate <- function(coefs, model, covar, type, eps = NULL, level = NULL,
                     newdata, alt = NULL){
    dta <- newdata
    if (is.null(alt)){
        selalt <- rep(TRUE, nrow(dta))
    } else {
        alts <- dta[[idx_name(model, 2, 1)]]
        selalt <- (alts == alt)
    }
    if (type == "dbl"){
#        .sd <- sd(dta[[covar]], na.rm = TRUE)
        .sd <- sd(get_data(model)[[covar]], na.rm  = TRUE)
        probs_init <- predict_est(coefs, model, dta)
        if (is.null(eps))  eps <- .sd * sqrt(.Machine$double.eps)
        dta[[covar]] <- dta[[covar]] + eps * selalt
    }
    if (type == "fct"){
        .levels <- levels(dta[[covar]])
        .reflevel <- .levels[1]
        dta[[covar]][] <- .reflevel
        probs_init <- predict_est(coefs, model, dta)
        dta[[covar]][] <- level
    }
    if (type == "dummy"){
        dta[[covar]][selalt] <- 0L
        probs_init <- predict_est(coefs, model, dta)
        dta[[covar]][selalt] <- 1L
    }
    if (type == "count"){
        probs_init <- predict_est(coefs, model, dta)
        dta[[covar]][selalt] <- dta[[covar]][selalt] + 1
    }
    probs_fin <- predict_est(coefs, model, dta)
    (probs_fin - probs_init) / ifelse(type == "dbl", eps, 1)
}

# get the standard deviations of the slopes using the delta
std.er <- function(model, covar, type, level = NULL, newdata, alt = NULL){
    .coefs <- model$coefficients
    jacob <- numDeriv::jacobian(estimate, .coefs, model = model,
                                covar = covar, level = level, type = type,
                                newdata = newdata, alt = alt)
    sqrt(apply(jacob, 1, function(x) x %*% vcov(model, subset = "all") %*% x))
}

#' @rdname micsr
#' @importFrom stats effects get_all_vars
#' @importFrom dfidx dfidx idx idx_name
#' @method effects micsr
#' @export
effects.micsr <- function(object, ..., newdata = NULL, covariates = NULL, se = TRUE){
    is_tibble <- inherits(model.frame(object), "tbl_df")
    is_tibble <- TRUE
    slps <- data.frame()
    MEM <- is.character(newdata) && length(newdata) == 1 && newdata == "mem"
    if (is.null(covariates)){
        covariates <- get_covariates(object)
        if (inherits(object, "mlogit")){
            covariates_1 <- get_covariates(object, rhs = 1)
            if (length(formula(object))[2] > 1){
                covariates_2 <- get_covariates(object, rhs = 2)
            } else {
                covariates_2 <- character(0)
            }
            covariates <- c(covariates_1, covariates_2)
        }
    }
    .coefs <- object$coefficients
    # newdata is by default the data argument of object
    if (inherits(object, "mlogit")) alts <- factor(unique(idx(object, 2, 1)),
                                                   levels = levels(idx(object, 2, 1)))
    if (is.null(newdata)){
        newdata <- get_data(object)
    } else {
        if (MEM){
            if (inherits(object, "mlogit")){
                nm_id <- idx_name(model.frame(object), 1, 1)
                nm_alt <- idx_name(model.frame(object), 2, 1)
                newdata <- mean(dfidx(get_data(object)))
                newdata <- cbind(id = 1, alt = factor(alts), newdata)
                names(newdata)[1:2] <- c(nm_id, nm_alt)
#                newdata[[get_response(object)]] <- c(F, T)
            } else {
                newdata <- mean(object)
            }
            if (inherits(object, "ordreg")){
                z <- levels(factor(get_data(object)[[get_response(object)]]))
                newdata[[get_response(object)]] <- factor(z[1], levels = z)
            } else {
                newdata[[get_response(object)]] <- 1
            }
        }
    }
    for (i in covariates){
        .series <- newdata[[i]]
        .type <- covar_type(.series)
        if (.type == "fct"){
            .levels <- levels(.series)
            for (l in .levels[-1]){
                aslp <- add_a_slope(coefs = .coefs, object = object, covar = i,
                                    level = l, type = .type, newdata = newdata, se = se)
                slps <- rbind(slps, aslp)
            }
        } else {
            if (inherits(object, "mlogit") && i %in% covariates_1){
                for (a in alts){
                    aslp <- add_a_slope(coefs = .coefs, object = object, covar = i,
                                        type = .type, newdata = newdata, alt = a, se = se)
                    slps <- rbind(slps, aslp)
                }
            } else {
                aslp <- add_a_slope(coefs = .coefs, object = object, covar = i,
                                    type = .type, newdata = newdata, se = se)
                slps <- rbind(slps, aslp)
            }
        }
    }
    if (se){
        slps$statistic <- slps$estimate / slps$std.er
        slps$p.value <- 2 * pnorm(abs(slps$statistic), lower.tail = FALSE)
    }
    if (is_tibble){
        class(slps) <- c("effects", "tbl_df", "tbl", "data.frame")
    } else {
        class(slps) <- c("effects", "data.frame")
    }        
    slps
}

# summary for the effects objects, compute the average slopes
#' @rdname micsr
#' @method summary effects
#' @export
summary.effects <- function(object, ...){
    se <- "std.error" %in% names(object)
    N <- nrow(object) / length(unique(object$term))
    is_tibble <- inherits(object, "tbl_df")
    is_tibble <- TRUE
    if ("alt" %in% names(object)){
        z <- tapply(object$estimate, list(object$term, object$alt), mean)
        z <- as.data.frame(z)
        names(z) <- paste("alt", names(z), sep = ".")
        z$term <- rownames(z)
        z <- reshape(z, varying = 1:4, direction = "long")[, 1:3]
        names(z) <- c("term", "alt", "estimate")
        rownames(z) <- NULL
        if (se){
            s <- tapply(object$std.error, list(object$term, object$alt),
                        function(x) sqrt(sum(x ^ 2) / nrow(object)))
            s <- as.data.frame(s)
            names(s) <- paste("alt", names(s), sep = ".")
            s$term <- rownames(s)
            s <- reshape(s, varying = 1:4, direction = "long")[, 1:3]
            names(s) <- c("term", "alt", "std.error")
            rownames(s) <- NULL
            result <- merge(z, s)
        } else {
            result <- z
        }
    } else {
        z <- tapply(object$estimate, object$term, mean)
        if (se){
            s <- tapply(object$std.error, object$term,
                        function(x) sqrt(mean(x ^ 2)))
        }
        .effect_term <- unique(object[c("term", "effect")])
        .effect <- .effect_term$effect
        names(.effect) <- .effect_term$term
        result <- data.frame(term = rownames(z),
                             effect = .effect[rownames(z)],
                             estimate = unname(z))
        if (se) result[["std.error"]] <- s
    }
    if (se){
        result$statistic <- result$estimate / result$std.error
        result$p.value <- 2 * pnorm(abs(result$statistic), lower.tail = FALSE)
    }
    if (is_tibble){
        class(result) <- c("tbl_df", "tbl", "data.frame")
    } else {
        class(result) <- "data.frame"
    }
    result
}

# summary for the predict objects, compute the average prediction and
# its standard deviations
#' @rdname micsr
#' @method summary predict
#' @export
summary.predict <- function(object, ...){
    is_tibble <- inherits(object, "tbl_df")
    is_tibble <- TRUE
    cls <- class(object)
    class(object) <- setdiff(cls, c("dfidx", "preds"))
    e <- mean(object$estimate)
    s <- sqrt(mean(object$std.error ^ 2) / nrow(object))
    result <- data.frame(estimate = unname(e), std.error = unname(s))
    if (is_tibble) class(result) <- c("tbl_df", "tbl", "data.frame")
    result
}

#' @rdname micsr
#' @method mean micsr
#' @export
mean.micsr <- function(x, ...){
    dta <- get_covariates(x)
    dta <- get_data(x)[dta]
    result <- data.frame(lapply(dta, function(x){
        if (is.numeric(x)){
            result <- mean(x, na.rm = TRUE)
        }
        if (is_dummy(x)) result <- ifelse(sum(x) > sum(1 - x), 1L, 0L)
        if (is.character(x)) x <- factor(x, levels = unique(x))
        if (is.factor(x)){
            result <- factor(names(which.max(table(x))),
                             levels = levels(x))
        }
        result
    }))
    result
}
