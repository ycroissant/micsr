#' Compute the mean of a variable for bins of another variable
#'
#' Plot average values of the outcome for bins of the forcing
#' variable, a common plot in regression discontinuity analysis
#'
#' @name binmeans
#' @param x,y either two numeric vector for the default method, or a
#'     formula and a data frame for the formula method
#' @param width the width of the bins
#' @param center the cuting point of the forcing variable
#' @param g a grouping variable
#' @param name_g internally used by the geom
#' @param ... further arguments
#' @param center the cuting value of the forcing variable
#' @param width the width of the bins
#' @importFrom dplyr bind_cols group_by left_join mutate n select summarise %>%
#' @importFrom rlang set_names `:=` .data
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#' @importFrom tibble as_tibble tibble
#' @importFrom ggplot2 ggproto Stat
#' @importFrom ggplot2 GeomPoint
#' @keywords plot 
#' @export
binmeans <- function(x, y, width = NULL, center = NULL, g = NULL, ...){
    UseMethod("binmeans")
}

#' @rdname binmeans
#' @export
binmeans.default <- function(x, y, width = NULL, center = NULL, g = NULL, ..., name_g = "colour"){
    .gps <- ! is.null(g)
    if (! .gps) g <- rep(1, length(x))
    .df <- tibble(x = x, y = y, g = g)
    .df <- compute_binmeans(.df, width, center)
    if (! .gps) .df <- .df %>% select( - 3)
    names(.df)[names(.df) == "g"] <- name_g
    .df
}

#' @rdname binmeans
#' @export
binmeans.formula <- function(x, y, width = NULL, center = NULL, ...){
    .formula <- Formula(x)
#    x <- update(x, . ~ . - 1)
    .df <- model.frame(.formula, y)
    .x <- model.part(.formula, .df, rhs = 1)
    .gps <- length(.formula)[2] == 2
    if (! .gps) .g <- data.frame(g = rep(1, nrow(.df)))
    else .g <- model.part(.formula, .df, rhs = 2)
    .y <- model.part(.formula, .df, lhs = 1)
    name_x <- names(.x)
    name_y <- names(.y)
    name_g <- names(.g)
    .df <- .x %>%
        bind_cols(.y) %>%
        bind_cols(.g) %>%
        as_tibble %>%
        set_names(c("x", "y", "g"))
    .df <- compute_binmeans(.df, width, center) %>% 
        set_names(c(name_x, name_y, name_g, "n"))
    if (! .gps) .df <- .df %>% select( - 3)
    .df
}

compute_binmeans <- function(x, width, center){
    .df <- x
    x_range <- range(.df$x)
    x_dist <- x_range[2] - x_range[1]
    if (is.null(width) | is.null(center)){
        default_width <- 10 ^ (floor(log10(x_dist)) - 1)
        if (is.null(center))
            center <- round(mean(.df$x) / default_width) * default_width
        if (is.null(width))
            width <- default_width
    }
    low <- ceiling((center - x_range[1]) / width)
    up <- ceiling((x_range[2] - center) / width)
    bks <- sort(c(center, center - (1:low) * width, center + (1:up) * width))
    bks <- round(bks, 5)
    pts <- (bks[-1] + bks[- length(bks)]) / 2
    .df <- .df %>%
        mutate(levs = cut(.data$x, breaks = bks))# %>%
    .levs <- tibble(levs = levels(.df$levs), pts)
#    left_join(tibble(levs = levels(.$levs), pts)) %>%
    .df <- .df %>% left_join(.levs, by = "levs") %>% 
        select(x = pts, "y", "g" ) %>%
        group_by(.data$g, .data$x) %>%
        summarise(y = mean(.data$y), n = n(), .groups = "drop") %>%
#    summarise("{{ y }}" := mean({{ y }}), "{{ n }}" = n()) %>% 
    select("x", "y", "g", "n")
    .df
}

#' @rdname binmeans
#' @export
StatBinmeans <-  ggplot2::ggproto("StatBinmeans", ggplot2::Stat,
                     ## setup_data = function(data, params) {
                     ##     if (anyDuplicated(data$group)) {
                     ##         data$group <- paste(data$group, seq_len(nrow(data)), sep = "-")
                     ##     }
                     ##     data
                     ## },
                         compute_panel = function(data,
                                                  scales,
                                                  width = NULL,
                                                  center = NULL) {
                             .shape <- ! is.null(data$shape)
                             .colour <- ! is.null(data$colour)
                             .fill <- ! is.null(data$fill)
                             .name <- case_when(! is.null(data$shape) ~ "shape",
                                                ! is.null(data$colour) ~ "colour",
                                                ! is.null(data$fill) ~ "fill")
                             binmeans(data$x, data$y, width, center, g = data[[.name]], name_g = .name)      
                         },
                         required_aes = c("x", "y"),
                         optional_aes = c("shape", "colour", "fill")
                         )

#' @rdname binmeans
#' @export
geom_binmeans <- function(mapping = NULL,
                          data = NULL, 
                          stat = "binmeans",
                          position = "identity", 
                          ..., 
                          center = NULL,
                          width = NULL, 
                          na.rm = FALSE, 
                          show.legend = NA, 
                          inherit.aes = TRUE
                          ) {
    ggplot2::layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = ggplot2::GeomPoint,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            center = center,
            width = width,
            na.rm = na.rm,
            ...
        )
    )
}



## va <- z %>%
##     filter(net_score >= -5, net_score <= 4) %>% 
##     ggplot(aes(net_score, participate, fill = gender)) +
##     geom_binmeans(aes(size = after_stat(n))) +
##     geom_smooth(aes(linetype = factor(above)), method = "lm", se = FALSE)
## print(va)

#bdh %>% mutate(eligible = ifelse(selben > 0, "no", "yes")) %>% ggplot(aes(selben, treated)) + geom_binmeans(width = .1, center = 0, shape = 21,  color = "black", fill = "red") + geom_smooth(aes(linetype = eligible), se = FALSE)


