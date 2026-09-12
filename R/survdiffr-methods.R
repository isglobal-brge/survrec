# Methods for grouped survdiffr results (class "survdiffr": a named list
# of boot objects, one per group).

# Replicates of one group on the extended real line: -1 codes "the
# survival estimate never fell below q", i.e. a quantile beyond the
# observed range, which is treated as +Inf rather than dropped.
survdiffr_reps <- function(b) {
  t <- as.vector(b$t)
  t[!is.na(t) & t < 0] <- Inf
  t
}

#' Summarize a grouped bootstrap of survival quantiles
#'
#' Computes, for each group, the observed quantile and a percentile
#' bootstrap confidence interval, and for each pair of groups the
#' difference of quantiles with its percentile interval and a two-sided
#' bootstrap p-value for the null hypothesis of no difference.
#'
#' Replicates in which the survival estimate never fell below the
#' requested level are treated as `+Inf` (the quantile exceeds the
#' observed follow-up); pairs where both replicates are infinite are
#' dropped from the difference. The p-value is
#' \eqn{2 \min(P(D \le 0), P(D \ge 0))} over the replicated differences
#' \eqn{D}.
#'
#' @param object a grouped result of [survdiffr()].
#' @param level confidence level of the percentile intervals.
#' @param ... ignored; kept for compatibility with the generic.
#'
#' @return A list of class `summary.survdiffr` with components `groups`
#'   (data frame with the observed quantile and its interval per group)
#'   and `contrasts` (data frame with the pairwise differences, intervals
#'   and p-values).
#'
#' @seealso [survdiffr()], [boot::boot.ci()]
#'
#' @examples
#' data(colon)
#' fit <- survdiffr(Survr(hc, time, event) ~ as.factor(dukes),
#'   data = colon, q = 0.5, seed = 1
#' )
#' summary(fit)
#' @export
summary.survdiffr <- function(object, level = 0.95, ...) {
  a <- (1 - level) / 2
  k <- names(object)

  groups <- do.call(rbind, lapply(k, function(g) {
    t <- survdiffr_reps(object[[g]])
    qq <- quantile(t, c(a, 1 - a), na.rm = TRUE, names = FALSE)
    data.frame(
      group = g, estimate = object[[g]]$t0, lower = qq[1],
      upper = qq[2], replicates = sum(!is.na(t))
    )
  }))

  pairs <- utils::combn(k, 2, simplify = FALSE)
  contrasts <- do.call(rbind, lapply(pairs, function(p) {
    d <- survdiffr_reps(object[[p[1]]]) - survdiffr_reps(object[[p[2]]])
    d <- d[!is.na(d) & !is.nan(d)] # NaN: both replicates were infinite
    qq <- quantile(d, c(a, 1 - a), names = FALSE)
    data.frame(
      contrast = paste(p[1], "-", p[2]),
      estimate = object[[p[1]]]$t0 - object[[p[2]]]$t0,
      lower = qq[1], upper = qq[2],
      p.value = min(2 * min(mean(d <= 0), mean(d >= 0)), 1),
      replicates = length(d)
    )
  }))

  ans <- list(groups = groups, contrasts = contrasts, level = level)
  oldClass(ans) <- "summary.survdiffr"
  ans
}

#' @rdname summary.survdiffr
#' @param x a `summary.survdiffr` object.
#' @param digits number of digits to print.
#' @export
print.summary.survdiffr <- function(x, digits = max(options()$digits - 4, 3),
                                    ...) {
  cat(
    "Bootstrap comparison of survival quantiles (",
    x$level * 100, "% percentile intervals)\n\n",
    sep = ""
  )
  print.data.frame(x$groups, row.names = FALSE, digits = digits)
  cat("\nPairwise differences:\n")
  print.data.frame(x$contrasts, row.names = FALSE, digits = digits)
  invisible(x)
}

#' Plot the bootstrap distributions of survival quantiles
#'
#' Draws the density of the bootstrap replicates of each group, with a
#' vertical line at the observed quantile. Infinite replicates (survival
#' estimate never below the requested level) are excluded from the
#' densities and reported in the caption.
#'
#' @param object a grouped result of [survdiffr()].
#' @param ... ignored; kept for compatibility with the generic.
#'
#' @return A `ggplot` object.
#'
#' @seealso [survdiffr()], [summary.survdiffr()], [theme_survrec()]
#'
#' @examples
#' data(colon)
#' fit <- survdiffr(Survr(hc, time, event) ~ as.factor(dukes),
#'   data = colon, q = 0.5, seed = 1
#' )
#' autoplot(fit)
#' @export
autoplot.survdiffr <- function(object, ...) {
  k <- names(object)
  df <- do.call(rbind, lapply(k, function(g) {
    data.frame(group = g, t = survdiffr_reps(object[[g]]))
  }))
  finite <- is.finite(df$t)
  dropped <- sum(!finite, na.rm = TRUE) + sum(is.na(df$t))
  obs <- data.frame(
    group = k,
    t0 = vapply(object, function(b) b$t0, numeric(1))
  )

  ggplot2::ggplot(
    df[finite & !is.na(df$t), ],
    ggplot2::aes(x = .data$t, colour = .data$group, fill = .data$group)
  ) +
    ggplot2::geom_density(alpha = 0.15, na.rm = TRUE) +
    ggplot2::geom_vline(
      data = obs,
      ggplot2::aes(xintercept = .data$t0, colour = .data$group),
      linetype = 2
    ) +
    ggplot2::scale_colour_manual(values = survrec_palette(length(k))) +
    ggplot2::scale_fill_manual(values = survrec_palette(length(k))) +
    ggplot2::labs(
      x = "Bootstrap quantile of the survival time", y = "Density",
      colour = "Group", fill = "Group",
      caption = if (dropped > 0) {
        paste(
          dropped,
          "replicates beyond the observed follow-up were excluded"
        )
      }
    ) +
    theme_survrec()
}
