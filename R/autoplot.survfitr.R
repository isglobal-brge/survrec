# ggplot2 graphics for survfitr objects. The base-graphics plot() method
# is kept for backwards compatibility; autoplot() is the modern interface.

# Long-format data of a survfitr fit (or one per group) with pointwise
# confidence limits on the survival scale, prepended with the (0, 1)
# origin so the curves start at time zero.
survfitr_frame <- function(x, level, conf.type) {
  s <- strata_fits(x)
  one <- function(fit, group) {
    se <- if (is.null(fit$std.error)) rep(NA_real_, length(fit$time)) else fit$std.error
    ci <- surv_ci(fit$survfunc, se, level = level, conf.type = conf.type)
    data.frame(
      group = group,
      time = c(0, fit$time),
      surv = c(1, fit$survfunc),
      lower = c(1, ci$lower),
      upper = c(1, ci$upper)
    )
  }
  if (is.null(s$groups)) {
    one(s$fits[[1]], group = "all")
  } else {
    do.call(rbind, Map(one, s$fits, s$groups))
  }
}

#' Plot survival curves with ggplot2
#'
#' Draws the estimated survival function of a `survfitr` object as a step
#' function, one colour per group, with an optional pointwise confidence
#' band. The cumulative distribution (`fun = "event"`) and the cumulative
#' hazard (`fun = "cumhaz"`) transformations are also available.
#'
#' The confidence band is computed from the standard errors of the fit
#' when the estimator provides them (Peña-Strawderman-Hollander and
#' Wang-Chang). With `conf.type = "log-log"` (default) the limits are
#' obtained on the \eqn{\log(-\log S)} scale, so they always stay inside
#' \eqn{[0, 1]}; `"plain"` gives the symmetric limits of the base
#' [plot.survfitr()] method.
#'
#' @param object an object of class `survfitr` (output of [survfitr()],
#'   [psh_fit()], [wc_fit()] or [mlefrailty_fit()]).
#' @param fun transformation of the survival curve: `"surv"` (default,
#'   the survival function), `"event"` (the cumulative probability of an
#'   event, 1 - S) or `"cumhaz"` (the cumulative hazard, -log S).
#' @param conf.int draw the pointwise confidence band? Default `TRUE`
#'   (ignored for estimators without standard errors).
#' @param level confidence level of the band.
#' @param conf.type transformation used for the band: `"log-log"`
#'   (default) or `"plain"`.
#' @param ... ignored; kept for compatibility with the generic.
#'
#' @return A `ggplot` object, which can be further styled with the usual
#'   `ggplot2` syntax.
#'
#' @seealso [plot.survfitr()] for the base-graphics version,
#'   [plotEstimators()] to compare the three estimators,
#'   [theme_survrec()]
#'
#' @examples
#' data(MMC)
#' fit <- survfitr(Survr(id, time, event) ~ group,
#'   data = MMC,
#'   type = "wang-chang"
#' )
#' autoplot(fit)
#' autoplot(fit, fun = "cumhaz", conf.int = FALSE)
#'
#' # it is a ggplot, so it can be restyled
#' autoplot(fit) + ggplot2::labs(title = "Migratory Motor Complex")
#' @export
autoplot.survfitr <- function(object, fun = c("surv", "event", "cumhaz"),
                              conf.int = TRUE, level = 0.95,
                              conf.type = c("log-log", "plain"), ...) {
  fun <- match.arg(fun)
  conf.type <- match.arg(conf.type)

  df <- survfitr_frame(object, level = level, conf.type = conf.type)
  # transform estimate and band; "event" flips the band order
  trans <- switch(fun,
    surv = identity,
    event = function(s) 1 - s,
    cumhaz = function(s) -log(s)
  )
  y <- trans(df$surv)
  lo <- trans(if (fun == "event") df$upper else df$lower)
  hi <- trans(if (fun == "event") df$lower else df$upper)
  df$y <- y
  df$lower <- lo
  df$upper <- hi
  ylab <- switch(fun,
    surv = "Survival probability",
    event = "Cumulative probability of event",
    cumhaz = "Cumulative hazard"
  )

  many <- length(unique(df$group)) > 1
  st <- do.call(rbind, lapply(split(df, df$group), step_frame))
  has.band <- conf.int && !all(is.na(st$lower))

  p <- ggplot2::ggplot(st, ggplot2::aes(x = .data$time))
  if (has.band) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$lower, ymax = .data$upper,
        fill = .data$group
      ),
      alpha = 0.15, na.rm = TRUE
    )
  }
  p <- p + ggplot2::geom_line(
    ggplot2::aes(y = .data$y, colour = .data$group),
    linewidth = 0.8
  )
  npal <- length(unique(df$group))
  p <- p +
    ggplot2::scale_colour_manual(values = survrec_palette(npal)) +
    ggplot2::scale_fill_manual(values = survrec_palette(npal)) +
    ggplot2::labs(
      x = "Time", y = ylab,
      colour = attr(object, "group"), fill = attr(object, "group"),
      caption = if (has.band) {
        paste0(level * 100, "% pointwise confidence band (", conf.type, ")")
      }
    ) +
    theme_survrec()
  if (!many) {
    p <- p + ggplot2::guides(colour = "none", fill = "none")
  }
  p
}

#' Compare the three survival estimators on the same data
#'
#' Fits the Peña-Strawderman-Hollander, the Wang-Chang and the gamma
#' frailty MLE estimators to the same sample and draws the three curves
#' on the same axes. Agreement between the PSH and the frailty estimates
#' suggests independent inter-occurrence times, while a Wang-Chang curve
#' separated from the PSH one points to within-subject correlation.
#'
#' @param x a survival recurrent event object, see [Survr()].
#' @param conf.int draw the pointwise confidence bands of the estimators
#'   that provide them? Default `FALSE`.
#' @param level confidence level of the bands.
#' @param alpha.console passed to [mlefrailty_fit()]; default `FALSE`.
#' @param ... additional arguments passed to the three fitting functions.
#'
#' @return A `ggplot` object.
#'
#' @seealso [autoplot.survfitr()], [psh_fit()], [wc_fit()],
#'   [mlefrailty_fit()]
#'
#' @examples
#' data(MMC)
#' plotEstimators(Survr(MMC$id, MMC$time, MMC$event))
#' @export
plotEstimators <- function(x, conf.int = FALSE, level = 0.95,
                           alpha.console = FALSE, ...) {
  if (!is.Survr(x)) {
    stop("\n x must be a Survr object")
  }
  fits <- list(
    "PSH" = psh_fit(x, ...),
    "Wang-Chang" = wc_fit(x, ...),
    "MLE frailty" = mlefrailty_fit(x, alpha.console = alpha.console, ...)
  )
  ans <- Map(function(fit, name) {
    se <- if (is.null(fit$std.error)) rep(NA_real_, length(fit$time)) else fit$std.error
    ci <- surv_ci(fit$survfunc, se, level = level)
    data.frame(
      group = name,
      time = c(0, fit$time),
      y = c(1, fit$survfunc),
      lower = c(1, ci$lower),
      upper = c(1, ci$upper)
    )
  }, fits, names(fits))
  df <- do.call(rbind, ans)
  df$group <- factor(df$group, levels = names(fits))
  st <- do.call(rbind, lapply(split(df, df$group), step_frame))

  p <- ggplot2::ggplot(st, ggplot2::aes(x = .data$time))
  if (conf.int) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$lower, ymax = .data$upper,
        fill = .data$group
      ),
      alpha = 0.15, na.rm = TRUE
    )
  }
  p +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$y, colour = .data$group),
      linewidth = 0.8
    ) +
    ggplot2::scale_colour_manual(values = survrec_palette(3)) +
    ggplot2::scale_fill_manual(values = survrec_palette(3)) +
    ggplot2::labs(
      x = "Time", y = "Survival probability",
      colour = "Estimator", fill = "Estimator"
    ) +
    theme_survrec()
}
