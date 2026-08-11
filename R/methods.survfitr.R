# S3 methods for survfitr objects: plot, lines, print, summary.

# Coordinates of the step function through (x, y), dropping repeated y
# values so consecutive segments merge.
dostep <- function(x, y) {
  if (is.na(x[1] + y[1])) {
    x <- x[-1]
    y <- y[-1]
  }
  n <- length(x)
  if (n > 2) {
    dupy <- c(TRUE, diff(y[-n]) != 0, TRUE)
    n2 <- sum(dupy)
    xrep <- rep(x[dupy], c(1, rep(2, n2 - 1)))
    yrep <- rep(y[dupy], c(rep(2, n2 - 1), 1))
    list(x = xrep, y = yrep)
  } else if (n == 1) {
    list(x = x, y = y)
  } else {
    list(x = x[c(1, 2, 2)], y = y[c(1, 1, 2)])
  }
}

# Draws one survival (or probability) curve with optional pointwise
# normal confidence bands.
plot_one_curve <- function(fit, prob, conf.int, col = 1, add = FALSE, ...) {
  y <- if (prob) 1 - fit$survfunc else fit$survfunc
  if (add) {
    lines(dostep(fit$time, y), col = col)
  } else {
    plot(dostep(fit$time, y), type = "l", col = col, ...)
  }
  if (conf.int && !is.null(fit$std.error)) {
    e <- qnorm(0.975) * fit$std.error
    lines(dostep(fit$time, y + e), col = col, lty = 2)
    lines(dostep(fit$time, y - e), col = col, lty = 2)
  }
  invisible()
}

#' Plot survival curves of recurrent event data
#'
#' Plots the estimated survival (or probability) function from a
#' `survfitr` object, one curve per group when the fit has strata.
#' Additional fits can be added to the same axes with the `lines` method.
#'
#' @param x an object of class `survfitr` (output of [survfitr()],
#'   [psh_fit()], [wc_fit()] or [mlefrailty_fit()]).
#' @param conf.int if `TRUE` (default), pointwise 95\% confidence bands
#'   are drawn for the estimators that provide standard errors.
#' @param prob if `TRUE` the probability function (1 - survival) is drawn
#'   instead of the survival function.
#' @param ... additional arguments passed to [plot()].
#'
#' @return No return value; called for its side effect.
#'
#' @seealso [psh_fit()], [wc_fit()], [mlefrailty_fit()]
#'
#' @examples
#' data(MMC)
#' fit <- survfitr(Survr(id, time, event) ~ group,
#'   data = MMC,
#'   type = "wang-chang"
#' )
#' plot(fit)
#' @keywords survival
#' @export
plot.survfitr <- function(x, conf.int = TRUE, prob = FALSE, ...) {
  y.lab <- ifelse(prob, "Probability Estimates",
    "Survivor Probability Estimates"
  )
  if (!is.null(attr(x, "strata"))) {
    y1 <- if (prob) 1 - x[[1]]$survfunc else x[[1]]$survfunc
    plot(dostep(x[[1]]$time, y1),
      type = "n", xlab = "Time",
      ylab = y.lab, ...
    )
    for (i in 1:attr(x, "strata")) {
      plot_one_curve(x[[i]], prob, conf.int, col = i, add = TRUE)
    }
  } else {
    plot_one_curve(x, prob, conf.int,
      xlab = "Time", ylab = y.lab,
      ylim = c(0, max(x$survfunc)), ...
    )
  }
  return(invisible())
}

#' @rdname plot.survfitr
#' @export
lines.survfitr <- function(x, prob = FALSE, ...) {
  if (!prob) {
    lines(dostep(x$time, x$survfunc), ...)
  } else {
    lines(dostep(x$time, 1 - x$survfunc), ...)
  }
  return(invisible())
}

#' Print a short summary of survival curves for recurrent event data
#'
#' Prints, for each curve, the number of subjects, the number of events,
#' the restricted mean survival and its standard error, the median
#' survival and the minimum, maximum and median number of recurrences per
#' subject.
#'
#' The restricted mean and its standard error are based on a truncated
#' estimator: if the survival curve does not reach zero, the reported
#' quantity is the mean survival restricted to the time before the last
#' censoring. The median is the time at which the survival curve crosses
#' 0.5.
#'
#' @param x the result of a call to [survfitr()], [psh_fit()], [wc_fit()]
#'   or [mlefrailty_fit()].
#' @param scale a numeric value to rescale the survival time, e.g. if the
#'   input data are in days, `scale = 365` scales the printout to years.
#' @param digits number of digits to print.
#' @param ... other unused arguments.
#'
#' @return `x`, invisibly.
#'
#' @seealso [summary.survfitr()], [survfitr()]
#'
#' @examples
#' data(MMC)
#' fit <- survfitr(Survr(id, time, event) ~ group, data = MMC)
#' print(fit)
#' @keywords survival
#' @export
print.survfitr <- function(x, scale = 1,
                           digits = max(options()$digits - 4, 3), ...) {
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  plab <- c(
    "n", "events", "mean", "se(mean)", "median",
    "recurrences: min", "max", "median"
  )

  # summary statistics (restricted mean, its s.e. and median) of one fit
  pfun <- function(x) {
    minmin <- function(y, xx) {
      if (any(!is.na(y) & y == 0.5)) {
        if (any(!is.na(y) & y < 0.5)) {
          0.5 * (min(xx[!is.na(y) & y == 0.5]) +
            min(xx[!is.na(y) & y < 0.5]))
        } else {
          0.5 * (min(xx[!is.na(y) & y == 0.5]) +
            max(xx[!is.na(y) & y == 0.5]))
        }
      } else {
        min(xx[!is.na(y) & y <= 0.5])
      }
    }

    stime <- x$time / scale
    n <- length(stime)
    if (is.matrix(x$AtRisk)) {
      n.risk <- colSums(x$AtRisk)
    } else {
      n.risk <- x$AtRisk
    }
    hh <- c(
      x$n.event[-n] / (n.risk[-n] * (n.risk[-n] - x$n.event[-n])),
      0
    )
    med <- minmin(x$survfunc, x$time)

    dif.time <- c(diff(c(0, stime)), 0)
    mean <- dif.time * c(1, x$survfunc)
    varmean <- sum(rev(cumsum(rev(mean))^2)[-1] * hh)

    c(
      x$n, sum(x$m), sum(mean), sqrt(varmean), med, min(x$m),
      max(x$m), median(x$m)
    )
  }

  if (is.null(attr(x, "strata"))) {
    x1 <- rbind(pfun(x))
    cat("Survival for recurrent event data")
    cat("\n")
    dimnames(x1) <- list(" ", plab)
    print(x1)
    cat("\n")
  } else {
    cat("Survival for recurrent event data. Group=", attr(x, "group"))
    cat("\n")
    x1 <- NULL
    for (i in 1:attr(x, "strata")) {
      x1 <- rbind(x1, pfun(x[[i]]))
    }
    dimnames(x1) <- list(names(x), plab)
    print(x1)
    cat("\n")
  }
  invisible(x)
}

#' Summary of survival curves for recurrent event data
#'
#' Returns a matrix with the distinct event times and, at each of them,
#' the number of events, the number of subjects at risk, the survival
#' estimate and (when available) its standard error. If the fit has
#' multiple curves, a list with one matrix per curve is returned.
#'
#' @param object output of a call to [survfitr()], [psh_fit()], [wc_fit()]
#'   or [mlefrailty_fit()].
#' @param ... other unused arguments.
#'
#' @return For a single curve, a matrix; for multiple curves, a list of
#'   class `summary.survfitr` with one matrix per curve.
#'
#' @seealso [survfitr()]
#'
#' @examples
#' data(MMC)
#' summary(survfitr(Survr(id, time, event) ~ group, data = MMC))
#' @keywords survival
#' @export
summary.survfitr <- function(object, ...) {
  x <- object
  if (!inherits(x, "survfitr")) {
    stop("Invalid data")
  }

  # one row per distinct time; std.error column only when available
  tabulate_one <- function(fit, has.se) {
    n.risk <- if (is.matrix(fit$AtRisk)) {
      colSums(fit$AtRisk)
    } else {
      fit$AtRisk
    }
    temp <- cbind(
      fit$time, fit$n.event, n.risk, fit$survfunc,
      if (has.se) fit$std.error
    )
    plab <- c(
      "time", "n.event", "n.risk", "surv",
      if (has.se) "std.error"
    )
    dimnames(temp) <- list(rep("", nrow(temp)), plab)
    temp
  }

  if (!is.null(attr(x, "strata"))) {
    has.se <- !is.null(x[[1]]$std.error)
    ans <- lapply(unclass(x)[1:attr(x, "strata")], tabulate_one,
      has.se = has.se
    )
    names(ans) <- names(x)
    oldClass(ans) <- "summary.survfitr"
    attr(ans, "strata") <- attr(x, "strata")
  } else {
    ans <- tabulate_one(x, has.se = !is.null(x$std.error))
  }
  ans
}

#' @rdname summary.survfitr
#' @param x a `summary.survfitr` object.
#' @param scale a numeric value to rescale the survival time.
#' @param digits number of digits to print.
#' @export
print.summary.survfitr <- function(x, scale = 1,
                                   digits = max(options()$digits - 4, 3),
                                   ...) {
  savedig <- options(digits = digits)
  on.exit(options(savedig))

  if (is.null(attr(x, "strata"))) {
    cat("\n")
    temp <- x
    oldClass(temp) <- NULL
    print(temp)
    cat("\n")
  } else {
    for (i in 1:attr(x, "strata")) {
      cat("\n      Group=", names(x)[i])
      cat("\n")
      print(x[[i]])
      cat("\n")
    }
  }
  invisible(x)
}
