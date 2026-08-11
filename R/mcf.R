#' Mean cumulative function of recurrent event data
#'
#' Nonparametric (Nelson-Aalen type) estimate of the mean cumulative
#' function, i.e. the expected number of events experienced by a subject
#' up to each calendar time, optionally by group.
#'
#' Calendar times are reconstructed as the within-subject cumulative sums
#' of the gap times of the [Survr()] response. At each distinct event
#' calendar time \eqn{t} the estimate increases by \eqn{dN(t)/Y(t)}, where
#' \eqn{dN(t)} is the number of events observed at \eqn{t} and \eqn{Y(t)}
#' the number of subjects still under observation (total observation time
#' \eqn{\ge t}).
#'
#' @param formula a formula object with a [Survr()] object as the response
#'   on the left of the `~` operator and a grouping term (or `1`) on the
#'   right.
#' @param data a data frame in which to interpret the variables named in
#'   the formula.
#'
#' @return An object of class `mcf`: a data frame with columns `group`,
#'   `time`, `n.event`, `n.risk` and `mcf`, with `print` and `autoplot`
#'   methods.
#'
#' @references
#' Lawless, J.F. and Nadeau, C. (1995). Some Simple Robust Methods for
#' the Analysis of Recurrent Events. *Technometrics* **37**, 158--168.
#'
#' @seealso [autoplot.mcf()], [survfitr()]
#'
#' @examples
#' data(MMC)
#' m <- mcf(Survr(id, time, event) ~ group, data = MMC)
#' m
#' autoplot(m)
#' @keywords survival
#' @export
mcf <- function(formula, data) {
  m <- match.call(expand.dots = FALSE)
  Terms <- terms(formula, "strata")
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  Y <- model.extract(m, "response")
  if (!is.Survr(Y)) {
    stop("Response must be a survival recurrent object")
  }
  ll <- attr(Terms, "term.labels")

  # MCF of one group from its Survr rows
  mcf_one <- function(Sr, group) {
    fid <- factor(Sr[, 1], levels = unique(Sr[, 1]))
    caltime <- ave(Sr[, 2], fid, FUN = cumsum)
    tau <- as.vector(rowsum(Sr[, 2], fid))

    ev <- sort(caltime[Sr[, 3] == 1])
    time <- unique(ev)
    n.event <- as.vector(table(factor(ev, levels = time)))
    # subjects still under observation at each event time
    n.risk <- length(tau) - findInterval(time, sort(tau), left.open = TRUE)
    data.frame(
      group = group, time = time, n.event = n.event,
      n.risk = n.risk, mcf = cumsum(n.event / n.risk)
    )
  }

  if (ncol(m) > 1) {
    group <- m[ll][, 1]
    k <- levels(group)
    ans <- do.call(rbind, lapply(k, function(g) {
      mcf_one(Y[group == g, , drop = FALSE], g)
    }))
    attr(ans, "group") <- ll
  } else {
    ans <- mcf_one(Y, "all")
  }
  oldClass(ans) <- c("mcf", "data.frame")
  ans
}

#' @rdname mcf
#' @param x an `mcf` object.
#' @param ... further arguments passed to [print.data.frame()].
#' @export
print.mcf <- function(x, ...) {
  cat("Mean cumulative function for recurrent event data\n\n")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}

#' Plot a mean cumulative function
#'
#' Draws the mean cumulative function estimated by [mcf()] as a step
#' function, one colour per group.
#'
#' @param object an object of class `mcf`.
#' @param ... ignored; kept for compatibility with the generic.
#'
#' @return A `ggplot` object.
#'
#' @seealso [mcf()], [theme_survrec()]
#'
#' @examples
#' data(colon)
#' autoplot(mcf(Survr(hc, time, event) ~ as.factor(dukes), data = colon))
#' @export
autoplot.mcf <- function(object, ...) {
  df <- as.data.frame(object)
  # start every curve at (0, 0)
  origin <- data.frame(
    group = unique(df$group), time = 0, n.event = 0,
    n.risk = NA, mcf = 0
  )
  df <- rbind(origin, df)
  df <- df[order(df$group, df$time), ]
  names(df)[names(df) == "mcf"] <- "y"
  st <- do.call(rbind, lapply(split(df, df$group), step_frame))

  many <- length(unique(df$group)) > 1
  p <- ggplot2::ggplot(
    st,
    ggplot2::aes(x = .data$time, y = .data$y, colour = .data$group)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(
      values = survrec_palette(length(unique(df$group)))
    ) +
    ggplot2::labs(
      x = "Calendar time", y = "Mean cumulative number of events",
      colour = attr(object, "group")
    ) +
    theme_survrec()
  if (!many) {
    p <- p + ggplot2::guides(colour = "none")
  }
  p
}
