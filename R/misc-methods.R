# Additional generic methods for survfitr objects.

#' Quantiles of the survival time
#'
#' Returns, for each requested probability `p`, the first time at which
#' the estimated survival function falls to `1 - p` or below (so
#' `probs = 0.5` is the median survival time). `NA` is returned when the
#' curve does not reach the requested level.
#'
#' @param x a `survfitr` object.
#' @param probs vector of probabilities.
#' @param ... ignored; kept for compatibility with the generic.
#'
#' @return A named vector of survival times, or a matrix with one row per
#'   group when the fit has strata.
#'
#' @seealso [q_search()], [survfitr()]
#'
#' @examples
#' data(MMC)
#' fit <- survfitr(Survr(id, time, event) ~ group,
#'   data = MMC,
#'   type = "wang-chang"
#' )
#' quantile(fit, probs = c(0.25, 0.5, 0.75))
#' @importFrom stats quantile
#' @export
quantile.survfitr <- function(x, probs = c(0.25, 0.5, 0.75), ...) {
  s <- strata_fits(x)
  qtime <- function(fit) {
    vapply(probs, function(p) {
      pos <- which(fit$survfunc <= 1 - p)
      if (length(pos)) fit$time[min(pos)] else NA_real_
    }, numeric(1))
  }
  ans <- t(vapply(s$fits, qtime, numeric(length(probs))))
  colnames(ans) <- paste0(format(100 * probs, trim = TRUE), "%")
  if (is.null(s$groups)) {
    ans <- drop(ans)
    if (length(probs) == 1) names(ans) <- paste0(100 * probs, "%")
  } else {
    rownames(ans) <- s$groups
  }
  ans
}

#' Coerce a survfitr object to a data frame
#'
#' Returns the estimated curves in long ("tidy") format, one row per
#' distinct event time (and group), ready for further processing or
#' plotting.
#'
#' @param x a `survfitr` object.
#' @param row.names,optional see [as.data.frame()]; ignored.
#' @param ... ignored.
#'
#' @return A data frame with columns `group` (only when the fit has
#'   strata), `time`, `n.event`, `n.risk`, `surv` and `std.error` (`NA`
#'   for estimators without standard errors).
#'
#' @examples
#' data(MMC)
#' fit <- survfitr(Survr(id, time, event) ~ group,
#'   data = MMC,
#'   type = "wang-chang"
#' )
#' head(as.data.frame(fit))
#' @export
as.data.frame.survfitr <- function(x, row.names = NULL, optional = FALSE,
                                   ...) {
  s <- strata_fits(x)
  one <- function(fit, group) {
    d <- data.frame(
      time = fit$time,
      n.event = as.numeric(fit$n.event),
      n.risk = if (is.matrix(fit$AtRisk)) {
        colSums(fit$AtRisk)
      } else {
        as.numeric(fit$AtRisk)
      },
      surv = fit$survfunc,
      std.error = if (is.null(fit$std.error)) {
        NA_real_
      } else {
        fit$std.error
      }
    )
    if (!is.null(group)) d <- cbind(group = group, d)
    d
  }
  if (is.null(s$groups)) {
    one(s$fits[[1]], group = NULL)
  } else {
    ans <- do.call(rbind, Map(one, s$fits, s$groups))
    rownames(ans) <- NULL
    ans
  }
}
