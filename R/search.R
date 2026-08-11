#' Survival estimate at selected times
#'
#' Evaluates a step survival function at arbitrary times. The estimate is
#' a decreasing piecewise-constant function with jumps at the event times,
#' so the value at any time is the one of the previous event time.
#'
#' @param tvals vector of times at which the survival function is
#'   evaluated.
#' @param time vector of ordered distinct event times.
#' @param surv vector of survival estimates at each event time.
#'
#' @return The survival estimate at each of the `tvals` times.
#'
#' @examples
#' time <- c(4, 7, 9, 15, 21, 67)
#' surv <- c(0.8, 0.7, 0.65, 0.55, 0.43, 0.22)
#'
#' # survival at times 1, 10, 32 and 74
#' surv_search(c(1, 10, 32, 74), time, surv)
#' @keywords survival
#' @export
surv_search <- function(tvals, time, surv) {
  # step-function lookup: estimate at the largest event time <= each tval
  time.c <- c(0, time)
  surv.c <- c(1, surv)
  return(surv.c[findInterval(tvals, time.c)])
}

#' Survival time of a selected quantile
#'
#' Given a `survfitr` object, returns the first time at which the
#' estimated survival function falls to the quantile `q` or below.
#'
#' @param f a `survfitr` object.
#' @param q quantile. Default is 0.5 (the median survival time).
#'
#' @return The survival time of the selected quantile, or `NA` (with a
#'   warning) when the survival estimates are not available.
#'
#' @examples
#' data(MMC)
#' fit <- survfitr(Survr(id, time, event) ~ 1, data = MMC)
#'
#' # 75th percentile of the survival function
#' q_search(fit, q = 0.75)
#' @keywords survival
#' @export
q_search <- function(f, q = 0.5) {
  tt <- c(0, f$time)
  ss <- c(1, f$survfunc)
  if (anyNA(ss)) {
    warning("survival estimates are not available; returning NA")
    return(NA_real_)
  }
  if (ss[length(ss)] > q) {
    stop(paste("\noverall survival estimate does not fall below ", q))
  }
  ans <- min(tt[ss <= q])
  return(ans)
}
