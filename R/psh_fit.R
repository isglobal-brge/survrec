#' Peña-Strawderman-Hollander survival estimator
#'
#' Estimates the survival function of the inter-occurrence times of
#' recurrent event data with the generalized product-limit estimator (PLE)
#' of Peña, Strawderman and Hollander (2001).
#'
#' The estimator computed by this function is the nonparametric estimator
#' of the inter-event time survivor function under the assumption of a
#' renewal or IID model. It generalizes the product-limit estimator to the
#' situation where the event is recurrent. For details and the theory
#' behind the estimator, please refer to Peña, Strawderman and Hollander
#' (2001).
#'
#' @param x a survival recurrent event object, see [Survr()].
#' @param tvals optional vector of times at which the survival function is
#'   also evaluated.
#'
#' @return A list with class `survfitr` containing:
#' \describe{
#'   \item{n}{number of units or subjects observed.}
#'   \item{m}{number of events of each subject.}
#'   \item{failed}{uncensored gap times, ordered by subject.}
#'   \item{censored}{censored gap time of each subject.}
#'   \item{time}{ordered distinct event times.}
#'   \item{n.event}{number of events at each distinct time.}
#'   \item{AtRisk}{matrix of units at risk at each distinct time (rows are
#'     subjects).}
#'   \item{survfunc}{survival estimate at each distinct time.}
#'   \item{std.error}{standard error of the survival estimate.}
#'   \item{tvals}{copy of the `tvals` argument.}
#'   \item{PSHpleAttvals}{survival estimate at the `tvals` times.}
#' }
#'
#' @references
#' Peña, E.A., Strawderman, R. and Hollander, M. (2001). Nonparametric
#' Estimation with Recurrent Event Data. *Journal of the American
#' Statistical Association* **96**, 1299--1315.
#'
#' @seealso [survfitr()], [Survr()]
#'
#' @examples
#' data(MMC)
#' fit <- psh_fit(Survr(MMC$id, MMC$time, MMC$event))
#' fit
#' plot(fit, conf.int = FALSE)
#'
#' # compare with the gamma frailty MLE
#' fit <- mlefrailty_fit(Survr(MMC$id, MMC$time, MMC$event))
#' lines(fit, lty = 2)
#'
#' # and with Wang-Chang
#' fit <- wc_fit(Survr(MMC$id, MMC$time, MMC$event))
#' lines(fit, lty = 3)
#' @keywords survival
#' @export
psh_fit <- function(x, tvals) {
  if (!is.Survr(x)) {
    stop("\n x must be a Survr object")
  }
  d <- survr_pieces(x)
  n <- d$n
  failed <- d$failed
  censored <- d$censored
  m <- d$m

  summ <- .distinct_failed(
    as.integer(m), as.double(failed),
    as.double(censored)
  )
  distinct <- summ$time
  numdeaths <- summ$n.event
  AtRisk <- summ$at.risk

  AtRiskTotals <- colSums(AtRisk)
  survfuncPSHple <- cumprod(1 - numdeaths / AtRiskTotals)

  # Nelson-Aalen-based standard error of the product-limit estimate
  se.NA <- cumsum(numdeaths / AtRiskTotals^2)
  se.PLE <- sqrt(se.NA * survfuncPSHple^2)

  if (!missing(tvals)) {
    PSHpleAttvals <- surv_search(sort(tvals), distinct, survfuncPSHple)
  } else {
    tvals <- NA
    PSHpleAttvals <- NA
  }
  ans <- list(
    n = n, m = m, failed = failed, censored = censored,
    time = distinct, n.event = numdeaths, AtRisk = AtRisk,
    survfunc = survfuncPSHple, std.error = se.PLE, tvals = tvals,
    PSHpleAttvals = PSHpleAttvals
  )
  oldClass(ans) <- "survfitr"
  ans
}
