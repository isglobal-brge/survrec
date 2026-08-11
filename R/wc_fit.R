#' Wang-Chang survival estimator
#'
#' Estimates the survival function of correlated or i.i.d. inter-occurrence
#' times of recurrent event data with the product-limit estimator of Wang
#' and Chang (1999).
#'
#' Wang and Chang (1999) proposed an estimator of the common marginal
#' survivor function for the case where within-unit inter-occurrence times
#' are correlated. The correlation structure they considered is quite
#' general and contains both the i.i.d. and the multiplicative (hence
#' gamma) frailty model as special cases.
#'
#' This estimator removes the bias noted for the product-limit estimator of
#' Peña, Strawderman and Hollander (2001) when inter-occurrence times are
#' correlated within units. However, when applied to i.i.d.
#' inter-occurrence times it is not expected to perform as well as the PSH
#' estimator, especially with regard to efficiency.
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
#'   \item{n.event}{weighted number of events at each distinct time.}
#'   \item{AtRisk}{weighted number of units at risk at each distinct
#'     time.}
#'   \item{survfunc}{survival estimate at each distinct time.}
#'   \item{std.error}{standard error of the survival estimate.}
#'   \item{tvals}{copy of the `tvals` argument.}
#'   \item{WCpleAttvals}{survival estimate at the `tvals` times.}
#' }
#'
#' @references
#' Wang, M.-C. and Chang, S.-H. (1999). Nonparametric Estimation of a
#' Recurrent Survival Function. *Journal of the American Statistical
#' Association* **94**, 146--153.
#'
#' @note The maintainers wish to thank Professors Chiung-Yu Huang and
#'   Shu-Hui Chang for providing the original Fortran code computing the
#'   standard errors of Wang and Chang's estimator.
#'
#' @seealso [survfitr()], [Survr()]
#'
#' @examples
#' data(MMC)
#' fit <- wc_fit(Survr(MMC$id, MMC$time, MMC$event))
#' fit
#' plot(fit, conf.int = FALSE)
#'
#' # compare with Peña-Strawderman-Hollander
#' fit <- psh_fit(Survr(MMC$id, MMC$time, MMC$event))
#' lines(fit, lty = 2)
#'
#' # and with the gamma frailty MLE
#' fit <- mlefrailty_fit(Survr(MMC$id, MMC$time, MMC$event))
#' lines(fit, lty = 3)
#' @keywords survival
#' @export
wc_fit <- function(x, tvals) {
  if (!is.Survr(x)) {
    stop("\n x must be a Survr object")
  }
  d <- survr_pieces(x)
  n <- d$n
  failed <- d$failed
  cen.gap <- d$censored
  m <- d$m
  distinct <- sort(unique(failed))

  summ <- .wc_fit(as.integer(m), as.double(failed), as.double(cen.gap),
    as.double(distinct),
    variance = TRUE
  )
  survfuncWCple <- summ$surv

  if (!missing(tvals)) {
    WCpleAttvals <- surv_search(sort(tvals), distinct, survfuncWCple)
  } else {
    tvals <- NA
    WCpleAttvals <- NA
  }
  ans <- list(
    n = n, m = m, failed = failed, censored = cen.gap,
    time = distinct, n.event = summ$n.event, AtRisk = summ$at.risk,
    survfunc = survfuncWCple, std.error = summ$std.error, tvals = tvals,
    WCpleAttvals = WCpleAttvals
  )
  oldClass(ans) <- "survfitr"
  ans
}
