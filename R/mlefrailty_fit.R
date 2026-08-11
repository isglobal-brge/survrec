#' Survival estimator under a gamma frailty model
#'
#' Maximum likelihood estimation of the survival function of correlated
#' inter-occurrence times of recurrent event data under a gamma frailty
#' model.
#'
#' The product-limit estimator of Peña, Strawderman and Hollander (2001)
#' is valid when the inter-occurrence times are an IID sample from some
#' underlying distribution F. This assumption is clearly restrictive in
#' biomedical applications, and one obvious generalization that allows
#' association between inter-occurrence times is a frailty model.
#'
#' A common and convenient choice of frailty distribution is a gamma with
#' shape and scale parameters set equal to an unknown parameter
#' \eqn{\alpha}. The common marginal survival function can be written as
#'
#' \deqn{\bar{F}(t) = \left[\frac{\alpha}{\alpha + \Lambda_0(t)}
#'   \right]^{\alpha}}{1 - F(t) = (alpha/(alpha + Lambda_0(t)))^alpha}
#'
#' The parameter \eqn{\alpha} controls the degree of association between
#' inter-occurrence times within a unit. Peña, Strawderman and Hollander
#' (2001) showed that \eqn{\alpha} and \eqn{\Lambda_0} can be estimated by
#' maximising the marginal likelihood function with an
#' expectation-maximisation (EM) algorithm.
#'
#' To achieve good convergence, a seed value for \eqn{\alpha} is estimated
#' first by maximising the profile likelihood for \eqn{\alpha} over the
#' interval (`alpha.min`, `alpha.max`) with Brent's method. If the EM
#' algorithm does not converge from that seed, it is restarted from up to
#' six additional seeds around it. If convergence still fails, the initial
#' value of alpha can be used as the `alpha.min` argument and the fit
#' recomputed.
#'
#' @param x a survival recurrent event object, see [Survr()].
#' @param tvals optional vector of times at which the survival function is
#'   also evaluated.
#' @param lambda optional vector of baseline hazard probabilities at the
#'   distinct event times. Defaults to the occurrence/exposure rates
#'   `n.event/colSums(AtRisk)`.
#' @param alpha optional shape and scale parameter of the frailty
#'   distribution. If unknown it is estimated via the EM algorithm,
#'   starting from a seed obtained by maximising the profile likelihood
#'   (see details).
#' @param alpha.min optional left bound of the interval used to find the
#'   seed of `alpha`. Default is 0.5.
#' @param alpha.max optional right bound of the interval used to find the
#'   seed of `alpha`. Default is the maximum distinct event time.
#' @param tol convergence tolerance of the EM algorithm. Default is 1e-7.
#' @param maxiter maximum number of EM iterations. Default is 500.
#' @param alpha.console if `TRUE`, prints the seed value and the final
#'   estimate of `alpha`.
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
#'   \item{status}{1 if the EM algorithm converged, -1 if it reached
#'     `maxiter` without converging (in that case no estimates are
#'     provided).}
#'   \item{alpha}{estimate of the gamma frailty parameter.}
#'   \item{lambda}{estimated hazard probabilities at the distinct event
#'     times.}
#'   \item{frailties}{posterior frailty estimate of each subject; values
#'     spread away from 1 indicate heterogeneity between subjects.}
#'   \item{survfunc}{survival estimate at each distinct time.}
#'   \item{tvals}{copy of the `tvals` argument.}
#'   \item{MLEAttvals}{survival estimate at the `tvals` times.}
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
#' fit <- mlefrailty_fit(Survr(MMC$id, MMC$time, MMC$event))
#' fit
#' plot(fit)
#'
#' # compare with Peña-Strawderman-Hollander
#' fit <- psh_fit(Survr(MMC$id, MMC$time, MMC$event))
#' lines(fit, lty = 2)
#'
#' # and with Wang-Chang
#' fit <- wc_fit(Survr(MMC$id, MMC$time, MMC$event))
#' lines(fit, lty = 3)
#' @keywords survival
#' @export
mlefrailty_fit <- function(x, tvals, lambda = NULL, alpha = NULL, alpha.min,
                           alpha.max, tol = 1e-07, maxiter = 500,
                           alpha.console = TRUE) {
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
  numdistinct <- length(summ$time)
  distinct <- summ$time
  numdeaths <- summ$n.event
  AtRisk <- summ$at.risk

  if (is.null(lambda[1])) {
    lambda <- numdeaths / colSums(AtRisk)
  }
  if (is.null(alpha)) {
    if (alpha.console) {
      cat("\nNeeds to Determine a Seed Value for Alpha")
    }
    if (missing(alpha.min)) {
      alpha.min <- 0.5
    }
    if (missing(alpha.max)) {
      alpha.max <- max(distinct)
    }
    if (alpha.min >= alpha.max) {
      stop("alpha.min must be smaller than alpha.max")
    }
    tol.max <- (alpha.max - alpha.min) / 50
    alpha <- .search_seed(
      as.integer(m), as.double(distinct),
      as.integer(numdeaths), AtRisk, as.double(lambda),
      as.double(alpha.min), as.double(alpha.max), as.double(tol.max)
    )
    if (alpha.console) {
      cat("\n Seed Alpha: ", alpha)
    }
  }

  # EM restarts: the seed first, then six values spread around it
  alphadel <- alpha / 4
  alphaseeds <- c(
    alpha, alpha - alphadel, alpha - 2 * alphadel,
    alpha - 3 * alphadel, alpha + alphadel, alpha + 2 * alphadel,
    alpha + 3 * alphadel
  )
  status <- 0
  ind <- 0
  while ((status != 1) && (ind < 7)) {
    ind <- ind + 1
    alpha <- alphaseeds[ind]
    Estimates <- .em_frailty(
      as.integer(m), as.double(distinct),
      as.integer(numdeaths), AtRisk, as.double(lambda),
      as.double(alpha), as.double(tol), as.integer(maxiter)
    )
    status <- Estimates$status
  }
  alpha <- Estimates$alpha
  if (alpha.console) {
    cat("\n ")
    cat("\n Alpha estimate=", alpha)
    cat("\n ")
  }
  lambda <- Estimates$lambda
  frailties <- Estimates$frailties

  if (!missing(tvals)) {
    tvalslen <- length(tvals)
  }
  if (!(status == 1)) {
    cat("\n\n WARNING: No estimates will be provided!")
    cat(
      "\n Value of (status,alpha) from iteration is ",
      c(status, alpha), "\n\n"
    )
    alpha <- NA
    survfuncMLE <- c(rep(NA, numdistinct))
    if (!missing(tvals)) {
      MLEAttvals <- c(rep(NA, tvalslen))
    } else {
      tvals <- NA
      MLEAttvals <- NA
    }
  } else {
    if (alpha >= 1e+05) {
      alpha <- 1e+05
    }
    survfuncMLE <- .mle_survival(as.double(lambda), as.double(alpha))
    if (!missing(tvals)) {
      MLEAttvals <- surv_search(sort(tvals), distinct, survfuncMLE)
    } else {
      tvals <- NA
      MLEAttvals <- NA
    }
  }
  ans <- list(
    n = n, m = m, failed = failed, censored = censored,
    time = distinct, n.event = numdeaths, AtRisk = AtRisk,
    status = status, alpha = alpha, lambda = lambda,
    frailties = frailties, survfunc = survfuncMLE, tvals = tvals,
    MLEAttvals = MLEAttvals
  )
  oldClass(ans) <- "survfitr"
  ans
}
