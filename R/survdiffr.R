#' Bootstrap survival quantiles for recurrent event data
#'
#' Obtains bootstrap replicates of a quantile (e.g. the median) of the
#' survival time for each group of subjects. Confidence intervals of the
#' replicates can be computed with [boot::boot.ci()].
#'
#' Three resampling schemes are available through `boot.F`: nonparametric
#' resampling of gap times from the Peña-Strawderman-Hollander (`"PSH"`)
#' or the Wang-Chang (`"WC"`) estimates of the inter-occurrence time
#' distribution, and semiparametric resampling under the fitted gamma
#' frailty model (`"semiparametric"`). With `boot.G = "empirical"` the
#' per-subject censoring times are additionally resampled from their
#' empirical distribution.
#'
#' All resampling uses R's random number generator, so results are
#' reproducible with [set.seed()] or the `seed` argument. The re-estimation
#' of the survival curve on each replicate can run in parallel, see
#' [survrecThreads()]. Some procedures (notably the semiparametric one,
#' which refits the frailty model on every replicate) can be slow.
#'
#' @param formula a formula object with a [Survr()] object as the response
#'   on the left of the `~` operator and a grouping term on the right.
#' @param data a data frame in which to interpret the variables named in
#'   the formula.
#' @param q quantile of the survival function whose time is bootstrapped.
#' @param B number of bootstrap samples.
#' @param boot.F a character string specifying the bootstrap procedure:
#'   `"PSH"` or `"WC"` (nonparametric) or `"semiparametric"`. The default
#'   is `"WC"`. Only the first letters are required, e.g. `"P"`, `"W"`,
#'   `"se"`.
#' @param boot.G a character string specifying whether the censoring times
#'   are also resampled from their empirical distribution: `"none"`
#'   (default) or `"empirical"`. Only the first letters are required.
#' @param seed optional random seed for the bootstrap resampling (passed
#'   to [set.seed()]). If `NULL` (default) the current state of the random
#'   number generator is used.
#' @param ... additional arguments passed to the estimator used for the
#'   observed (non-bootstrap) quantile.
#'
#' @return An object of class `boot` (or a named list with one per group)
#'   whose `t0` is the observed quantile and `t` the bootstrap replicates.
#'   A replicate is `-1` when the survival estimate does not fall below
#'   `q`, and `NA` when the resampled data had no events. For the
#'   semiparametric procedure the component `alpha` contains the frailty
#'   parameter estimated on each replicate.
#'
#' @references
#' González, J.R. and Peña, E.A. (2003). Bootstrapping median survival
#' with recurrent event data. *IX Conferencia Española de Biometría*.
#'
#' @seealso [survfitr()], [boot::boot.ci()], [survrecThreads()]
#'
#' @examples
#' data(colon)
#'
#' # compare the median survival time between the three Dukes stages
#' fit <- survdiffr(Survr(hc, time, event) ~ as.factor(dukes),
#'   data = colon, q = 0.5, seed = 1
#' )
#' boot::boot.ci(fit$"1", type = c("norm", "basic", "perc"))
#' boot::boot.ci(fit$"2", type = c("norm", "basic", "perc"))
#' boot::boot.ci(fit$"3", type = c("norm", "basic", "perc"))
#'
#' # 75th quantile with percentile confidence intervals
#' fit <- survdiffr(Survr(hc, time, event) ~ as.factor(dukes),
#'   data = colon, q = 0.75, seed = 1
#' )
#' quantile(fit$"1"$t, c(0.025, 0.975))
#'
#' # resampling from the PSH estimate instead
#' fit <- survdiffr(Survr(hc, time, event) ~ as.factor(dukes),
#'   data = colon, q = 0.5, boot.F = "PSH", seed = 1
#' )
#' @keywords survival
#' @export
survdiffr <- function(formula, data, q, B = 500, boot.F = "WC",
                      boot.G = "none", seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  method.F <- charmatch(boot.F, c("PSH", "WC", "semiparametric"),
    nomatch = 0
  )
  if (method.F == 0) {
    stop("bootstrap from F must be PSH, WC or semiparametric")
  }
  method.G <- charmatch(boot.G, c("none", "empirical"), nomatch = 0)
  if (method.G == 0) {
    stop("bootstrap from G must be none or empirical")
  }

  # plans 2-7 of the resampling engine; `type` selects the estimator
  # used for the observed quantile t0
  if (method.F == 1) {
    type.boot <- 2 + method.G - 1
    type <- "p"
  }
  if (method.F == 2) {
    type.boot <- 4 + method.G - 1
    type <- "w"
  }
  if (method.F == 3) {
    type.boot <- 6 + method.G - 1
    type <- "M"
  }

  call <- match.call()
  if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("Survr")) ||
    inherits(formula, "Survr")) {
    stop("formula.default(object): invalid formula")
  }

  m <- match.call(expand.dots = FALSE)
  m$q <- m$B <- m$boot.F <- m$boot.G <- m$seed <- m$... <- NULL
  Terms <- terms(formula, "strata")
  ord <- attr(Terms, "order")
  if (length(ord) & any(ord != 1)) {
    stop("Interaction terms are not valid for this function")
  }
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  Y <- model.extract(m, "response")
  if (!is.Survr(Y)) {
    stop("Response must be a survival recurrent object")
  }
  ll <- attr(Terms, "term.labels")
  group <- m[ll][, 1]

  # builds one boot-class object from the replicates of one group
  boot_one <- function(Sr) {
    d <- survr_pieces(Sr)
    summ <- .boot_quantile(
      as.integer(d$m), as.double(d$failed),
      as.double(d$censored), as.double(d$tau), as.integer(B),
      as.integer(type.boot), as.double(q), resolve_threads()
    )

    ans <- list()
    ans$t0 <- q_search(survfitr(Sr ~ 1, type = type), q = q)
    ans$t <- cbind(summ$t)
    if (!is.null(summ$alpha)) ans$alpha <- summ$alpha
    ans$R <- B
    ans$data <- unclass(Sr)
    ans$seed <- get(".Random.seed", envir = globalenv())
    ans$statistic <- NULL
    ans$sim <- c("ordinary")
    ans$call <- call
    ans$stype <- c("i")
    ans$strata <- rep(1, nrow(Sr))
    ans$weights <- rep(1 / nrow(Sr), nrow(Sr))
    oldClass(ans) <- "boot"
    ans
  }

  if (!is.null(group)) {
    k <- levels(group)
    ans <- vector("list", length(k))
    for (i in 1:length(k)) {
      temp <- Y[group == k[i], ]
      ans[[i]] <- boot_one(Survr(temp[, 1], temp[, 2], temp[, 3]))
    }
    names(ans) <- k
    # the per-group elements keep class "boot"; the container gets its
    # own class so that summary()/autoplot() can compare the groups
    oldClass(ans) <- "survdiffr"
  } else {
    ans <- boot_one(Survr(Y[, 1], Y[, 2], Y[, 3]))
  }
  ans
}
