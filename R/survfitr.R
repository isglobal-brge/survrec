#' Survival curves for recurrent event data
#'
#' Computes an estimate of the survival curve of inter-occurrence times
#' for recurrent event data, optionally by group, using the
#' Peña-Strawderman-Hollander, the Wang-Chang or the gamma frailty MLE
#' estimators, together with their standard errors.
#'
#' See [psh_fit()], [wc_fit()] and [mlefrailty_fit()] for details on each
#' estimator; additional arguments in `...` are passed to the chosen
#' fitting function.
#'
#' @param formula a formula object with a [Survr()] object as the response
#'   on the left of the `~` operator and a grouping term on the right.
#'   For a single survival curve use `~ 1`.
#' @param data a data frame in which to interpret the variables named in
#'   the formula.
#' @param type a character string specifying the estimator:
#'   `"pena-strawderman-hollander"`, `"wang-chang"` or `"MLEfrailty"`
#'   (default). Only the first letters are required, e.g. `"pe"`, `"wa"`,
#'   `"ML"`.
#' @param ... additional arguments passed to the fitting function of the
#'   chosen estimator.
#'
#' @return A `survfitr` object: the fit of the chosen estimator, or a
#'   named list with one fit per group when the formula defines groups.
#'   Methods are provided for `print`, `plot`, `lines` and `summary`.
#'
#' @references
#' Peña, E.A., Strawderman, R. and Hollander, M. (2001). Nonparametric
#' Estimation with Recurrent Event Data. *Journal of the American
#' Statistical Association* **96**, 1299--1315.
#'
#' Wang, M.-C. and Chang, S.-H. (1999). Nonparametric Estimation of a
#' Recurrent Survival Function. *Journal of the American Statistical
#' Association* **94**, 146--153.
#'
#' @seealso [print.survfitr()], [plot.survfitr()], [summary.survfitr()],
#'   [Survr()], [psh_fit()], [wc_fit()], [mlefrailty_fit()]
#'
#' @examples
#' data(colon)
#' # fit a Peña-Strawderman-Hollander estimator by Dukes stage
#' fit <- survfitr(Survr(hc, time, event) ~ as.factor(dukes),
#'   data = colon, type = "pena"
#' )
#' plot(fit, ylim = c(0, 1), xlim = c(0, 2000))
#' fit
#' summary(fit)
#' @keywords survival
#' @export
survfitr <- function(formula, data, type = "MLEfrailty", ...) {
  method <- charmatch(type, c(
    "pena-strawderman-hollander",
    "wang-chang", "MLEfrailty"
  ), nomatch = 0)
  if (method == 0) {
    stop("estimator must be pena-strawderman-hollander wang-chang or MLEfrailty")
  }

  call <- match.call()
  if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("Survr")) ||
    inherits(formula, "Survr")) {
    stop("formula.default(object): invalid formula")
  }

  m <- match.call(expand.dots = FALSE)
  m$type <- m$... <- NULL
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

  if (method == 1) FUN <- psh_fit
  if (method == 2) FUN <- wc_fit
  if (method == 3) FUN <- mlefrailty_fit

  if (ncol(m) > 1) {
    group <- m[ll][, 1]
    k <- levels(group)
    ans <- NULL
    for (i in 1:length(k)) {
      temp <- Y[group == k[i], ]
      temp1 <- Survr(temp[, 1], temp[, 2], temp[, 3])
      ans[[i]] <- FUN(temp1, ...)
    }
    names(ans) <- k
    oldClass(ans) <- "survfitr"
    attr(ans, "strata") <- length(k)
    attr(ans, "group") <- ll
  } else {
    temp <- Survr(Y[, 1], Y[, 2], Y[, 3])
    ans <- FUN(temp, ...)
  }
  ans
}
