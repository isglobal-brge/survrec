#' Deprecated functions in survrec
#'
#' The dot-separated names of survrec 1.x were renamed to snake_case in
#' version 2.0.0. The old names still work but emit a deprecation warning
#' and will be removed in a future release.
#'
#' \describe{
#'   \item{`psh.fit()`}{use [psh_fit()] instead.}
#'   \item{`wc.fit()`}{use [wc_fit()] instead.}
#'   \item{`mlefrailty.fit()`}{use [mlefrailty_fit()] instead.}
#'   \item{`surv.search()`}{use [surv_search()] instead.}
#'   \item{`q.search()`}{use [q_search()] instead.}
#' }
#'
#' @param x,tvals,lambda,alpha,alpha.min,alpha.max,tol,maxiter,alpha.console
#'   see the replacement functions.
#' @param time,surv,f,q see the replacement functions.
#'
#' @return Each deprecated function emits a deprecation warning and then
#'   returns exactly what its replacement returns: `psh.fit()`, `wc.fit()`
#'   and `mlefrailty.fit()` return a `survfitr` object, `surv.search()`
#'   returns the interpolated survival probabilities as a numeric vector,
#'   and `q.search()` returns the estimated quantile as a single number.
#'   See the replacement functions for the full description of the value.
#' @name survrec-deprecated
#' @keywords internal
NULL

#' @rdname survrec-deprecated
#' @export
psh.fit <- function(x, tvals) {
  .Deprecated("psh_fit")
  if (missing(tvals)) psh_fit(x) else psh_fit(x, tvals)
}

#' @rdname survrec-deprecated
#' @export
wc.fit <- function(x, tvals) {
  .Deprecated("wc_fit")
  if (missing(tvals)) wc_fit(x) else wc_fit(x, tvals)
}

#' @rdname survrec-deprecated
#' @export
mlefrailty.fit <- function(x, tvals, lambda = NULL, alpha = NULL, alpha.min,
                           alpha.max, tol = 1e-07, maxiter = 500,
                           alpha.console = TRUE) {
  .Deprecated("mlefrailty_fit")
  cl <- match.call()
  cl[[1]] <- quote(survrec::mlefrailty_fit)
  eval.parent(cl)
}

#' @rdname survrec-deprecated
#' @export
surv.search <- function(tvals, time, surv) {
  .Deprecated("surv_search")
  surv_search(tvals, time, surv)
}

#' @rdname survrec-deprecated
#' @export
q.search <- function(f, q = 0.5) {
  .Deprecated("q_search")
  q_search(f, q)
}
