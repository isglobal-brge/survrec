#' Create a survival recurrent object
#'
#' Creates a survival recurrent object, usually used as the response
#' variable in a model formula.
#'
#' Each subject must contribute one censored time (its last, incomplete
#' gap time), so `id` must have exactly as many distinct values as there
#' are zeros in `event`. Rows are assumed to be grouped by subject and in
#' chronological order within subject.
#'
#' @param id identifier of each subject; the same value for all the
#'   recurrent times of a subject.
#' @param time gap times of recurrence. The last time of each subject is
#'   censored.
#' @param event status indicator: 1 = recurrence, 0 = censored. Only these
#'   values are accepted.
#' @param x any R object.
#'
#' @return An object of class `Survr`, implemented as a matrix with
#'   columns `id`, `time` and `event`.
#'
#'   `is.Survr()` returns `TRUE` if `x` inherits from class `Survr` and
#'   `FALSE` otherwise.
#'
#' @seealso [survfitr()], [psh_fit()], [wc_fit()], [mlefrailty_fit()]
#'
#' @examples
#' data(MMC)
#' x <- Survr(MMC$id, MMC$time, MMC$event)
#' is.Survr(x)
#' @keywords survival
#' @export
Survr <- function(id, time, event) {
  if (length(unique(event)) > 2 || max(event) != 1 || min(event) != 0) {
    stop("event must be 0-1")
  }
  if (any(tapply(event, id, function(e) sum(e == 0)) != 1)) {
    stop("Data doesn't match. Every subject must have a censored time")
  }

  ans <- cbind(id, time, event)
  oldClass(ans) <- "Survr"
  invisible(ans)
}

#' @rdname Survr
#' @export
is.Survr <- function(x) inherits(x, "Survr")
