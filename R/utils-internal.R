# Internal helpers shared by the fitting and plotting functions.

# Splits a Survr object into the subject-grouped pieces expected by the
# C++ core. Rows are regrouped by subject in order of first appearance
# (with a stable sort, so the within-subject chronological order is kept),
# which makes the result well defined even when the input rows are not
# contiguous by subject. This replaces the table()-based extraction of
# survrec <= 1.2-5, which silently misaligned subjects and gap times when
# the ids did not appear in sorted order.
survr_pieces <- function(x) {
  fid <- factor(x[, 1], levels = unique(x[, 1]))
  ord <- order(fid)
  time <- x[ord, 2]
  event <- x[ord, 3]
  fid <- fid[ord]

  m <- tabulate(fid) - 1L
  names(m) <- levels(fid)
  list(
    n = nlevels(fid),
    m = m,
    failed = time[event == 1],
    censored = time[event == 0],
    # total observation time per subject
    tau = as.vector(rowsum(time, fid))
  )
}

# Duplicates the interior points of (time, ...) so that geom_ribbon and
# geom_line draw proper step functions: (t1,y1),(t2,y1),(t2,y2),...
step_frame <- function(d) {
  n <- nrow(d)
  if (n < 2) {
    return(d)
  }
  out <- d[rep(seq_len(n), each = 2)[-(2 * n)], , drop = FALSE]
  out$time <- rep(d$time, each = 2)[-1]
  rownames(out) <- NULL
  out
}

# Pointwise confidence limits of a survival estimate. The log-log
# transformation keeps the limits inside [0, 1] (the default of
# survival::survfit); "plain" reproduces the symmetric bands of the base
# plot method.
surv_ci <- function(surv, se, level = 0.95,
                    conf.type = c("log-log", "plain")) {
  conf.type <- match.arg(conf.type)
  z <- qnorm(1 - (1 - level) / 2)
  if (conf.type == "plain") {
    return(list(
      lower = pmax(surv - z * se, 0),
      upper = pmin(surv + z * se, 1)
    ))
  }
  # se of log(-log S) by the delta method; degenerate at S = 0 or 1
  ok <- surv > 0 & surv < 1 & se > 0
  se.ll <- ifelse(ok, se / (surv * abs(log(surv))), NA_real_)
  lower <- ifelse(ok, surv^exp(z * se.ll), surv)
  upper <- ifelse(ok, surv^exp(-z * se.ll), surv)
  list(lower = lower, upper = upper)
}

# survfitr objects are either a single fit or a list of per-group fits;
# this returns a uniform named list of fits.
strata_fits <- function(x) {
  if (is.null(attr(x, "strata"))) {
    list(fits = list(x), groups = NULL)
  } else {
    list(
      fits = unclass(x)[seq_len(attr(x, "strata"))],
      groups = names(x)
    )
  }
}
