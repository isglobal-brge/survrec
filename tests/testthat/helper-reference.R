# Shared helpers for the regression tests.
#
# The reference ("gold standard") was generated with the original Fortran
# implementation (v1.2-5) via tools/build-reference-fits.R. The C++ port
# must reproduce these values, except where a documented Fortran bug was
# fixed (see tools/build-reference-fits.R for the list):
#
#  * mlefrailty_fit: the Fortran EM silently relaxed the convergence
#    tolerance to 1e-4 after the first iteration (estalpha overwrote the
#    aliased `tol` argument) and computed part of the convergence distance
#    in single precision. The port converges to the user tolerance in
#    double precision, so frailty results are compared with TOL_MLE.
#
#  * wc_fit: the Fortran read an undefined loop index in the event count of
#    subjects whose only gap is the censored one. The fixed code gives
#    identical results whenever no censored gap of such subjects ties a
#    distinct failure time (true for the reference datasets).

reference_fits <- readRDS(test_path("fixtures", "fortran-reference.rds"))

TOL_EXACT <- 1e-8  # deterministic estimators, identical algorithm
TOL_MLE <- 1e-3    # frailty model: convergence-tolerance bug fixed

# Components whose values must match between a fresh survfitr-type fit and
# the stored reference. `m` is compared without names/class (the reference
# stores `table` objects).
expect_matches_reference <- function(fit, ref, tolerance = TOL_EXACT,
                                     skip = character()) {
  numeric_parts <- c("n", "failed", "censored", "time", "n.event",
                     "survfunc", "std.error", "tvals")
  for (part in setdiff(intersect(numeric_parts, names(ref)), skip)) {
    expect_equal(fit[[part]], ref[[part]], tolerance = tolerance,
                 info = part, ignore_attr = TRUE)
  }
  if (!("m" %in% skip)) {
    expect_equal(as.integer(fit$m), as.integer(ref$m), info = "m")
  }
  if (!("AtRisk" %in% skip) && !is.null(ref$AtRisk)) {
    expect_equal(unname(as.matrix(fit$AtRisk)), unname(as.matrix(ref$AtRisk)),
                 tolerance = tolerance, info = "AtRisk", ignore_attr = TRUE)
  }
  invisible(fit)
}
