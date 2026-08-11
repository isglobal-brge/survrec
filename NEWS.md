# survrec 2.0.0

## Major changes

* The dot-separated function names of survrec 1.x were renamed to
  snake_case: `psh_fit()`, `wc_fit()`, `mlefrailty_fit()`,
  `surv_search()` and `q_search()`. The old names (`psh.fit()`,
  `wc.fit()`, `mlefrailty.fit()`, `surv.search()`, `q.search()`) still
  work but emit a deprecation warning and will be removed in a future
  release. S3 methods (`plot.survfitr()`, ...) and the argument names
  shared with the survival ecosystem (`conf.int`, `boot.F`, ...) are
  unchanged.

* The computational core has been rewritten from Fortran 77 to modern C++
  (via Rcpp). The public R API (`Survr()`, `survfitr()`, `psh_fit()`,
  `wc_fit()`, `mlefrailty_fit()`, `survdiffr()`, `q_search()`,
  `surv_search()` and the S3 methods) is unchanged.
* New maintainer: Dolors Pelegri-Siso.
* The bootstrap (`survdiffr()`) now uses R's random number generator, so
  results are reproducible with `set.seed()`. The previous implementation
  used the non-standard GNU `rand()` generator, which R could not seed. The
  `seed` argument now defaults to `NULL` (use the current RNG state) and is
  honoured when supplied; in 1.2-5 it was silently ignored.
* The semiparametric (gamma frailty) bootstrap of `survdiffr()`
  (`boot.F = "semiparametric"`) is now available; it existed in the Fortran
  code but was disabled in the R layer.
* `mlefrailty_fit()` now returns the posterior frailty estimates
  (`$frailties`), previously computed internally and discarded.
* The proprietary IMSL routine `ZXGSP` (golden-section search) and the
  copied Netlib Brent/`D1MACH` code have been replaced by a clean
  implementation of Brent's minimization algorithm.
* All fixed-size internal buffers (800 elements) have been removed; data
  size is now unlimited.

## Bug fixes

* `summary()` on a fit without strata no longer fails with "$ operator is
  invalid for atomic vectors".
* `wc_fit()` no longer reads an undefined loop variable when counting
  events of subjects whose only gap time is censored (undefined behaviour
  in the Fortran `wc2` routine).
* The EM algorithm of `mlefrailty_fit()` now honours the requested
  convergence tolerance `tol`. In 1.2-5 the inner alpha search overwrote it
  with 1e-4 after the first iteration and computed part of the convergence
  criterion in single precision, so estimates may differ slightly (they are
  now more accurate).
* The 7-seed retry loop of `mlefrailty_fit()` is now really executed when
  the EM does not converge (a status-flag bug made it dead code).
* `survdiffr()` no longer fails when `boot.G` or `seed` are supplied
  together with a grouping formula (they leaked into the model frame).
* `q_search()` returns `NA` with a warning (instead of failing) when the
  survival estimates are unavailable.
* `wc_fit()` now reports `m` as the number of recurrences per subject,
  consistently with `psh_fit()`; it previously counted the censored record
  as well.
* Several `plot`, `summary` and `q_search` internals accessed `$surv`
  relying on partial matching; they now use the actual component name.

## New features

* New `ggplot2` graphics sharing the visual style of the companion
  package gcmrec: `autoplot()` methods for survival fits (with
  `fun = "surv"/"event"/"cumhaz"` and pointwise log-log confidence
  bands), `plotEstimators()` to overlay the three estimators on the same
  data, `autoplot()` for `mcf` and grouped `survdiffr` results, and the
  exported `theme_survrec()`. The base-graphics `plot()`/`lines()`
  methods are kept for backwards compatibility.
* New `mcf()`: Nelson-Aalen type estimate of the mean cumulative function
  (expected number of events by calendar time), optionally by group.
* New `quantile()` and `as.data.frame()` methods for `survfitr` objects.
* Grouped `survdiffr()` results now have class `survdiffr` with a
  `summary()` method that reports percentile bootstrap intervals per
  group and pairwise differences of quantiles with bootstrap p-values
  (the per-group elements keep class `boot`, so `boot::boot.ci()` keeps
  working).
* Fits no longer depend on the ordering of the input rows: subjects are
  regrouped by order of first appearance. In 1.2-5 the estimators
  silently misaligned subjects and gap times when the ids did not appear
  in sorted, contiguous blocks.
* `Survr()` now checks that every subject has exactly one censored time
  (the old check could be fooled by compensating subjects).

* The bootstrap re-estimation of `survdiffr()` runs in parallel with
  OpenMP when available. The number of threads is controlled with the new
  `survrecThreads()` function (default: let OpenMP decide, honouring
  `OMP_NUM_THREADS`). Resampling stays serial on R's RNG, so results are
  identical for any number of threads.
* `plot.survfitr()` no longer draws meaningless confidence bands when the
  fit provides no standard errors or `conf.int = FALSE` is given.

## Infrastructure

* Registered native routines via Rcpp; `DESCRIPTION` uses `Authors@R`;
  documentation migrated from hand-written latin1 `.Rd` files to roxygen2
  (UTF-8), with the `NAMESPACE` generated by roxygen2.
* R sources reorganized into one file per exported function; interpreted
  code vectorized (`colSums()`, `rowsum()`, `findInterval()` replace
  explicit loops and `outer()`/`apply()` constructions).
* Regression test suite (testthat) validating the C++ port against
  reference values generated with the original Fortran implementation
  (see `tools/build-reference-fits.R`).
