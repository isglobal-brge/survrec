# Regression tests of the C++ port against the original Fortran results
# (fixtures generated with tools/build-reference-fits.R, survrec 1.2-5).

surv_mmc <- Survr(MMC$id, MMC$time, MMC$event)
surv_colon <- Survr(colon$hc, colon$time, colon$event)

test_that("psh_fit reproduces the Fortran results", {
  expect_matches_reference(psh_fit(surv_mmc), reference_fits$psh_mmc)
  expect_matches_reference(psh_fit(surv_colon), reference_fits$psh_colon)
})

test_that("psh_fit with tvals reproduces the Fortran results", {
  tvals <- attr(reference_fits, "meta")$tvals_mmc
  fit <- expect_matches_reference(psh_fit(surv_mmc, tvals = tvals),
                                  reference_fits$psh_mmc_tvals)
  expect_equal(fit$PSHpleAttvals, reference_fits$psh_mmc_tvals$PSHpleAttvals)
})

test_that("wc_fit reproduces the Fortran results", {
  # the reference `m` counted the censored record (see helper); everything
  # numeric must match exactly
  fit <- wc_fit(surv_mmc)
  expect_matches_reference(fit, reference_fits$wc_mmc, skip = c("m", "AtRisk"))
  expect_equal(as.integer(fit$m), as.integer(reference_fits$wc_mmc$m) - 1L)
  # reference stored the weighted at-risk vector under AtRisk
  expect_equal(fit$AtRisk, reference_fits$wc_mmc$AtRisk,
               tolerance = TOL_EXACT, ignore_attr = TRUE)

  fit <- wc_fit(surv_colon)
  expect_matches_reference(fit, reference_fits$wc_colon,
                           skip = c("m", "AtRisk"))
  expect_equal(fit$AtRisk, reference_fits$wc_colon$AtRisk,
               tolerance = TOL_EXACT, ignore_attr = TRUE)
})

test_that("wc_fit with tvals reproduces the Fortran results", {
  tvals <- attr(reference_fits, "meta")$tvals_mmc
  fit <- wc_fit(surv_mmc, tvals = tvals)
  expect_equal(fit$WCpleAttvals, reference_fits$wc_mmc_tvals$WCpleAttvals,
               tolerance = TOL_EXACT)
})

test_that("mlefrailty_fit reproduces the Fortran results", {
  fit <- mlefrailty_fit(surv_mmc, alpha.console = FALSE)
  ref <- reference_fits$mle_mmc
  expect_equal(fit$status, ref$status)
  expect_equal(fit$alpha, ref$alpha, tolerance = TOL_MLE)
  expect_equal(fit$lambda, ref$lambda, tolerance = TOL_MLE)
  expect_equal(fit$survfunc, ref$survfunc, tolerance = TOL_MLE)
  expect_matches_reference(fit, ref, tolerance = TOL_MLE,
                           skip = c("survfunc"))
  # new in 2.0.0: posterior frailties are returned
  expect_length(fit$frailties, fit$n)
  expect_true(all(fit$frailties > 0))

  fit <- mlefrailty_fit(surv_colon, alpha.console = FALSE)
  ref <- reference_fits$mle_colon
  expect_equal(fit$status, ref$status)
  expect_equal(fit$alpha, ref$alpha, tolerance = TOL_MLE)
  expect_equal(fit$survfunc, ref$survfunc, tolerance = TOL_MLE)
})

test_that("mlefrailty_fit honours alpha, alpha.min/alpha.max and tvals", {
  fit <- mlefrailty_fit(surv_mmc, alpha = 5, alpha.console = FALSE)
  ref <- reference_fits$mle_mmc_alpha_fixed
  expect_equal(fit$alpha, ref$alpha, tolerance = TOL_MLE)
  expect_equal(fit$survfunc, ref$survfunc, tolerance = TOL_MLE)

  fit <- mlefrailty_fit(surv_mmc, alpha.min = 1, alpha.max = 100,
                        alpha.console = FALSE)
  ref <- reference_fits$mle_mmc_alpha_range
  expect_equal(fit$alpha, ref$alpha, tolerance = TOL_MLE)

  tvals <- attr(reference_fits, "meta")$tvals_mmc
  fit <- mlefrailty_fit(surv_mmc, tvals = tvals, alpha.console = FALSE)
  expect_equal(fit$MLEAttvals, reference_fits$mle_mmc_tvals$MLEAttvals,
               tolerance = TOL_MLE)
})

test_that("survfitr formula interface with strata matches the reference", {
  fit <- survfitr(Survr(id, time, event) ~ group, data = MMC,
                  type = "pena-strawderman-hollander")
  ref <- reference_fits$survfitr_mmc_group_psh
  expect_equal(attr(fit, "strata"), attr(ref, "strata"))
  expect_equal(names(fit), names(ref))
  for (g in names(ref)) {
    expect_matches_reference(fit[[g]], ref[[g]])
  }

  fit <- survfitr(Survr(id, time, event) ~ group, data = MMC,
                  type = "wang-chang")
  ref <- reference_fits$survfitr_mmc_group_wc
  for (g in names(ref)) {
    expect_matches_reference(fit[[g]], ref[[g]], skip = c("m", "AtRisk"))
    expect_equal(fit[[g]]$AtRisk, ref[[g]]$AtRisk,
                 tolerance = TOL_EXACT, ignore_attr = TRUE)
  }

  suppressWarnings({
    fit <- survfitr(Survr(hc, time, event) ~ as.factor(dukes),
                    data = colon, type = "wang-chang")
  })
  ref <- reference_fits$survfitr_colon_dukes_wc
  expect_equal(attr(fit, "strata"), attr(ref, "strata"))
  for (g in names(ref)) {
    expect_matches_reference(fit[[g]], ref[[g]], skip = c("m", "AtRisk"))
  }
})

test_that("q_search and surv_search match the reference", {
  expect_equal(q_search(psh_fit(surv_mmc), q = 0.5),
               reference_fits$qsearch_psh_mmc_q50)
  expect_equal(q_search(wc_fit(surv_mmc), q = 0.5),
               reference_fits$qsearch_wc_mmc_q50)
  expect_equal(q_search(mlefrailty_fit(surv_mmc, alpha.console = FALSE),
                        q = 0.5),
               reference_fits$qsearch_mle_mmc_q50)

  tvals <- attr(reference_fits, "meta")$tvals_mmc
  fit <- psh_fit(surv_mmc)
  expect_equal(surv_search(tvals, fit$time, fit$survfunc),
               reference_fits$survsearch_psh_mmc)
})

# The stored summary_* fixtures are not usable as reference: in 1.2-5
# summary() on a single (non-strata) fit stopped with "$ operator is invalid
# for atomic vectors" (it read x[[1]]$std.error, i.e. n), a bug fixed in
# 2.0.0. The table contents are checked directly against the fit instead.
test_that("summary.survfitr returns the estimate table", {
  fit <- psh_fit(surv_mmc)
  s <- summary(fit)
  expect_equal(colnames(s), c("time", "n.event", "n.risk", "surv",
                              "std.error"))
  expect_equal(unname(s[, "time"]), fit$time)
  expect_equal(unname(s[, "n.event"]), as.numeric(fit$n.event))
  expect_equal(unname(s[, "n.risk"]), unname(apply(fit$AtRisk, 2, sum)))
  expect_equal(unname(s[, "surv"]), fit$survfunc)
  expect_equal(unname(s[, "std.error"]), fit$std.error)

  fit <- survfitr(Survr(id, time, event) ~ group, data = MMC,
                  type = "wang-chang")
  s <- summary(fit)
  expect_equal(names(s), names(fit))
  expect_equal(unname(s$Males[, "surv"]), fit$Males$survfunc)
})
