# survdiffr uses R's RNG in 2.0.0 (the Fortran version used the GNU rand()
# extension, which set.seed could not control), so exact reproduction of the
# old bootstrap replicates is impossible by design. We test reproducibility,
# the deterministic pieces, and distributional sanity against the reference.

test_that("survdiffr is reproducible with a seed", {
  b1 <- survdiffr(Survr(id, time, event) ~ group, data = MMC,
                  q = 0.5, B = 25, boot.F = "WC", seed = 42)
  b2 <- survdiffr(Survr(id, time, event) ~ group, data = MMC,
                  q = 0.5, B = 25, boot.F = "WC", seed = 42)
  expect_identical(b1$Males$t, b2$Males$t)
  expect_identical(b1$Females$t, b2$Females$t)

  set.seed(7)
  b3 <- survdiffr(Survr(id, time, event) ~ group, data = MMC,
                  q = 0.5, B = 25, boot.F = "WC")
  set.seed(7)
  b4 <- survdiffr(Survr(id, time, event) ~ group, data = MMC,
                  q = 0.5, B = 25, boot.F = "WC")
  expect_identical(b3$Males$t, b4$Males$t)
})

test_that("results are independent of the number of threads", {
  # resampling is serial (R RNG); only the deterministic re-estimation is
  # parallel, so the thread count must not change the results
  old <- survrecThreads(1)
  on.exit(survrecThreads(old))
  b1 <- survdiffr(Survr(id, time, event) ~ group, data = MMC,
                  q = 0.5, B = 25, boot.F = "WC", seed = 5)
  survrecThreads(2)
  b2 <- survdiffr(Survr(id, time, event) ~ group, data = MMC,
                  q = 0.5, B = 25, boot.F = "WC", seed = 5)
  expect_identical(b1$Males$t, b2$Males$t)
  expect_identical(b1$Females$t, b2$Females$t)
})

test_that("survdiffr point estimates match the Fortran reference", {
  b <- survdiffr(Survr(id, time, event) ~ group, data = MMC,
                 q = 0.5, B = 10, boot.F = "WC", seed = 1)
  expect_equal(b$Males$t0, reference_fits$survdiffr_mmc_wc$Males$t0)
  expect_equal(b$Females$t0, reference_fits$survdiffr_mmc_wc$Females$t0)
  expect_s3_class(b$Males, "boot")
  expect_equal(b$Males$R, 10)

  b <- survdiffr(Survr(id, time, event) ~ group, data = MMC,
                 q = 0.5, B = 10, boot.F = "PSH", seed = 1)
  expect_equal(b$Males$t0, reference_fits$survdiffr_mmc_psh$Males$t0)
  expect_equal(b$Females$t0, reference_fits$survdiffr_mmc_psh$Females$t0)
})

test_that("bootstrap replicates are statistically compatible with Fortran", {
  b <- survdiffr(Survr(id, time, event) ~ group, data = MMC,
                 q = 0.5, B = 200, boot.F = "WC", seed = 123)
  for (g in c("Males", "Females")) {
    t_new <- b[[g]]$t[b[[g]]$t > 0]
    ref_q <- reference_fits$survdiffr_mmc_wc[[g]]$t_quantiles
    # medians of the bootstrap distributions should be in the same range
    expect_gt(median(t_new), ref_q[["10%"]])
    expect_lt(median(t_new), ref_q[["90%"]])
    # replicates are valid: -1 (quantile not reached), NA, or observed times
    expect_true(all(is.na(b[[g]]$t) | b[[g]]$t == -1 | b[[g]]$t > 0))
  }
})

test_that("semiparametric bootstrap works in 2.0.0", {
  # this plan existed in the Fortran code but was disabled in the R layer;
  # colon/chemoter is used because the frailty model converges in both groups
  b <- suppressWarnings(
    survdiffr(Survr(hc, time, event) ~ as.factor(chemoter), data = colon,
              q = 0.5, B = 5, boot.F = "semiparametric", seed = 99)
  )
  for (g in names(b)) {
    expect_s3_class(b[[g]], "boot")
    expect_true(is.finite(b[[g]]$t0))
    expect_length(b[[g]]$alpha, 5)
    expect_true(all(b[[g]]$alpha[!is.na(b[[g]]$alpha)] > 0))
  }
})
