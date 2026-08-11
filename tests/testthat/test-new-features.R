# Features introduced in 2.0.0: robust subject marshalling, quantile and
# as.data.frame methods, MCF, and the ggplot2 graphics.

surv_mmc <- Survr(MMC$id, MMC$time, MMC$event)

test_that("fits do not depend on the ordering of the input rows", {
  ref <- psh_fit(surv_mmc)

  # subjects in reversed block order
  ids <- unique(MMC$id)
  rev_idx <- unlist(lapply(rev(ids), function(i) which(MMC$id == i)))
  rev_fit <- psh_fit(Survr(
    MMC$id[rev_idx], MMC$time[rev_idx],
    MMC$event[rev_idx]
  ))
  expect_equal(rev_fit$survfunc, ref$survfunc)
  expect_equal(rev_fit$std.error, ref$std.error)
  expect_equal(sort(as.integer(rev_fit$m)), sort(as.integer(ref$m)))

  # subjects interleaved row by row (ids no longer contiguous)
  int_idx <- order(ave(seq_along(MMC$id), MMC$id, FUN = seq_along))
  int_fit <- wc_fit(Survr(
    MMC$id[int_idx], MMC$time[int_idx],
    MMC$event[int_idx]
  ))
  expect_equal(int_fit$survfunc, wc_fit(surv_mmc)$survfunc)

  mle_ref <- mlefrailty_fit(surv_mmc, alpha.console = FALSE)
  mle_rev <- mlefrailty_fit(
    Survr(MMC$id[rev_idx], MMC$time[rev_idx], MMC$event[rev_idx]),
    alpha.console = FALSE
  )
  # summing subjects in a different order changes floating-point rounding
  expect_equal(mle_rev$alpha, mle_ref$alpha, tolerance = 1e-6)
})

test_that("quantile method matches q_search", {
  fit <- wc_fit(surv_mmc)
  expect_equal(
    unname(quantile(fit, probs = 0.5)),
    q_search(fit, q = 0.5)
  )
  q <- quantile(fit, probs = c(0.25, 0.5, 0.75))
  expect_named(q, c("25%", "50%", "75%"))

  fit <- survfitr(Survr(id, time, event) ~ group,
    data = MMC,
    type = "wang-chang"
  )
  q <- quantile(fit, probs = c(0.25, 0.5))
  expect_equal(rownames(q), c("Males", "Females"))
  expect_equal(unname(q["Males", "50%"]), q_search(fit$Males, 0.5))
})

test_that("as.data.frame returns tidy curves", {
  fit <- wc_fit(surv_mmc)
  df <- as.data.frame(fit)
  expect_named(df, c("time", "n.event", "n.risk", "surv", "std.error"))
  expect_equal(df$surv, fit$survfunc)

  fit <- survfitr(Survr(id, time, event) ~ group,
    data = MMC,
    type = "pena"
  )
  df <- as.data.frame(fit)
  expect_named(
    df,
    c("group", "time", "n.event", "n.risk", "surv", "std.error")
  )
  expect_equal(unique(df$group), c("Males", "Females"))
  expect_equal(df$surv[df$group == "Males"], fit$Males$survfunc)
})

test_that("mcf estimates the mean cumulative function", {
  m <- mcf(Survr(id, time, event) ~ 1, data = MMC)
  expect_s3_class(m, "mcf")
  expect_true(all(diff(m$mcf) > 0))
  expect_equal(m$n.risk[1], length(unique(MMC$id)))
  # total expected events per subject <= mean events per subject observed
  expect_gt(max(m$mcf), 1)

  mg <- mcf(Survr(id, time, event) ~ group, data = MMC)
  expect_equal(sort(unique(mg$group)), sort(levels(MMC$group)))
})

test_that("ggplot2 graphics build without errors", {
  fit <- survfitr(Survr(id, time, event) ~ group,
    data = MMC,
    type = "wang-chang"
  )
  p <- autoplot(fit)
  expect_s3_class(p, "ggplot")
  expect_s3_class(autoplot(fit, fun = "cumhaz", conf.int = FALSE), "ggplot")
  expect_s3_class(autoplot(fit, conf.type = "plain"), "ggplot")
  expect_s3_class(autoplot(wc_fit(surv_mmc)), "ggplot")
  expect_s3_class(plotEstimators(surv_mmc), "ggplot")
  expect_s3_class(
    autoplot(mcf(Survr(id, time, event) ~ group, data = MMC)),
    "ggplot"
  )
  # they must render without errors
  expect_silent(invisible(ggplot2::ggplot_build(p)))
})

test_that("grouped survdiffr supports summary and autoplot", {
  b <- survdiffr(Survr(id, time, event) ~ group,
    data = MMC,
    q = 0.5, B = 50, boot.F = "WC", seed = 11
  )
  expect_s3_class(b, "survdiffr")
  expect_s3_class(b$Males, "boot") # boot.ci compatibility preserved

  s <- summary(b)
  expect_s3_class(s, "summary.survdiffr")
  expect_equal(s$groups$estimate, c(b$Males$t0, b$Females$t0))
  expect_equal(nrow(s$contrasts), 1)
  expect_true(s$contrasts$p.value >= 0 && s$contrasts$p.value <= 1)
  expect_true(s$contrasts$lower <= s$contrasts$upper)
  expect_output(print(s), "Pairwise differences")

  expect_s3_class(autoplot(b), "ggplot")
})

test_that("survrec 1.x names are deprecated but functional", {
  expect_warning(fit <- psh.fit(surv_mmc), "deprecated")
  expect_equal(fit$survfunc, psh_fit(surv_mmc)$survfunc)
  expect_warning(wc.fit(surv_mmc), "deprecated")
  expect_warning(
    fit <- mlefrailty.fit(surv_mmc, alpha.console = FALSE),
    "deprecated"
  )
  expect_equal(fit$alpha, mlefrailty_fit(surv_mmc, alpha.console = FALSE)$alpha)
  expect_warning(q.search(psh_fit(surv_mmc)), "deprecated")
  expect_warning(surv.search(10, c(4, 7), c(0.9, 0.8)), "deprecated")
})

test_that("log-log confidence bands stay inside [0, 1]", {
  fit <- wc_fit(surv_mmc)
  ci <- survrec:::surv_ci(fit$survfunc, fit$std.error)
  expect_true(all(ci$lower >= 0 & ci$lower <= 1, na.rm = TRUE))
  expect_true(all(ci$upper >= 0 & ci$upper <= 1, na.rm = TRUE))
  expect_true(all(ci$lower <= fit$survfunc & fit$survfunc <= ci$upper))
})
