# Builds the reference ("gold standard") fixture for the regression tests.
#
# It must be run with the ORIGINAL Fortran implementation (v1.2-5) installed,
# BEFORE the C++ rewrite replaces src/. The resulting file is committed at
# tests/testthat/fixtures/fortran-reference.rds and any reimplementation of
# the numerical core must reproduce these values (up to documented bug fixes).
#
# Usage:
#   R CMD INSTALL -l <libdir> .          # from a tree still containing the Fortran code
#   Rscript tools/build-reference-fits.R <libdir>
#
# Notes:
#  * The deterministic estimators (psh.fit, wc.fit, mlefrailty.fit) are exact
#    references, comparable to tight numerical tolerance.
#  * survdiffr() is driven by the non-standard GNU rand() generator inside the
#    Fortran code, which ignores R's RNG. Its bootstrap replicates are NOT
#    reproducible references; only t0 (deterministic point estimate) and loose
#    distribution summaries are stored, for sanity checks.
#  * Known Fortran bugs (fixed in the C++ port and documented in NEWS.md):
#      - wc2 reads the loop variable `j` outside its defining loop when
#        count(i) <= 1 (undefined behavior affecting d(k) for subjects with a
#        single time).
#      - estalpha overwrites the `tol` argument (aliased with emalgo's tol),
#        silently changing the EM convergence criterion after the first
#        iteration.
#      - emalgo uses implicitly-typed single precision REALs (dalpha1, dalpha2)
#        in the convergence distance.
#    Reference values therefore embed those behaviors; tests against them may
#    need relaxed tolerances or per-component exclusions where the fixes have
#    a numerical effect.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) .libPaths(c(args[[1]], .libPaths()))

library(survrec)

stopifnot(packageVersion("survrec") < "2.0.0")

data(MMC)
data(colon)

ref <- list()

# strip environments/calls so the fixture is stable across sessions
clean <- function(x) {
  x[setdiff(names(x), c("call"))]
}

# ---------------------------------------------------------------------------
# Deterministic estimators, whole samples
# ---------------------------------------------------------------------------

ref$psh_mmc <- clean(psh.fit(Survr(MMC$id, MMC$time, MMC$event)))
ref$wc_mmc  <- clean(wc.fit(Survr(MMC$id, MMC$time, MMC$event)))
ref$mle_mmc <- clean(mlefrailty.fit(Survr(MMC$id, MMC$time, MMC$event),
                                    alpha.console = FALSE))

ref$psh_colon <- clean(psh.fit(Survr(colon$hc, colon$time, colon$event)))
ref$wc_colon  <- clean(wc.fit(Survr(colon$hc, colon$time, colon$event)))
ref$mle_colon <- clean(mlefrailty.fit(Survr(colon$hc, colon$time, colon$event),
                                      alpha.console = FALSE))

# with tvals (exercises surv.search through the fit functions)
tvals_mmc <- c(20, 50, 100, 150, 250)
ref$psh_mmc_tvals <- clean(psh.fit(Survr(MMC$id, MMC$time, MMC$event),
                                   tvals = tvals_mmc))
ref$wc_mmc_tvals  <- clean(wc.fit(Survr(MMC$id, MMC$time, MMC$event),
                                  tvals = tvals_mmc))
ref$mle_mmc_tvals <- clean(mlefrailty.fit(Survr(MMC$id, MMC$time, MMC$event),
                                          tvals = tvals_mmc,
                                          alpha.console = FALSE))

# mlefrailty.fit with fixed alpha and with explicit search interval
ref$mle_mmc_alpha_fixed <- clean(mlefrailty.fit(
  Survr(MMC$id, MMC$time, MMC$event), alpha = 5, alpha.console = FALSE))
ref$mle_mmc_alpha_range <- clean(mlefrailty.fit(
  Survr(MMC$id, MMC$time, MMC$event),
  alpha.min = 1, alpha.max = 100, alpha.console = FALSE))

# ---------------------------------------------------------------------------
# Formula interface with strata
# ---------------------------------------------------------------------------

fit_grouped <- function(type) {
  fit <- survfitr(Survr(id, time, event) ~ group, data = MMC, type = type)
  out <- lapply(fit, clean)
  attributes(out) <- attributes(fit)[c("names", "strata", "group")]
  out
}

ref$survfitr_mmc_group_psh <- fit_grouped("pena-strawderman-hollander")
ref$survfitr_mmc_group_wc  <- fit_grouped("wang-chang")

suppressWarnings({
  fit <- survfitr(Survr(hc, time, event) ~ as.factor(dukes),
                  data = colon, type = "wang-chang")
})
out <- lapply(fit, clean)
attributes(out) <- attributes(fit)[c("names", "strata", "group")]
ref$survfitr_colon_dukes_wc <- out

# ---------------------------------------------------------------------------
# Helper functions on top of the fits
# ---------------------------------------------------------------------------

ref$qsearch_psh_mmc_q50 <- q.search(ref$psh_mmc, q = 0.5)
ref$qsearch_wc_mmc_q50  <- q.search(ref$wc_mmc,  q = 0.5)
ref$qsearch_mle_mmc_q50 <- q.search(ref$mle_mmc, q = 0.5)

ref$survsearch_psh_mmc <- surv.search(tvals_mmc, ref$psh_mmc$time,
                                      ref$psh_mmc$survfunc)

# summary/print numerical cores
ref$summary_psh_mmc <- unclass(summary(ref$psh_mmc))
ref$summary_wc_mmc  <- unclass(summary(ref$wc_mmc))

# ---------------------------------------------------------------------------
# survdiffr: only deterministic pieces + loose distribution summaries
# (Fortran uses GNU rand(); replicates are not a reproducible reference)
# ---------------------------------------------------------------------------

set.seed(20260806)
boot_wc <- survdiffr(Survr(id, time, event) ~ group, data = MMC,
                     q = 0.5, B = 500, boot.F = "WC")
boot_psh <- survdiffr(Survr(id, time, event) ~ group, data = MMC,
                      q = 0.5, B = 500, boot.F = "PSH")

boot_summary <- function(b) {
  lapply(b, function(g) list(
    t0 = g$t0,
    R = g$R,
    t_quantiles = quantile(g$t[g$t > 0], probs = c(.1, .25, .5, .75, .9)),
    prop_undefined = mean(g$t < 0)   # -1 codes "quantile not reached"
  ))
}

ref$survdiffr_mmc_wc  <- boot_summary(boot_wc)
ref$survdiffr_mmc_psh <- boot_summary(boot_psh)

# ---------------------------------------------------------------------------

attr(ref, "meta") <- list(
  survrec_version = as.character(packageVersion("survrec")),
  r_version = R.version.string,
  platform = R.version$platform,
  date = format(Sys.Date()),
  tvals_mmc = tvals_mmc
)

out_dir <- file.path("tests", "testthat", "fixtures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(ref, file.path(out_dir, "fortran-reference.rds"), version = 2)

cat("Written", file.path(out_dir, "fortran-reference.rds"),
    "with", length(ref), "reference cases\n")
