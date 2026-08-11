# survrec

<!-- badges: start -->
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
<!-- badges: end -->

Survival analysis for recurrent event data.

`survrec` estimates the survival function of the time between occurrences of
a recurrent event — repeated hospitalizations, tumour relapses, successive
failures of a machine — from censored gap-time data. It implements three
estimators: the generalized product-limit estimator of Peña, Strawderman and
Hollander (2001), valid for independent gap times; the estimator of Wang and
Chang (1999), which remains consistent when gap times are correlated within
subjects; and maximum likelihood estimation under a gamma frailty model.
Survival quantiles can be compared between groups through several bootstrap
schemes.

Version 2.0.0 is a complete rewrite of the package: the Fortran 77 core has
been replaced by C++ (validated against the original implementation), the
bootstrap is reproducible with `set.seed()` and runs in parallel, and the
package gains `ggplot2` graphics, the mean cumulative function and formal
bootstrap comparisons of survival quantiles. The last Fortran-based version
is preserved as the [v1.2-5
release](https://github.com/isglobal-brge/survrec/releases/tag/v1.2-5).

## Installation

Development version:

```r
# install.packages("remotes")
remotes::install_github("isglobal-brge/survrec")
```

## Usage

```r
library(survrec)
data(colon)

# survival of the rehospitalization gap times, by Dukes stage
fit <- survfitr(Survr(hc, time, event) ~ as.factor(dukes),
                data = colon, type = "wang-chang")

fit                                    # summary per group
quantile(fit, probs = c(0.25, 0.5))    # survival quantiles
autoplot(fit)                          # survival curves with confidence bands
```

The three estimators can be compared on the same data — the natural check
for within-subject correlation:

```r
x <- Survr(colon$hc, colon$time, colon$event)
plotEstimators(x)
mlefrailty_fit(x)$alpha   # small alpha = strong association between gap times
```

Beyond estimation, the package provides the mean cumulative function of the
events (`mcf()`), and bootstrap comparisons of survival quantiles between
groups with percentile intervals and p-values:

```r
b <- survdiffr(Survr(hc, time, event) ~ as.factor(dukes),
               data = colon, q = 0.5, seed = 1)
summary(b)
autoplot(b)
```

The vignette walks through a complete analysis:

```r
vignette("survrec")
```

## References

Peña, E.A., Strawderman, R. and Hollander, M. (2001). Nonparametric
estimation with recurrent event data. *Journal of the American Statistical
Association*, 96(456), 1299–1315.
[doi:10.1198/016214501753381922](https://doi.org/10.1198/016214501753381922)

Wang, M.-C. and Chang, S.-H. (1999). Nonparametric estimation of a recurrent
survival function. *Journal of the American Statistical Association*,
94(445), 146–153.
[doi:10.1080/01621459.1999.10473831](https://doi.org/10.1080/01621459.1999.10473831)

González, J.R. and Peña, E.A. (2003). Bootstrapping median survival with
recurrent event data. *IX Conferencia Española de Biometría*, A Coruña,
Spain.

## See also

[`gcmrec`](https://github.com/isglobal-brge/gcmrec): regression models for
recurrent event data with effective age, from the same research group.

## License

GPL (>= 2)
