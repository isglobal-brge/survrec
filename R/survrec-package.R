#' survrec: Survival Analysis for Recurrent Event Data
#'
#' Estimation of the survival function of inter-occurrence times for
#' recurrent event data. Three estimators are available: the generalized
#' product-limit estimator of Peña, Strawderman and Hollander (2001), the
#' estimator of Wang and Chang (1999) for correlated inter-occurrence
#' times, and maximum likelihood estimation under a gamma frailty model.
#' Survival quantiles can be compared between groups through several
#' bootstrap schemes.
#'
#' The main entry points are [Survr()] to build the response object,
#' [survfitr()] to estimate survival curves (optionally by group) and
#' [survdiffr()] to bootstrap survival quantiles.
#'
#' @references
#' Peña, E.A., Strawderman, R. and Hollander, M. (2001). Nonparametric
#' Estimation with Recurrent Event Data. *Journal of the American
#' Statistical Association* **96**, 1299--1315.
#'
#' Wang, M.-C. and Chang, S.-H. (1999). Nonparametric Estimation of a
#' Recurrent Survival Function. *Journal of the American Statistical
#' Association* **94**, 146--153.
#'
#' @useDynLib survrec, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom boot boot.ci
#' @importFrom graphics lines plot
#' @importFrom stats median model.extract qnorm terms ave
#' @importFrom rlang .data
#' @importFrom utils combn
#' @importFrom grDevices colorRampPalette
"_PACKAGE"

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
