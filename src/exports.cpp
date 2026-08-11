/// @file exports.cpp
/// @brief Rcpp interface: conversion of the data coming from R into the
///   internal structures and functions exported to the package.
///
/// All exported functions are internal to the package ("." prefix);
/// they are not part of the public API.

#include <Rcpp.h>

#include "at_risk.h"
#include "bootstrap.h"
#include "frailty_em.h"
#include "survrec_data.h"
#include "wang_chang.h"

using namespace Rcpp;

namespace {

survrec::RecurrentData make_data(const IntegerVector& m,
                                 const NumericVector& failed,
                                 const NumericVector& censored) {
  survrec::RecurrentData data;
  data.m.assign(m.begin(), m.end());
  data.failed.assign(failed.begin(), failed.end());
  data.censored.assign(censored.begin(), censored.end());
  if (static_cast<int>(data.censored.size()) != data.n())
    stop("every subject must have exactly one censored time");
  if (data.total_events() != static_cast<int>(data.failed.size()))
    stop("sum(m) must equal length(failed)");
  return data;
}

survrec::DistinctFailures make_distinct(const NumericVector& time,
                                        const IntegerVector& n_event,
                                        const IntegerMatrix& at_risk) {
  survrec::DistinctFailures df;
  df.n = at_risk.nrow();
  df.time.assign(time.begin(), time.end());
  df.n_event.assign(n_event.begin(), n_event.end());
  df.at_risk.assign(at_risk.begin(), at_risk.end());
  if (at_risk.ncol() != df.ndistinct())
    stop("ncol(at.risk) must equal length(time)");
  return df;
}

}  // namespace

// [[Rcpp::export(".distinct_failed")]]
List distinct_failed(IntegerVector m, NumericVector failed,
                     NumericVector censored) {
  const survrec::RecurrentData data = make_data(m, failed, censored);
  const survrec::DistinctFailures df = survrec::distinct_failures(data);

  IntegerMatrix at_risk(df.n, df.ndistinct());
  std::copy(df.at_risk.begin(), df.at_risk.end(), at_risk.begin());

  return List::create(_["time"] = wrap(df.time),
                      _["n.event"] = wrap(df.n_event),
                      _["at.risk"] = at_risk);
}

// [[Rcpp::export(".wc_fit")]]
List wc_fit(IntegerVector m, NumericVector failed, NumericVector censored,
            NumericVector distinct, bool variance = true) {
  const survrec::RecurrentData data = make_data(m, failed, censored);
  const std::vector<double> dist(distinct.begin(), distinct.end());
  const survrec::WangChangFit fit =
      survrec::wang_chang(data, dist, variance);

  List out = List::create(_["time"] = wrap(fit.time),
                          _["at.risk"] = wrap(fit.at_risk),
                          _["n.event"] = wrap(fit.n_event),
                          _["surv"] = wrap(fit.surv));
  if (variance) out["std.error"] = wrap(fit.std_error);
  return out;
}

// [[Rcpp::export(".search_seed")]]
double search_seed(IntegerVector m, NumericVector time, IntegerVector n_event,
                   IntegerMatrix at_risk, NumericVector lambda,
                   double alpha_min, double alpha_max, double tol) {
  const survrec::DistinctFailures df = make_distinct(time, n_event, at_risk);
  const std::vector<int> mm(m.begin(), m.end());
  const std::vector<double> lam(lambda.begin(), lambda.end());
  return survrec::optimize_alpha(df, mm, lam, alpha_min, alpha_max, tol);
}

// [[Rcpp::export(".em_frailty")]]
List em_frailty(IntegerVector m, NumericVector time, IntegerVector n_event,
                IntegerMatrix at_risk, NumericVector lambda, double alpha,
                double tol, int maxiter) {
  const survrec::DistinctFailures df = make_distinct(time, n_event, at_risk);
  const std::vector<int> mm(m.begin(), m.end());
  const std::vector<double> lambda0(lambda.begin(), lambda.end());
  const survrec::EmResult res =
      survrec::em_frailty(df, mm, lambda0, alpha, tol, maxiter);

  return List::create(_["lambda"] = wrap(res.lambda),
                      _["alpha"] = res.alpha,
                      _["frailties"] = wrap(res.frailties),
                      _["status"] = res.status,
                      _["iter"] = res.iter);
}

// [[Rcpp::export(".mle_survival")]]
NumericVector mle_survival_r(NumericVector lambda, double alpha) {
  const std::vector<double> lam(lambda.begin(), lambda.end());
  return wrap(survrec::mle_survival(lam, alpha));
}

// [[Rcpp::export(".boot_quantile")]]
List boot_quantile(IntegerVector m, NumericVector failed,
                   NumericVector censored, NumericVector tau, int nboot,
                   int plan, double percentile, int nthreads = 0) {
  if (plan < 1 || plan > 7) stop("plan must be between 1 and 7");
  const survrec::RecurrentData data = make_data(m, failed, censored);
  const std::vector<double> tt(tau.begin(), tau.end());
  const survrec::BootResult res =
      survrec::boot_quantile(data, tt, nboot, plan, percentile, nthreads);

  List out = List::create(_["t"] = wrap(res.stat));
  if (!res.alpha.empty()) out["alpha"] = wrap(res.alpha);
  return out;
}
