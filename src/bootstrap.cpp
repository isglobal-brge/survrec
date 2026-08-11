#include "bootstrap.h"

#include <R.h>
#include <Rmath.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>

#include "at_risk.h"
#include "frailty_em.h"
#include "wang_chang.h"

namespace survrec {

namespace {

/// Replicates with no events are redrawn up to this many times (the original
/// Fortran looped forever on plan 1 and up to 100 times on plans 3/5/7).
const int kMaxRedraws = 100;

/// Number of threads to use: `nthreads` if positive; otherwise whatever
/// OpenMP decides (which honors OMP_NUM_THREADS).
int resolve_threads(int nthreads) {
#ifdef _OPENMP
  const int available = omp_get_max_threads();
  return (nthreads > 0) ? std::min(nthreads, available) : available;
#else
  (void)nthreads;
  return 1;
#endif
}

int sample_int(int n) {
  int j = static_cast<int>(unif_rand() * n);
  return j >= n ? n - 1 : j;
}

/// Inverse-CDF draw: index of the first entry of `cdf` >= u.
int sample_cdf(const std::vector<double>& cdf) {
  const double u = unif_rand();
  std::vector<double>::const_iterator it =
      std::lower_bound(cdf.begin(), cdf.end(), u);
  if (it == cdf.end()) return static_cast<int>(cdf.size()) - 1;
  return static_cast<int>(it - cdf.begin());
}

/// Estimated CDF with an appended sentinel atom at +Inf carrying the mass
/// beyond the last observed time (the Fortran used a 1e30 sentinel); a draw
/// of the sentinel always exceeds the censoring time and closes the subject.
struct GapDistribution {
  std::vector<double> time;
  std::vector<double> cdf;

  GapDistribution() {}

  GapDistribution(const std::vector<double>& t, const std::vector<double>& F)
      : time(t), cdf(F) {
    if (cdf.empty() || cdf.back() < 1.0) {
      time.push_back(std::numeric_limits<double>::infinity());
      cdf.push_back(1.0);
    }
  }
};

/// Draws one subject: gaps iid from `dist` until the accumulated time would
/// exceed tau0; the remainder becomes the censored gap.
void draw_subject(const GapDistribution& dist, double tau0,
                  RecurrentData& out) {
  int count = 0;
  double cum = 0.0;
  for (;;) {
    const double g = dist.time[sample_cdf(dist.cdf)];
    if (cum + g > tau0) break;
    out.failed.push_back(g);
    cum += g;
    ++count;
  }
  out.m.push_back(count);
  out.censored.push_back(tau0 - cum);
}

/// Plan 1: resample whole subject trajectories with replacement.
RecurrentData resample_subjects(const RecurrentData& data) {
  const std::vector<int> off = data.offsets();
  RecurrentData out;
  for (int i = 0; i < data.n(); ++i) {
    const int j = sample_int(data.n());
    out.m.push_back(data.m[j]);
    out.failed.insert(out.failed.end(), data.failed.begin() + off[j],
                      data.failed.begin() + off[j + 1]);
    out.censored.push_back(data.censored[j]);
  }
  return out;
}

/// Plans 2-5: gaps from a nonparametric estimate of F; censoring times are
/// the subjects' own (resample_tau = false) or resampled (true).
RecurrentData resample_from_f(const GapDistribution& dist,
                              const std::vector<double>& tau,
                              bool resample_tau) {
  const int n = static_cast<int>(tau.size());
  RecurrentData out;
  for (int i = 0; i < n; ++i) {
    const double tau0 = resample_tau ? tau[sample_int(n)] : tau[i];
    draw_subject(dist, tau0, out);
  }
  return out;
}

/// Plans 6-7: semiparametric resampling under the gamma frailty model with
/// baseline hazard increments `lambda0` at times `time`.
RecurrentData resample_semiparametric(const std::vector<double>& time,
                                      const std::vector<double>& lambda0,
                                      double alpha,
                                      const std::vector<double>& tau,
                                      bool resample_tau) {
  const int n = static_cast<int>(tau.size());
  RecurrentData out;
  std::vector<double> cdf(time.size());
  for (int i = 0; i < n; ++i) {
    const double tau0 = resample_tau ? tau[sample_int(n)] : tau[i];
    // frailty Z ~ Gamma(shape = alpha, rate = alpha), mean 1
    const double z = rgamma(alpha, 1.0 / alpha);
    double cum = 0.0;
    for (std::size_t k = 0; k < time.size(); ++k) {
      cum += lambda0[k];
      cdf[k] = 1.0 - std::exp(-z * cum);
    }
    draw_subject(GapDistribution(time, cdf), tau0, out);
  }
  return out;
}

/// First time at which `surv` falls to `q` or below; -1 if it never does.
double quantile_at(const std::vector<double>& time,
                   const std::vector<double>& surv, double q) {
  if (surv.empty() || surv.back() > q) return -1.0;
  for (std::size_t k = 0; k < surv.size(); ++k)
    if (surv[k] <= q) return time[k];
  return -1.0;
}

std::vector<double> psh_survival(const DistinctFailures& df) {
  const std::vector<double> total = df.risk_totals();
  std::vector<double> surv(df.ndistinct());
  double s = 1.0;
  for (int j = 0; j < df.ndistinct(); ++j) {
    s *= 1.0 - df.n_event[j] / total[j];
    surv[j] = s;
  }
  return surv;
}

std::vector<double> psh_cdf(const DistinctFailures& df) {
  std::vector<double> f = psh_survival(df);
  for (std::size_t k = 0; k < f.size(); ++k) f[k] = 1.0 - f[k];
  return f;
}

std::vector<double> wc_cdf(const RecurrentData& data,
                           const DistinctFailures& df) {
  const WangChangFit fit = wang_chang(data, df.time, /*variance=*/false);
  std::vector<double> f(fit.surv.size());
  for (std::size_t k = 0; k < f.size(); ++k) f[k] = 1.0 - fit.surv[k];
  return f;
}

/// Full gamma-frailty fit as the Fortran mleALL: initial Nelson-Aalen-type
/// increments, Brent seed search, then EM restarted from up to 7 seeds
/// around the initial alpha until it converges.
EmResult mle_fit(const DistinctFailures& df, const std::vector<int>& m) {
  std::vector<double> lambda0(df.ndistinct());
  const std::vector<double> total = df.risk_totals();
  for (int j = 0; j < df.ndistinct(); ++j) {
    lambda0[j] = df.n_event[j] / total[j];
    if (lambda0[j] == 0.0) lambda0[j] = 1e-6;
  }

  const double alpha_seed =
      optimize_alpha(df, m, lambda0, 0.5, df.time.back(), 1e-6);
  const double del = alpha_seed / 4.0;
  const double seeds[7] = {alpha_seed,           alpha_seed - del,
                           alpha_seed - 2 * del, alpha_seed - 3 * del,
                           alpha_seed + del,     alpha_seed + 2 * del,
                           alpha_seed + 3 * del};

  EmResult res;
  for (int k = 0; k < 7; ++k) {
    res = em_frailty(df, m, lambda0, seeds[k], 1e-5, 500);
    if (res.status == 1) break;
  }
  return res;
}

}  // namespace

// The computation is split in two phases so that the expensive one can run
// in parallel while complying with R's constraints:
//  * resampling consumes R's RNG, which is not thread-safe -> serial;
//  * re-estimation per replicate is deterministic and free of R API
//    calls -> OpenMP parallel loop.
BootResult boot_quantile(const RecurrentData& data,
                         const std::vector<double>& tau, int nboot, int plan,
                         double percentile, int nthreads) {
  const DistinctFailures df0 = distinct_failures(data);
  const bool resample_tau = (plan == 3 || plan == 5 || plan == 7);
  const bool semiparametric = (plan == 6 || plan == 7);

  // gap distribution / frailty fit estimated once from the original sample
  GapDistribution fhat;
  if (plan >= 2 && plan <= 5)
    fhat = GapDistribution(df0.time, plan == 4 || plan == 5
                                         ? wc_cdf(data, df0)
                                         : psh_cdf(df0));
  EmResult mle0;
  if (semiparametric) mle0 = mle_fit(df0, data.m);

  // phase 1 (serial): draw all bootstrap samples
  std::vector<RecurrentData> replicates(nboot);
  for (int ss = 0; ss < nboot; ++ss) {
    for (int attempt = 0; attempt < kMaxRedraws; ++attempt) {
      switch (plan) {
        case 1:
          replicates[ss] = resample_subjects(data);
          break;
        case 2:
        case 3:
        case 4:
        case 5:
          replicates[ss] = resample_from_f(fhat, tau, resample_tau);
          break;
        default:
          replicates[ss] = resample_semiparametric(
              df0.time, mle0.lambda, mle0.alpha, tau, resample_tau);
          break;
      }
      if (replicates[ss].total_events() > 0) break;
    }
  }

  BootResult out;
  out.stat.assign(nboot, NA_REAL);
  if (semiparametric) out.alpha.assign(nboot, NA_REAL);

  // phase 2 (parallel): re-estimate the survival curve on each replicate
  const int threads = resolve_threads(nthreads);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
  for (int ss = 0; ss < nboot; ++ss) {
    const RecurrentData& rep = replicates[ss];
    if (rep.total_events() == 0) continue;  // stat stays NA

    const DistinctFailures df = distinct_failures(rep);
    std::vector<double> surv;
    switch (plan) {
      case 1:
      case 2:
      case 3:
        surv = psh_survival(df);
        break;
      case 4:
      case 5:
        surv = wang_chang(rep, df.time, /*variance=*/false).surv;
        break;
      default: {
        const EmResult fit = mle_fit(df, rep.m);
        surv = mle_survival(fit.lambda, fit.alpha);
        out.alpha[ss] = fit.alpha;
        break;
      }
    }
    out.stat[ss] = quantile_at(df.time, surv, percentile);
  }

  return out;
}

}  // namespace survrec
