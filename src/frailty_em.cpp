#include "frailty_em.h"

#include <algorithm>
#include <cfloat>
#include <cmath>

namespace survrec {

namespace {

// tolerance of the inner Brent search for alpha within each EM iteration
const double kAlphaTol = 1e-8;

/// Brent's method (golden-section search with successive parabolic
/// interpolation) for the minimum of f on [ax, bx]. Endpoints are never
/// evaluated. Written from the published algorithm (Brent, 1973, ch. 5);
/// same scheme R's optimize() uses. Brent_fmin from R's C code is not part
/// of the exported API, hence this implementation.
template <typename F>
double brent_minimize(double ax, double bx, F f, double tol) {
  const double golden = 0.5 * (3.0 - std::sqrt(5.0));
  const double eps = std::sqrt(DBL_EPSILON);

  double a = ax, b = bx;
  double x = a + golden * (b - a), w = x, v = x;
  double fx = f(x), fw = fx, fv = fx;
  double d = 0.0, e = 0.0;

  for (;;) {
    const double xm = 0.5 * (a + b);
    const double tol1 = eps * std::fabs(x) + tol / 3.0;
    const double tol2 = 2.0 * tol1;
    if (std::fabs(x - xm) <= tol2 - 0.5 * (b - a)) break;

    bool use_golden = true;
    if (std::fabs(e) > tol1) {
      // try a parabola through x, v, w
      const double r = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0) p = -p;
      q = std::fabs(q);
      const double e_old = e;
      e = d;
      if (std::fabs(p) < std::fabs(0.5 * q * e_old) && p > q * (a - x) &&
          p < q * (b - x)) {
        d = p / q;
        double u = x + d;
        if (u - a < tol2 || b - u < tol2) d = x < xm ? tol1 : -tol1;
        use_golden = false;
      }
    }
    if (use_golden) {
      e = (x < xm ? b : a) - x;
      d = golden * e;
    }

    const double u = x + (std::fabs(d) >= tol1 ? d : (d > 0.0 ? tol1 : -tol1));
    const double fu = f(u);

    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }
  return x;
}

}  // namespace

double profile_loglik(const DistinctFailures& df, const std::vector<int>& m,
                      const std::vector<double>& lambda, double alpha) {
  if (alpha <= 0.0) return -DBL_MAX;
  const int n = df.n;
  const int nd = df.ndistinct();

  double c3 = 0.0, q1 = 0.0, q2 = 0.0;
  for (int i = 0; i < n; ++i)
    for (int j = 1; j <= m[i]; ++j) c3 += std::log(alpha + j - 1.0);

  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int j = 0; j < nd; ++j) s += df.risk(i, j) * lambda[j];
    q2 += (alpha + m[i]) * std::log1p(s / alpha);
  }

  for (int j = 0; j < nd; ++j)
    q1 += df.n_event[j] * std::log(lambda[j] / alpha);

  return q1 - q2 + c3;
}

double optimize_alpha(const DistinctFailures& df, const std::vector<int>& m,
                      const std::vector<double>& lambda, double alpha_min,
                      double alpha_max, double tol) {
  return brent_minimize(alpha_min, alpha_max,
                        [&](double alpha) {
                          return -profile_loglik(df, m, lambda, alpha);
                        },
                        tol);
}

EmResult em_frailty(const DistinctFailures& df, const std::vector<int>& m,
                    const std::vector<double>& lambda0, double alpha0,
                    double tol, int maxiter) {
  const int n = df.n;
  const int nd = df.ndistinct();

  EmResult res;
  res.lambda = lambda0;
  res.alpha = alpha0;
  res.frailties.assign(n, 0.0);

  std::vector<double> s(n), zold(n), lold(nd);

  for (;;) {
    ++res.iter;

    for (int i = 0; i < n; ++i) {
      s[i] = 0.0;
      for (int j = 0; j < nd; ++j) s[i] += df.risk(i, j) * res.lambda[j];
    }

    zold = res.frailties;
    for (int i = 0; i < n; ++i)
      res.frailties[i] =
          (1.0 + m[i] / res.alpha) / (1.0 + s[i] / res.alpha);

    lold = res.lambda;
    for (int j = 0; j < nd; ++j) {
      double denom = 0.0;
      for (int i = 0; i < n; ++i) denom += df.risk(i, j) * res.frailties[i];
      res.lambda[j] = df.n_event[j] / denom;
    }

    const double alpha_old = res.alpha;
    res.alpha = optimize_alpha(df, m, res.lambda,
                               std::max(res.alpha - 50.0, 0.0),
                               res.alpha + 50.0, kAlphaTol);

    double dz = 0.0, dl = 0.0;
    for (int i = 0; i < n; ++i)
      dz += (res.frailties[i] - zold[i]) * (res.frailties[i] - zold[i]);
    for (int j = 0; j < nd; ++j)
      dl += (res.lambda[j] - lold[j]) * (res.lambda[j] - lold[j]);
    const double da = std::min(std::fabs(alpha_old - res.alpha),
                               std::fabs(1.0 / alpha_old - 1.0 / res.alpha));
    const double dist = std::max(std::sqrt(dz), std::max(std::sqrt(dl), da));

    if (dist <= tol) {
      res.status = 1;
      break;
    }
    if (res.iter > maxiter) {
      res.status = -1;
      break;
    }
  }

  return res;
}

std::vector<double> mle_survival(const std::vector<double>& lambda,
                                 double alpha) {
  std::vector<double> surv(lambda.size());
  double cum = 0.0;
  for (std::size_t j = 0; j < lambda.size(); ++j) {
    cum += lambda[j];
    surv[j] = std::pow(alpha / (alpha + cum), alpha);
  }
  return surv;
}

}  // namespace survrec
