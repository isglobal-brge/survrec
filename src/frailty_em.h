/// @file frailty_em.h
/// @brief Gamma-frailty MLE via the EM algorithm (replaces the Fortran
///   subroutines emalgo, estalpha, loglik, mlevalue and SearchForSeed;
///   the IMSL ZXGSP golden-section search and the copied Brent routine are
///   replaced by Brent_fmin from R's C API).

#ifndef SURVREC_FRAILTY_EM_H
#define SURVREC_FRAILTY_EM_H

#include <vector>

#include "at_risk.h"

namespace survrec {

struct EmResult {
  std::vector<double> lambda;
  double alpha = 0.0;
  std::vector<double> frailties;  ///< posterior frailty estimates Zhat_i
  int status = 0;                 ///< 1 converged, -1 maxiter exceeded
  int iter = 0;
};

/// Profile log-likelihood in alpha given the hazard increments lambda.
double profile_loglik(const DistinctFailures& df, const std::vector<int>& m,
                      const std::vector<double>& lambda, double alpha);

/// Maximizer of profile_loglik over [alpha_min, alpha_max] (Brent).
double optimize_alpha(const DistinctFailures& df, const std::vector<int>& m,
                      const std::vector<double>& lambda, double alpha_min,
                      double alpha_max, double tol);

/// EM algorithm for (lambda, alpha) started at (lambda0, alpha0).
/// Convergence when max(||dZhat||, ||dlambda||, min(|dalpha|, |d(1/alpha)|))
/// <= tol. Unlike the Fortran version, `tol` is never modified: the inner
/// alpha search uses its own fixed tolerance.
EmResult em_frailty(const DistinctFailures& df, const std::vector<int>& m,
                    const std::vector<double>& lambda0, double alpha0,
                    double tol, int maxiter);

/// Marginal survival estimate S(t_j) = (alpha / (alpha + Lambda_j))^alpha.
std::vector<double> mle_survival(const std::vector<double>& lambda,
                                 double alpha);

}  // namespace survrec

#endif  // SURVREC_FRAILTY_EM_H
