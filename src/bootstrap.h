/// @file bootstrap.h
/// @brief Bootstrap of survival quantiles for recurrent event data
///   (replaces the Fortran bootMedian and the boot1/boot6-boot9 plans).
///
/// All randomness comes from R's RNG (unif_rand / rgamma), so results are
/// reproducible with set.seed() -- the original Fortran used the GNU rand()
/// extension, which R could not seed.

#ifndef SURVREC_BOOTSTRAP_H
#define SURVREC_BOOTSTRAP_H

#include <vector>

#include "survrec_data.h"

namespace survrec {

/// Bootstrap plans, numbered as in the original implementation:
///   1: resample subjects (whole trajectories) from the empirical data
///   2: gaps from the PSH estimate of F, subject's own censoring time
///   3: gaps from the PSH estimate of F, censoring times resampled
///   4: gaps from the WC estimate of F, subject's own censoring time
///   5: gaps from the WC estimate of F, censoring times resampled
///   6: semiparametric (gamma frailty MLE), subject's own censoring time
///   7: semiparametric (gamma frailty MLE), censoring times resampled
struct BootResult {
  std::vector<double> stat;   ///< bootstrap quantile per replicate (NA if
                              ///< the replicate had no events; -1 if the
                              ///< survival estimate never falls below q)
  std::vector<double> alpha;  ///< frailty alpha per replicate (plans 6-7)
};

/// tau[i] is subject i's total observation time (sum of all gaps).
/// `percentile` is the survival level q whose associated time is extracted
/// from each replicate (q = 0.5 gives the median). `nthreads` caps the
/// OpenMP threads of the re-estimation phase (0 = let OpenMP decide);
/// resampling itself is serial because it consumes R's RNG.
BootResult boot_quantile(const RecurrentData& data,
                         const std::vector<double>& tau, int nboot, int plan,
                         double percentile, int nthreads);

}  // namespace survrec

#endif  // SURVREC_BOOTSTRAP_H
