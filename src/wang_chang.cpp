#include "wang_chang.h"

#include <algorithm>
#include <cmath>

namespace survrec {

// Port of the Fortran wc2 with two changes:
//  * the undefined read cen(i,j) outside its defining loop (subjects whose
//    only gap is the censored one never contribute to d(k) -- their event
//    indicator is 0 by construction) is removed;
//  * r/d and the variance are accumulated in one sweep over the ordered
//    distinct times with per-subject pointers instead of recomputing every
//    inner sum from scratch (same sums, O(n * ndistinct) instead of
//    O(n * ndistinct^2 * m)).
WangChangFit wang_chang(const RecurrentData& data,
                        const std::vector<double>& distinct,
                        bool variance) {
  const int n = data.n();
  const int nd = static_cast<int>(distinct.size());
  const std::vector<int> off = data.offsets();

  WangChangFit fit;
  fit.time = distinct;
  fit.at_risk.assign(nd, 0.0);
  fit.n_event.assign(nd, 0.0);
  fit.surv.assign(nd, 0.0);

  // per-subject sorted gaps and weight 1/m_i (subjects with no uncensored
  // gap enter through their censored time with weight 1)
  std::vector<std::vector<double> > gaps(n);
  std::vector<double> wt(n);
  for (int i = 0; i < n; ++i) {
    gaps[i].assign(data.failed.begin() + off[i],
                   data.failed.begin() + off[i + 1]);
    std::sort(gaps[i].begin(), gaps[i].end());
    wt[i] = data.m[i] > 0 ? 1.0 / data.m[i] : 1.0;
  }

  // risk_i(k) for the current k, advanced with a pointer per subject
  std::vector<std::size_t> less(n, 0);
  auto risk_at = [&](int i, int k) -> double {
    if (data.m[i] == 0) return data.censored[i] >= distinct[k] ? 1.0 : 0.0;
    while (less[i] < gaps[i].size() && gaps[i][less[i]] < distinct[k])
      ++less[i];
    return (data.m[i] - static_cast<double>(less[i])) * wt[i];
  };

  for (int k = 0; k < nd; ++k) {
    for (int i = 0; i < n; ++i) {
      fit.at_risk[k] += risk_at(i, k);
      if (data.m[i] > 0) {
        // ties: number of subject i's gaps equal to distinct[k]
        std::size_t p = less[i];
        while (p < gaps[i].size() && gaps[i][p] == distinct[k]) ++p;
        fit.n_event[k] += (p - less[i]) * wt[i];
      }
    }
  }

  double s = 1.0;
  for (int k = 0; k < nd; ++k) {
    if (fit.at_risk[k] > 0.0) s *= 1.0 - fit.n_event[k] / fit.at_risk[k];
    fit.surv[k] = s;
  }

  if (!variance) return fit;

  // var(k) = sqrt( sum_i (w_i(k) - fai_i(k))^2 ) * S(k)  with
  //   w_i(k)   = sum_{k2 <= k} risk_i(k2) d(k2) / r(k2)^2
  //   fai_i(k) = sum_{gaps g_ij < t_k} 1 / (m_i r(pos(g_ij)))
  // both accumulated incrementally in k.
  fit.std_error.assign(nd, 0.0);
  std::vector<double> w(n, 0.0), fai(n, 0.0);
  std::vector<std::size_t> gap_ptr(n, 0);
  std::fill(less.begin(), less.end(), 0);

  for (int k = 0; k < nd; ++k) {
    // fai uses gaps strictly below t_k: advance before adding w's k-th term
    for (int i = 0; i < n; ++i) {
      while (gap_ptr[i] < gaps[i].size() && gaps[i][gap_ptr[i]] < distinct[k]) {
        std::size_t pos = std::lower_bound(distinct.begin(), distinct.end(),
                                           gaps[i][gap_ptr[i]]) -
                          distinct.begin();
        fai[i] += wt[i] / fit.at_risk[pos];
        ++gap_ptr[i];
      }
    }
    const double dk_over_rk2 =
        fit.n_event[k] / (fit.at_risk[k] * fit.at_risk[k]);
    double phi = 0.0;
    for (int i = 0; i < n; ++i) {
      w[i] += risk_at(i, k) * dk_over_rk2;
      phi += (w[i] - fai[i]) * (w[i] - fai[i]);
    }
    fit.std_error[k] = std::sqrt(phi) * fit.surv[k];
  }

  return fit;
}

}  // namespace survrec
