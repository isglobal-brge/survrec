#include "at_risk.h"

#include <algorithm>

namespace survrec {

std::vector<double> DistinctFailures::risk_totals() const {
  std::vector<double> tot(ndistinct(), 0.0);
  for (int j = 0; j < ndistinct(); ++j)
    for (int i = 0; i < n; ++i) tot[j] += risk(i, j);
  return tot;
}

DistinctFailures distinct_failures(const RecurrentData& data) {
  DistinctFailures out;
  out.n = data.n();

  std::vector<double> sorted(data.failed);
  std::sort(sorted.begin(), sorted.end());

  for (double t : sorted) {
    if (out.time.empty() || t != out.time.back()) {
      out.time.push_back(t);
      out.n_event.push_back(1);
    } else {
      ++out.n_event.back();
    }
  }

  const int nd = out.ndistinct();
  const std::vector<int> off = data.offsets();
  out.at_risk.assign(static_cast<std::size_t>(out.n) * nd, 0);

  for (int i = 0; i < out.n; ++i) {
    // subject's gaps sorted ascending: #{gaps >= t_j} = m_i - #{gaps < t_j}
    std::vector<double> gaps(data.failed.begin() + off[i],
                             data.failed.begin() + off[i + 1]);
    std::sort(gaps.begin(), gaps.end());
    std::size_t less = 0;
    for (int j = 0; j < nd; ++j) {
      while (less < gaps.size() && gaps[less] < out.time[j]) ++less;
      int risk = data.m[i] - static_cast<int>(less);
      if (data.censored[i] >= out.time[j]) ++risk;
      out.at_risk[static_cast<std::size_t>(i) +
                  static_cast<std::size_t>(out.n) * j] = risk;
    }
  }

  return out;
}

}  // namespace survrec
