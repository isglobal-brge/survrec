/// @file at_risk.h
/// @brief Distinct failure times, death counts and at-risk matrix
///   (replaces the Fortran subroutine DistinctFailed).

#ifndef SURVREC_AT_RISK_H
#define SURVREC_AT_RISK_H

#include <cstddef>
#include <vector>

#include "survrec_data.h"

namespace survrec {

/// Ordered distinct failure times with death counts and the subject-by-time
/// at-risk matrix used by the PSH and frailty estimators:
///   risk(i, j) = #{gaps of subject i >= time[j]} + 1{censored[i] >= time[j]}
/// The matrix is stored column-major with n rows.
struct DistinctFailures {
  int n = 0;
  std::vector<double> time;
  std::vector<int> n_event;
  std::vector<int> at_risk;

  int ndistinct() const { return static_cast<int>(time.size()); }

  int risk(int i, int j) const {
    return at_risk[static_cast<std::size_t>(i) +
                   static_cast<std::size_t>(n) * j];
  }

  /// Column sums of the at-risk matrix (total at risk at each time)
  std::vector<double> risk_totals() const;
};

DistinctFailures distinct_failures(const RecurrentData& data);

}  // namespace survrec

#endif  // SURVREC_AT_RISK_H
