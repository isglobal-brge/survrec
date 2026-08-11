/// @file wang_chang.h
/// @brief Wang-Chang estimator of the common gap-time survival function
///   (replaces the Fortran subroutines wc2, wcPLE, wc and wcFn).

#ifndef SURVREC_WANG_CHANG_H
#define SURVREC_WANG_CHANG_H

#include <vector>

#include "survrec_data.h"

namespace survrec {

/// Wang-Chang product-limit fit evaluated at `time` (the ordered distinct
/// failure times). `at_risk` and `n_event` are the weighted r(k) and d(k)
/// of Wang & Chang (1999); `std_error` is empty unless requested.
struct WangChangFit {
  std::vector<double> time;
  std::vector<double> at_risk;
  std::vector<double> n_event;
  std::vector<double> surv;
  std::vector<double> std_error;
};

/// `distinct` must be the sorted distinct uncensored gap times of `data`.
/// Subjects with no uncensored gaps contribute their censored gap to the
/// risk set with weight 1; subjects with m_i >= 1 contribute each
/// uncensored gap with weight 1/m_i.
WangChangFit wang_chang(const RecurrentData& data,
                        const std::vector<double>& distinct,
                        bool variance);

}  // namespace survrec

#endif  // SURVREC_WANG_CHANG_H
