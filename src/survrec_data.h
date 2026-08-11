/// @file survrec_data.h
/// @brief Internal representation of a recurrent-event sample.

#ifndef SURVREC_DATA_H
#define SURVREC_DATA_H

#include <cstddef>
#include <vector>

namespace survrec {

/// Recurrent-event sample in subject-grouped form: subject i has m[i]
/// uncensored gap times stored consecutively in `failed` (subjects in the
/// same order as `m`), plus exactly one censored gap in `censored[i]`.
struct RecurrentData {
  std::vector<int> m;
  std::vector<double> failed;
  std::vector<double> censored;

  int n() const { return static_cast<int>(m.size()); }

  int total_events() const {
    int tot = 0;
    for (int mi : m) tot += mi;
    return tot;
  }

  /// offsets[i] is the first position of subject i's gaps inside `failed`
  std::vector<int> offsets() const {
    std::vector<int> off(m.size() + 1, 0);
    for (std::size_t i = 0; i < m.size(); ++i) off[i + 1] = off[i] + m[i];
    return off;
  }
};

}  // namespace survrec

#endif  // SURVREC_DATA_H
