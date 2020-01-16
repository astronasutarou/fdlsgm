/**
 * @file baseline.cc
 * @brief Base Line Class
 * @author Ryou Ohsawa
 * @year 2020
 */
#include "fdlsgm.h"


namespace fdlsgm {

  const std::vector<baseline>
  find_segments(const std::vector<dls>& pool,
                const size_t& size_limit,
                const parameter& param_insert,
                const parameter& param_reallocate,
                const parameter& param_coalesce)
  {
    accumulator<360> ac;
    for (auto& dls: pool) ac.insert(dls);
    ac.reallocate();
    ac.coalesce();
    std::vector<baseline> ret;
    const size_t n = ac.count_baseline();
    for (size_t i=0; i<n; i++) {
      const auto& b = ac[i];
      if (b.size()>=size_limit) ret.push_back(b);
    }
    return ret;
  }

}
