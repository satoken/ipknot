#ifndef IPKNOT_NMR_CANDIDATE_FILTER_H
#define IPKNOT_NMR_CANDIDATE_FILTER_H

#include <algorithm>
#include <utility>
#include <vector>

// Necessary conditions from the existing endpoint-support, base-uniqueness,
// and coaxial topology rows. Storage is O(sequence length + candidate pairs),
// not a dense pair-by-pair compatibility matrix. Passing does not prove that a
// complete RNA structure exists; failing proves the witness cannot be selected.
class NMRCandidateFilter {
public:
  using Pair = std::pair<int, int>;
  struct CoaxialTopology { Pair closing, first_child, second_child; };

  NMRCandidateFilter(int length, const std::vector<Pair>& pairs, int distance)
      : pairs_(pairs), by_left_(length), by_right_(length), distance_(distance) {
    for (const auto& pair : pairs) {
      by_left_.at(pair.first).push_back(pair);
      by_right_.at(pair.second).push_back(pair);
    }
  }

  bool contains(const Pair& pair) const {
    return std::binary_search(pairs_.begin(), pairs_.end(), pair);
  }

  static bool compatible(const Pair& candidate,
                         const std::vector<Pair>& required,
                         const CoaxialTopology* topology = nullptr) {
    for (const auto& pair : required) {
      if (candidate != pair && shares_base(candidate, pair)) return false;
      if (topology && crosses(candidate, pair)) return false;
    }
    if (topology && encloses(topology->closing, candidate) &&
        (encloses(candidate, topology->first_child) ||
         encloses(candidate, topology->second_child))) return false;
    return true;
  }

  bool has_endpoint_support(const std::vector<Pair>& required,
                            const CoaxialTopology* topology = nullptr) const {
    for (const auto& pair : required) {
      if (!contains(pair) || !compatible(pair, required, topology)) return false;
      // --allow-isolated omits the ordinary endpoint-support inequalities.
      if (distance_ == 0) continue;
      if (!has_neighbor(pair.first, by_left_, required, topology) ||
          !has_neighbor(pair.second, by_right_, required, topology)) return false;
    }
    return true;
  }

private:
  static bool shares_base(const Pair& a, const Pair& b) {
    return a.first == b.first || a.first == b.second ||
           a.second == b.first || a.second == b.second;
  }
  static bool crosses(const Pair& a, const Pair& b) {
    return (a.first < b.first && b.first < a.second && a.second < b.second) ||
           (b.first < a.first && a.first < b.second && b.second < a.second);
  }
  static bool encloses(const Pair& outer, const Pair& inner) {
    return outer.first < inner.first && inner.second < outer.second;
  }
  bool has_neighbor(int position, const std::vector<std::vector<Pair>>& rows,
                    const std::vector<Pair>& required,
                    const CoaxialTopology* topology) const {
    for (int distance = 1; distance <= distance_; ++distance) {
      for (int neighbor : {position - distance, position + distance}) {
        if (neighbor < 0 || neighbor >= static_cast<int>(rows.size())) continue;
        for (const auto& pair : rows[neighbor]) {
          if (compatible(pair, required, topology)) return true;
        }
      }
    }
    return false;
  }

  const std::vector<Pair>& pairs_;
  std::vector<std::vector<Pair>> by_left_, by_right_;
  int distance_;
};

#endif
