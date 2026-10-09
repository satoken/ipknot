#include "nmr_candidate_filter.h"
#include <cassert>
#include <random>

using Pair = NMRCandidateFilter::Pair;
using Topology = NMRCandidateFilter::CoaxialTopology;

static bool encloses(Pair a, Pair b) {
  return a.first < b.first && b.second < a.second;
}
static bool crosses(Pair a, Pair b) {
  return (a.first < b.first && b.first < a.second && a.second < b.second) ||
         (b.first < a.first && a.first < b.second && b.second < a.second);
}
static bool allowed(Pair p, const std::vector<Pair>& required, const Topology* t) {
  for (Pair q : required) {
    if (p != q && (p.first == q.first || p.first == q.second ||
                   p.second == q.first || p.second == q.second)) return false;
    if (t && crosses(p, q)) return false;
  }
  return !t || !encloses(t->closing, p) ||
      (!encloses(p, t->first_child) && !encloses(p, t->second_child));
}
static bool brute(const std::vector<Pair>& pairs, const std::vector<Pair>& required,
                  int distance, const Topology* t) {
  for (Pair p : required) {
    if (std::find(pairs.begin(), pairs.end(), p) == pairs.end() ||
        !allowed(p, required, t)) return false;
    if (!distance) continue;
    bool left = false, right = false;
    for (Pair q : pairs) {
      if (!allowed(q, required, t)) continue;
      const int dl = std::abs(p.first - q.first), dr = std::abs(p.second - q.second);
      left |= dl > 0 && dl <= distance;
      right |= dr > 0 && dr <= distance;
    }
    if (!left || !right) return false;
  }
  return true;
}

int main() {
  const std::vector<Pair> required{{1,23}, {2,9}, {10,20}};
  const Topology topology{required[0], required[1], required[2]};
  std::vector<Pair> pairs = required;
  pairs.insert(pairs.end(), {{0,24}, {3,8}, {11,19}});
  std::sort(pairs.begin(), pairs.end());
  assert(NMRCandidateFilter(26, pairs, 1).has_endpoint_support(required, &topology));
  pairs.erase(std::find(pairs.begin(), pairs.end(), Pair{0,24}));
  assert(!NMRCandidateFilter(26, pairs, 1).has_endpoint_support(required, &topology));
  assert(NMRCandidateFilter(26, pairs, 0).has_endpoint_support(required, &topology));

  std::mt19937 random(20261004);
  for (int length = 0; length <= 40; ++length) {
    std::vector<Pair> sparse;
    for (int i = 0; i < length; ++i)
      for (int j = i + 4; j < length; ++j)
        if (random() % 4 == 0) sparse.emplace_back(i,j);
    for (int distance = 0; distance <= 2; ++distance) {
      NMRCandidateFilter filter(length, sparse, distance);
      for (int trial = 0; trial < 300 && !sparse.empty(); ++trial) {
        std::vector<Pair> wanted;
        for (int n = 0; n < 2 + int(random() % 3); ++n)
          wanted.push_back(sparse[random() % sparse.size()]);
        const Topology t{wanted[0],wanted[1],wanted.back()};
        assert(filter.has_endpoint_support(wanted) == brute(sparse,wanted,distance,nullptr));
        assert(filter.has_endpoint_support(wanted,&t) == brute(sparse,wanted,distance,&t));
      }
    }
  }

  // Exhaustively check feasible subsets: no witness admitted by the original
  // ordinary-stacking rows may be rejected by the necessary-condition filter.
  for (int trial = 0; trial < 100; ++trial) {
    std::vector<Pair> sparse = required;
    for (int n = 0; n < 9; ++n) {
      int i = random() % 21, j = i + 4 + random() % (26-i-4);
      sparse.emplace_back(i,j);
    }
    std::sort(sparse.begin(),sparse.end());
    sparse.erase(std::unique(sparse.begin(),sparse.end()),sparse.end());
    for (int distance = 1; distance <= 2; ++distance) {
      NMRCandidateFilter filter(26,sparse,distance);
      if (filter.has_endpoint_support(required,&topology)) continue;
      for (unsigned mask = 0; mask < (1u << sparse.size()); ++mask) {
        std::vector<Pair> selected;
        for (size_t i = 0; i < sparse.size(); ++i)
          if (mask & (1u << i)) selected.push_back(sparse[i]);
        if (!std::all_of(required.begin(),required.end(),[&](Pair p) {
          return std::find(selected.begin(),selected.end(),p) != selected.end();
        })) continue;
        const bool compatible_subset = std::all_of(selected.begin(),selected.end(),[&](Pair p) {
          return allowed(p,selected,&topology);
        });
        if (!compatible_subset) continue;
        assert(!brute(selected,selected,distance,nullptr));
      }
    }
  }
}
