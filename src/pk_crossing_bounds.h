#ifndef IPKNOT_PK_CROSSING_BOUNDS_H
#define IPKNOT_PK_CROSSING_BOUNDS_H

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <vector>

struct PKCrossingBoundPair {
  int left, right;
  double weight; // nonnegative magnitude for one sign of the crossing score
};

// Maximum weight of a noncrossing, base-disjoint subset of pairs crossing
// (left,right). Each side forms a strictly nested chain. A left-side chain
// and a right-side chain coexist only if max(left.right) < min(right.left).
// Weighted dominance queries give O(m log m) time and O(m) memory.
inline double pk_crossing_chain_bound(int left, int right,
                                     std::vector<PKCrossingBoundPair> pairs) {
  if (pairs.empty()) return 0;
  std::sort(pairs.begin(), pairs.end(), [](const auto& a, const auto& b) {
    return a.left > b.left || (a.left == b.left && a.right < b.right);
  });
  std::vector<int> coordinates;
  double total = 0;
  for (const auto& pair : pairs) {
    coordinates.push_back(pair.right);
    total += pair.weight;
  }
  std::sort(coordinates.begin(), coordinates.end());
  coordinates.erase(std::unique(coordinates.begin(), coordinates.end()), coordinates.end());
  std::vector<double> tree(coordinates.size() + 1), chain(pairs.size());
  auto query = [&](size_t index) {
    double value = 0;
    for (; index; index -= index & -index) value = std::max(value, tree[index]);
    return value;
  };
  for (size_t first = 0; first < pairs.size();) {
    size_t last = first + 1;
    while (last < pairs.size() && pairs[last].left == pairs[first].left) ++last;
    // Batch equal left endpoints so they cannot be used together.
    for (size_t i = first; i < last; ++i) {
      const size_t index = std::lower_bound(coordinates.begin(), coordinates.end(), pairs[i].right)
          - coordinates.begin();
      chain[i] = pairs[i].weight + query(index); // strictly smaller inner right
    }
    for (size_t i = first; i < last; ++i) {
      size_t index = std::lower_bound(coordinates.begin(), coordinates.end(), pairs[i].right)
          - coordinates.begin() + 1;
      for (; index < tree.size(); index += index & -index)
        tree[index] = std::max(tree[index], chain[i]);
    }
    first = last;
  }
  std::vector<std::pair<int, double>> prefix;
  for (size_t i = 0; i < pairs.size(); ++i)
    if (pairs[i].left < left && pairs[i].right < right)
      prefix.emplace_back(pairs[i].right, chain[i]);
  std::sort(prefix.begin(), prefix.end());
  double bound = 0;
  for (auto& entry : prefix) {
    bound = std::max(bound, entry.second);
    entry.second = bound;
  }
  for (size_t i = 0; i < pairs.size(); ++i) {
    if (pairs[i].left <= left) continue;
    const auto pos = std::lower_bound(prefix.begin(), prefix.end(), pairs[i].left,
        [](const auto& entry, int coordinate) { return entry.first < coordinate; });
    const double before = pos == prefix.begin() ? 0 : std::prev(pos)->second;
    bound = std::max(bound, before + chain[i]);
  }
  // Protect a valid integer bound against differing floating-point sum orders.
  const double roundoff = 16 * std::numeric_limits<double>::epsilon() * (pairs.size() + 1) * total;
  return std::nextafter(bound + roundoff, std::numeric_limits<double>::infinity());
}

#endif
