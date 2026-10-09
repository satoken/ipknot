#ifndef IPKNOT_NMR_PAIR_INDEX_H
#define IPKNOT_NMR_PAIR_INDEX_H

#include <algorithm>
#include <array>
#include <initializer_list>
#include <limits>
#include <utility>
#include <vector>

// A sparse, range-pruning index over lexicographically sorted base pairs.
// Storage is O(number of pairs), never sequence_length squared. Queries emit
// indices in the original order, preserving the solver's row/column ordering.
class NMRPairRangeIndex {
public:
  using Pair = std::pair<int, int>;
  static constexpr int low = std::numeric_limits<int>::min();
  static constexpr int high = std::numeric_limits<int>::max();
  struct Rectangle { int left_low, left_high, right_low, right_high; };

  explicit NMRPairRangeIndex(const std::vector<Pair>& pairs) : pairs_(pairs) {
    if (!pairs.empty()) {
      bounds_.resize(4 * pairs.size());
      build(1, 0, pairs.size());
    }
  }

  // Both coordinate ranges are open intervals.
  void append(int left_low, int left_high, int right_low, int right_high,
              std::vector<size_t>& output, size_t* visits = nullptr) const {
    if (pairs_.empty() || left_low >= left_high || right_low >= right_high) return;
    const auto begin = std::upper_bound(pairs_.begin(), pairs_.end(),
        Pair(left_low, high)) - pairs_.begin();
    const auto end = std::lower_bound(pairs_.begin(), pairs_.end(),
        Pair(left_high, low)) - pairs_.begin();
    if (begin < end) query(1, 0, pairs_.size(), begin, end,
                           right_low, right_high, output, visits);
  }

  void append_crossing(const Pair& pair, std::vector<size_t>& output) const {
    append(low, pair.first, pair.first, pair.second, output);
    append(pair.first, pair.second, pair.second, high, output);
  }

  void append_between(const Pair& outer, const Pair& inner,
                      std::vector<size_t>& output) const {
    append(outer.first, inner.first, inner.second, outer.second, output);
  }

  // Query a small union without sorting, duplicate output, or an O(L*L) mask.
  // Blocker predicates use at most six rectangles. Traversal is in pair order.
  void append_union(std::initializer_list<Rectangle> rectangles,
                    std::vector<size_t>& output) const {
    if (pairs_.empty()) return;
    std::array<Query, 6> queries;
    size_t count = 0;
    for (const auto& r : rectangles) {
      if (r.left_low >= r.left_high || r.right_low >= r.right_high) continue;
      const size_t begin = std::upper_bound(pairs_.begin(), pairs_.end(),
          Pair(r.left_low, high)) - pairs_.begin();
      const size_t end = std::lower_bound(pairs_.begin(), pairs_.end(),
          Pair(r.left_high, low)) - pairs_.begin();
      if (begin < end) queries.at(count++) = {begin, end, r.right_low, r.right_high};
    }
    if (count) query_union(1, 0, pairs_.size(), queries, count, output);
  }

private:
  struct Query { size_t begin, end; int right_low, right_high; };

  void build(size_t node, size_t begin, size_t end) {
    if (end - begin == 1) {
      bounds_[node] = {pairs_[begin].second, pairs_[begin].second};
      return;
    }
    const size_t middle = begin + (end - begin) / 2;
    build(node * 2, begin, middle);
    build(node * 2 + 1, middle, end);
    bounds_[node] = {std::min(bounds_[node * 2].first, bounds_[node * 2 + 1].first),
                     std::max(bounds_[node * 2].second, bounds_[node * 2 + 1].second)};
  }

  void query(size_t node, size_t begin, size_t end, size_t query_begin,
             size_t query_end, int right_low, int right_high,
             std::vector<size_t>& output, size_t* visits) const {
    if (visits) ++*visits;
    if (query_end <= begin || end <= query_begin ||
        bounds_[node].second <= right_low || bounds_[node].first >= right_high) return;
    if (query_begin <= begin && end <= query_end &&
        right_low < bounds_[node].first && bounds_[node].second < right_high) {
      for (size_t i = begin; i < end; ++i) output.push_back(i);
      return;
    }
    if (end - begin == 1) return;
    const size_t middle = begin + (end - begin) / 2;
    query(node * 2, begin, middle, query_begin, query_end,
          right_low, right_high, output, visits);
    query(node * 2 + 1, middle, end, query_begin, query_end,
          right_low, right_high, output, visits);
  }

  void query_union(size_t node, size_t begin, size_t end,
                   const std::array<Query, 6>& queries, size_t count,
                   std::vector<size_t>& output) const {
    bool intersects = false;
    for (size_t q = 0; q < count; ++q) {
      const auto& query = queries[q];
      if (query.end <= begin || end <= query.begin ||
          bounds_[node].second <= query.right_low ||
          bounds_[node].first >= query.right_high) continue;
      if (query.begin <= begin && end <= query.end &&
          query.right_low < bounds_[node].first &&
          bounds_[node].second < query.right_high) {
        for (size_t i = begin; i < end; ++i) output.push_back(i);
        return;
      }
      intersects = true;
    }
    if (!intersects || end - begin == 1) return;
    const size_t middle = begin + (end - begin) / 2;
    query_union(node * 2, begin, middle, queries, count, output);
    query_union(node * 2 + 1, middle, end, queries, count, output);
  }

  const std::vector<Pair>& pairs_;
  std::vector<Pair> bounds_;
};

#endif
