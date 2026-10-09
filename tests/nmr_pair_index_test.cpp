#include "nmr_pair_index.h"
#include <cassert>
#include <random>

int main() {
  std::mt19937 random(7);
  for (int length = 0; length < 100; ++length) {
    std::vector<NMRPairRangeIndex::Pair> pairs;
    for (int i = 0; i < length; ++i)
      for (int j = i + 1; j < length; ++j)
        if (random() % 5 == 0) pairs.emplace_back(i, j);
    NMRPairRangeIndex index(pairs);
    for (int test = 0; test < 100; ++test) {
      int a = random() % 110 - 5, b = random() % 110 - 5;
      int c = random() % 110 - 5, d = random() % 110 - 5;
      std::vector<size_t> actual, expected;
      index.append(a, b, c, d, actual);
      for (size_t k = 0; k < pairs.size(); ++k)
        if (a < pairs[k].first && pairs[k].first < b &&
            c < pairs[k].second && pairs[k].second < d) expected.push_back(k);
      assert(actual == expected);
      actual.clear(); expected.clear();
      if (a > b) std::swap(a, b);
      index.append_crossing({a, b}, actual);
      for (size_t k = 0; k < pairs.size(); ++k)
        if ((a < pairs[k].first && pairs[k].first < b && b < pairs[k].second) ||
            (pairs[k].first < a && a < pairs[k].second && pairs[k].second < b))
          expected.push_back(k);
      assert(actual == expected);
      actual.clear(); expected.clear();
      const NMRPairRangeIndex::Rectangle rectangles[] = {{a, b, c, d}, {c, d, a, b},
          {-5, a, b, 110}, {b, 110, -5, c}};
      index.append_union({rectangles[0], rectangles[1], rectangles[2], rectangles[3]}, actual);
      for (size_t k = 0; k < pairs.size(); ++k) {
        bool contained = false;
        for (const auto& r : rectangles)
          contained = contained || (r.left_low < pairs[k].first && pairs[k].first < r.left_high &&
                                   r.right_low < pairs[k].second && pairs[k].second < r.right_high);
        if (contained) expected.push_back(k);
      }
      assert(actual == expected);
      actual.clear(); expected.clear();
      index.append_between({a, d}, {b, c}, actual);
      for (size_t k = 0; k < pairs.size(); ++k)
        if (a < pairs[k].first && pairs[k].first < b &&
            c < pairs[k].second && pairs[k].second < d) expected.push_back(k);
      assert(actual == expected);
    }
  }
  // An empty-result query on a long sparse input must prune the entire tree.
  std::vector<NMRPairRangeIndex::Pair> pairs;
  for (int i = 0; i < 100000; ++i) pairs.emplace_back(i, i + 100001);
  NMRPairRangeIndex index(pairs);
  std::vector<size_t> result;
  size_t visits = 0;
  index.append(-1, 100000, -1, 100001, result, &visits);
  assert(result.empty() && visits == 1);
}
