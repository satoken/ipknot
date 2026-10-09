#include "ip.h"
#include "pk_score.h"
#include "pk_crossing_bounds.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>

struct Pair { int left, right, level; double weight; };
using Scores = std::map<std::pair<int, int>, double>;

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

bool crosses(const Pair& a, const Pair& b) {
  return (a.left < b.left && b.left < a.right && a.right < b.right) ||
         (b.left < a.left && a.left < b.right && b.right < a.right);
}

bool feasible(const std::vector<Pair>& pairs, unsigned mask) {
  for (size_t a = 0; a < pairs.size(); ++a) if (mask & (1u << a)) {
    for (size_t b = a + 1; b < pairs.size(); ++b) if (mask & (1u << b)) {
      const auto& p = pairs[a];
      const auto& q = pairs[b];
      if (p.left == q.left || p.left == q.right || p.right == q.left || p.right == q.right ||
          (p.level == q.level && crosses(p, q))) return false;
    }
    for (int below = 0; below < pairs[a].level; ++below) {
      bool witness = false;
      for (size_t b = 0; b < pairs.size(); ++b)
        if ((mask & (1u << b)) && pairs[b].level == below && crosses(pairs[a], pairs[b])) witness = true;
      if (!witness) return false;
    }
  }
  return true;
}

// Direct pair products are an independent score oracle. No shape enumerator,
// aggregate variable or big-M bound is used to compute expected values.
double oracle(const Scores& scores, unsigned mask) {
  double result = 0;
  for (auto [edge, weight] : scores)
    if ((mask & (1u << edge.first)) && (mask & (1u << edge.second))) result += weight;
  return result;
}

void check(const std::vector<Pair>& pairs, const Scores& scores, int fixed_mask,
           bool simplify = false, int mode = 0) {
  double expected = -std::numeric_limits<double>::infinity();
  for (unsigned mask = 0; mask < (1u << pairs.size()); ++mask) {
    if ((fixed_mask >= 0 && mask != static_cast<unsigned>(fixed_mask)) || !feasible(pairs, mask)) continue;
    double value = oracle(scores, mask);
    for (size_t i = 0; i < pairs.size(); ++i)
      if (mask & (1u << i)) value += pairs[i].weight;
    expected = std::max(expected, value);
  }
  require(std::isfinite(expected), "invalid oracle fixture");
  IP ip(IP::MAX, 1);
  int levels_count = 0, length = 0;
  for (auto pair : pairs) {
    levels_count = std::max(levels_count, pair.level + 1);
    length = std::max(length, pair.right + 1);
  }
  PKLevelPairs levels(levels_count, std::vector<std::vector<std::pair<unsigned int, int>>>(length));
  std::vector<int> variables;
  for (size_t i = 0; i < pairs.size(); ++i) {
    int var = ip.make_variable(pairs[i].weight);
    variables.push_back(var);
    levels[pairs[i].level][pairs[i].left].emplace_back(pairs[i].right, var);
    if (fixed_mask >= 0) {
      const int selected = (fixed_mask >> i) & 1;
      int row = ip.make_constraint(IP::FX, selected, selected);
      ip.add_constraint(row, var, 1);
    }
  }
  for (int position = 0; position < length; ++position) {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (size_t i = 0; i < pairs.size(); ++i)
      if (pairs[i].left == position || pairs[i].right == position) ip.add_constraint(row, variables[i], 1);
  }
  for (size_t a = 0; a < pairs.size(); ++a) for (size_t b = a + 1; b < pairs.size(); ++b)
    if (pairs[a].level == pairs[b].level && crosses(pairs[a], pairs[b])) {
      int row = ip.make_constraint(IP::UP, 0, 1);
      ip.add_constraint(row, variables[a], 1);
      ip.add_constraint(row, variables[b], 1);
    }
  PKScoreModel model;
  for (auto [edge, weight] : scores) model.crossing_scores[variables[edge.first]][variables[edge.second]] = weight;
  PKScoreOptions options;
  options.crossing = true;
  options.simplify_crossing = simplify;
  options.normalize_crossing = mode & 1;
  options.tight_crossing_bounds = mode & 2;
  options.crossing_hypograph = mode & 4;
  add_pk_crossing_constraints(ip, levels, options, model);
  size_t original_rows = 0, original_nonzeros = 0;
  for (const auto& p : pairs) for (int below = 0; below < p.level; ++below) {
    ++original_rows;
    ++original_nonzeros;
    for (const auto& q : pairs) if (q.level == below && crosses(p, q)) ++original_nonzeros;
  }
  require(model.crossing_rows == original_rows && model.crossing_nonzeros == original_nonzeros &&
          model.tightened_rows == 0 && model.removed_crossings == 0,
          "aggregate scoring changed the original crossing constraints");
  require(model.continuous_variables <= original_rows && model.rows <= 4 * model.continuous_variables,
          "aggregate auxiliaries grew faster than the original support rows");
  ip.update();
  require(std::abs(ip.solve() - expected) < 1e-7,
          "aggregate crossing optimum differs from direct pair-product oracle");
  unsigned selected = 0;
  for (size_t i = 0; i < pairs.size(); ++i) if (ip.get_value(variables[i]) > 0.5) selected |= 1u << i;
  require(std::abs(model.value(ip) - oracle(scores, selected)) < 1e-7,
          "aggregate score evaded a negative charge or rewarded an absent upper pair");
}

void check_growth(int alternatives, bool simplify) {
  IP ip(IP::MAX, 1);
  PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(alternatives + 10));
  int upper = ip.make_variable(0);
  levels[1][0].emplace_back(8, upper);
  PKScoreModel model;
  int capacity = ip.make_constraint(IP::UP, 0, 1);
  for (int i = 0; i < alternatives; ++i) {
    int other = ip.make_variable(0);
    levels[0][4].emplace_back(9 + i, other);
    ip.add_constraint(capacity, other, 1);
    model.crossing_scores[upper][other] = i % 2 ? -0.2 : 0.3;
  }
  PKScoreOptions options;
  options.crossing = true;
  options.simplify_crossing = simplify;
  add_pk_crossing_constraints(ip, levels, options, model);
  require(model.continuous_variables == (simplify && alternatives == 1 ? 0u : 1u) &&
          model.weighted_crossings == static_cast<size_t>(alternatives),
          "one auxiliary was created for each crossing partner instead of each row");
  ip.update();
  require(std::abs(ip.solve() - 0.3) < 1e-7, "endpoint bounds excluded a valid crossing score");
}

void check_equal_weights(double weight, bool unscored) {
  IP ip(IP::MAX, 1);
  PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(16));
  int upper = ip.make_variable(0);
  levels[1][0].emplace_back(8, upper);
  int fixed = ip.make_constraint(IP::FX, 1, 1);
  ip.add_constraint(fixed, upper, 1);
  int capacity = ip.make_constraint(IP::UP, 0, 1);
  PKScoreModel model;
  for (int i = 0; i < 4; ++i) {
    int other = ip.make_variable(0);
    levels[0][4].emplace_back(10 + i, other);
    ip.add_constraint(capacity, other, 1);
    if (!unscored || i < 3) model.crossing_scores[upper][other] = weight;
  }
  PKScoreOptions options;
  options.crossing = options.simplify_crossing = true;
  add_pk_crossing_constraints(ip, levels, options, model);
  require(model.direct_crossing_rows == (unscored ? 0u : 1u) &&
          model.continuous_variables == (unscored ? 1u : 0u),
          "exclusive constant witnesses were not reduced exactly");
  ip.update();
  const double expected = unscored ? std::max(0., weight) : weight;
  require(std::abs(ip.solve() - expected) < 1e-7 &&
          std::abs(model.value(ip) - expected) < 1e-7,
          "direct crossing term changed a signed score or ignored an unscored witness");
}

void check_shape_selection(double ac_score, bool select_upper) {
  IP ip(IP::MAX, 1);
  PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(14));
  for (int group = 0; group < 3; ++group) {
    int left = group == 0 ? 0 : (group == 1 ? 4 : 2);
    int right = group == 0 ? 9 : (group == 1 ? 13 : 11);
    for (int d = 0; d < 2; ++d) {
      int var = ip.make_variable(0);
      levels[group == 0 ? 1 : 0][left + d].emplace_back(right - d, var);
      const int selected = group == 0 ? select_upper : group == 2;
      int row = ip.make_constraint(IP::FX, selected, selected);
      ip.add_constraint(row, var, 1);
    }
  }
  PKScoreOptions options;
  options.crossing = true;
  options.table[{2, 2, 2, 2, 2}] = 1;
  options.table[{2, 2, 0, 4, 0}] = ac_score;
  auto model = add_pk_h_score(ip, levels, options);
  add_pk_crossing_constraints(ip, levels, options, model);
  ip.update();
  require(std::abs(ip.solve() - (select_upper ? ac_score : 0)) < 1e-7,
          "shape projection rewarded absent A-B instead of the realized A-C crossing");
}

void check_asymmetric_levels(bool swapped) {
  IP ip(IP::MAX, 1);
  PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(14));
  for (int group = 0; group < 2; ++group)
    for (int d = 0; d < (group == 0 ? 3 : 2); ++d) {
      int left = (group == 0 ? 0 : 4) + d;
      int right = (group == 0 ? 9 : 13) - d;
      int var = ip.make_variable(0);
      levels[(group + swapped) % 2][left].emplace_back(right, var);
      int row = ip.make_constraint(IP::FX, 1, 1);
      ip.add_constraint(row, var, 1);
    }
  PKScoreOptions options;
  options.crossing = true;
  options.table[{3, 2, 1, 1, 2}] = 0.9;
  auto model = add_pk_h_score(ip, levels, options);
  add_pk_crossing_constraints(ip, levels, options, model);
  ip.update();
  require(std::abs(ip.solve() - 0.9) < 1e-7,
          "crossing-score normalization depended on which stem occupied the upper level");
}

void check_chain_bounds() {
  std::mt19937 randomizer(20261004);
  for (int trial = 0; trial < 160; ++trial) {
    const Pair upper{7, 14, 1, 0};
    std::vector<Pair> candidates;
    while (candidates.size() < 10) {
      int i = randomizer() % 21, j = randomizer() % 21;
      if (i > j) std::swap(i, j);
      Pair pair{i, j, 0, (1 + randomizer() % 29) / 23.0};
      if (!crosses(upper, pair) || std::any_of(candidates.begin(), candidates.end(),
          [&](const Pair& old) { return old.left == i && old.right == j; })) continue;
      candidates.push_back(pair);
    }
    double expected = 0;
    for (unsigned mask = 0; mask < (1u << candidates.size()); ++mask) {
      if (!feasible(candidates, mask)) continue;
      double value = 0;
      for (size_t i = 0; i < candidates.size(); ++i)
        if (mask & (1u << i)) value += candidates[i].weight;
      expected = std::max(expected, value);
    }
    std::vector<PKCrossingBoundPair> pairs;
    for (const auto& pair : candidates) pairs.push_back({pair.left, pair.right, pair.weight});
    const double actual = pk_crossing_chain_bound(upper.left, upper.right, pairs);
    require(actual >= expected - 1e-13 && actual - expected < 1e-10,
            "chain bound differs from exhaustive noncrossing-subset optimum");
  }
}

int main() {
  try {
    const std::vector<Pair> two{{0, 8, 1, 0.03}, {1, 7, 1, -0.02},
                              {3, 11, 0, 0.08}, {4, 10, 0, 0.02},
                              {3, 12, 0, 0.07}, {5, 9, 0, -0.1}};
    const std::vector<Pair> three{{0, 9, 2, 0.03}, {1, 8, 2, -0.02},
                                {2, 11, 1, 0.08}, {3, 10, 1, 0.02},
                                {4, 13, 0, 0.07}, {5, 12, 0, -0.1}};
    int checked = 0;
    for (int mode = 0; mode < 8; ++mode)
    for (bool simplify : {false, true})
    for (const auto& pairs : {two, three}) for (int pattern = 0; pattern < 4; ++pattern) {
      Scores scores;
      for (size_t a = 0; a < pairs.size(); ++a) for (size_t b = 0; b < pairs.size(); ++b)
        if (pairs[a].level > pairs[b].level && crosses(pairs[a], pairs[b])) {
          double weight = 0.1 * (1 + (a + b) % 4);
          if (pattern == 1 || (pattern == 2 && (a + b) % 2)) weight = -weight;
          if (pattern == 3) weight = 0;
          scores[{a, b}] = weight;
        }
      for (unsigned mask = 0; mask < (1u << pairs.size()); ++mask) if (feasible(pairs, mask)) {
        check(pairs, scores, mask, simplify, mode);
        ++checked;
      }
      check(pairs, scores, -1, simplify, mode);
      ++checked;
    }
    for (bool simplify : {false, true})
      for (int alternatives : {1, 2, 8, 32}) check_growth(alternatives, simplify);
    for (double weight : {-0.2, 0.3})
      for (bool unscored : {false, true}) check_equal_weights(weight, unscored);
    check_shape_selection(0.2, true);
    check_shape_selection(-0.2, true);
    check_shape_selection(0.2, false);
    check_asymmetric_levels(false);
    check_asymmetric_levels(true);
    check_chain_bounds();
    std::cout << "Verified " << checked << " fixed assignments/optima and linear auxiliary growth\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
