#include "ip.h"
#include "pk_score.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

struct Candidate {
  int left, right, level, group;
  double weight;
};

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

bool crosses(const Candidate& a, const Candidate& b) {
  return (a.left < b.left && b.left < a.right && a.right < b.right) ||
         (b.left < a.left && a.left < b.right && b.right < a.right);
}

// Three candidate helices A=(0,9),(1,8), B=(4,13),(5,12),
// C=(2,11),(3,10). A-B has geometry (2,2,2,2,2), whereas A-C
// and C-B have geometry (2,2,0,4,0). The expected compatibility
// and projected coefficients below are specified from this fixture,
// independently of the production motif enumerator and row builder.
void check(double ab, double ac, double scale, double tolerance,
           bool three_levels, bool incomplete_c, bool prefer_c) {
  std::vector<Candidate> candidates{
      {0, 9, three_levels ? 2 : 1, 0, 0.08},
      {1, 8, three_levels ? 2 : 1, 0, 0.03},
      {4, 13, 0, 1, prefer_c ? 0.02 : 0.21},
      {5, 12, 0, 1, prefer_c ? 0.04 : 0.23},
      {2, 11, three_levels ? 1 : 0, 2, prefer_c ? 0.15 : -0.05},
      {3, 10, three_levels ? 1 : 0, 2, prefer_c ? 0.17 : -0.07}};
  if (incomplete_c) candidates.pop_back();
  const double ab_pair = scale * ab / 2;
  const double ac_pair = incomplete_c ? 0 : scale * ac / 2;
  const double a_bonus = three_levels ? std::min(ab_pair, ac_pair)
                                     : std::max(ab_pair, ac_pair);
  auto bonus = [&](const Candidate& candidate) {
    if (candidate.group == 0) return a_bonus;
    return three_levels && candidate.group == 2 ? ac_pair : 0.0;
  };
  auto quality = [&](const Candidate& a, const Candidate& b) {
    return a.group == 0 && b.group == 1 ? ab_pair : ac_pair;
  };

  double expected = -std::numeric_limits<double>::infinity();
  for (unsigned mask = 0; mask < (1u << candidates.size()); ++mask) {
    bool feasible = true;
    double value = 0;
    for (size_t a = 0; a < candidates.size(); ++a) if (mask & (1u << a)) {
      const auto& selected = candidates[a];
      value += selected.weight + bonus(selected);
      for (size_t b = a + 1; b < candidates.size(); ++b) if (mask & (1u << b)) {
        const auto& other = candidates[b];
        if (selected.left == other.left || selected.left == other.right ||
            selected.right == other.left || selected.right == other.right ||
            (selected.level == other.level && crosses(selected, other))) feasible = false;
      }
      for (int below = 0; below < selected.level; ++below) {
        bool witness = false;
        for (size_t b = 0; b < candidates.size(); ++b) if (mask & (1u << b)) {
          const auto& other = candidates[b];
          if (other.level == below && crosses(selected, other) &&
              (bonus(selected) == 0 || quality(selected, other) >= bonus(selected) - tolerance))
            witness = true;
        }
        feasible &= witness;
      }
    }
    if (feasible) expected = std::max(expected, value);
  }

  IP ip(IP::MAX, 1);
  PKLevelPairs levels(three_levels ? 3 : 2,
      std::vector<std::vector<std::pair<unsigned int, int>>>(14));
  std::vector<int> variables;
  for (const auto& candidate : candidates) {
    int var = ip.make_variable(candidate.weight);
    variables.push_back(var);
    levels[candidate.level][candidate.left].emplace_back(candidate.right, var);
  }
  for (size_t a = 0; a < candidates.size(); ++a)
    for (size_t b = a + 1; b < candidates.size(); ++b)
      if (candidates[a].level == candidates[b].level && crosses(candidates[a], candidates[b])) {
        int row = ip.make_constraint(IP::UP, 0, 1);
        ip.add_constraint(row, variables[a], 1);
        ip.add_constraint(row, variables[b], 1);
      }
  PKScoreOptions options;
  options.supported = true;
  options.h_weight = scale;
  options.support_tolerance = tolerance;
  options.table[{2, 2, 2, 2, 2}] = ab;
  options.table[{2, 2, 0, 4, 0}] = ac;
  auto model = add_pk_h_score(ip, levels, options);
  add_pk_crossing_constraints(ip, levels, options, model);
  require(model.continuous_variables == 0 && model.rows == 0,
          "supported projection introduced extra variables or rows");
  size_t original_rows = 0, original_nonzeros = 0;
  for (const auto& a : candidates) for (int below = 0; below < a.level; ++below) {
    ++original_rows;
    ++original_nonzeros;
    for (const auto& b : candidates)
      if (b.level == below && crosses(a, b)) ++original_nonzeros;
  }
  require(model.crossing_rows == original_rows &&
          model.crossing_nonzeros + model.removed_crossings == original_nonzeros,
          "existing crossing rows or their nonzero budget changed");
  ip.update();
  const double actual = ip.solve();
  require(std::abs(actual - expected) < 1e-7,
          "supported decoder optimum disagrees with exhaustive structural oracle");
  double realized = 0;
  for (size_t a = 0; a < candidates.size(); ++a)
    realized += bonus(candidates[a]) * ip.get_value(variables[a]);
  require(std::abs(model.value(ip) - realized) < 1e-7,
          "supported selection score differs from its objective contribution");
}

void check_neutral_alternative() {
  IP ip(IP::MAX, 1);
  PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(14));
  for (auto pair : {std::pair<int, int>{0, 9}, {1, 8}, {2, 7}})
    levels[1][pair.first].emplace_back(pair.second, ip.make_variable(0.1));
  for (auto pair : {std::pair<int, int>{4, 13}, {5, 12}})
    levels[0][pair.first].emplace_back(pair.second, ip.make_variable(0.1));
  PKScoreOptions options;
  options.supported = true;
  options.table[{3, 2, 1, 1, 2}] = -1;
  // Each crossing edge also admits a shorter, neutral candidate geometry.
  // Its compatibility must be 0, rather than the negative longer-stem score.
  auto model = add_pk_h_score(ip, levels, options);
  add_pk_crossing_constraints(ip, levels, options, model);
  require(model.terms.empty() && model.tightened_rows == 0,
          "a neutral geometry was omitted in favor of a negative score");
  ip.update();
  require(std::abs(ip.solve() - 0.5) < 1e-7,
          "neutral alternatives changed the original optimum");
}

void check_split_candidate(bool complete) {
  IP ip(IP::MAX, 1);
  PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(14));
  std::vector<int> upper_a;
  for (auto pair : {std::pair<int, int>{0, 9}, {1, 8}}) {
    int var = ip.make_variable(0);
    upper_a.push_back(var);
    levels[1][pair.first].emplace_back(pair.second, var);
  }
  for (auto pair : {std::pair<int, int>{4, 13}, {5, 12}, {2, 11}})
    levels[0][pair.first].emplace_back(pair.second, ip.make_variable(0));
  // C exists in the union of candidates, but is split across the two levels.
  levels[1][3].emplace_back(10, ip.make_variable(0));
  PKScoreOptions options;
  options.supported = true;
  options.support_complete_stems = complete;
  options.table[{2, 2, 2, 2, 2}] = 1;
  options.table[{2, 2, 0, 4, 0}] = 2;
  auto model = add_pk_h_score(ip, levels, options);
  add_pk_crossing_constraints(ip, levels, options, model);
  for (int var : upper_a) {
    auto term = std::find_if(model.terms.begin(), model.terms.end(),
                            [&](auto value) { return value.first == var; });
    require(term != model.terms.end() && std::abs(term->second - (complete ? 0.5 : 1.0)) < 1e-7,
            "candidate level availability was handled incorrectly");
  }
}

void check_pinned_structure(bool supported, double tolerance) {
  IP ip(IP::MAX, 1);
  PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(14));
  for (int group = 0; group < 3; ++group) {
    const int left = group == 0 ? 0 : (group == 1 ? 4 : 2);
    const int right = group == 0 ? 9 : (group == 1 ? 13 : 11);
    for (int d = 0; d < 2; ++d) {
      int var = ip.make_variable(0);
      levels[group == 0 ? 1 : 0][left + d].emplace_back(right - d, var);
      int row = ip.make_constraint(IP::FX, group != 1, group != 1);
      ip.add_constraint(row, var, 1);
    }
  }
  PKScoreOptions options;
  options.supported = supported;
  options.projected = !supported;
  options.support_tolerance = tolerance;
  options.table[{2, 2, 2, 2, 2}] = 1;
  options.table[{2, 2, 0, 4, 0}] = 0.2;
  auto model = add_pk_h_score(ip, levels, options);
  add_pk_crossing_constraints(ip, levels, options, model);
  ip.update();
  bool feasible = true;
  try { ip.solve(); } catch (const std::runtime_error&) { feasible = false; }
  require(feasible == (!supported || tolerance >= 0.4),
          "pinned incompatible structure bypassed support restrictions");
}

int main() {
  try {
    int checked = 0;
    for (auto scores : {std::pair<double, double>{1, 0.2}, {-1, -0.2},
                        {1, -0.2}, {0, 0}, {0.2, 1}})
      for (bool three_levels : {false, true})
        for (bool incomplete_c : {false, true})
          for (double tolerance : {0.0, 0.4, 1.0})
            for (bool prefer_c : {false, true}) {
              check(scores.first, scores.second, 1, tolerance,
                    three_levels, incomplete_c, prefer_c);
              ++checked;
            }
    for (double scale : {0.0, 0.25}) {
      check(1, 0.2, scale, 0, false, false, true);
      ++checked;
    }
    check_neutral_alternative();
    check_split_candidate(false);
    check_split_candidate(true);
    check_pinned_structure(false, 0);
    check_pinned_structure(true, 0);
    check_pinned_structure(true, 0.4);
    std::cout << "Verified " << checked
              << " optima against exhaustive crossing-witness oracles\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
