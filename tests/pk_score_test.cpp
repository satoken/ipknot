#include "ip.h"
#include "pk_score.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>

using Pair = std::pair<int, int>;
using Pairs = std::set<Pair>;

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

// Independent oracle: find maximal selected helices, then test whether two
// helices occupy a complete H region with no other paired loop nucleotides.
double oracle(const Pairs& selected, const PKScoreOptions& options) {
  std::vector<std::array<int, 3>> stems;
  std::set<int> paired;
  for (const auto& [i, j] : selected) {
    paired.insert(i);
    paired.insert(j);
    if (selected.count({i - 1, j + 1})) continue;
    int length = 1;
    while (selected.count({i + length, j - length})) ++length;
    if (length >= 2 && length <= options.max_stem) stems.push_back({i, j, length});
  }
  double score = 0;
  for (auto a : stems) for (auto b : stems) {
    if (!(a[0] < b[0] && b[0] < a[1] && a[1] < b[1])) continue;
    std::array<int, 5> g{a[2], b[2], b[0] - a[0] - a[2],
        a[1] - a[2] - b[0] - b[2] + 1, b[1] - b[2] - a[1]};
    bool valid = true;
    for (int k = 2; k < 5; ++k) valid &= g[k] >= 0 && g[k] <= options.max_loop;
    if (!valid) continue;
    if (options.unpaired_loops) for (int position : paired) {
      if ((a[0] + a[2] <= position && position < b[0]) ||
          (b[0] + b[2] <= position && position <= a[1] - a[2]) ||
          (a[1] < position && position <= b[1] - b[2])) valid = false;
    }
    if (valid) score += options.score(g);
  }
  return score;
}

void check_assignment(const Pairs& candidates, const Pairs& selected,
                      PKScoreOptions options, bool swap_levels) {
  IP ip(IP::MAX, 1);
  PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(16));
  for (auto [i, j] : candidates) {
    int level = (i + (swap_levels ? 1 : 0)) % 2;
    int var = ip.make_variable(0);
    levels[level][i].emplace_back(j, var);
    int row = ip.make_constraint(IP::FX, selected.count({i, j}), selected.count({i, j}));
    ip.add_constraint(row, var, 1);
  }
  auto model = add_pk_h_score(ip, levels, options);
  ip.update();
  ip.solve();
  require(std::abs(model.value(ip) - oracle(selected, options)) < 1e-7,
          "H-motif score disagrees with maximal-helix oracle");
  std::vector<int> bpseq(16, -1);
  for (auto [i, j] : selected) bpseq[i] = j, bpseq[j] = i;
  require(std::abs(score_pk_h_structure(bpseq, options) - oracle(selected, options)) < 1e-7,
          "candidate reranking disagrees with maximal-helix oracle");
  for (auto [var, coefficient] : model.terms) {
    double value = ip.get_value(var);
    require(std::abs(value - std::round(value)) < 1e-7,
            "continuous motif auxiliary was not integral at an integer pair assignment");
  }
}

void check_optimization(double coefficient, bool expect_all) {
  IP ip(IP::MAX, 1);
  PKLevelPairs levels(1, std::vector<std::vector<std::pair<unsigned int, int>>>(14));
  std::vector<int> vars;
  for (auto [i, j] : Pairs{{0, 9}, {1, 8}, {4, 13}, {5, 12}}) {
    int v = ip.make_variable(0.1);
    vars.push_back(v);
    levels[0][i].emplace_back(j, v);
  }
  PKScoreOptions options;
  options.intercept = coefficient;
  auto model = add_pk_h_score(ip, levels, options);
  ip.update();
  ip.solve();
  int count = 0;
  for (int v : vars) count += ip.get_value(v) > 0.5;
  require((count == 4) == expect_all, "motif score did not change the optimum as expected");
  require(std::abs(model.value(ip) - (expect_all ? coefficient : 0)) < 1e-7,
          "negative motif score was evaded or positive motif score was charged without its pairs");
}

int main() {
  try {
    std::vector<Pair> pool{{0, 9}, {1, 8}, {2, 7}, {4, 13}, {5, 12}, {3, 14}, {2, 6}, {10, 11}};
    Pairs candidates(pool.begin(), pool.end());
    int checked = 0;
    for (unsigned mask = 0; mask < (1u << pool.size()); ++mask) {
      Pairs selected;
      std::set<int> endpoints;
      bool valid = true;
      for (unsigned k = 0; k < pool.size(); ++k) if (mask & (1u << k)) {
        selected.insert(pool[k]);
        valid &= endpoints.insert(pool[k].first).second;
        valid &= endpoints.insert(pool[k].second).second;
      }
      if (!valid) continue;
      for (double sign : {-1.0, 1.0}) {
        PKScoreOptions options;
        options.intercept = sign * 0.4;
        options.loop_penalty = sign * 0.03;
        options.stem_reward = sign * 0.02;
        options.coax_bonus = sign * 0.1;
        options.table[{2, 2, 2, 2, 2}] = sign * 0.75;
        for (bool unpaired : {false, true}) {
          options.unpaired_loops = unpaired;
          check_assignment(candidates, selected, options, sign > 0);
          ++checked;
        }
      }
    }
    check_optimization(1.0, true);
    check_optimization(-1.0, false);
    // Empty middle/outer loops, asymmetric stems, and finite score domain.
    for (int a = 2; a <= 4; ++a) for (int b = 2; b <= 4; ++b)
      for (int l1 = 0; l1 <= 1; ++l1) for (int l2 = 0; l2 <= 1; ++l2)
        for (int l3 = 0; l3 <= 1; ++l3) {
          int i = 0, k = a + l1, j = k + b + l2 + a - 1, l = j + l3 + b;
          // This fixture needs more than 16 positions in some combinations.
          if (l >= 16) continue;
          Pairs selected;
          for (int d = 0; d < a; ++d) selected.insert({i + d, j - d});
          for (int d = 0; d < b; ++d) selected.insert({k + d, l - d});
          PKScoreOptions options;
          options.intercept = -0.5;
          options.max_stem = 3;
          options.max_loop = 1;
          check_assignment(selected, selected, options, false);
          ++checked;
        }
    IP ip(IP::MAX, 1);
    PKLevelPairs levels(1, std::vector<std::vector<std::pair<unsigned int, int>>>(14));
    for (auto [i, j] : Pairs{{0, 9}, {1, 8}, {4, 13}, {5, 12}})
      levels[0][i].emplace_back(j, ip.make_variable(0));
    auto disabled = add_pk_h_score(ip, levels, PKScoreOptions());
    require(disabled.continuous_variables == 0 && disabled.rows == 0,
            "disabled scoring changed the model");
    std::cout << "Verified " << checked << " fixed assignments plus positive/negative optimization and disabled mode\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
