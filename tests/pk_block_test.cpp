#include "ip.h"
#include "pk_score.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
struct Candidate { int left, right, level, block; };

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

bool crosses(const Candidate& a, const Candidate& b) {
  return (a.left < b.left && b.left < a.right && a.right < b.right) ||
         (b.left < a.left && a.left < b.right && b.right < a.right);
}

bool feasible(const std::vector<Candidate>& pairs, unsigned mask) {
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

// The known fixture contains two distinct H block pairs A-B and A-C.
// B and C overlap and do not form a complete H, although individual contacts
// can cross. Count selected contacts by block independently of the allocator.
double score_oracle(const std::vector<Candidate>& pairs, unsigned mask,
                    double ab, double ac) {
  int selected[3] = {0, 0, 0};
  for (size_t i = 0; i < pairs.size(); ++i)
    if (mask & (1u << i)) ++selected[pairs[i].block];
  return ab * selected[0] * selected[1] / 6.0 +
         ac * selected[0] * selected[2] / 6.0;
}

void check_assignment(const std::vector<Candidate>& candidates, unsigned mask,
                      double ab, double ac, bool enabled) {
  IP ip(IP::MAX, 1);
  int nlevels = 0, length = 0;
  for (const auto& c : candidates) {
    nlevels = std::max(nlevels, c.level + 1);
    length = std::max(length, c.right + 1);
  }
  PKLevelPairs levels(nlevels, std::vector<std::vector<std::pair<unsigned int, int>>>(length));
  std::vector<int> variables;
  for (size_t i = 0; i < candidates.size(); ++i) {
    const auto& c = candidates[i];
    const int v = ip.make_variable(0);
    variables.push_back(v);
    levels[c.level][c.left].emplace_back(c.right, v);
    const int selected = (mask >> i) & 1;
    int row = ip.make_constraint(IP::FX, selected, selected);
    ip.add_constraint(row, v, 1);
  }
  for (int position = 0; position < length; ++position) {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (size_t i = 0; i < candidates.size(); ++i)
      if (candidates[i].left == position || candidates[i].right == position)
        ip.add_constraint(row, variables[i], 1);
  }
  for (size_t a = 0; a < candidates.size(); ++a)
    for (size_t b = a + 1; b < candidates.size(); ++b)
      if (candidates[a].level == candidates[b].level && crosses(candidates[a], candidates[b])) {
        int row = ip.make_constraint(IP::UP, 0, 1);
        ip.add_constraint(row, variables[a], 1);
        ip.add_constraint(row, variables[b], 1);
      }
  PKScoreOptions options;
  options.crossing = true;
  options.fixed_blocks = true;
  options.table[{3, 2, 1, 1, 2}] = ab;
  options.table[{3, 2, 0, 2, 0}] = ac;
  // Optimistic proper sub-stems deliberately have inconsistent scores.
  // Fixed allocation must ignore them in favor of the full blocks above.
  options.table[{2, 2, 2, 2, 2}] = 20;
  options.h_weight = enabled ? 1 : 0;
  auto model = add_pk_h_score(ip, levels, options);
  if (enabled) {
    require(model.candidate_blocks == 3 && model.eligible_block_pairs == 2,
            "full candidate blocks were replaced by overlapping sub-stems");
    require(model.scored_block_pairs == static_cast<size_t>((ab != 0) + (ac != 0)),
            "nonzero fixed block-pair accounting is incorrect");
    require(model.covered_crossing_pairs == static_cast<size_t>(6 * ((ab != 0) + (ac != 0))),
            "fixed block contacts were duplicated or omitted");
  } else {
    require(model.candidate_blocks == 0 && model.crossing_scores.empty(),
            "disabled scoring allocated candidate scores");
  }
  add_pk_crossing_constraints(ip, levels, options, model);
  size_t rows = 0, nonzeros = 0;
  for (const auto& a : candidates) for (int below = 0; below < a.level; ++below) {
    ++rows;
    ++nonzeros;
    for (const auto& b : candidates)
      if (b.level == below && crosses(a, b)) ++nonzeros;
  }
  require(model.crossing_rows == rows && model.crossing_nonzeros == nonzeros &&
          model.tightened_rows == 0 && model.removed_crossings == 0,
          "fixed allocation altered the original witness constraints");
  require(model.continuous_variables <= rows && model.rows <= 4 * model.continuous_variables,
          "fixed allocation created more than one auxiliary per support row");
  ip.update();
  const double expected = enabled ? score_oracle(candidates, mask, ab, ac) : 0;
  require(std::abs(ip.solve() - expected) < 1e-7,
          "fixed block score differs from partial/full block activation oracle");
  require(std::abs(model.value(ip) - expected) < 1e-7,
          "fixed block score did not charge the selected contacts");
}

void check_domain_and_union() {
  IP ip(IP::MAX, 1);
  PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(14));
  // A has three physical pairs. The middle pair appears in both levels;
  // candidate multiplicity must not split or duplicate its maximal chain.
  for (int d = 0; d < 3; ++d)
    levels[1][d].emplace_back(9 - d, ip.make_variable(0));
  levels[0][1].emplace_back(8, ip.make_variable(0));
  for (int d = 0; d < 2; ++d)
    levels[0][4 + d].emplace_back(13 - d, ip.make_variable(0));
  PKScoreOptions options;
  options.crossing = true;
  options.fixed_blocks = true;
  options.intercept = 1;
  auto model = add_pk_h_score(ip, levels, options);
  require(model.candidate_blocks == 2 && model.eligible_block_pairs == 1 &&
          model.covered_crossing_pairs == 6,
          "level union duplicated physical candidate blocks");
  options.max_stem = 2;
  model = add_pk_h_score(ip, levels, options);
  require(model.candidate_blocks == 2 && model.out_of_domain_blocks == 1 &&
          model.eligible_block_pairs == 0 && model.crossing_scores.empty(),
          "out-of-domain full blocks were silently trimmed");
  options.max_stem = 12;
  options.max_loop = 1;
  model = add_pk_h_score(ip, levels, options);
  require(model.eligible_block_pairs == 0 && model.skipped_block_pairs == 1,
          "out-of-domain outer loop received a block score");
  options.max_loop = 30;
  options.max_motifs = 1;
  levels[0][3].emplace_back(11, ip.make_variable(0));
  levels[0][4].emplace_back(10, ip.make_variable(0));
  bool failed = false;
  try { add_pk_h_score(ip, levels, options); }
  catch (const std::runtime_error&) { failed = true; }
  require(failed, "block-pair budget silently dropped a physical motif");
}

void check_energy_integration() {
  for (bool swapped : {false, true}) for (double intercept : {0.0, 6.0}) {
    IP ip(IP::MAX, 1);
    PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(14));
    for (int block = 0; block < 2; ++block)
      for (int d = 0; d < (block == 0 ? 3 : 2); ++d) {
        int v = ip.make_variable(0);
        levels[(block + swapped) % 2][(block == 0 ? 0 : 4) + d].emplace_back(
            (block == 0 ? 9 : 13) - d, v);
        int row = ip.make_constraint(IP::FX, 1, 1);
        ip.add_constraint(row, v, 1);
      }
    PKScoreOptions options;
    options.crossing = true;
    options.fixed_blocks = true;
    options.energy.model = PKLoopEnergyModel::DP;
    options.energy_scale = 0.2;
    options.energy_intercept = intercept;
    options.h_weight = 0.5;
    auto model = add_pk_h_score(ip, levels, options);
    require(model.energy_dp_motifs == 1 && model.energy_cc06_motifs == 0 &&
            model.energy_cc09_motifs == 0 && model.energy_fallback_motifs == 0,
            "DP block source accounting is incorrect");
    add_pk_crossing_constraints(ip, levels, options, model);
    ip.update();
    // Pure H loop U=1+1+2. Use the documented local DP parameters directly,
    // independently of PKLoopEnergy::evaluate and PKScoreOptions::score.
    const double rt = 0.00198720425864083 * (273.15 + 37);
    const double expected = 0.5 * (intercept - 0.2 * (9.8 + 0.1 * 4) / rt);
    require(std::abs(ip.solve() - expected) < 1e-7,
            "DP block allocation duplicated or evaded its signed physical correction");
  }
}
} // namespace

int main() {
  try {
    size_t checked = 0;
    for (bool three_levels : {false, true}) {
      const int upper = three_levels ? 2 : 1;
      const std::vector<Candidate> candidates{{0, 9, upper, 0}, {1, 8, upper, 0},
          {2, 7, upper, 0}, {4, 13, three_levels ? 1 : 0, 1},
          {5, 12, three_levels ? 1 : 0, 1}, {3, 11, 0, 2}, {4, 10, 0, 2}};
      for (const auto& weights : {std::pair<double, double>{2.4, 0.6},
                                  {-2.4, -0.6}, {2.4, -0.6}, {0, 0}})
        for (unsigned mask = 0; mask < (1u << candidates.size()); ++mask)
          if (feasible(candidates, mask)) {
            check_assignment(candidates, mask, weights.first, weights.second, true);
            ++checked;
          }
      for (unsigned mask = 0; mask < (1u << candidates.size()); ++mask)
        if (feasible(candidates, mask)) {
          check_assignment(candidates, mask, 2.4, -0.6, false);
          ++checked;
        }
    }
    check_domain_and_union();
    check_energy_integration();
    std::cout << "Verified " << checked << " fixed-block assignments and domain/union checks\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
