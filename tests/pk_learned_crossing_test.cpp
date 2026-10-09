#include "ip.h"
#include "pk_score.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
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
        witness |= (mask & (1u << b)) && pairs[b].level == below && crosses(pairs[a], pairs[b]);
      if (!witness) return false;
    }
  }
  return true;
}

struct Files {
  std::filesystem::path model;
  std::filesystem::path features;
  Files() {
    const auto token = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto stem = std::filesystem::temp_directory_path() /
        ("ipknot-learned-crossing-" + std::to_string(token));
    model = stem.string() + ".model";
    features = stem.string() + ".tsv";
    std::ofstream output(model);
    output << "IPKNOT_PK_LINEAR_V1\n";
    for (const char* name : pk_learned_feature_names())
      output << name << ' ' << (std::string(name) == "bias" ? 1 : 0) << '\n';
    require(static_cast<bool>(output), "could not write test learned model");
  }
  ~Files() {
    std::error_code error;
    std::filesystem::remove(model, error);
    std::filesystem::remove(features, error);
  }
};

void check_assignment(const std::vector<Candidate>& candidates, unsigned mask,
                      const PKLearnedModel& imported, double delta,
                      bool second_anchor, bool dump_only, const Files& files,
                      double shape = 0, bool simplify = false, bool projected = false) {
  IP ip(IP::MAX, 1);
  int nlevels = 0, length = 0;
  for (const auto& pair : candidates) {
    nlevels = std::max(nlevels, pair.level + 1);
    length = std::max(length, pair.right + 1);
  }
  PKLevelPairs levels(nlevels, std::vector<std::vector<std::pair<unsigned int, int>>>(length));
  PKPosteriorPairs posterior(length + 1);
  std::vector<int> variables;
  for (size_t i = 0; i < candidates.size(); ++i) {
    const auto& pair = candidates[i];
    const int variable = ip.make_variable(0);
    variables.push_back(variable);
    levels[pair.level][pair.left].emplace_back(pair.right, variable);
    const float evidence = pair.block == 0 ? (second_anchor ? 0.4f : 0.9f) :
        pair.block == 1 ? (second_anchor ? 0.95f : 0.15f) : 0.05f;
    posterior[pair.left + 1].emplace_back(pair.right + 1, evidence);
    int row = ip.make_constraint(IP::FX, (mask >> i) & 1, (mask >> i) & 1);
    ip.add_constraint(row, variable, 1);
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
  PKPosteriorContext context(posterior);
  PKScoreOptions options;
  options.crossing = !projected;
  options.projected = projected;
  options.fixed_blocks = true;
  options.learned = imported;
  options.learned.weights[0] = delta;
  options.learned_scale = 0.5;
  options.hybrid_shape = shape != 0;
  options.intercept = shape;
  options.simplify_crossing = simplify;
  if (dump_only) {
    options.learned_scale = 0;
    options.feature_output = files.features.string();
  }
  if (delta == 0 && !dump_only && shape == 0) options.max_motifs = 1;
  auto model = add_pk_h_score(ip, levels, options, &context);
  if (delta != 0 || dump_only || shape != 0)
    require(model.candidate_blocks == 3 && model.eligible_block_pairs == 2,
            "learned allocation changed maximal blocks or full H spans");
  else
    require(model.candidate_blocks == 0 && !options.needs_posterior_context() && !options.enabled(),
            "all-zero imported model still analyzed blocks or hit the score budget");
  require(model.missing_posterior_block_pairs == 0 && model.unsupported_block_pairs == 0,
          "complete available learned interactions were excluded");
  if (dump_only) {
    require(model.block_features.size() == 2 && model.crossing_scores.empty() && model.terms.empty(),
            "dump-only mode added weights or lost native-supported records");
    require(!options.enabled(), "dump-only metadata activated candidate reranking");
  }
  add_pk_crossing_constraints(ip, levels, options, model);
  size_t original_rows = 0, original_nonzeros = 0;
  for (const auto& a : candidates) for (int below = 0; below < a.level; ++below) {
    ++original_rows;
    ++original_nonzeros;
    for (const auto& b : candidates)
      original_nonzeros += b.level == below && crosses(a, b);
  }
  require(model.crossing_rows == original_rows && model.crossing_nonzeros == original_nonzeros &&
          model.tightened_rows == 0 && model.removed_crossings == 0,
          "learned score changed original crossing support conditions");
  require(model.continuous_variables <= original_rows && model.rows <= 4 * model.continuous_variables,
          "learned score exceeded one continuous auxiliary per support row");
  if (projected)
    require(model.continuous_variables == 0 && model.rows == 0 && model.crossing_scores.empty(),
            "projection retained partner products or added variables/constraints");
  if ((dump_only || delta == 0) && shape == 0)
    require(model.continuous_variables == 0 && model.rows == 0,
            "neutral learned score created ILP auxiliaries");
  int selected[3] = {};
  for (size_t i = 0; i < candidates.size(); ++i)
    if (mask & (1u << i)) ++selected[candidates[i].block];
  // The bias-only residual is delta*scale per weak-stem pair. A full anchor
  // contributes once to each weak pair; a partial anchor contributes its
  // selected fraction. B has the higher evidence only in second_anchor mode.
  double expected = dump_only ? 0 : 0.5 * delta * selected[0] *
      (selected[1] / (second_anchor ? 2.0 : 3.0) + selected[2] / 3.0)
      + shape * selected[0] / 3.0 * (selected[1] / 2.0 + selected[2] / 2.0);
  if (projected) {
    // Analytic full-partner look-ahead for the two H candidates A/B and A/C.
    // In the same lower level these compete by signed maximum, even if that
    // best partner is unselected. In different lower levels their scores add.
    const int level[] = {candidates[0].level, candidates[3].level, candidates[5].level};
    const int size[] = {3, 2, 2};
    const double contact[] = {0.5 * delta / (second_anchor ? 2 : 3) + shape / 6,
                              0.5 * delta / 3 + shape / 6};
    expected = 0;
    for (int block = 0; block < 3; ++block) for (int below = 0; below < level[block]; ++below) {
      std::vector<double> potential;
      for (int other = 0; other < 3; ++other)
        if ((block == 0) != (other == 0) && level[other] == below) {
          const int interaction = block == 1 || other == 1 ? 0 : 1;
          if (contact[interaction] != 0)
            potential.push_back(contact[interaction] * size[other]);
        }
      if (!potential.empty())
        expected += selected[block] * *std::max_element(potential.begin(), potential.end());
    }
  }
  ip.update();
  require(std::abs(ip.solve() - expected) < 1e-7,
          "learned signed interaction differs from independent partial-anchor oracle");
  require(std::abs(model.value(ip) - expected) < 1e-7,
          "learned score value differs from selected crossing contact products");
}

void check_export_and_missing(const PKLearnedModel& imported, const Files& files) {
  for (bool missing : {false, true}) {
    IP ip(IP::MAX, 1);
    PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(14));
    PKPosteriorPairs posterior(15);
    for (int block = 0; block < 2; ++block)
      for (int offset = 0; offset < (block == 0 ? 3 : 2); ++offset) {
        const int left = (block == 0 ? 0 : 4) + offset;
        const int right = (block == 0 ? 9 : 13) - offset;
        int variable = ip.make_variable(0);
        levels[block == 0 ? 1 : 0][left].emplace_back(right, variable);
        const int chosen = block == 0 ? offset < 2 : offset < 1;
        int row = ip.make_constraint(IP::FX, chosen, chosen);
        ip.add_constraint(row, variable, 1);
        if (!missing || block != 0 || offset != 1)
          posterior[left + 1].emplace_back(right + 1, block == 0 ? .8f : .2f);
      }
    PKPosteriorContext context(posterior);
    PKScoreOptions options;
    options.crossing = options.fixed_blocks = true;
    options.feature_output = files.features.string();
    auto model = add_pk_h_score(ip, levels, options, &context);
    require(model.missing_posterior_block_pairs == static_cast<size_t>(missing),
            "forced posterior-absent pair was not counted explicitly");
    require(model.block_features.size() == static_cast<size_t>(!missing),
            "export did not abstain from undefined block confidence");
    add_pk_crossing_constraints(ip, levels, options, model);
    ip.update();
    require(std::abs(ip.solve()) < 1e-7, "feature export changed selected pair feasibility/objective");
    if (!missing) {
      model.write_features(files.features.string(), "GGGAACCCUUAGGU", {.125f, .0625f}, ip);
      std::ifstream input(files.features);
      std::string line, block_line;
      while (std::getline(input, line)) if (line.rfind("B\t", 0) == 0) block_line = line;
      std::istringstream fields(block_line);
      std::vector<std::string> columns;
      for (std::string field; std::getline(fields, field, '\t');) columns.push_back(field);
      require(columns.size() == 28 && columns[1] == "0" && columns[4] == "4" &&
              columns[12] == "1" && columns[13] == "0" &&
              columns[26] == "2" && columns[27] == "1",
              "zero-based coordinates, anchor, or selected physical counts are wrong in TSV");
    }
  }
  PKScoreOptions options;
  options.crossing = options.fixed_blocks = true;
  options.learned = imported;
  options.learned_scale = 0;
  require(!options.has_h_score() && !options.needs_posterior_context() && !options.enabled(),
          "zero learned scale did not preserve the disabled formulation");
  options.learned_scale = 1;
  options.learned.weights[0] = std::numeric_limits<double>::quiet_NaN();
  bool rejected = false;
  try { options.validate(); } catch (const std::invalid_argument&) { rejected = true; }
  require(rejected, "manually altered nonfinite learned weight was accepted");
}
void check_hybrid_missing(const PKLearnedModel& imported, bool missing, bool projected = false) {
  IP ip(IP::MAX, 1);
  PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(14));
  PKPosteriorPairs posterior(15);
  for (int block = 0; block < 2; ++block)
    for (int d = 0; d < (block == 0 ? 3 : 2); ++d) {
      const int left = (block == 0 ? 0 : 4) + d;
      const int right = (block == 0 ? 9 : 13) - d;
      int var = ip.make_variable(0);
      levels[block == 0 ? 1 : 0][left].emplace_back(right, var);
      const int chosen = block == 0 ? d < 2 : d < 1;
      int row = ip.make_constraint(IP::FX, chosen, chosen);
      ip.add_constraint(row, var, 1);
      if (!missing || block != 0 || d != 1)
        posterior[left + 1].emplace_back(right + 1, block == 0 ? .8f : .2f);
    }
  PKPosteriorContext context(posterior);
  PKScoreOptions options;
  options.crossing = !projected;
  options.projected = projected;
  options.fixed_blocks = options.hybrid_shape = true;
  options.learned = imported;
  options.learned_scale = .5;
  options.intercept = 1;
  auto model = add_pk_h_score(ip, levels, options, &context);
  add_pk_crossing_constraints(ip, levels, options, model);
  ip.update();
  const double expected = (missing ? 1./3 : 2./3) * (projected ? 2 : 1);
  require(std::abs(ip.solve() - expected) < 1e-7 &&
          std::abs(model.value(ip) - expected) < 1e-7,
          "undefined posterior removed the well-defined hybrid geometry component");
}

void check_sparse_projected_partner(const PKLearnedModel& imported) {
  // B has three physical candidate pairs but only two exist in the lower
  // level. Only one is selected; projection must count two, not one or three.
  for (double delta : {-1., 1.}) {
    IP ip(IP::MAX, 1);
    PKLevelPairs levels(2, std::vector<std::vector<std::pair<unsigned int, int>>>(15));
    PKPosteriorPairs posterior(16);
    for (int block = 0; block < 2; ++block) for (int d = 0; d < 3; ++d) {
      const int left = (block == 0 ? 0 : 4) + d;
      const int right = (block == 0 ? 9 : 14) - d;
      const int level = block == 0 || d == 2 ? 1 : 0;
      const int variable = ip.make_variable(0);
      levels[level][left].emplace_back(right, variable);
      posterior[left + 1].emplace_back(right + 1, block == 0 ? .9f : .2f);
      const int chosen = d == 0;
      int row = ip.make_constraint(IP::FX, chosen, chosen);
      ip.add_constraint(row, variable, 1);
    }
    PKPosteriorContext context(posterior);
    PKScoreOptions options;
    options.projected = options.fixed_blocks = options.hybrid_shape = true;
    options.learned = imported;
    options.learned.weights[0] = delta;
    options.learned_scale = .5;
    options.intercept = .6;
    auto model = add_pk_h_score(ip, levels, options, &context);
    add_pk_crossing_constraints(ip, levels, options, model);
    require(model.candidate_blocks == 2 && model.eligible_block_pairs == 1 &&
            model.terms.size() == 3 && model.crossing_scores.empty() &&
            model.continuous_variables == 0 && model.rows == 0,
            "sparse full-partner projection changed allocation or added auxiliaries");
    ip.update();
    const double expected = 2 * (.5 * delta / 3 + .6 / 9);
    require(std::abs(ip.solve() - expected) < 1e-7 &&
            std::abs(model.value(ip) - expected) < 1e-7,
            "projection used selected or physical partner size instead of level candidates");
  }
}
} // namespace

int main() {
  try {
    Files files;
    PKLearnedModel model;
    model.load(files.model.string());
    size_t checked = 0;
    for (bool three_levels : {false, true}) for (bool swapped : {false, true}) {
      const int A = swapped ? 0 : three_levels ? 2 : 1;
      const int B = swapped ? (three_levels ? 2 : 1) : three_levels ? 1 : 0;
      const int C = swapped ? 1 : 0;
      const std::vector<Candidate> candidates{{0, 9, A, 0}, {1, 8, A, 0}, {2, 7, A, 0},
          {4, 13, B, 1}, {5, 12, B, 1}, {3, 11, C, 2}, {4, 10, C, 2}};
      for (bool second_anchor : {false, true})
        for (double delta : {-1.2, 0.0, 1.2})
          for (unsigned mask = 0; mask < (1u << candidates.size()); ++mask)
            if (feasible(candidates, mask)) {
              check_assignment(candidates, mask, model, delta, second_anchor, false, files);
              ++checked;
            }
      for (bool simplify : {false, true})
        for (double shape : {-0.6, 0.6})
          for (double delta : {-1.2, 0., 1.2})
            for (unsigned mask = 0; mask < (1u << candidates.size()); ++mask)
              if (feasible(candidates, mask)) {
                check_assignment(candidates, mask, model, delta, true, false, files, shape, simplify);
                ++checked;
                if (!simplify) {
                  check_assignment(candidates, mask, model, delta, true, false, files, shape, false, true);
                  ++checked;
                }
              }
      for (unsigned mask = 0; mask < (1u << candidates.size()); ++mask)
        if (feasible(candidates, mask)) {
          check_assignment(candidates, mask, model, 1.2, false, true, files);
          ++checked;
        }
    }
    check_export_and_missing(model, files);
    check_hybrid_missing(model, false);
    check_hybrid_missing(model, true);
    check_hybrid_missing(model, false, true);
    check_hybrid_missing(model, true, true);
    check_sparse_projected_partner(model);
    std::cout << "Verified " << checked << " learned crossing assignments and TSV/missing-evidence checks\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
