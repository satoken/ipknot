#include "pk_score.h"
#include "pk_crossing_bounds.h"
#include "ip.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

bool PKScoreOptions::has_h_score() const {
  if (needs_posterior_context()) return true;
  if (h_weight == 0) return false;
  if (energy.model != PKLoopEnergyModel::None)
    return energy_scale != 0 || energy_intercept != 0;
  return intercept != 0 || loop_penalty != 0 || stem_reward != 0 ||
         coax_bonus != 0 || !table.empty();
}

bool PKScoreOptions::needs_posterior_context() const {
  return (learned.active() && learned_scale != 0) || !feature_output.empty();
}

bool PKScoreOptions::enabled() const {
  if (hybrid_shape) return has_h_score() || level_penalty != 0;
  if (learned.loaded() || !feature_output.empty())
    return (learned.active() && learned_scale != 0) || level_penalty != 0;
  return has_h_score() || level_penalty != 0;
}

double PKScoreOptions::score(const std::array<int, 5>& g) const {
  if (energy.model != PKLoopEnergyModel::None)
    return h_weight * (energy_intercept - energy_scale * energy.evaluate(g).q);
  auto it = table.find(g);
  if (it != table.end()) return h_weight * it->second;
  return h_weight * (intercept + stem_reward * std::min(g[0], g[1])
      - loop_penalty * (std::log1p(g[2]) + std::log1p(g[3]) + std::log1p(g[4]))
      + (g[3] == 0 ? coax_bonus : 0.0));
}

void PKScoreOptions::validate() const {
  if(!std::isfinite(rank_scale) || rank_scale<0 || (rank_scale && !ranker.loaded()))
    throw std::invalid_argument("Nonnegative PK rank scale requires a ranking model");
  if(ensemble && (ranker.loaded() || !rank_output.empty()))
    throw std::invalid_argument("PK ranking and finite ensemble selection are separate modes");
  if (core_width != 0 && core_width != 2 && core_width != 3)
    throw std::invalid_argument("PK core width must be 0, 2 or 3");
  if (core_width && (!fixed_blocks || (!crossing && !projected)))
    throw std::invalid_argument("PK short cores require crossing or projected block allocation");
  if (best_partner && (!crossing || !fixed_blocks))
    throw std::invalid_argument("PK best partner requires crossing block allocation");
  if (!std::isfinite(ensemble_scale) || ensemble_scale < 0 ||
      !std::isfinite(ensemble_intercept) || !std::isfinite(ensemble_temperature) ||
      ensemble_temperature <= 0 || !std::isfinite(ensemble_threshold) ||
      ensemble_threshold < 0 || ensemble_threshold > 1)
    throw std::invalid_argument("Invalid PK ensemble parameters");
  if (!ensemble && (!ensemble_output.empty() || ensemble_scale != 0 || ensemble_intercept != 0))
    throw std::invalid_argument("PK ensemble parameters require --pk-ensemble");
  if (ensemble && (enabled() || learned.loaded() || !feature_output.empty() || candidate_threshold != -1))
    throw std::invalid_argument("PK ensemble requires unmodified baseline candidate generation");
  for (double w : {level_penalty, intercept, loop_penalty, stem_reward,
                   coax_bonus, h_weight, candidate_threshold, selection_weight,
                   support_tolerance, energy_scale, energy_intercept}) {
    if (!std::isfinite(w)) throw std::invalid_argument("PK score weights must be finite");
  }
  if (!std::isfinite(learned_scale) || learned_scale < 0)
    throw std::invalid_argument("PK learned scale must be finite and nonnegative");
  for (double weight : learned.weights)
    if (!std::isfinite(weight))
      throw std::invalid_argument("PK learned weights must be finite");
  if (hybrid_shape && !learned.loaded())
    throw std::invalid_argument("PK hybrid shape requires a learned model");
  if ((learned.bounded() || learned.dp_conversion() || learned.exclusion()) && (hybrid_shape || !feature_output.empty()))
    throw std::invalid_argument("Bounded PK models cannot use hybrid shape or legacy feature export");
  if ((simplify_crossing || normalize_crossing || tight_crossing_bounds || crossing_hypograph) && !crossing)
    throw std::invalid_argument("PK crossing simplification requires crossing mode");
  if (learned.loaded() || !feature_output.empty()) {
    if ((!crossing && !projected) || !fixed_blocks)
      throw std::invalid_argument("PK learned scoring requires crossing or projected fixed-block allocation");
    if (!feature_output.empty() && !crossing)
      throw std::invalid_argument("PK feature export requires crossing fixed-block allocation");
    if (energy.model != PKLoopEnergyModel::None || energy_scale != 0 || energy_intercept != 0 ||
        level_penalty != 0 || !table.empty() ||
        (!hybrid_shape && (intercept != 0 || loop_penalty != 0 || stem_reward != 0 ||
                           coax_bonus != 0 || h_weight != 1)))
      throw std::invalid_argument("PK learned scoring and feature export cannot be mixed with other PK scores");
  }
  if (h_weight < 0) throw std::invalid_argument("PK H weight must be nonnegative");
  if (support_tolerance < 0) throw std::invalid_argument("PK support tolerance must be nonnegative");
  if (energy_scale < 0) throw std::invalid_argument("PK energy scale must be nonnegative");
  energy.validate();
  if (energy.model != PKLoopEnergyModel::None &&
      (intercept != 0 || loop_penalty != 0 || stem_reward != 0 ||
       coax_bonus != 0 || !table.empty()))
    throw std::invalid_argument("PK energy models cannot be mixed with feature weights or a geometry score table");
  if (energy.model == PKLoopEnergyModel::None &&
      (energy_scale != 0 || energy_intercept != 0))
    throw std::invalid_argument("PK energy weights require an energy model");
  if (candidate_threshold != -1.0 &&
      (candidate_threshold < 0 || candidate_threshold > 1)) {
    throw std::invalid_argument("PK candidate threshold must be -1 or in [0,1]");
  }
  if (max_stem < 2 || max_loop < 0 || max_motifs < 1) {
    throw std::invalid_argument("PK score domain requires max-stem >= 2, max-loop >= 0 and max-motifs >= 1");
  }
  if (fixed_blocks && !crossing && !projected)
    throw std::invalid_argument("PK fixed-block allocation requires crossing or projected mode");
  for (const auto& [g, w] : table) {
    if (g[0] < 2 || g[1] < 2 || g[2] < 0 || g[3] < 0 || g[4] < 0 || !std::isfinite(w)) {
      throw std::invalid_argument("Invalid PK score table geometry or coefficient");
    }
  }
}

void PKScoreOptions::load_table(const std::string& filename) {
  std::ifstream input(filename);
  if (!input) throw std::invalid_argument("Cannot read PK score table: " + filename);
  std::string line;
  int number = 0;
  while (std::getline(input, line)) {
    ++number;
    line.erase(line.find('#') == std::string::npos ? line.size() : line.find('#'));
    if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
    std::istringstream row(line);
    std::array<int, 5> geometry;
    double coefficient;
    std::string extra;
    if (!(row >> geometry[0] >> geometry[1] >> geometry[2] >> geometry[3] >> geometry[4] >> coefficient)
        || (row >> extra) || !table.emplace(geometry, coefficient).second) {
      throw std::invalid_argument("Invalid or duplicate PK score table row " + std::to_string(number));
    }
  }
  if (table.empty()) throw std::invalid_argument("PK score table is empty");
  validate();
}

double PKScoreModel::value(const IP& ip) const {
  double result = 0;
  for (const auto& [var, coefficient] : terms) result += coefficient * ip.get_value(var);
  return result;
}

void PKScoreModel::write_features(const std::string& filename,
                                 const std::string& sequence,
                                 const std::vector<float>& thresholds,
                                 const IP& ip) const {
  std::ofstream output(filename, std::ios::app);
  if (!output) throw std::runtime_error("Cannot append PK features: " + filename);
  output << std::setprecision(17);
  output << "#IPKNOT_PK_FEATURES_V1\nG\t" << sequence << '\t';
  for (size_t i = 0; i < thresholds.size(); ++i) {
    if (i) output << ',';
    output << thresholds[i];
  }
  output << '\t' << candidate_blocks << '\t' << out_of_domain_blocks
         << '\t' << eligible_block_pairs << '\t' << skipped_block_pairs
         << '\t' << missing_posterior_block_pairs << '\t' << unsupported_block_pairs << '\n';
  auto selected = [&](const std::vector<std::vector<int>>& pair_variables) {
    int count = 0;
    for (const auto& variables : pair_variables) {
      double value = 0;
      for (int variable : variables) value += ip.get_value(variable);
      count += value > 0.5;
    }
    return count;
  };
  for (const auto& block : block_features) {
    output << 'B';
    for (int value : block.stems) output << '\t' << value;
    for (int value : block.geometry) output << '\t' << value;
    output << '\t' << block.features.anchor_first << '\t' << block.delta;
    for (double value : block.features.values) output << '\t' << value;
    output << '\t' << selected(block.first_pair_variables)
           << '\t' << selected(block.second_pair_variables) << '\n';
  }
  output << "E\t" << block_features.size() << '\n';
  output.flush();
  if (!output) throw std::runtime_error("Failed to append PK features: " + filename);
}

namespace {
using Pair = std::pair<int, int>;
struct Stem {
  int left;
  int right;
  int length;
};
struct Motif {
  int stem1;
  int stem2;
  double coefficient;
};

class Builder {
public:
  IP& ip;
  PKScoreModel model;
  std::map<Pair, std::vector<int>> pair_vars;
  std::vector<std::vector<int>> incident;
  std::vector<Stem> stems;
  std::map<int, int> stem_vars;
  std::map<Pair, int> interval_vars;

  Builder(IP& ip, const PKLevelPairs& levels) : ip(ip) {
    if (levels.empty()) return;
    incident.resize(levels.front().size());
    for (const auto& level : levels) {
      for (int i = 0; i < static_cast<int>(level.size()); ++i) {
        for (const auto& [j, var] : level[i]) {
          pair_vars[{i, j}].push_back(var);
          incident[i].push_back(var);
          incident[j].push_back(var);
        }
      }
    }
  }

  int row(IP::BoundType bound, double rhs) {
    ++model.rows;
    return ip.make_constraint(bound, rhs, rhs);
  }
  void add(int row, int var, double coefficient) {
    ip.add_constraint(row, var, coefficient);
    ++model.nonzeros;
  }
  int variable(double coefficient = 0) {
    ++model.continuous_variables;
    return ip.make_continuous_variable(coefficient);
  }

  // Presence of every base pair in one candidate stem. No level-specific
  // stem variable is introduced; each pair sum is binary by uniqueness.
  int stem_var(int index) {
    auto found = stem_vars.find(index);
    if (found != stem_vars.end()) return found->second;
    const auto& s = stems[index];
    int h = variable();
    stem_vars[index] = h;
    int lower = row(IP::LO, 1 - s.length);
    add(lower, h, 1);
    for (int d = 0; d < s.length; ++d) {
      const auto& vars = pair_vars.at({s.left + d, s.right - d});
      int upper = row(IP::UP, 0);
      add(upper, h, 1);
      for (int v : vars) {
        add(upper, v, -1);
        add(lower, v, -1);
      }
    }
    // Charge the selected maximal stem exactly once, including span mode
    // where loop emptiness cannot exclude shorter sub-stems for us.
    for (Pair neighbor : {Pair{s.left - 1, s.right + 1},
                          Pair{s.left + s.length, s.right - s.length}}) {
      auto found = pair_vars.find(neighbor);
      if (found == pair_vars.end()) continue;
      int upper = row(IP::UP, 1);
      add(upper, h, 1);
      for (int v : found->second) {
        add(upper, v, 1);
        add(lower, v, 1);
      }
    }
    return h;
  }

  // Empty intervals are shared between motifs. Each nucleotide's paired
  // sum is already <= 1, so its complement is an affine binary literal.
  int interval_var(int left, int right) {
    if (left > right) return -1; // the empty interval is the constant 1
    Pair key{left, right};
    auto found = interval_vars.find(key);
    if (found != interval_vars.end()) return found->second;
    bool has_candidates = false;
    for (int i = left; i <= right; ++i) has_candidates |= !incident[i].empty();
    if (!has_candidates) return -1;
    int u = variable();
    interval_vars[key] = u;
    int lower = row(IP::LO, 1);
    add(lower, u, 1);
    // A pair with both endpoints in the interval occurs twice in the lower
    // row. Consolidate it for GLPK and other sparse matrix backends.
    std::map<int, int> counts;
    for (int i = left; i <= right; ++i) {
      if (incident[i].empty()) continue;
      int upper = row(IP::UP, 1);
      add(upper, u, 1);
      for (int v : incident[i]) {
        add(upper, v, 1);
        ++counts[v];
      }
    }
    for (const auto& [v, count] : counts) add(lower, v, count);
    return u;
  }

  void add_motif(const Motif& motif) {
    const auto& a = stems[motif.stem1];
    const auto& b = stems[motif.stem2];
    std::vector<int> literals{stem_var(motif.stem1), stem_var(motif.stem2)};
    if (unpaired_loops) for (Pair interval : {Pair{a.left + a.length, b.left - 1},
                          Pair{b.left + b.length, a.right - a.length},
                          Pair{a.right + 1, b.right - b.length}}) {
      int var = interval_var(interval.first, interval.second);
      if (var >= 0) literals.push_back(var);
    }
    int z = variable(motif.coefficient);
    model.terms.emplace_back(z, motif.coefficient);
    int lower = row(IP::LO, 1 - static_cast<int>(literals.size()));
    add(lower, z, 1);
    for (int v : literals) {
      int upper = row(IP::UP, 0);
      add(upper, z, 1);
      add(upper, v, -1);
      add(lower, v, -1);
    }
  }
  bool unpaired_loops = true;
};

// A candidate base pair belongs to exactly one maximal antidiagonal chain in
// the union of the level-specific candidate sets. Unlike the sub-stem maximum,
// this allocation uses one fixed geometry for every contact of a block pair.
// It is still a partial-activation score: selecting a of A and b of B pairs
// contributes Phi*a*b/(A*B), without requiring the complete stems or empty gaps.
PKScoreModel fixed_block_scores(Builder& builder, const PKLevelPairs& levels,
                               const PKScoreOptions& options,
                               const PKPosteriorContext* posterior) {
  const int length = levels.front().size();
  const int max_stem = std::min(options.max_stem, length);
  const int max_loop = std::min(options.max_loop, length);
  const bool learned_features = options.needs_posterior_context();
  if (learned_features && (!posterior || posterior->length() != length))
    throw std::invalid_argument("PK learned scoring requires posterior context with matching length");
  std::vector<std::vector<int>> by_left(length);
  std::vector<bool> stem_has_posterior;
  for (const auto& [pair, variables] : builder.pair_vars) {
    if (builder.pair_vars.count({pair.first - 1, pair.second + 1})) continue;
    int run = 1;
    while (pair.first + run < pair.second - run &&
           builder.pair_vars.count({pair.first + run, pair.second - run})) ++run;
    for (const auto [offset, size] : pk_core_segments(run, options.core_width)) {
    ++builder.model.candidate_blocks;
    if (size < 2 || size > max_stem) {
      ++builder.model.out_of_domain_blocks;
      continue;
    }
    by_left[pair.first + offset].push_back(builder.stems.size());
    builder.stems.push_back({pair.first + offset, pair.second - offset, size});
    if (learned_features) {
      bool complete_evidence = true;
      for (int d = 0; d < size; ++d)
        complete_evidence &= posterior->contains(pair.first + offset + d, pair.second - offset - d);
      stem_has_posterior.push_back(complete_evidence);
    }
    }
  }
  std::map<Pair, std::vector<int>> level_vars;
  // The best complete candidate partner in each required lower level. Keep
  // signed maxima: an all-negative set remains a penalty, as in the original
  // sub-stem projection. Zero/unscored contacts do not propose a geometry.
  std::map<Pair, double> projected_rows;
  for (size_t level = 0; level < levels.size(); ++level)
    for (int i = 0; i < length; ++i)
      for (auto [j, var] : levels[level][i]) {
        auto [it, inserted] = level_vars.try_emplace({i, j}, levels.size(), -1);
        it->second[level] = var;
      }
  for (size_t a_index = 0; a_index < builder.stems.size(); ++a_index) {
    const auto& a = builder.stems[a_index];
    const int first = a.left + a.length;
    const int last = std::min({length - 1, first + max_loop,
                              a.right - a.length - 1});
    // The index restricts the outer left loop before visiting any partners;
    // no unrestricted Cartesian product of all stem blocks is materialized.
    for (int k = first; k <= last; ++k) {
      const auto& list = by_left[k];
      auto begin = std::lower_bound(list.begin(), list.end(), a.right + 2,
          [&](int index, int right) { return builder.stems[index].right < right; });
      for (auto it = begin; it != list.end(); ++it) {
        const auto& b = builder.stems[*it];
        if (b.right > a.right + max_stem + max_loop) break;
        const std::array<int, 5> geometry{a.length, b.length,
            b.left - a.left - a.length,
            a.right - a.length - b.left - b.length + 1,
            b.right - b.length - a.right};
        if (geometry[3] < 0 || geometry[3] > options.max_loop ||
            geometry[4] < 0 || geometry[4] > options.max_loop) {
          ++builder.model.skipped_block_pairs;
          continue;
        }
        if (builder.model.eligible_block_pairs >= static_cast<size_t>(options.max_motifs))
          throw std::runtime_error("PK block-pair budget exceeded; reduce the scored domain or explicitly raise --pk-h-max-motifs (no candidates were silently dropped)");
        ++builder.model.eligible_block_pairs;
        double coefficient;
        double shape_coefficient = 0;
        PKLearnedFeatures features;
        bool have_features = false;
        if (learned_features) {
          // Forced NMR pairs can be absent from posterior probabilities.
          // Keep every original variable/constraint and abstain from scoring
          // blocks containing such pairs instead of inventing confidence.
          if (!stem_has_posterior[a_index] || !stem_has_posterior[*it]) {
            ++builder.model.missing_posterior_block_pairs;
            if (!options.hybrid_shape) continue;
            // Missing posterior confidence abstains only from the learned
            // component. The geometric prior remains well-defined.
            coefficient = 0;
          } else {
            features = options.learned.features(*posterior, {a.left, a.right, a.length},
                                                {b.left, b.right, b.length}, geometry);
            have_features = true;
            coefficient = options.learned.loaded()
                ? options.learned_scale * options.learned.score(features) : 0;
          }
          if (options.hybrid_shape) shape_coefficient = options.score(geometry);
        } else if (options.energy.model == PKLoopEnergyModel::None) {
          coefficient = options.score(geometry);
        } else {
          const auto energy = options.energy.evaluate(geometry);
          switch (energy.source) {
            case PKLoopEnergySource::DP: ++builder.model.energy_dp_motifs; break;
            case PKLoopEnergySource::CC06: ++builder.model.energy_cc06_motifs; break;
            case PKLoopEnergySource::CC09: ++builder.model.energy_cc09_motifs; break;
            case PKLoopEnergySource::DPFallback: ++builder.model.energy_fallback_motifs; break;
            case PKLoopEnergySource::None: break;
          }
          coefficient = options.h_weight * (options.energy_intercept - options.energy_scale * energy.q);
        }
        if (!std::isfinite(coefficient) || !std::isfinite(shape_coefficient))
          throw std::invalid_argument("PK block score overflow");
        if (coefficient == 0 && !learned_features) continue;
        // A learned residual is per pair of the weaker stem. A full selected
        // anchor realizes its complete confidence, including when the anchor
        // is assigned to the upper level. Physical geometry scores keep the
        // existing symmetric whole-motif normalization.
        const double contact_score = (learned_features
            ? (options.learned.block_normalized() ? coefficient / a.length / b.length
               : coefficient / (features.anchor_first ? a.length : b.length))
            : coefficient / a.length / b.length)
            + shape_coefficient / a.length / b.length;
        if (!std::isfinite(contact_score)) throw std::invalid_argument("PK contact score overflow");
        if (contact_score != 0) ++builder.model.scored_block_pairs;
        size_t supported_contacts = 0;
        std::map<Pair, double> block_projection;
        for (int da = 0; da < a.length; ++da)
          for (int db = 0; db < b.length; ++db) {
            const auto& a_vars = level_vars.at({a.left + da, a.right - da});
            const auto& b_vars = level_vars.at({b.left + db, b.right - db});
            bool covered = false;
            for (size_t upper = 1; upper < levels.size(); ++upper)
              for (size_t lower = 0; lower < upper; ++lower) {
                if (a_vars[upper] >= 0 && b_vars[lower] >= 0) {
                  if (contact_score != 0) {
                    if (options.projected)
                      block_projection[{a_vars[upper], static_cast<int>(lower)}] += contact_score;
                    else
                      builder.model.crossing_scores[a_vars[upper]].emplace(b_vars[lower], contact_score);
                  }
                  covered = true;
                }
                if (b_vars[upper] >= 0 && a_vars[lower] >= 0) {
                  if (contact_score != 0) {
                    if (options.projected)
                      block_projection[{b_vars[upper], static_cast<int>(lower)}] += contact_score;
                    else
                      builder.model.crossing_scores[b_vars[upper]].emplace(a_vars[lower], contact_score);
                  }
                  covered = true;
                }
              }
            if (covered) {
              ++supported_contacts;
              if (contact_score != 0) ++builder.model.covered_crossing_pairs;
            }
          }
        for (const auto& [row, score] : block_projection) {
          auto [it, inserted] = projected_rows.emplace(row, score);
          if (!inserted) it->second = std::max(it->second, score);
        }
        if (learned_features && supported_contacts == 0)
          ++builder.model.unsupported_block_pairs;
        if (!options.feature_output.empty() && have_features && supported_contacts > 0) {
          PKBlockFeatureRecord record;
          record.stems = {a.left, a.right, a.length, b.left, b.right, b.length};
          record.geometry = geometry;
          record.features = features;
          record.delta = coefficient;
          for (int offset = 0; offset < a.length; ++offset)
            record.first_pair_variables.push_back(builder.pair_vars.at({a.left + offset, a.right - offset}));
          for (int offset = 0; offset < b.length; ++offset)
            record.second_pair_variables.push_back(builder.pair_vars.at({b.left + offset, b.right - offset}));
          builder.model.block_features.push_back(std::move(record));
        }
      }
    }
  }
  // Replace selected-partner contact products by the best potential full
  // partner score. Existing crossing-support constraints are built unchanged;
  // they may be satisfied by a different or only partly selected stem.
  std::map<int, double> projected_coefficients;
  for (const auto& [row, score] : projected_rows)
    projected_coefficients[row.first] += score;
  for (const auto& [variable, coefficient] : projected_coefficients) {
    if (!std::isfinite(coefficient)) throw std::invalid_argument("PK projected score overflow");
    if (coefficient == 0) continue;
    builder.ip.add_objective_coefficient(variable, coefficient);
    builder.model.terms.emplace_back(variable, coefficient);
  }
  builder.model.motifs = builder.model.eligible_block_pairs;
  return std::move(builder.model);
}
} // namespace

PKScoreModel add_pk_h_score(IP& ip, const PKLevelPairs& pairs,
                           const PKScoreOptions& options,
                           const PKPosteriorContext* posterior) {
  if (!options.has_h_score() || pairs.empty() || options.rerank) return {};
  options.validate();
  Builder builder(ip, pairs);
  if (options.fixed_blocks && (options.crossing || options.projected))
    return fixed_block_scores(builder, pairs, options, posterior);
  builder.unpaired_loops = options.unpaired_loops;
  const int length = pairs.front().size();
  const int max_stem = std::min(options.max_stem, length);
  const int max_loop = std::min(options.max_loop, length);
  // At most (max_stem-1) candidates per sparse base-pair candidate. Enumerate
  // every sub-stem in this domain, not just maximal posterior chains.
  std::vector<std::vector<int>> by_left(length);
  for (const auto& [pair, vars] : builder.pair_vars) {
    int run = 1;
    while (run < max_stem && pair.first + run < pair.second - run &&
           builder.pair_vars.count({pair.first + run, pair.second - run})) {
      ++run;
      by_left[pair.first].push_back(builder.stems.size());
      builder.stems.push_back({pair.first, pair.second, run});
    }
  }
  // For each stem only inspect second stems whose left endpoint can leave
  // a loop within the scored domain; no global H x H candidate cross product.
  // Lists are sorted by right endpoint because pair_vars is lexicographic.
  std::vector<Motif> motifs;
  for (int ai = 0; ai < static_cast<int>(builder.stems.size()); ++ai) {
    const auto& a = builder.stems[ai];
    int first = a.left + a.length;
    int last = std::min({length - 1, first + max_loop,
                         a.right - a.length - 1});
    for (int k = first; k <= last; ++k) {
      const auto& list = by_left[k];
      auto start = std::lower_bound(list.begin(), list.end(), a.right + 2,
          [&](int index, int right) { return builder.stems[index].right < right; });
      for (auto it = start; it != list.end(); ++it) {
        const auto& b = builder.stems[*it];
        if (b.right > a.right + max_stem + max_loop) break;
        std::array<int, 5> geometry{a.length, b.length,
            b.left - a.left - a.length,
            a.right - a.length - b.left - b.length + 1,
            b.right - b.length - a.right};
        if (geometry[3] < 0 || geometry[3] > options.max_loop ||
            geometry[4] < 0 || geometry[4] > options.max_loop) continue;
        double coefficient = options.score(geometry);
        if (!std::isfinite(coefficient)) throw std::invalid_argument("PK motif score overflow");
        // Neutral shapes remain alternatives to a negative compatibility in
        // partner-scored modes; omitting them invents a negative best score.
        if (coefficient == 0 && !options.supported && !options.crossing) continue;
        if (motifs.size() >= static_cast<size_t>(options.max_motifs)) {
          throw std::runtime_error("PK motif budget exceeded; reduce the scored domain or explicitly raise --pk-h-max-motifs (no candidates were silently dropped)");
        }
        motifs.push_back({ai, *it, coefficient});
      }
    }
  }
  if (options.supported || options.crossing) {
    // A shape supports a level-specific crossing edge when both variables
    // exist. Optionally require the complete candidate stems in those levels.
    // Actual selection of full stems and loop states remains an approximation.
    std::map<Pair, std::vector<int>> level_vars;
    for (size_t level = 0; level < pairs.size(); ++level)
      for (int i = 0; i < length; ++i)
        for (auto [j, var] : pairs[level][i]) {
          auto [it, inserted] = level_vars.try_emplace({i, j}, pairs.size(), -1);
          it->second[level] = var;
        }
    std::vector<std::vector<std::vector<int>>> stem_levels(builder.stems.size(),
        std::vector<std::vector<int>>(pairs.size()));
    for (size_t index = 0; index < builder.stems.size(); ++index) {
      const auto& stem = builder.stems[index];
      for (size_t level = 0; level < pairs.size(); ++level) {
        auto& vars = stem_levels[index][level];
        for (int d = 0; d < stem.length; ++d) {
          int var = level_vars.at({stem.left + d, stem.right - d})[level];
          if (var >= 0) vars.push_back(var);
          else if (options.supported && options.support_complete_stems) { vars.clear(); break; }
        }
      }
    }
    for (const auto& motif : motifs)
      for (auto [upper, lower] : {Pair{motif.stem1, motif.stem2},
                                 Pair{motif.stem2, motif.stem1}}) {
        // Sum over a*b realized crossings contributes one candidate motif
        // score when the two complete stems are its only scored geometry.
        const double score = motif.coefficient / builder.stems[upper].length /
            (options.crossing ? builder.stems[lower].length : 1);
        for (size_t level = 1; level < pairs.size(); ++level) {
          const auto& upper_vars = stem_levels[upper][level];
          if (upper_vars.empty()) continue;
          for (size_t below = 0; below < level; ++below) {
            const auto& lower_vars = stem_levels[lower][below];
            for (int v : upper_vars) for (int w : lower_vars) {
              auto& scores = builder.model.crossing_scores[v];
              auto [it, inserted] = scores.emplace(w, score);
              if (!inserted) it->second = std::max(it->second, score);
            }
          }
        }
      }
  } else if (options.projected) {
    // A cheap look-ahead surrogate: distribute each H-core score over either
    // potential upper stem and retain the best score per pair. It does NOT
    // require that the proposed crossing partner or gap state be selected.
    // The existing levelwise constraints still require an actual crossing.
    std::map<Pair, double> coefficients;
    for (const auto& motif : motifs) for (int index : {motif.stem1, motif.stem2}) {
      const auto& s = builder.stems[index];
      const double coefficient = motif.coefficient / s.length;
      for (int d = 0; d < s.length; ++d) {
        Pair pair{s.left + d, s.right - d};
        auto [it, inserted] = coefficients.emplace(pair, coefficient);
        if (!inserted) it->second = std::max(it->second, coefficient);
      }
    }
    for (size_t level = 1; level < pairs.size(); ++level)
      for (int i = 0; i < static_cast<int>(pairs[level].size()); ++i)
        for (auto [j, var] : pairs[level][i]) {
          auto found = coefficients.find({i, j});
          if (found == coefficients.end()) continue;
          ip.add_objective_coefficient(var, found->second);
          builder.model.terms.emplace_back(var, found->second);
        }
  } else {
    // The model is unchanged if the score has no compatible motif.
    for (const auto& motif : motifs) builder.add_motif(motif);
  }
  builder.model.motifs = motifs.size();
  builder.model.stems = builder.stem_vars.size();
  builder.model.intervals = builder.interval_vars.size();
  return std::move(builder.model);
}

void add_pk_crossing_constraints(IP& ip, const PKLevelPairs& pairs,
                                const PKScoreOptions& options,
                                PKScoreModel& model) {
  for (size_t level = 1; level < pairs.size(); ++level)
    for (int k = 0; k < static_cast<int>(pairs[level].size()); ++k)
      for (auto [l, var] : pairs[level][k]) {
        auto visit_witnesses = [&](size_t below, auto visit) {
          for (int i = 0; i < k; ++i)
            for (auto [j, other] : pairs[below][i])
              if (k < j && j < l) visit(other, i, j);
          for (int i = k + 1; i < static_cast<int>(l); ++i)
            for (auto [j, other] : pairs[below][i])
              if (l < j) visit(other, i, j);
        };
        auto found = model.crossing_scores.find(var);
        if (options.crossing) {
          // Preserve the original OR constraint. Score its selected partners
          // with one product auxiliary, rather than one AND per crossing edge.
          for (size_t below = 0; below < level; ++below) {
            int support = ip.make_constraint(IP::LO, 0, 0);
            ip.add_constraint(support, var, -1);
            ++model.crossing_rows;
            ++model.crossing_nonzeros;
            std::vector<std::pair<int, double>> weights;
            std::map<int, std::array<double, 2>> positive, negative;
            std::vector<PKCrossingBoundPair> positive_pairs, negative_pairs;
            double positive_sum = 0, negative_sum = 0;
            size_t witness_count = 0;
            std::array<int, 2> common_bases{-1, -1};
            visit_witnesses(below, [&](int other, int i, int j) {
              ip.add_constraint(support, other, 1);
              ++model.crossing_nonzeros;
              if (witness_count++ == 0) common_bases = {i, j};
              else for (int& base : common_bases)
                if (base != i && base != j) base = -1;
              if (found == model.crossing_scores.end()) return;
              auto edge = found->second.find(other);
              if (edge == found->second.end() || edge->second == 0) return;
              const double weight = edge->second;
              weights.emplace_back(other, weight);
              if (options.tight_crossing_bounds)
                (weight > 0 ? positive_pairs : negative_pairs).push_back({i, j, std::abs(weight)});
              auto& endpoint_max = weight > 0 ? positive : negative;
              endpoint_max[i][0] = std::max(endpoint_max[i][0], std::abs(weight));
              endpoint_max[j][1] = std::max(endpoint_max[j][1], std::abs(weight));
              if (weight > 0) positive_sum += weight;
              else negative_sum -= weight;
            });
            if (weights.empty()) continue;
            model.weighted_crossings += weights.size();
            // The original support row gives x <= sum(y). If all witnesses
            // share a base, base uniqueness gives sum(y) <= 1 at integral
            // assignments. With equal weights c (and no unscored witness),
            // x*sum(c*y) = c*x, including negative c. Singleton rows are a
            // special case. Keep exact equality, never a weight tolerance.
            if (options.simplify_crossing && weights.size() == witness_count &&
                (common_bases[0] >= 0 || common_bases[1] >= 0) &&
                std::all_of(weights.begin(), weights.end(), [&](const auto& edge) {
                  return edge.second == weights.front().second;
                })) {
              const double coefficient = weights.front().second;
              ip.add_objective_coefficient(var, coefficient);
              model.terms.emplace_back(var, coefficient);
              ++model.direct_crossing_rows;
              model.direct_crossing_edges += weights.size();
              continue;
            }
            // Each selected pair uses a unique left and a unique right base.
            // Either side's sum of maxima bounds the weighted sum, including
            // fractional feasible pair values. Keep the smaller bound.
            auto endpoint_bound = [](const auto& maxima) {
              std::array<double, 2> sum{0, 0};
              for (const auto& [position, value] : maxima)
                for (int side = 0; side < 2; ++side) sum[side] += value[side];
              return std::min(sum[0], sum[1]);
            };
            double upper = std::min(positive_sum, endpoint_bound(positive));
            double lower = -std::min(negative_sum, endpoint_bound(negative));
            if (options.tight_crossing_bounds) {
              const double before = upper - lower;
              const double old_upper = upper, old_lower = lower;
              upper = std::min(upper, pk_crossing_chain_bound(k, l, std::move(positive_pairs)));
              lower = -std::min(-lower, pk_crossing_chain_bound(k, l, std::move(negative_pairs)));
              model.tightened_bound_rows += upper < old_upper || lower > old_lower;
              model.crossing_bound_width_before += before;
              model.crossing_bound_width_after += upper - lower;
            }
            if (!std::isfinite(upper) || !std::isfinite(lower))
              throw std::invalid_argument("PK crossing-score bound overflow");
            double scale = 1;
            if (options.normalize_crossing) {
              scale = 0;
              for (const auto& edge : weights) scale = std::max(scale, std::abs(edge.second));
              lower /= scale;
              upper /= scale;
              for (auto& edge : weights) edge.second /= scale;
              if (!std::isfinite(lower) || !std::isfinite(upper))
                throw std::invalid_argument("PK normalized crossing-score bound overflow");
            }
            int score = ip.make_continuous_variable(scale, lower, upper);
            ++model.continuous_variables;
            model.terms.emplace_back(score, scale);
            auto add = [&](int row, int column, double coefficient) {
              if (coefficient == 0) return;
              ip.add_constraint(row, column, coefficient);
              ++model.nonzeros;
            };
            // z=x*r, r=sum(w*y), L<=r<=U. Full bounds enforce equality.
            // With a positive objective coefficient on z, maximization
            // needs only z<=U*x and z<=r-L*(1-x): their maximum is exactly
            // 0 at x=0 and r at x=1. Lower links can then be omitted.
            if (upper != 0) {
              int row = ip.make_constraint(IP::UP, 0, 0);
              ++model.rows;
              add(row, score, 1);
              add(row, var, -upper);
            }
            if (lower != 0 && !options.crossing_hypograph) {
              int row = ip.make_constraint(IP::LO, 0, 0);
              ++model.rows;
              add(row, score, 1);
              add(row, var, -lower);
            }
            int lo = -1;
            if (!options.crossing_hypograph) {
              lo = ip.make_constraint(IP::LO, -upper, -upper);
              ++model.rows;
              add(lo, score, 1);
              add(lo, var, -upper);
            }
            int hi = ip.make_constraint(IP::UP, -lower, -lower);
            ++model.rows;
            add(hi, score, 1);
            add(hi, var, -lower);
            for (auto [other, weight] : weights) {
              if (lo >= 0) add(lo, other, -weight);
              add(hi, other, -weight);
            }
          }
          continue;
        }
        const bool scored = options.supported && found != model.crossing_scores.end();
        if (!scored) {
          // Preserve the original streaming construction when this variable
          // has no shape score, including the disabled/default formulation.
          for (size_t below = 0; below < level; ++below) {
            int row = ip.make_constraint(IP::LO, 0, 0);
            ip.add_constraint(row, var, -1);
            ++model.crossing_rows;
            ++model.crossing_nonzeros;
            visit_witnesses(below, [&](int other, int, int) {
              ip.add_constraint(row, other, 1);
              ++model.crossing_nonzeros;
            });
          }
          continue;
        }
        std::vector<std::vector<int>> witnesses(level);
        for (size_t below = 0; below < level; ++below)
          visit_witnesses(below, [&](int other, int, int) { witnesses[below].push_back(other); });
        auto quality = [&](int other) {
          auto edge = found->second.find(other);
          return edge == found->second.end() ? 0.0 : edge->second;
        };
        double coefficient = 0;
        if (scored) {
          // A witness is required in every lower level. The weakest level's
          // best score bounds the bonus that all those rows can support.
          coefficient = std::numeric_limits<double>::infinity();
          for (const auto& row : witnesses) {
            if (row.empty()) { coefficient = 0; break; }
            double best = -std::numeric_limits<double>::infinity();
            for (int other : row) best = std::max(best, quality(other));
            coefficient = std::min(coefficient, best);
          }
        }
        if (coefficient != 0) {
          ip.add_objective_coefficient(var, coefficient);
          model.terms.emplace_back(var, coefficient);
        }
        for (const auto& witnesses_in_level : witnesses) {
          int row = ip.make_constraint(IP::LO, 0, 0);
          ip.add_constraint(row, var, -1);
          ++model.crossing_rows;
          ++model.crossing_nonzeros;
          bool tightened = false;
          for (int other : witnesses_in_level) {
            if (coefficient != 0 && quality(other) < coefficient - options.support_tolerance) {
              ++model.removed_crossings;
              tightened = true;
              continue;
            }
            ip.add_constraint(row, other, 1);
            ++model.crossing_nonzeros;
          }
          if (tightened) ++model.tightened_rows;
        }
      }
  model.crossing_scores.clear();
}

double score_pk_h_structure(const std::vector<int>& bpseq,
                            const PKScoreOptions& options) {
  if (!options.has_h_score()) return 0;
  options.validate();
  std::vector<Stem> stems;
  const int length = bpseq.size();
  for (int i = 0; i < length; ++i) {
    int j = bpseq[i];
    if (j <= i || (i > 0 && bpseq[i - 1] == j + 1)) continue;
    int run = 1;
    while (i + run < j - run && bpseq[i + run] == j - run) ++run;
    if (run >= 2 && run <= options.max_stem) stems.push_back({i, j, run});
  }
  double result = 0;
  for (const auto& a : stems) {
    int first = a.left + a.length;
    int last = first + std::min(options.max_loop, length - first);
    auto begin = std::lower_bound(stems.begin(), stems.end(), first,
        [](const Stem& stem, int position) { return stem.left < position; });
    for (auto it = begin; it != stems.end() && it->left <= last; ++it) {
      const auto& b = *it;
      if (b.right <= a.right) continue;
      std::array<int, 5> geometry{a.length, b.length,
          b.left - a.left - a.length,
          a.right - a.length - b.left - b.length + 1,
          b.right - b.length - a.right};
      bool valid = true;
      for (int k = 2; k < 5; ++k)
        valid &= geometry[k] >= 0 && geometry[k] <= options.max_loop;
      if (!valid) continue;
      if (options.unpaired_loops) {
        for (Pair interval : {Pair{a.left + a.length, b.left - 1},
                              Pair{b.left + b.length, a.right - a.length},
                              Pair{a.right + 1, b.right - b.length}})
          for (int i = interval.first; i <= interval.second; ++i)
            if (bpseq[i] >= 0) valid = false;
      }
      if (valid) result += options.score(geometry);
    }
  }
  if (!std::isfinite(result)) throw std::invalid_argument("PK structure score overflow");
  return result;
}
