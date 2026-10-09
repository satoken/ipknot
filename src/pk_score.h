#ifndef IPKNOT_PK_SCORE_H
#define IPKNOT_PK_SCORE_H

#include "pk_energy.h"
#include "pk_learned.h"
#include "pk_ranker.h"

#include <array>
#include <algorithm>
#include <cstddef>
#include <map>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

class IP;
using PKLevelPairs = std::vector<std::vector<
    std::vector<std::pair<unsigned int, int>>>>;

// Experimental decoder scores, in IPknot objective units, not kcal/mol.
// All zero weights preserve the original formulation and candidate set.
struct PKScoreOptions {
  bool ensemble = false;
  int core_width = 0; // 0: maximal stems; 2/3: disjoint short cores
  bool best_partner = false; // DD crossing: max over actually selected partner groups
  PKRankModel ranker;
  double rank_scale = 0;
  std::string rank_output;
  double ensemble_scale = 0, ensemble_intercept = 0;
  double ensemble_temperature = 1, ensemble_threshold = .5;
  std::string ensemble_output;
  double level_penalty = 0.0; // cheap, decomposition-dependent comparator
  double intercept = 0.0;
  double loop_penalty = 0.0;
  double stem_reward = 0.0;
  double coax_bonus = 0.0; // potential flush junction, not a full energy model
  double h_weight = 1.0; // scale a fixed feature model or imported score table
  double candidate_threshold = -1.0; // negative means use the original cut
  double selection_weight = 1.0;
  double support_tolerance = 0.0; // allowed per-pair score gap for crossing witnesses
  PKLoopEnergy energy;
  double energy_scale = 0.0; // coefficient of dimensionless PK loop cost G/(RT)
  double energy_intercept = 0.0; // independently calibrated PK shape prior
  PKLearnedModel learned;
  double learned_scale = 1.0; // independent residual scale; does not use h_weight
  bool hybrid_shape = false; // add geometry/(A*B) to the learned per-pair correction
  bool simplify_crossing = false; // exact direct terms for constant exclusive witnesses
  bool normalize_crossing = false; // rescale auxiliaries, preserving the same objective
  bool tight_crossing_bounds = false; // noncrossing-chain bounds, no new rows/variables
  bool crossing_hypograph = false; // only upper product bounds; exact when maximizing
  std::string feature_output; // optional append-only selected-block TSV diagnostics
  int max_stem = 12;
  int max_loop = 30;
  int max_motifs = 10000; // fail rather than silently truncate the score
  bool unpaired_loops = true; // false scores gap spans, permitting nested loops
  bool projected = false; // best potential partner per upper pair; no activation
  bool supported = false; // project scores while restricting existing crossing rows
  bool crossing = false; // one continuous score auxiliary per crossing-support row
  bool fixed_blocks = false; // score each maximal posterior stem block pair once
  bool support_complete_stems = false; // require whole candidate stems within witness levels
  bool rerank = false; // score realized motifs only, leaving each ILP unchanged
  std::map<std::array<int, 5>, double> table;

  bool has_h_score() const;
  bool needs_posterior_context() const;
  bool enabled() const;
  double score(const std::array<int, 5>& geometry) const;
  void validate() const;
  void load_table(const std::string& filename);
};

// Disjoint cores avoid counting every overlapping window as another motif.
// Merge a one-pair tail into the preceding core; its size is at most width+1.
inline std::vector<std::pair<int, int>> pk_core_segments(int length, int width) {
  if (width != 0 && width != 2 && width != 3)
    throw std::invalid_argument("PK core width must be 0, 2 or 3");
  std::vector<std::pair<int, int>> result;
  for (int offset = 0; offset < length;) {
    int size = width ? std::min(width, length - offset) : length;
    if (width && length - offset - size == 1) ++size;
    result.push_back({offset, size}); offset += size;
  }
  return result;
}

struct PKBlockFeatureRecord {
  std::array<int, 6> stems; // first left/right/length, second left/right/length
  std::array<int, 5> geometry;
  PKLearnedFeatures features;
  double delta = 0;
  // All level-specific variables for each physical base pair. Metadata is
  // retained only when feature export is requested, rather than for scoring.
  std::vector<std::vector<int>> first_pair_variables;
  std::vector<std::vector<int>> second_pair_variables;
};

struct PKScoreModel {
  std::vector<std::pair<int, double>> terms;
  size_t motifs = 0;
  size_t stems = 0;
  size_t intervals = 0;
  size_t continuous_variables = 0;
  size_t rows = 0;
  size_t nonzeros = 0;
  size_t crossing_rows = 0;
  size_t crossing_nonzeros = 0;
  size_t tightened_rows = 0;
  size_t removed_crossings = 0;
  size_t weighted_crossings = 0;
  size_t direct_crossing_rows = 0;
  size_t direct_crossing_edges = 0;
  size_t tightened_bound_rows = 0;
  double crossing_bound_width_before = 0;
  double crossing_bound_width_after = 0;
  size_t candidate_blocks = 0;
  size_t out_of_domain_blocks = 0;
  size_t eligible_block_pairs = 0;
  size_t skipped_block_pairs = 0;
  size_t scored_block_pairs = 0;
  size_t covered_crossing_pairs = 0;
  size_t energy_dp_motifs = 0;
  size_t energy_cc06_motifs = 0;
  size_t energy_cc09_motifs = 0;
  size_t energy_fallback_motifs = 0;
  size_t missing_posterior_block_pairs = 0;
  size_t unsupported_block_pairs = 0;
  std::vector<PKBlockFeatureRecord> block_features;
  // Potential shape compatibility, not realized motif activation. The key
  // variables already encode their levels; unlisted crossing edges score 0.
  std::map<int, std::map<int, double>> crossing_scores;

  double value(const IP& ip) const;
  void write_features(const std::string& filename, const std::string& sequence,
                      const std::vector<float>& thresholds, const IP& ip) const;
};

// Scores pairs of maximal contiguous crossing stems. In unpaired mode these
// must occupy a complete H motif; span mode permits other pairs in the gaps.
// Shared stem/interval ANDs are continuous: integral pair assignments force
// every auxiliary to 0 or 1, for either sign of the objective coefficient.
PKScoreModel add_pk_h_score(IP& ip, const PKLevelPairs& pairs,
                           const PKScoreOptions& options,
                           const PKPosteriorContext* posterior = nullptr);

// Builds IPknot's existing requirement of a crossing witness in EVERY lower
// level. Supported narrows its witnesses. Crossing keeps the original rows
// and scores selected partners with <= 1 continuous auxiliary per row.
void add_pk_crossing_constraints(IP& ip, const PKLevelPairs& pairs,
                                const PKScoreOptions& options,
                                PKScoreModel& model);

// Layer-independent score on an already decoded structure. Used to choose
// among the candidates already generated by automatic threshold search.
double score_pk_h_structure(const std::vector<int>& bpseq,
                            const PKScoreOptions& options);

#endif
