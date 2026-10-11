#ifndef IPKNOT_DUAL_DECOMPOSITION_H
#define IPKNOT_DUAL_DECOMPOSITION_H

#include <string>
#include <utility>
#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

struct DDOptions {
  bool enabled = false;
  int max_iterations = 50;
  std::size_t constraint_states = 2048; // total visited repair states per solve; 0: unlimited
  int constraint_recovery_every = 0; // 0: final-only; finite budget shared by all calls
  bool linear_constraints = true;
  bool noe_ilp = false; // joint NOE factor; RNA pairs remain in the DP
  int constraint_passes = 8;
  int nmr_pair_beam = 64, nmr_witnesses = 64, nmr_pattern_beam = 32;
  int nmr_bulge_mode = 2; // none=0, fallback=1, all=2
  int beam = 100;
  int threads = 1; // independent level oracles; 0: OpenMP runtime default
  int dp_threads = 1; // beam/wavefront workers; 0: automatic; overrides level threads
  bool nussinov_dp = false; // exact Nussinov ignores the configured beam
  bool wavefront_dp = false; // exact interval wavefront, quadratic storage
  bool improved_beam = false; // suffix dominance before beam selection
  int dp_beam() const { return nussinov_dp || wavefront_dp ? 0 : beam; }
  int crossing_beam = 100;
  int witnesses = 16;
  int patience = 0;
  double step = 1.5;
  bool projected_norm = true;
  bool unpruned_bound = false;
  bool diminishing_step = false;
  int bound_block = 0;
  int bound_every = 10;
  bool bound_shift = false, bound_strict_stack = false;
  bool global_bound = false;
  int joint_bound_width = 0;
  bool joint_bound_clusters = false;
  bool joint_bound_shift = false, joint_bound_matching = false;
  std::size_t joint_bound_states = 20000; // zero: unlimited diagnostic search
  int exchange_width = 0, exchange_passes = 1, exchange_every = 0;
  std::size_t exchange_states = 20000; // zero: unlimited diagnostic search
  bool recovery_cache = true, recovery_share = true;
  int recovery_every = 0;
  bool recovery_target_best = false; // use extra proposals in the Polyak target
  unsigned recovery_mode = 7; // original=1, averaged=2, upper-aware=4

  // Streaming diagnostics; disabled by default. State dumps can be large.
  std::string trace_file;
  bool trace_state = false;
  void validate() const;
};

struct DDPair {
  int left, right, level;
  double weight;
};

// A sparse, left-to-right Nussinov chart with reusable buffers. beam=0
// disables beam pruning, for exact per-level decoding on the retained graph.
// Exact charts omit dominated starts and decomposable closed split candidates.
// wavefront=true selects an exact quadratic interval chart (requires beam=0).
// dp_threads controls wavefront/finite-beam teams; sparse exact stays serial.
class DDNussinov {
public:
  DDNussinov(int length, const std::vector<DDPair>& pairs,
             int level, int beam, bool no_lonely_pairs, bool improved_beam = false,
             int dp_threads = 1, bool wavefront = false);
  ~DDNussinov();
  DDNussinov(const DDNussinov&) = delete;
  DDNussinov& operator=(const DDNussinov&) = delete;
  std::size_t pruned_states() const;
  std::size_t sparsified_candidates() const;
  int dp_threads() const;
  int used_workers() const; // most recent decode (including serial fallback)
  std::size_t parallel_steps() const; // processed columns, or wavefront diagonals
  double decode(const std::vector<double>& weights,
                const std::vector<unsigned char>& allowed,
                std::vector<int>& selected);
private:
  struct Workspace;
  std::unique_ptr<Workspace> workspace_;
};

struct DDCrossingEvidence {
  int left1, right1, left2, right2;
  double product;
};
// Fixed-width physical crossing sample for automatic threshold selection.
std::vector<DDCrossingEvidence> dd_crossing_evidence(
    const std::vector<std::vector<std::pair<unsigned int, float>>>& posterior,
    int crossing_beam);

struct DDBoundedRow {
  int upper, lower_level;
  std::vector<int> contacts;
};
struct DDBoundedGraph {
  std::vector<DDBoundedRow> rows;
  std::size_t crossing_drops = 0, witness_drops = 0;
};
DDBoundedGraph dd_bounded_graph(int length, const std::vector<DDPair>& pairs,
    int levels, const DDOptions& options);

struct DDResult {
  std::vector<int> bpseq, levels;
  double objective = 0;
  // Bound on the bounded witness graph, NOT the unrestricted ILP model.
  double upper_bound = std::numeric_limits<double>::infinity();
  int iterations = 0;
  int bound_evaluations = 0, recovery_proposals = 0;
  std::size_t joint_windows = 0, joint_fallbacks = 0, joint_states = 0;
  std::size_t exchange_windows = 0, exchange_improvements = 0;
  std::size_t exchange_budget_windows = 0, exchange_states = 0;
  std::string stop_reason = "no_candidates";
  std::size_t pairs = 0, support_rows = 0, contacts = 0;
  std::size_t crossing_beam_drops = 0, witness_drops = 0;
};

DDResult solve_dual_decomposition(int length, const std::vector<DDPair>& pairs,
    int levels, bool no_lonely_pairs, const DDOptions& options);

#endif
