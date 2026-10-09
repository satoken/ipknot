// Sparse beam Nussinov follows DAFS src/nussinov.cpp (Kengo Sato, GPLv3).
// The stacked-pair states, bounded crossing graph and DD factors are IPknot
// extensions. See docs/dual-decomposition.md for the model and complexity.
#include "dual_decomposition.h"
#include "dd_bounds.h"
#include "dd_recovery.h"
#include "dd_joint_bound.h"
#include "dd_exchange.h"
#include "pk_best_partner.h"

#include <algorithm>
#include <atomic>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace {
std::uint64_t key(int i, int j) {
  return (std::uint64_t(static_cast<unsigned>(i)) << 32) | static_cast<unsigned>(j);
}
bool crosses(const DDPair& a, const DDPair& b) {
  return (a.left < b.left && b.left < a.right && a.right < b.right) ||
         (b.left < a.left && a.left < b.right && b.right < a.right);
}
struct Trace { int pair = -1, first = -1, second = -1; };
struct Cell { int start = 0, trace = -1; double score = 0; };
}

void DDOptions::validate() const {
  if (constraint_passes < 1 || nmr_pair_beam < 1 || nmr_witnesses < 1 || nmr_pattern_beam < 1)
    throw std::invalid_argument("DD constraint passes and NMR beams must be positive");
  if (constraint_recovery_every < 0)
    throw std::invalid_argument("DD constraint recovery interval must be nonnegative");
  if (bound_block < 0 || bound_block > 64 || bound_every < 1 || recovery_every < 0 ||
      recovery_mode < 1 || recovery_mode > 7)
    throw std::invalid_argument("DD requires bound-block in 0..64, bound-every >= 1, recovery-every >= 0 and valid recovery mode");
  if (joint_bound_width < 0 || joint_bound_width > 12 || exchange_width < 0 ||
      exchange_width > 12 || exchange_passes < 1 || exchange_passes > 4 || exchange_every < 0)
    throw std::invalid_argument("DD requires joint/exchange widths in 0..12, exchange passes in 1..4 and a nonnegative interval");
  if (joint_bound_clusters && joint_bound_shift)
    throw std::invalid_argument("DD contact clusters and shifted sequence windows are mutually exclusive");
  if (trace_state && trace_file.empty())
    throw std::invalid_argument("DD state tracing requires --dd-trace FILE");
  if (max_iterations < 1 || beam < 0 || crossing_beam < 0 || witnesses < 0 ||
      patience < 0 || !std::isfinite(step) || step <= 0 || step >= 2)
    throw std::invalid_argument("DD requires max-iter >= 1, nonnegative beams/witnesses/patience, and 0 < step < 2");
}

struct DDNussinov::Workspace {
  int length, beam;
  std::size_t pruned = 0;
  bool no_lonely, improved;
  const std::vector<DDPair>& pairs;
  std::vector<std::vector<int>> by_right;
  std::vector<int> inner, paired_trace;
  std::vector<double> paired_score;
  std::vector<std::vector<Cell>> chart;
  std::vector<Cell> candidates;
  std::vector<Trace> candidate_traces;
  std::vector<int> origin, scratch;
  std::vector<int> generation, starts, suffix_best, pending;
  std::vector<Trace> traces;
  std::vector<double> prefix;
  Workspace(int n, const std::vector<DDPair>& p, int level, int b, bool lonely, bool improve)
      : length(n), beam(b), no_lonely(lonely), improved(improve), pairs(p), by_right(n),
        inner(p.size(), -1), paired_trace(p.size(), -1), paired_score(p.size()),
        chart(n), candidates(n), candidate_traces(n), origin(n), generation(n, -1), prefix(n) {
    std::unordered_map<std::uint64_t, int> lookup;
    for (int id = 0; id < static_cast<int>(p.size()); ++id) if (p[id].level == level) {
      by_right[p[id].right].push_back(id);
      lookup.emplace(key(p[id].left, p[id].right), id);
    }
    for (const auto& row : by_right) for (int id : row) {
      auto found = lookup.find(key(p[id].left + 1, p[id].right - 1));
      if (found != lookup.end()) inner[id] = found->second;
    }
  }
};

DDNussinov::DDNussinov(int n, const std::vector<DDPair>& pairs, int level,
                     int beam, bool no_lonely_pairs, bool improved_beam) {
  if (n < 0 || level < 0 || beam < 0) throw std::invalid_argument("Invalid Nussinov dimensions");
  for (const auto& p : pairs)
    if (p.left < 0 || p.right >= n || p.left >= p.right || p.level < 0 || !std::isfinite(p.weight))
      throw std::invalid_argument("Invalid Nussinov pair");
  workspace_ = std::make_unique<Workspace>(n, pairs, level, beam, no_lonely_pairs, improved_beam);
}
DDNussinov::~DDNussinov() = default;
std::size_t DDNussinov::pruned_states() const { return workspace_->pruned; }

double DDNussinov::decode(const std::vector<double>& weights,
                         const std::vector<unsigned char>& allowed,
                         std::vector<int>& selected) {
  auto& w = *workspace_;
  if (weights.size() != w.pairs.size() || allowed.size() != weights.size())
    throw std::invalid_argument("Nussinov weight/mask size mismatch");
  selected.clear();
  w.pruned = 0;
  if (!w.length) return 0;
  w.traces.clear();
  std::fill(w.generation.begin(), w.generation.end(), -1);
  std::fill(w.paired_score.begin(), w.paired_score.end(), -std::numeric_limits<double>::infinity());
  const auto trace = [&](int pair, int first, int second) {
    w.traces.push_back({pair, first, second});
    return static_cast<int>(w.traces.size() - 1);
  };
  for (int right = 0; right < w.length; ++right) {
    w.starts.clear();
    const auto offer = [&](int start, double score, int pair, int first, int second) {
      if (w.improved) {
        const bool fresh = w.generation[start] != right;
        const int derivation = pair < 0 ? -1 : pair;
        if (fresh || score > w.candidates[start].score ||
            (score == w.candidates[start].score && derivation < w.origin[start])) {
          w.origin[start] = derivation;
          if (fresh) { w.generation[start] = right; w.starts.push_back(start); }
          w.candidates[start] = {start, pair < 0 ? first : -1, score};
          w.candidate_traces[start] = {pair, first, second};
        }
        return;
      }
      if (w.generation[start] != right) {
        w.generation[start] = right;
        w.starts.push_back(start);
        w.candidates[start] = {start, first, score};
        w.candidate_traces[start] = {pair, first, second};
      } else if (score > w.candidates[start].score) {
        w.candidates[start] = {start, first, score};
        w.candidate_traces[start] = {pair, first, second};
      }
    };
    offer(0, 0, -1, -1, -1);
    offer(right, 0, -1, -1, -1);
    if (right) {
      const auto& previous = w.chart[right - 1];
      w.suffix_best.resize(previous.size());
      int best = -1;
      for (int k = static_cast<int>(previous.size()); k-- > 0;) {
        if (best < 0 || previous[k].score > previous[best].score) best = k;
        w.suffix_best[k] = best;
      }
      for (const auto& cell : previous) offer(cell.start, cell.score, -1, cell.trace, -1);
    }
    for (int id : w.by_right[right]) {
      if (!allowed[id]) continue;
      const auto& p = w.pairs[id];
      double inside = 0;
      int inside_trace = -1;
      if (right) {
        const auto& previous = w.chart[right - 1];
        auto it = std::lower_bound(previous.begin(), previous.end(), p.left + 1,
            [](const Cell& cell, int start) { return cell.start < start; });
        if (it != previous.end()) {
          const auto& cell = previous[w.suffix_best[it - previous.begin()]];
          inside = cell.score; inside_trace = cell.trace;
        }
      }
      // P(i,j) may borrow stacking support from its parent. F intervals
      // contain only complete helices; a pair exposed to F must have a child.
      const int child = w.inner[id];
      if (w.no_lonely && child >= 0 && w.paired_score[child] > inside) {
        inside = w.paired_score[child]; inside_trace = w.paired_trace[child];
      }
      w.paired_score[id] = weights[id] + inside;
      w.paired_trace[id] = trace(id, inside_trace, -1);
      double complete = w.paired_score[id];
      int complete_inside = inside_trace;
      if (w.no_lonely) {
        if (child < 0 || !std::isfinite(w.paired_score[child])) continue;
        complete = weights[id] + w.paired_score[child];
        complete_inside = w.paired_trace[child];
      }
      // Negative pairs remain available inside a profitable stack. Only a
      // complete negative component can safely be discarded.
      if (!(complete > 0)) continue;
      offer(p.left, complete, id, complete_inside, -1);
      if (p.left) for (const auto& prefix : w.chart[p.left - 1])
        offer(prefix.start, prefix.score + complete, id, complete_inside, prefix.trace);
    }
    if (w.improved) {
      // Match the fixed-beam dominance/lazy-trace experiment, with heap and
      // adaptive widths disabled. Sorting large sparse columns by four fixed
      // radix passes keeps work linear before applying the fixed beam.
      if (w.starts.size() < 1024) std::sort(w.starts.begin(), w.starts.end());
      else {
        w.scratch.resize(w.starts.size());
        for (unsigned shift = 0; shift < 32; shift += 8) {
          std::array<std::size_t, 256> offsets{};
          for (int start : w.starts) ++offsets[(static_cast<unsigned>(start) >> shift) & 255];
          std::size_t offset = 0;
          for (auto& count : offsets) { const auto next = count; count = offset; offset += next; }
          for (int start : w.starts) w.scratch[offsets[(static_cast<unsigned>(start) >> shift) & 255]++] = start;
          w.starts.swap(w.scratch);
        }
      }
      // An earlier non-root interval with no greater score merely adds unused
      // leading bases. Its later suffix answers every future query at least
      // as well, including under signed pair weights and stacking support.
      double suffix_score = -std::numeric_limits<double>::infinity();
      w.scratch.clear();
      for (std::size_t k = w.starts.size(); k-- > 0;) {
        const int start = w.starts[k];
        const double score = w.candidates[start].score;
        if (start == 0 || score > suffix_score) {
          w.scratch.push_back(start);
          suffix_score = std::max(suffix_score, score);
        }
      }
      w.starts.swap(w.scratch);
    }
    auto better = [&](int a, int b) {
      const double sa = w.candidates[a].score + (a ? w.prefix[a - 1] : 0);
      const double sb = w.candidates[b].score + (b ? w.prefix[b - 1] : 0);
      return sa != sb ? sa > sb : a < b;
    };
    if (w.beam && w.starts.size() > static_cast<std::size_t>(w.beam)) {
      w.pruned += w.starts.size() - w.beam;
      std::nth_element(w.starts.begin(), w.starts.begin() + w.beam, w.starts.end(), better);
      w.starts.resize(w.beam);
      if (std::find(w.starts.begin(), w.starts.end(), 0) == w.starts.end()) {
        auto worst = std::max_element(w.starts.begin(), w.starts.end(), better);
        *worst = 0;
      }
    }
    std::sort(w.starts.begin(), w.starts.end()); // interval-start order for suffix queries
    auto& column = w.chart[right];
    column.clear();
    // Materialize only the winning traceback for each retained cell.
    // Eager traces for every improving split can use cubic memory in exact
    // mode; retained chart cells plus one paired trace per candidate need
    // O(n*min(beam,n)+M) records (O(n^2+M) when beam=0).
    for (int start : w.starts) {
      auto cell = w.candidates[start];
      const auto& proposal = w.candidate_traces[start];
      if (proposal.pair >= 0) cell.trace = trace(proposal.pair, proposal.first, proposal.second);
      column.push_back(cell);
    }
    w.prefix[right] = w.candidates[0].score;
  }
  w.pending.clear();
  w.pending.push_back(w.chart.back().front().trace);
  while (!w.pending.empty()) {
    int id = w.pending.back(); w.pending.pop_back();
    if (id < 0) continue;
    const auto& t = w.traces[id];
    if (t.pair >= 0) selected.push_back(t.pair);
    w.pending.push_back(t.first); w.pending.push_back(t.second);
  }
  return w.prefix.back();
}

namespace {
struct Contact { int lower; double score = 0, first_q = 0, second_q = 0; int first_z = 0, second_z = 0; int group = -1; };
struct SupportRow { int upper, lower_level; double multiplier = 0; std::vector<Contact> contacts;
  bool best_partner = false; double upper_q = 0; int upper_z = 0, groups = 0;
  std::vector<PKPartnerTerm> partner_terms;
  PKPartnerWorkspace partner_workspace;
};

double row_score(const SupportRow& row, const std::vector<unsigned char>& selected) {
  if (!row.best_partner) {
    double result=0; for(const auto& c:row.contacts) if(selected[c.lower]) result+=c.score;
    return result;
  }
  std::vector<double> sums(row.groups); std::vector<int> counts(row.groups);
  for(const auto& c:row.contacts) if(selected[c.lower]) { sums[c.group]+=c.score; ++counts[c.group]; }
  double result=-std::numeric_limits<double>::infinity();
  for(int g=0;g<row.groups;++g) if(counts[g]) result=std::max(result,sums[g]);
  return result;
}
struct Block { int left, right, length; };

class ContactScores {
  const PKScoreOptions& pk;
  const PKPosteriorContext* posterior;
  std::vector<Block> blocks;
  std::vector<int> pair_block;
  std::unordered_map<std::uint64_t, double> cache;
  std::size_t eligible = 0;
  const std::vector<DDPair>& pairs;
  struct Maximum { double value = 0; bool seen = false; };
  std::vector<std::vector<int>> level_counts;
  std::vector<std::vector<Maximum>> projected_best;
  std::vector<unsigned char> has_posterior;
  std::size_t projected_blocks = 0;
public:
  ContactScores(const std::vector<DDPair>& p, int levels, const PKScoreOptions& options,
                const PKPosteriorContext* context) : pk(options), posterior(context), pair_block(p.size(), -1), pairs(p) {
    if (!pk.has_h_score()) return;
    const bool learned = pk.needs_posterior_context();
    std::unordered_map<std::uint64_t, std::vector<int>> physical;
    for (int id = 0; id < static_cast<int>(p.size()); ++id)
      physical[key(p[id].left, p[id].right)].push_back(id);
    // Visit physical coordinates in deterministic input order, never hash order.
    for (int id = 0; id < static_cast<int>(p.size()); ++id) {
      const auto& pair = p[id];
      if (pair_block[id] >= 0 || physical.count(key(pair.left - 1, pair.right + 1))) continue;
      int run = 0;
      while (pair.left + run < pair.right - run && physical.count(key(pair.left + run, pair.right - run))) ++run;
      for (const auto [offset, size] : pk_core_segments(run, pk.core_width)) {
      const int block = blocks.size();
      blocks.push_back({pair.left + offset, pair.right - offset, size});
      bool complete = true;
      if (pk.projected) {
        level_counts.emplace_back(levels, 0);
        projected_best.emplace_back(levels - 1);
      }
      for (int d = offset; d < offset + size; ++d) {
        if (learned)
          complete &= posterior && posterior->contains(pair.left + d, pair.right - d);
        for (int var : physical.at(key(pair.left + d, pair.right - d))) {
          pair_block[var] = block;
          if (pk.projected) ++level_counts[block][p[var].level];
        }
      }
      has_posterior.push_back(complete);
      }
    }
  }
  int block_id(int id) const { return pair_block[id]; }
  double score(int a_id, int b_id) {
    if (!pk.has_h_score()) return 0;
    int a = pair_block[a_id], b = pair_block[b_id];
    if (blocks[a].left > blocks[b].left) std::swap(a, b);
    const auto token = key(a, b);
    auto cached = cache.find(token);
    if (cached != cache.end()) return cached->second;
    const auto& s = blocks[a]; const auto& t = blocks[b];
    const std::array<int, 5> g{s.length, t.length, t.left - s.left - s.length,
        s.right - s.length - t.left - t.length + 1, t.right - t.length - s.right};
    double result = 0;
    if (s.length >= 2 && t.length >= 2 && s.length <= pk.max_stem && t.length <= pk.max_stem &&
        g[2] >= 0 && g[3] >= 0 && g[4] >= 0 && g[2] <= pk.max_loop && g[3] <= pk.max_loop && g[4] <= pk.max_loop) {
      if (eligible++ >= static_cast<std::size_t>(pk.max_motifs))
        throw std::runtime_error("DD PK block-pair budget exceeded; raise --pk-h-max-motifs");
      if (pk.needs_posterior_context()) {
        if (!posterior) throw std::invalid_argument("DD learned PK scoring needs posterior context");
        // Forced pairs may be missing from BPP. Match native block scoring:
        // abstain from learned confidence while retaining the geometric prior.
        if (has_posterior[a] && has_posterior[b]) {
          auto features = pk.learned.features(*posterior, {s.left, s.right, s.length}, {t.left, t.right, t.length}, g);
          result = pk.learned_scale * pk.learned.score(features);
          if (pk.learned.block_normalized()) result = result / s.length / t.length;
          else result /= features.anchor_first ? s.length : t.length;
        }
        if (pk.hybrid_shape) result += pk.score(g) / s.length / t.length;
      } else result = pk.score(g) / s.length / t.length;
    }
    if (!std::isfinite(result)) throw std::invalid_argument("DD PK contact score overflow");
    cache.emplace(token, result);
    if (pk.projected && result != 0) {
      ++projected_blocks;
      for (const auto orientation : {std::pair<int,int>{a,b}, {b,a}}) {
        auto& maxima = projected_best[orientation.first];
        const auto& counts = level_counts[orientation.second];
        for (std::size_t lower = 0; lower < maxima.size(); ++lower) if (counts[lower]) {
          // Repeated addition matches the ILP full-partner coefficient,
          // including roundoff, without materializing A*B contact products.
          double value = 0;
          for (int n = 0; n < counts[lower]; ++n) value += result;
          if (!std::isfinite(value)) throw std::invalid_argument("DD PK projected score overflow");
          auto& best = maxima[lower];
          if (!best.seen || value > best.value) { best.value = value; best.seen = true; }
        }
      }
    }
    return result;
  }
  std::vector<double> projection() const {
    if (!pk.projected || !pk.has_h_score()) return {};
    std::vector<double> result(pairs.size());
    for (std::size_t id = 0; id < pairs.size(); ++id) {
      for (int lower = 0; lower < pairs[id].level; ++lower) {
        const auto& best = projected_best[pair_block[id]][lower];
        if (best.seen) result[id] += best.value;
      }
      if (!std::isfinite(result[id])) throw std::invalid_argument("DD PK projected score overflow");
    }
    return result;
  }
  std::size_t projection_blocks() const { return projected_blocks; }
};
}

DDBoundedGraph dd_bounded_graph(int length, const std::vector<DDPair>& pairs,
    int levels, const DDOptions& options, const PKScoreOptions& pk,
    const PKPosteriorContext* posterior) {
  DDBoundedGraph graph;
  ContactScores scores(pairs, levels, pk, posterior);
  std::vector<std::vector<int>> by_left(length), row_ids(pairs.size()), active(levels);
  for (int id = 0; id < static_cast<int>(pairs.size()); ++id) {
    by_left[pairs[id].left].push_back(id);
    for (int lower = 0; lower < pairs[id].level; ++lower) {
      row_ids[id].push_back(graph.rows.size()); graph.rows.push_back({id, lower, {}});
    }
  }
  const auto offer = [&](int upper, int lower) {
    auto& contacts = graph.rows[row_ids[upper][pairs[lower].level]].contacts;
    const double score = scores.score(upper, lower);
    contacts.push_back({lower, pk.projected ? 0 : score});
    if (options.witnesses && contacts.size() > static_cast<std::size_t>(options.witnesses)) {
      auto better = [&](const auto& a, const auto& b) {
        const double sa = pairs[a.first].weight + a.second, sb = pairs[b.first].weight + b.second;
        return sa != sb ? sa > sb : a.first < b.first;
      };
      contacts.erase(std::max_element(contacts.begin(), contacts.end(), better)); ++graph.witness_drops;
    }
  };
  for (int left = 0; left < length; ++left) {
    for (auto& beam : active) beam.erase(std::remove_if(beam.begin(), beam.end(),
        [&](int id) { return pairs[id].right <= left; }), beam.end());
    for (int id : by_left[left]) for (int level = 0; level < levels; ++level)
      if (level != pairs[id].level) for (int other : active[level]) if (crosses(pairs[id], pairs[other])) {
        if (pairs[id].level > level) offer(id, other); else offer(other, id);
      }
    for (int id : by_left[left]) {
      auto& beam = active[pairs[id].level]; beam.push_back(id);
      if (options.crossing_beam && beam.size() > static_cast<std::size_t>(options.crossing_beam)) {
        auto better = [&](int a, int b) {
          return pairs[a].weight != pairs[b].weight ? pairs[a].weight > pairs[b].weight : a < b;
        };
        beam.erase(std::max_element(beam.begin(), beam.end(), better)); ++graph.crossing_drops;
      }
    }
  }
  graph.projected_coefficients = scores.projection();
  graph.projected_blocks = scores.projection_blocks();
  return graph;
}

DDResult solve_dual_decomposition(int length, const std::vector<DDPair>& input_pairs,
    int levels, bool no_lonely, const DDOptions& options,
    const PKScoreOptions& pk, const PKPosteriorContext* posterior) {
  options.validate(); pk.validate();
  if (length < 0 || levels < 1) throw std::invalid_argument("Invalid DD dimensions");
  if (pk.has_h_score() && ((!pk.crossing && !pk.projected) || !pk.fixed_blocks))
    throw std::invalid_argument("DD PK scoring requires crossing or projected with --pk-h-allocation blocks");
  if (!pk.feature_output.empty() || pk.rerank || pk.supported)
    throw std::invalid_argument("DD does not support PK feature export, rerank or supported scoring; use --decoder ilp");
  if (pk.best_partner && (options.global_bound || options.joint_bound_width ||
      options.exchange_width || options.recovery_every))
    throw std::invalid_argument("DD best-partner scoring cannot use additive-contact bound/recovery helpers");
  std::vector<DDPair> projected_pairs;
  if (pk.projected && pk.has_h_score()) projected_pairs = input_pairs;
  const auto& pairs = projected_pairs.empty() ? input_pairs : projected_pairs;
  DDResult result;
  result.bpseq.assign(length, -1); result.levels.assign(length, -1); result.pairs = pairs.size();
  const int count = pairs.size();
  std::vector<std::vector<int>> by_left(length), row_ids(count), by_level(levels);
  std::vector<std::unordered_map<std::uint64_t, int>> unique(levels);
  for (int id = 0; id < count; ++id) {
    const auto& p = pairs[id];
    if (p.left < 0 || p.right >= length || p.left >= p.right || p.level < 0 || p.level >= levels || !std::isfinite(p.weight))
      throw std::invalid_argument("Invalid DD pair");
    // Separate physical pairs may occur in several levels, never twice in one.
    const auto token = key(p.left, p.right);
    if (!unique[p.level].emplace(token, id).second)
      throw std::invalid_argument("Duplicate DD pair in one level");
    by_left[p.left].push_back(id); by_level[p.level].push_back(id);
  }
  std::ofstream trace;
  unsigned long long solve_id = 0;
  if (!options.trace_file.empty()) {
    trace.open(options.trace_file, std::ios::app);
    if (!trace) throw std::runtime_error("Cannot open DD trace: " + options.trace_file);
    trace.exceptions(std::ios::badbit | std::ios::failbit);
    trace << std::setprecision(17);
    static std::atomic<unsigned long long> next_id{0};
    solve_id = ++next_id;
  }
  const auto emit_summary = [&]() {
    if (!trace.is_open()) return;
    trace << "{\"event\":\"summary\",\"solve_id\":" << solve_id
          << ",\"bound_evaluations\":" << result.bound_evaluations
          << ",\"recovery_proposals\":" << result.recovery_proposals
          << ",\"exchange_windows\":" << result.exchange_windows
          << ",\"exchange_improvements\":" << result.exchange_improvements
          << ",\"exchange_budget_windows\":" << result.exchange_budget_windows
          << ",\"exchange_states\":" << result.exchange_states
          << ",\"iterations\":" << result.iterations << ",\"lower_bound\":" << result.objective
          << ",\"upper_bound\":" << result.upper_bound << ",\"stop\":\"" << result.stop_reason << "\"}\n";
  };
  if (!count) { result.upper_bound = 0; emit_summary(); return result; }
  ContactScores scores(pairs, levels, pk, posterior);
  std::vector<SupportRow> rows;
  for (int id = 0; id < count; ++id) for (int lower = 0; lower < pairs[id].level; ++lower) {
    row_ids[id].push_back(rows.size()); rows.push_back({id, lower, 0, {}});
  }
  std::vector<std::vector<int>> active(levels);
  const auto offer_contact = [&](int upper, int lower) {
    auto& row = rows[row_ids[upper][pairs[lower].level]];
    const double potential = scores.score(upper, lower);
    const double score = pk.projected ? 0 : potential;
    row.contacts.push_back({lower, score});
    if (options.witnesses && row.contacts.size() > static_cast<std::size_t>(options.witnesses)) {
      // Prefer posterior/objective evidence, then PK compatibility, then ID.
      auto better = [&](const Contact& a, const Contact& b) {
        const double sa = pairs[a.lower].weight + a.score, sb = pairs[b.lower].weight + b.score;
        return sa != sb ? sa > sb : a.lower < b.lower;
      };
      auto worst = std::max_element(row.contacts.begin(), row.contacts.end(), better);
      row.contacts.erase(worst); ++result.witness_drops;
    }
  };
  for (int left = 0; left < length; ++left) {
    for (auto& beam : active) beam.erase(std::remove_if(beam.begin(), beam.end(),
        [&](int id) { return pairs[id].right <= left; }), beam.end());
    // Query before inserting ANY pair at this left endpoint. Shared endpoints
    // are never crossings, and insertion order cannot create extra witnesses.
    for (int id : by_left[left]) for (int level = 0; level < levels; ++level) if (level != pairs[id].level)
      for (int other : active[level]) if (crosses(pairs[id], pairs[other])) {
        if (pairs[id].level > level) offer_contact(id, other);
        else offer_contact(other, id);
      }
    for (int id : by_left[left]) {
      auto& beam = active[pairs[id].level]; beam.push_back(id);
      if (options.crossing_beam && beam.size() > static_cast<std::size_t>(options.crossing_beam)) {
        auto better = [&](int a, int b) {
          return pairs[a].weight != pairs[b].weight ? pairs[a].weight > pairs[b].weight : a < b;
        };
        beam.erase(std::max_element(beam.begin(), beam.end(), better)); ++result.crossing_beam_drops;
      }
    }
  }
  // Build support with the original weights. Projection changes only unary
  // objectives, including for 3+ levels, not witness ranking or feasibility.
  const auto projected_coefficients = scores.projection();
  result.projected_blocks = scores.projection_blocks();
  for (std::size_t id = 0; id < projected_coefficients.size(); ++id) {
    projected_pairs[id].weight += projected_coefficients[id];
    if (!std::isfinite(projected_pairs[id].weight))
      throw std::invalid_argument("DD PK projected pair weight overflow");
    result.projected_pairs += projected_coefficients[id] != 0;
  }
  // Propagate impossible witnesses and stacking support in O(M+E). This
  // removes factors known to be zero before the first subgradient iteration.
  std::vector<unsigned char> allowed(count, 1), selected(count), recovered(count);
  std::vector<std::array<int, 2>> neighbors(count, {-1, -1});
  std::vector<int> stack_degree(count), support_count(rows.size()), pending;
  std::vector<std::vector<int>> incoming(count);
  if (no_lonely) for (int id = 0; id < count; ++id) {
    const auto& p = pairs[id];
    int side = 0;
    for (const auto token : {key(p.left + 1, p.right - 1), key(p.left - 1, p.right + 1)}) {
      auto it = unique[p.level].find(token);
      if (it != unique[p.level].end()) { neighbors[id][side] = it->second; ++stack_degree[id]; }
      ++side;
    }
    if (!stack_degree[id]) pending.push_back(id);
  }
  for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
    support_count[r] = rows[r].contacts.size();
    if (!support_count[r]) pending.push_back(rows[r].upper);
    for (const auto& c : rows[r].contacts) incoming[c.lower].push_back(r);
  }
  while (!pending.empty()) {
    int id = pending.back(); pending.pop_back();
    if (!allowed[id]) continue;
    allowed[id] = 0;
    if (no_lonely) for (int other : neighbors[id])
      if (other >= 0 && --stack_degree[other] == 0) pending.push_back(other);
    for (int r : incoming[id]) if (--support_count[r] == 0) pending.push_back(rows[r].upper);
  }
  result.support_rows = rows.size();
  for (auto& row : rows) {
    if (!allowed[row.upper]) row.contacts.clear();
    else row.contacts.erase(std::remove_if(row.contacts.begin(), row.contacts.end(),
        [&](const Contact& c) { return !allowed[c.lower]; }), row.contacts.end());
    result.contacts += row.contacts.size();
    for (const auto& c : row.contacts) result.scored_contacts += c.score != 0;
    row.best_partner = pk.best_partner && std::any_of(row.contacts.begin(), row.contacts.end(),
        [](const Contact& c) { return c.score != 0; });
    if (row.best_partner) {
      std::unordered_map<int,int> groups;
      for(auto& c:row.contacts) {
        const int block=scores.block_id(c.lower);
        auto found=groups.find(block);
        if(found==groups.end()) found=groups.emplace(block,groups.size()).first;
        c.group=found->second;
        row.partner_terms.push_back({c.group,c.score,0});
      }
      row.groups=groups.size();
    }
  }
  double static_bound = std::numeric_limits<double>::infinity();
  if (options.global_bound) {
    std::vector<DDScoredContact> contacts;
    for (const auto& row : rows) for (const auto& c : row.contacts)
      if (c.score != 0) contacts.push_back({row.upper,c.lower,c.score});
    static_bound = dd_global_bound(length,pairs,contacts,allowed);
  }
  std::vector<DDRecoveryRow> model_rows;
  if (options.recovery_every || options.joint_bound_width || options.exchange_width)
    for (const auto& row : rows) {
      DDRecoveryRow copied; copied.upper=row.upper;
      for (const auto& c : row.contacts) copied.contacts.push_back({c.lower,c.score});
      model_rows.push_back(std::move(copied));
    }
  if (options.joint_bound_width) {
    for (int offset : {0, options.joint_bound_width / 2}) {
      if (offset && (!options.joint_bound_shift || options.joint_bound_clusters)) break;
      // Width one has only one distinct partition.
      if (!offset && result.joint_windows) break;
      DDJointBound bound(length,pairs,levels,no_lonely,allowed,model_rows,
          options.joint_bound_width,offset,options.joint_bound_matching,options.joint_bound_states,options.joint_bound_clusters);
      const auto evaluated=bound.evaluate();
      static_bound=std::min(static_bound,evaluated.upper_bound);
      result.joint_windows+=evaluated.exact_windows+evaluated.fallback_windows;
      result.joint_fallbacks+=evaluated.fallback_windows;
      result.joint_states+=evaluated.states;
    }
  }
  std::vector<std::unique_ptr<DDBlockBound>> block_bounds, shifted_bounds;
  if (options.bound_block && options.dp_beam())
    for (int level=0;level<levels;++level)
      block_bounds.push_back(std::make_unique<DDBlockBound>(length,pairs,level,options.bound_block,0,no_lonely,options.bound_strict_stack));
  if (!block_bounds.empty() && options.bound_shift && options.bound_block > 1)
    for (int level=0;level<levels;++level)
      shifted_bounds.push_back(std::make_unique<DDBlockBound>(length,pairs,level,options.bound_block,options.bound_block/2,no_lonely,options.bound_strict_stack));
  // Recovery can reuse these buffers after the dual state and gradients have
  // been saved. Declare the helpers after the decoders to preserve lifetime.
  std::vector<std::unique_ptr<DDNussinov>> decoders;
  std::vector<DDNussinov*> shared_decoders;
  for (int level=0;level<levels;++level) {
    decoders.push_back(std::make_unique<DDNussinov>(length,pairs,level,options.dp_beam(),no_lonely,options.improved_beam));
    if (options.recovery_share) shared_decoders.push_back(decoders.back().get());
  }
  std::unique_ptr<DDPrimalRecovery> recovery;
  if (options.recovery_every)
    recovery=std::make_unique<DDPrimalRecovery>(length,pairs,levels,options.dp_beam(),no_lonely,
        allowed,model_rows,8,options.recovery_mode,options.recovery_cache,shared_decoders,options.improved_beam);
  std::unique_ptr<DDLocalExchange> exchange;
  if (options.exchange_width)
    exchange=std::make_unique<DDLocalExchange>(length,pairs,levels,no_lonely,allowed,model_rows);
  auto evaluate_block = [&](int level, const std::vector<double>& coefficients, int& evaluated) {
    double value=block_bounds[level]->evaluate(coefficients,allowed);
    ++evaluated; ++result.bound_evaluations;
    if (!shifted_bounds.empty()) {
      value=std::min(value,shifted_bounds[level]->evaluate(coefficients,allowed));
      ++evaluated; ++result.bound_evaluations;
    }
    return value;
  };
  if (trace.is_open()) {
    trace << "{\"event\":\"problem\",\"solve_id\":" << solve_id
          << ",\"length\":" << length << ",\"levels\":" << levels
          << ",\"no_lonely\":" << int(no_lonely) << ",\"beam\":" << options.dp_beam()
          << ",\"improved_beam\":" << int(options.improved_beam)
          << ",\"crossing_beam\":" << options.crossing_beam << ",\"witnesses\":" << options.witnesses
          << ",\"max_iterations\":" << options.max_iterations << ",\"patience\":" << options.patience
          << ",\"best_partner\":" << int(pk.best_partner) << ",\"core_width\":" << pk.core_width
          << ",\"bound_block\":" << options.bound_block << ",\"bound_every\":" << options.bound_every
          << ",\"global_bound\":" << int(options.global_bound) << ",\"recovery_every\":" << options.recovery_every
          << ",\"recovery_target_best\":" << int(options.recovery_target_best)
          << ",\"recovery_mode\":" << options.recovery_mode
          << ",\"bound_shift\":" << int(options.bound_shift)
          << ",\"bound_strict_stack\":" << int(options.bound_strict_stack)
          << ",\"joint_bound_width\":" << options.joint_bound_width
          << ",\"joint_bound_clusters\":" << int(options.joint_bound_clusters)
          << ",\"joint_bound_shift\":" << int(options.joint_bound_shift)
          << ",\"joint_bound_matching\":" << int(options.joint_bound_matching)
          << ",\"joint_bound_states\":" << options.joint_bound_states
          << ",\"joint_windows\":" << result.joint_windows
          << ",\"joint_fallbacks\":" << result.joint_fallbacks
          << ",\"joint_states\":" << result.joint_states
          << ",\"exchange_width\":" << options.exchange_width
          << ",\"exchange_passes\":" << options.exchange_passes
          << ",\"exchange_every\":" << options.exchange_every
          << ",\"exchange_state_budget\":" << options.exchange_states
          << ",\"recovery_cache\":" << int(options.recovery_cache)
          << ",\"recovery_share\":" << int(options.recovery_share)
          << ",\"unpruned_bound\":" << int(options.unpruned_bound)
          << ",\"diminishing_step\":" << int(options.diminishing_step)
          << ",\"step\":" << options.step << ",\"projected_norm\":" << int(options.projected_norm)
          << ",\"crossing_beam_drops\":" << result.crossing_beam_drops
          << ",\"witness_drops\":" << result.witness_drops;
    if (options.global_bound || options.joint_bound_width) trace << ",\"static_certificate\":" << static_bound;
    if (options.trace_state) {
      trace << ",\"projected_coefficients\":[";
      for (std::size_t id = 0; id < projected_coefficients.size(); ++id) {
        if (id) trace << ',';
        trace << projected_coefficients[id];
      }
      trace << "],\"partner_groups\":[";
      for(std::size_t r=0;r<rows.size();++r) {
        if(r) trace<<','; trace<<'[';
        for(std::size_t k=0;k<rows[r].contacts.size();++k) {if(k)trace<<',';trace<<rows[r].contacts[k].group;}
        trace<<']';
      }
      trace<<"],\"best_partner_rows\":[";
      for(std::size_t r=0;r<rows.size();++r) {if(r)trace<<',';trace<<int(rows[r].best_partner);}
      trace<<']';
      trace << ",\"pairs\":[";
      for (int id = 0; id < count; ++id) {
        const auto& p = pairs[id];
        if (id) trace << ',';
        trace << '[' << p.left << ',' << p.right << ',' << p.level << ',' << p.weight << ',' << int(allowed[id]) << ']';
      }
      trace << "],\"rows\":[";
      for (std::size_t r = 0; r < rows.size(); ++r) {
        if (r) trace << ',';
        trace << '[' << rows[r].upper << ",[";
        for (std::size_t k = 0; k < rows[r].contacts.size(); ++k) {
          if (k) trace << ',';
          const auto& c = rows[r].contacts[k];
          trace << '[' << c.lower << ',' << c.score << ']';
        }
        trace << "]]";
      }
      trace << ']';
    }
    trace << "}\n";
  }
  std::vector<double> lambda(length), weights(count), degree(length), gradient(length), maxima(length);
  std::vector<int> decoded, used(length), best;
  std::vector<double> level_bounds(levels);
  std::vector<unsigned char> exact_levels(levels);
  double best_bound = static_bound, best_beam = std::numeric_limits<double>::infinity();
  double baseline_best = 0;
  double scale = 0;
  for (const auto& p : pairs) scale = std::max(scale, std::abs(p.weight));
  for (const auto& row : rows) for (const auto& c : row.contacts) scale = std::max(scale, std::abs(c.score));
  scale = std::max(scale, 1e-3);
  result.stop_reason = "max_iterations";
  int stagnant = 0;
  double relaxation = options.step;
  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    std::fill(selected.begin(), selected.end(), 0);
    std::fill(degree.begin(), degree.end(), 0);
    for (int id = 0; id < count; ++id) weights[id] = pairs[id].weight - lambda[pairs[id].left] - lambda[pairs[id].right];
    double factor_value = 0;
    for (auto& row : rows) {
      if(row.best_partner) {
        weights[row.upper]-=row.upper_q;
        for(std::size_t k=0;k<row.contacts.size();++k) {
          auto& c=row.contacts[k]; weights[c.lower]-=c.second_q;
          row.partner_terms[k].multiplier=c.second_q;
        }
        const auto& state=pk_best_partner_factor(row.upper_q,row.partner_terms,row.groups,row.partner_workspace);
        factor_value+=state.value; row.upper_z=state.upper;
        for(std::size_t k=0;k<row.contacts.size();++k) row.contacts[k].second_z=state.lower[k];
        continue;
      }
      weights[row.upper] -= row.multiplier;
      for (auto& c : row.contacts) {
        weights[c.lower] += row.multiplier;
        if (c.score == 0) continue;
        weights[row.upper] -= c.first_q; weights[c.lower] -= c.second_q;
        c.first_z = c.second_z = 0;
        double value = 0;
        for (int a = 0; a <= 1; ++a) for (int b = 0; b <= 1; ++b) {
          const double candidate = c.score * a * b + c.first_q * a + c.second_q * b;
          if (candidate > value) { value = candidate; c.first_z = a; c.second_z = b; }
        }
        factor_value += value;
      }
    }
    const double constant = std::accumulate(lambda.begin(), lambda.end(), factor_value);
    double beam_value = constant, upper_bound = constant, endpoint_certificate = constant;
    int bound_evaluations = 0;
    const bool block_due = iteration == 0 || (iteration + 1) % options.bound_every == 0 || iteration + 1 == options.max_iterations;
    std::size_t pruned_states = 0;
    for (int level = 0; level < levels; ++level) {
      const double score = decoders[level]->decode(weights, allowed, decoded);
      beam_value += score;
      pruned_states += decoders[level]->pruned_states();
      if (!options.dp_beam() || (options.unpruned_bound && !decoders[level]->pruned_states())) {
        upper_bound += score; endpoint_certificate += score;
        level_bounds[level] = score; exact_levels[level] = 1;
      } else {
        std::fill(maxima.begin(), maxima.end(), 0);
        for (int id : by_level[level]) if (allowed[id]) {
          maxima[pairs[id].left] = std::max(maxima[pairs[id].left], weights[id]);
          maxima[pairs[id].right] = std::max(maxima[pairs[id].right], weights[id]);
        }
        const double endpoint_bound = std::accumulate(maxima.begin(), maxima.end(), 0.0) / 2;
        // Endpoint maxima count both ends, hence /2. This remains a bound
        // with signed weights, beam pruning, and the stricter stacking rule.
        endpoint_certificate += endpoint_bound;
        double level_bound = endpoint_bound; exact_levels[level] = 0;
        if (!block_bounds.empty() && block_due) {
          level_bound = std::min(level_bound,evaluate_block(level,weights,bound_evaluations));
        }
        upper_bound += level_bound; level_bounds[level] = level_bound;
      }
      for (int id : decoded) { selected[id] = 1; ++degree[pairs[id].left]; ++degree[pairs[id].right]; }
    }
    if (!std::isfinite(beam_value))
      throw std::overflow_error("DD dual value overflow; reduce the score scale");
    best_bound = std::min(best_bound, upper_bound);
    if (recovery) recovery->observe_adjusted(weights);
    // Save the oracle coefficients before primal recovery overwrites weights.
    std::vector<double> trace_weights;
    if ((trace.is_open() && options.trace_state) || !block_bounds.empty()) trace_weights = weights;
    // Recover a feasible primal in level order. Each DP sees only unused
    // bases and partners actually selected in EVERY required lower level.
    std::fill(recovered.begin(), recovered.end(), 0); std::fill(used.begin(), used.end(), 0);
    for (int level = 0; level < levels; ++level) {
      auto mask = allowed;
      for (int id : by_level[level]) {
        mask[id] = allowed[id] && !used[pairs[id].left] && !used[pairs[id].right];
        weights[id] = pairs[id].weight;
        for (int row_id : row_ids[id]) {
          bool supported = false;
          for (const auto& c : rows[row_id].contacts) if (recovered[c.lower]) supported = true;
          if(supported) {
            if(rows[row_id].best_partner) weights[id] += row_score(rows[row_id],recovered);
            else for(const auto& c:rows[row_id].contacts) if(recovered[c.lower]) weights[id]+=c.score;
          }
          mask[id] &= supported;
        }
      }
      if (level == 0) {
        decoded.clear();
        for (int id : by_level[0]) if (selected[id]) decoded.push_back(id);
      } else decoders[level]->decode(weights, mask, decoded);
      for (int id : decoded) { recovered[id] = 1; used[pairs[id].left] = used[pairs[id].right] = 1; }
    }
    double primal = 0, pk_value = 0;
    for (int id = 0; id < count; ++id) if (recovered[id]) {
      primal += pairs[id].weight;
      if (!projected_coefficients.empty()) pk_value += projected_coefficients[id];
      if (pairs[id].level > 0) pk_value -= pk.level_penalty;
    }
    for (const auto& row : rows) if (recovered[row.upper]) {
      if(row.best_partner) { const double value=row_score(row,recovered); primal+=value; pk_value+=value; }
      else for(const auto& c:row.contacts) if(recovered[c.lower]) { primal+=c.score; pk_value+=c.score; }
    }
    if (!std::isfinite(primal) || !std::isfinite(pk_value))
      throw std::overflow_error("DD primal objective overflow; reduce the score scale");
    const double baseline_primal = primal;
    const bool baseline_improved = baseline_primal > baseline_best;
    baseline_best = std::max(baseline_best,baseline_primal);
    result.iterations = iteration + 1;
    double full_norm = 0, projected_norm = 0;
    int capacity_violations = 0, support_violations = 0, copy_violations = 0;
    double complementarity = 0;
    const auto add_norm = [&](double g, double multiplier, bool inequality) {
      full_norm += g * g;
      if (!inequality || multiplier > 0 || g < 0) projected_norm += g * g;
    };
    for (int base = 0; base < length; ++base) {
      gradient[base] = 1 - degree[base];
      add_norm(gradient[base], lambda[base], true);
      capacity_violations += gradient[base] < 0;
      complementarity += std::abs(lambda[base] * gradient[base]);
    }
    for (const auto& row : rows) {
      if(row.best_partner) {
        add_norm(row.upper_z-selected[row.upper],0,false);
        copy_violations+=row.upper_z!=selected[row.upper];
        int support=0;
        for(const auto& c:row.contacts) {
          support+=selected[c.lower]; add_norm(c.second_z-selected[c.lower],0,false);
          copy_violations+=c.second_z!=selected[c.lower];
        }
        support_violations+=selected[row.upper] && !support;
        continue;
      }
      double g = -selected[row.upper];
      for (const auto& c : row.contacts) {
        g += selected[c.lower];
        if (c.score != 0) {
          add_norm(c.first_z - selected[row.upper], 0, false);
          add_norm(c.second_z - selected[c.lower], 0, false);
          copy_violations += c.first_z != selected[row.upper];
          copy_violations += c.second_z != selected[c.lower];
        }
      }
      add_norm(g, row.multiplier, true);
      support_violations += g < 0;
      complementarity += std::abs(row.multiplier * g);
    }
    const double norm = options.projected_norm ? projected_norm : full_norm;
    int recovery_proposals = 0;
    std::string recovery_method;
    const bool possible_stop = norm == 0 || iteration + 1 == options.max_iterations ||
        (options.patience && !baseline_improved && stagnant + 1 >= options.patience);
    // Before an off-cadence terminal iteration, tighten its certificate as
    // well. These are the saved dual coefficients, not recovery weights.
    if (!block_bounds.empty() && possible_stop && !block_due) {
      double terminal_bound = constant;
      for (int level=0;level<levels;++level) {
        if (!exact_levels[level]) {
          level_bounds[level] = std::min(level_bounds[level],evaluate_block(level,trace_weights,bound_evaluations));
        }
        terminal_bound += level_bounds[level];
      }
      upper_bound = std::min(upper_bound,terminal_bound);
      best_bound = std::min(best_bound,upper_bound);
    }
    if (recovery && (iteration == 0 || (iteration + 1) % options.recovery_every == 0 || possible_stop) &&
        (!std::isfinite(best_bound) ||
         best_bound - std::max(result.objective,primal) > 1e-9 * std::max(1.0,std::abs(best_bound)))) {
      auto proposal = recovery->propose(selected);
      recovery_proposals = proposal.proposals;
      result.recovery_proposals += recovery_proposals;
      if (proposal.objective > primal) {
        recovered = std::move(proposal.selected); primal = 0; pk_value = 0;
        for (int id=0;id<count;++id) if (recovered[id]) {
          primal += pairs[id].weight;
          if (!projected_coefficients.empty()) pk_value += projected_coefficients[id];
          if (pairs[id].level > 0) pk_value -= pk.level_penalty;
        }
        for (const auto& row : rows) if (recovered[row.upper]) for (const auto& c : row.contacts)
          if (recovered[c.lower]) { primal += c.score; pk_value += c.score; }
        recovery_method = proposal.method;
      }
    }
    if (exchange && (iteration == 0 || possible_stop ||
        (options.exchange_every && (iteration+1)%options.exchange_every == 0)) &&
        (!std::isfinite(best_bound) ||
         best_bound-std::max(result.objective,primal) > 1e-9*std::max(1.0,std::abs(best_bound)))) {
      auto seed=recovered;
      if (result.objective > primal) {
        std::fill(seed.begin(),seed.end(),0);
        for (int id : best) seed[id]=1;
      }
      auto proposal=exchange->improve(seed,options.exchange_width,options.exchange_passes,options.exchange_states);
      result.exchange_windows+=proposal.windows;
      result.exchange_improvements+=proposal.improved_windows;
      result.exchange_budget_windows+=proposal.budget_windows;
      result.exchange_states+=proposal.states_visited;
      // Compare using the same double accumulation order as the baseline.
      // The helper's extended-precision score may round differently.
      double candidate=0, candidate_pk=0;
      for(int id=0;id<count;++id) if(proposal.selected[id]) {
        candidate+=pairs[id].weight;
        if (!projected_coefficients.empty()) candidate_pk+=projected_coefficients[id];
        if(pairs[id].level>0) candidate_pk-=pk.level_penalty;
      }
      for(const auto& row:rows) if(proposal.selected[row.upper]) for(const auto& c:row.contacts)
        if(proposal.selected[c.lower]) { candidate+=c.score; candidate_pk+=c.score; }
      if (candidate > primal) {
        recovered=std::move(proposal.selected); primal=candidate; pk_value=candidate_pk;
        recovery_method="joint_exchange";
      }
    }
    if (!std::isfinite(primal) || !std::isfinite(pk_value))
      throw std::overflow_error("DD primal objective overflow; reduce the score scale");
    if (primal > result.objective) {
      result.objective = primal; result.pk_score = pk_value; best.clear();
      for (int id=0;id<count;++id) if (recovered[id]) best.push_back(id);
      stagnant = 0;
    } else if (baseline_improved) stagnant = 0;
    else ++stagnant;
    std::string stop;
    if (std::isfinite(best_bound) && std::isfinite(result.objective) && best_bound - result.objective <= 1e-9 * std::max(1.0, std::abs(best_bound))) stop = "bound_gap";
    else if (norm == 0) stop = options.dp_beam() ? "beam_stationary" : "exact_stationary";
    else if (options.patience && stagnant >= options.patience) stop = "patience";
    else if (iteration + 1 == options.max_iterations) stop = "max_iterations";
    if (beam_value < best_beam - 1e-9) best_beam = beam_value;
    else if (iteration && iteration % 10 == 0) relaxation = std::max(0.1, relaxation * 0.7);
    // With a beam, the Lagrangian value is not a dual bound. A diminishing
    // fallback avoids zero/reversed steps when a primal exceeds that value.
    // Additional proposals improve the returned incumbent. By default they
    // leave the original Polyak target unchanged, retaining its trajectory
    // until a valid bound proves an earlier stop. Best-target is optional.
    const double step_target = options.recovery_target_best ? result.objective : baseline_best;
    const double eta = !stop.empty() ? 0 : options.diminishing_step
        ? options.step * scale / (std::pow(iteration + 1.0, 0.75) * std::sqrt(norm))
        : beam_value > step_target + 1e-9
        ? relaxation * (beam_value - step_target) / norm
        : scale / std::sqrt((iteration + 1.0) * norm);
    if (trace.is_open()) {
      trace << "{\"event\":\"iteration\",\"solve_id\":" << solve_id
            << ",\"endpoint_certificate\":" << endpoint_certificate
            << ",\"bound_evaluations\":" << bound_evaluations
            << ",\"baseline_primal\":" << baseline_primal
            << ",\"baseline_lower_bound\":" << baseline_best
            << ",\"step_target\":" << step_target
            << ",\"recovery_proposals\":" << recovery_proposals
            << ",\"recovery_method\":\"" << recovery_method << "\""
            << ",\"pruned_states\":" << pruned_states
            << ",\"iteration\":" << iteration + 1 << ",\"primal\":" << primal
            << ",\"lower_bound\":" << result.objective << ",\"oracle_value\":" << beam_value
            << ",\"certificate\":" << upper_bound << ",\"upper_bound\":" << best_bound
            << ",\"full_norm_squared\":" << full_norm << ",\"projected_norm_squared\":" << projected_norm
            << ",\"capacity_violations\":" << capacity_violations << ",\"support_violations\":" << support_violations
            << ",\"copy_violations\":" << copy_violations << ",\"complementarity\":" << complementarity
            << ",\"stagnant\":" << stagnant << ",\"eta\":" << eta << ",\"relaxation\":" << relaxation
            << ",\"stop\":\"" << stop << "\"";
      if (options.trace_state) {
        trace << ",\"constant\":" << constant << ",\"weights\":[";
        for (int id = 0; id < count; ++id) { if (id) trace << ','; trace << trace_weights[id]; }
        trace << "],\"selected\":[";
        bool comma = false;
        for (int id = 0; id < count; ++id) if (selected[id]) { if (comma) trace << ','; trace << id; comma = true; }
        trace << "],\"recovered\":["; comma = false;
        for (int id = 0; id < count; ++id) if (recovered[id]) { if (comma) trace << ','; trace << id; comma = true; }
        trace << ']';
      }
      trace << "}\n";
    }
    if (!stop.empty()) { result.stop_reason = stop; break; }
    for (int base = 0; base < length; ++base) lambda[base] = std::max(0.0, lambda[base] - eta * gradient[base]);
    for (auto& row : rows) {
      if(row.best_partner) {
        row.upper_q-=eta*(row.upper_z-selected[row.upper]);
        for(auto& c:row.contacts) c.second_q-=eta*(c.second_z-selected[c.lower]);
        continue;
      }
      double g = -selected[row.upper];
      for (auto& c : row.contacts) {
        g += selected[c.lower];
        if (c.score != 0) {
          c.first_q -= eta * (c.first_z - selected[row.upper]);
          c.second_q -= eta * (c.second_z - selected[c.lower]);
        }
      }
      row.multiplier = std::max(0.0, row.multiplier - eta * g);
    }
  }
  result.upper_bound = best_bound;
  for (int id : best) {
    const auto& p = pairs[id];
    result.bpseq[p.left] = p.right; result.bpseq[p.right] = p.left;
    result.levels[p.left] = result.levels[p.right] = p.level;
  }
  emit_summary();
  return result;
}

std::vector<DDCrossingEvidence> dd_crossing_evidence(
    const PKPosteriorPairs& posterior, int crossing_beam) {
  if (posterior.empty() || crossing_beam < 0)
    throw std::invalid_argument("Invalid DD crossing evidence dimensions");
  struct Evidence { int left, right; double probability; };
  std::vector<Evidence> active;
  std::vector<DDCrossingEvidence> contacts;
  for (int left = 1; left < static_cast<int>(posterior.size()); ++left) {
    active.erase(std::remove_if(active.begin(), active.end(),
        [&](const Evidence& e) { return e.right <= left - 1; }), active.end());
    for (const auto& [right, value] : posterior[left]) {
      if (right == 0 || right >= posterior.size() || right == static_cast<unsigned>(left) ||
          !std::isfinite(value) || value < 0)
        throw std::invalid_argument("Invalid DD crossing evidence pair");
      if (right <= static_cast<unsigned>(left)) continue;
      for (const auto& e : active) if (e.right < static_cast<int>(right - 1))
        contacts.push_back({e.left, e.right, left - 1, static_cast<int>(right - 1), e.probability * value});
    }
    for (const auto& [right, value] : posterior[left]) if (right > static_cast<unsigned>(left)) {
      active.push_back({left - 1, static_cast<int>(right - 1), value});
      if (crossing_beam && active.size() > static_cast<std::size_t>(crossing_beam)) {
        auto better = [](const Evidence& a, const Evidence& b) {
          if (a.probability != b.probability) return a.probability > b.probability;
          return a.left != b.left ? a.left < b.left : a.right < b.right;
        };
        active.erase(std::max_element(active.begin(), active.end(), better));
      }
    }
  }
  return contacts;
}
