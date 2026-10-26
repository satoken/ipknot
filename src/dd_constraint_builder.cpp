#include "dd_constraint_builder.h"
#include "bpseq.h"
#include "dd_constrained.h"
#include "ip.h"
#include "spdlog/spdlog.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

namespace {
using Pair = std::pair<int, int>;
std::uint64_t key(Pair p) {
  return (std::uint64_t(unsigned(p.first)) << 32) | unsigned(p.second);
}
struct Physical {
  Pair pair;
  double probability;
  std::vector<int> columns;
};
struct Builder {
  const std::string &seq;
  const VSVF &posterior;
  const VF &th;
  const VF &alpha;
  int levels, length;
  const DDOptions &options;
  const VI &fixed;
  bool constrained;
  IPModel model;
  IP ip;
  std::vector<Physical> physical;
  std::unordered_map<std::uint64_t, int> lookup;
  std::unordered_map<std::uint64_t, double> probabilities;

  Builder(const std::string &s, const VSVF &bp, const VF &t, const VF &a, int k,
          const DDOptions &o, const VI &f, bool c)
      : seq(s), posterior(bp), th(t), alpha(a), levels(k), length(s.size()),
        options(o), fixed(f), constrained(c), ip(model) {
    for (int i = 1; i <= length; ++i)
      for (const auto &[j, p] : posterior[i]) {
        if (j == 0 || j >= posterior.size() || j == static_cast<unsigned>(i) ||
            !std::isfinite(p) || p < 0)
          throw std::invalid_argument("Invalid sparse posterior pair for DD");
        if (i < static_cast<int>(j))
          probabilities.emplace(key({i - 1, int(j) - 1}), p);
      }
  }
  bool compatible(Pair p) const {
    if (p.first < 0 || p.second >= length || p.first >= p.second)
      return false;
    if (!constrained)
      return true;
    for (auto end : {p.first, p.second}) {
      const int other = end == p.first ? p.second : p.first, value = fixed[end];
      if (value == BPSEQ::U || (value >= 0 && value != other) ||
          (value == BPSEQ::L && other < end) ||
          (value == BPSEQ::R && other > end))
        return false;
    }
    return true;
  }
  double probability(Pair p) const {
    auto it = probabilities.find(key(p));
    return it == probabilities.end() ? 0 : it->second;
  }
  int add(Pair p, double probability, bool all) {
    if (!compatible(p))
      return -1;
    auto inserted = lookup.emplace(key(p), physical.size());
    if (inserted.second)
      physical.push_back({p, probability, std::vector<int>(levels, -1)});
    auto &pair = physical[inserted.first->second];
    for (int lv = 0; lv < levels; ++lv) {
      if (pair.columns[lv] < 0 && (all || probability > th[lv]))
        pair.columns[lv] = ip.make_variable((probability - th[lv]) * alpha[lv]);
    }
    return inserted.first->second;
  }
  void add_pair_sum(int row, int id, double coefficient) {
    for (int col : physical[id].columns)
      if (col >= 0)
        ip.add_constraint(row, col, coefficient);
  }
};
} // namespace

void decode_linear_constraints(
    const std::string &sequence, const VSVF &posterior, const VF &thresholds,
    const VF &alpha, int levels, bool stacking,
    const DDOptions &options, VI &bpseq, VI &plevel, bool fixed) {
  options.validate();
  const int length = sequence.size();
  if (posterior.size() != sequence.size() + 1 ||
      thresholds.size() != static_cast<std::size_t>(levels) ||
      alpha.size() != thresholds.size() ||
      (fixed && bpseq.size() != sequence.size()))
    throw std::invalid_argument("Invalid constrained DD input dimensions");
  if (options.crossing_beam == 0 ||
      options.witnesses == 0 || options.constraint_states == 0)
    throw std::invalid_argument(
        "Linear constrained DD needs positive crossing beams and repair budget; use "
        "--dd-constraints full for unlimited diagnostics");
  if (fixed)
    for (int i = 0; i < length; ++i)
      if (bpseq[i] >= length || bpseq[i] == i || bpseq[i] < BPSEQ::LR ||
          (bpseq[i] >= 0 && bpseq[bpseq[i]] != i))
        throw std::invalid_argument("Invalid fixed structure for DD");
  const VI specification = fixed ? bpseq : VI(length, BPSEQ::DOT);
  Builder b(sequence, posterior, thresholds, alpha, levels, options,
            specification, fixed);
  for (int i = 1; i <= length; ++i)
    for (const auto &[j, p] : posterior[i])
      if (i < static_cast<int>(j)) {
        bool eligible = false;
        for (int lv = 0; lv < levels; ++lv)
          eligible |= p > thresholds[lv];
        const Pair pair{i - 1, int(j) - 1};
        if (eligible || (fixed && specification[i - 1] == int(j) - 1))
          b.add(pair, p, fixed && specification[i - 1] == int(j) - 1);
      }
  if (fixed)
    for (int i = 0; i < length; ++i)
      if (specification[i] > i)
        b.add({i, specification[i]}, b.probability({i, specification[i]}), true);
  std::vector<DDPair> pairs;
  std::vector<int> columns;
  std::vector<std::vector<int>> left(levels * length), right(levels * length),
      endpoints(length);
  for (const auto &physical : b.physical)
    for (int lv = 0; lv < levels; ++lv)
      if (physical.columns[lv] >= 0) {
        const int col = physical.columns[lv], i = physical.pair.first,
                  j = physical.pair.second;
        pairs.push_back({i, j, lv, b.model.variables[col].coefficient});
        columns.push_back(col);
        left[lv * length + i].push_back(col);
        right[lv * length + j].push_back(col);
        endpoints[i].push_back(col);
        endpoints[j].push_back(col);
      }
  for (int i = 0; i < length; ++i) {
    const int row =
        b.ip.make_constraint(IP::UP, 0, specification[i] == BPSEQ::U ? 0 : 1);
    for (int col : endpoints[i])
      b.ip.add_constraint(row, col, 1);
    const int value = specification[i];
    if (value == BPSEQ::L || value == BPSEQ::R || value == BPSEQ::LR) {
      if (endpoints[i].empty())
        throw DDInfeasible("Fixed paired-base constraint has no candidate "
                           "partner in retained DD graph");
      const int required = b.ip.make_constraint(IP::FX, 1, 1);
      for (int col : endpoints[i])
        b.ip.add_constraint(required, col, 1);
    } else if (value > i) {
      const int required = b.ip.make_constraint(IP::FX, 1, 1);
      b.add_pair_sum(required, b.lookup.at(key({i, value})), 1);
    }
  }
  if (stacking) {
    const int distance = 1;
    spdlog::debug("Stacking neighbor distance: {}", distance);
    for (int lv = 0; lv < levels; ++lv)
      for (int i = 0; i < length; ++i)
        for (auto *side : {&left, &right}) {
          const int row = b.ip.make_constraint(IP::LO, 0, 0);
          for (int col : (*side)[lv * length + i])
            b.ip.add_constraint(row, col, -1);
          for (int d = 1; d <= distance; ++d)
            for (int neighbor : {i - d, i + d})
              if (neighbor >= 0 && neighbor < length)
                for (int col : (*side)[lv * length + neighbor])
                  b.ip.add_constraint(row, col, 1);
        }
  }
  const auto graph = dd_bounded_graph(length, pairs, levels, options);
  for (const auto &support : graph.rows) {
    const int upper = columns[support.upper],
              row = b.ip.make_constraint(IP::LO, 0, 0);
    b.ip.add_constraint(row, upper, -1);
    for (int id : support.contacts)
      b.ip.add_constraint(row, columns[id], 1);
  }
  const auto result =
      solve_constrained_dd(length, pairs, columns, levels, b.model, options);
  bpseq.assign(length, -1);
  plevel.assign(length, -1);
  for (std::size_t id = 0; id < pairs.size(); ++id)
    if (b.model.solution[columns[id]] > .5) {
      const auto &p = pairs[id];
      bpseq[p.left] = p.right;
      bpseq[p.right] = p.left;
      plevel[p.left] = plevel[p.right] = p.level;
    }
  spdlog::info("DD: iterations={}, pairs={}, constraint_rows={}, nonzeros={}, "
               "objective={:.12g}, graph_upper_bound={:.12g}, "
               "repair_states={}, propagation_work={}, structural_work={}, "
               "crossing_beam_drops={}, witness_drops={}, linear=1, stop={}",
               result.iterations, pairs.size(), b.model.rows.size(),
               result.nonzeros, result.objective, result.upper_bound,
               result.repair_states, result.propagation_work, result.structural_work,
               graph.crossing_drops, graph.witness_drops, result.stop_reason);
}
