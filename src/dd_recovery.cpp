#include "dd_recovery.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

struct DDPrimalRecovery::Workspace {
  int length, levels, history_width;
  unsigned proposal_mask;
  bool cache_static;
  const std::vector<DDPair>& pairs;
  std::vector<unsigned char> allowed;
  std::vector<DDRecoveryRow> rows;
  std::vector<std::vector<int>> by_level, row_ids;
  std::vector<std::unique_ptr<DDNussinov>> owned_decoders;
  std::vector<DDNussinov*> decoders;
  // Only level zero changes between seeds. Keep this fixed support compact.
  std::vector<std::vector<double>> history;
  std::size_t history_next = 0, history_count = 0;
  std::vector<double> weights;
  std::vector<unsigned char> selected, mask;
  std::vector<int> decoded, used;
  DDRecoveryProposal original;
  std::vector<unsigned char> original_seed;
  bool original_ready = false, all_ready = false;
  DDRecoveryProposal all_upper;
  std::vector<unsigned char> all_seed;

  Workspace(int n, const std::vector<DDPair>& p, int l, int beam, bool lonely,
            const std::vector<unsigned char>& a,
            const std::vector<DDRecoveryRow>& r, int h, unsigned methods, bool cache,
            const std::vector<DDNussinov*>& shared, bool improved)
      : length(n), levels(l), history_width(h), proposal_mask(methods), cache_static(cache), pairs(p), allowed(a), rows(r),
        by_level(l), row_ids(p.size()), history(h), weights(p.size()),
        selected(p.size()), mask(p.size()), used(n) {
    for (int id = 0; id < static_cast<int>(p.size()); ++id)
      by_level[p[id].level].push_back(id);
    for (int id = 0; id < static_cast<int>(rows.size()); ++id)
      row_ids[rows[id].upper].push_back(id);
    if (!shared.empty()) decoders = shared;
    else for (int level = 0; level < l; ++level) {
      owned_decoders.push_back(std::make_unique<DDNussinov>(n, p, level, beam, lonely, improved));
      decoders.push_back(owned_decoders.back().get());
    }
  }

  void original_weights() {
    for (std::size_t id = 0; id < pairs.size(); ++id) weights[id] = pairs[id].weight;
  }

  // Complete a level-zero seed in level order. When level k is optimized,
  // every contact with an already selected lower level is charged exactly;
  // all required witness rows are checked independently.
  DDRecoveryProposal complete(const std::vector<unsigned char>& seed) {
    selected = seed;
    std::fill(used.begin(), used.end(), 0);
    for (int id : by_level[0]) if (selected[id])
      used[pairs[id].left] = used[pairs[id].right] = 1;
    for (int level = 1; level < levels; ++level) {
      mask = allowed;
      for (int id : by_level[level]) {
        const auto& p = pairs[id];
        mask[id] = allowed[id] && !used[p.left] && !used[p.right];
        weights[id] = p.weight;
        for (int r : row_ids[id]) {
          bool supported = false;
          for (const auto& c : rows[r].contacts) if (selected[c.first]) {
            supported = true;
            weights[id] += c.second;
          }
          mask[id] &= supported;
        }
      }
      decoders[level]->decode(weights, mask, decoded);
      for (int id : decoded) {
        selected[id] = 1;
        used[pairs[id].left] = used[pairs[id].right] = 1;
      }
    }
    DDRecoveryProposal result;
    result.selected = selected;
    for (std::size_t id = 0; id < pairs.size(); ++id)
      if (selected[id]) result.objective += pairs[id].weight;
    for (const auto& row : rows) if (selected[row.upper])
      for (const auto& c : row.contacts) if (selected[c.first]) result.objective += c.second;
    // A seed may have a negative original score despite a profitable proposal
    // coefficient. The empty structure is a feasible zero lower bound.
    if (!(result.objective > 0)) {
      result.objective = 0;
      std::fill(result.selected.begin(), result.selected.end(), 0);
    }
    return result;
  }

  std::vector<unsigned char> seed() {
    decoders[0]->decode(weights, allowed, decoded);
    std::vector<unsigned char> result(pairs.size(), 0);
    for (int id : decoded) result[id] = 1;
    return result;
  }

  void guide(const std::vector<unsigned char>* upper_selected, bool opportunity_cost = false) {
    original_weights();
    if (opportunity_cost && upper_selected) {
      // A level-zero pair using an endpoint of a proposed upper pair prevents
      // that upper pair from surviving completion. Per-base positive reward
      // maxima estimate this opportunity cost in O(n+M). The /2 convention
      // keeps two endpoint penalties from charging the same reward twice.
      // This heuristic is only a seed score, never an upper-bound formula.
      std::vector<double> lost(length);
      for (std::size_t id = 0; id < pairs.size(); ++id)
        if (pairs[id].level > 0 && allowed[id] && (*upper_selected)[id]) {
          const double reward = std::max(0.0, pairs[id].weight);
          lost[pairs[id].left] = std::max(lost[pairs[id].left], reward);
          lost[pairs[id].right] = std::max(lost[pairs[id].right], reward);
        }
      for (int id : by_level[0]) {
        weights[id] -= (lost[pairs[id].left] + lost[pairs[id].right]) / 2;
      }
    }
    for (const auto& row : rows) {
      const auto& upper = pairs[row.upper];
      if (!allowed[row.upper] || row.contacts.empty() ||
          (upper_selected && !(*upper_selected)[row.upper])) continue;
      // Only rows into level zero influence the seed. An upper pair's reward
      // is shared across its required rows and alternative witnesses. This
      // is a proposal heuristic, so multiple overlapping upper rewards may
      // still be counted; completion and evaluation remove that optimism.
      const double divisor = static_cast<double>(upper.level) * row.contacts.size();
      for (const auto& c : row.contacts) if (pairs[c.first].level == 0 && allowed[c.first])
        weights[c.first] += std::max(0.0, upper.weight + c.second) / divisor;
    }
  }
};

DDPrimalRecovery::DDPrimalRecovery(int n, const std::vector<DDPair>& pairs,
    int levels, int beam, bool no_lonely,
    const std::vector<unsigned char>& allowed,
    const std::vector<DDRecoveryRow>& rows, int history, unsigned proposal_mask,
    bool cache_static, const std::vector<DDNussinov*>& shared_decoders, bool improved_beam) {
  if (n < 0 || levels < 1 || beam < 0 || history < 1 || history > 32 ||
      allowed.size() != pairs.size() || !proposal_mask || proposal_mask > 7)
    throw std::invalid_argument("Invalid DD recovery dimensions/history");
  if (!shared_decoders.empty() && (shared_decoders.size() != static_cast<std::size_t>(levels) ||
      std::find(shared_decoders.begin(),shared_decoders.end(),nullptr) != shared_decoders.end()))
    throw std::invalid_argument("Invalid shared DD recovery decoders");
  for (const auto& p : pairs) if (p.left < 0 || p.right >= n || p.left >= p.right ||
      p.level < 0 || p.level >= levels || !std::isfinite(p.weight))
    throw std::invalid_argument("Invalid DD recovery pair");
  std::vector<std::vector<unsigned char>> lower_rows(pairs.size());
  for (std::size_t id = 0; id < pairs.size(); ++id)
    lower_rows[id].resize(pairs[id].level, 0);
  for (const auto& row : rows) {
    if (row.upper < 0 || row.upper >= static_cast<int>(pairs.size()) ||
        pairs[row.upper].level == 0)
      throw std::invalid_argument("Invalid DD recovery witness row");
    int lower = -1;
    for (const auto& c : row.contacts) {
      if (c.first < 0 || c.first >= static_cast<int>(pairs.size()) ||
          pairs[c.first].level >= pairs[row.upper].level || !std::isfinite(c.second))
        throw std::invalid_argument("Invalid DD recovery contact");
      const auto& u = pairs[row.upper];
      const auto& l = pairs[c.first];
      if (!((u.left < l.left && l.left < u.right && u.right < l.right) ||
            (l.left < u.left && u.left < l.right && l.right < u.right)))
        throw std::invalid_argument("Noncrossing DD recovery contact");
      if (lower < 0) lower = pairs[c.first].level;
      else if (lower != pairs[c.first].level)
        throw std::invalid_argument("Mixed lower levels in DD recovery row");
    }
    // Empty rows are legitimate only for candidates already eliminated by
    // support propagation. They cannot reveal which lower level was required.
    if (lower < 0) {
      if (allowed[row.upper])
        throw std::invalid_argument("Allowed DD recovery pair has an empty witness row");
    } else if (lower_rows[row.upper][lower]++) {
      throw std::invalid_argument("Duplicate DD recovery lower-level row");
    }
  }
  for (std::size_t id = 0; id < pairs.size(); ++id) if (allowed[id])
    for (unsigned char count : lower_rows[id]) if (!count)
      throw std::invalid_argument("Missing DD recovery lower-level witness row");
  workspace_ = std::make_unique<Workspace>(n, pairs, levels, beam, no_lonely, allowed, rows, history, proposal_mask, cache_static, shared_decoders, improved_beam);
}

DDPrimalRecovery::~DDPrimalRecovery() = default;

void DDPrimalRecovery::observe_adjusted(const std::vector<double>& weights) {
  auto& w = *workspace_;
  if (weights.size() != w.pairs.size())
    throw std::invalid_argument("DD recovery adjusted weight size mismatch");
  auto& slot = w.history[w.history_next];
  slot.resize(w.by_level[0].size());
  for (std::size_t k = 0; k < slot.size(); ++k) {
    slot[k] = weights[w.by_level[0][k]];
    if (!std::isfinite(slot[k])) throw std::invalid_argument("Nonfinite DD recovery adjusted weight");
  }
  w.history_next = (w.history_next + 1) % w.history_width;
  w.history_count = std::min(w.history_count + 1, static_cast<std::size_t>(w.history_width));
}

DDRecoveryProposal DDPrimalRecovery::propose(const std::vector<unsigned char>& dual_selected) {
  auto& w = *workspace_;
  if (dual_selected.size() != w.pairs.size())
    throw std::invalid_argument("DD recovery selected size mismatch");
  DDRecoveryProposal result;
  result.selected.assign(w.pairs.size(), 0);
  std::vector<std::vector<unsigned char>> seen;
  int proposals = 0;
  // Retain completed seeds only while initializing the static all-upper
  // cache, so duplicate seeds reuse exactly the same deterministic completion.
  const bool cache_cold = w.cache_static && (w.proposal_mask & 4) && !w.all_ready;
  std::vector<DDRecoveryProposal> cold_completed;
  if ((w.proposal_mask & 1) && !w.original_ready) {
    w.original_weights();
    auto seed = w.seed();
    seen.push_back(seed);
    w.original = w.complete(seed);
    w.original.method = "original";
    w.original_seed = seed;
    w.original_ready = true;
    ++proposals;
  }
  if (w.proposal_mask & 1) {
    result = w.original;
    if (seen.empty()) seen.push_back(w.original_seed);
    if (cache_cold) cold_completed.push_back(w.original);
  }
  const auto consider = [&](const char* method, bool static_all = false) {
    auto seed = w.seed();
    auto found = std::find(seen.begin(), seen.end(), seed);
    if (found != seen.end()) {
      if (static_all && cache_cold) {
        w.all_upper = cold_completed[static_cast<std::size_t>(found-seen.begin())];
        w.all_upper.method = method; w.all_seed=seed; w.all_ready=true;
      }
      return;
    }
    seen.push_back(seed);
    auto candidate = w.complete(seed);
    candidate.method = method;
    ++proposals;
    if (cache_cold) cold_completed.push_back(candidate);
    if (static_all && cache_cold) {
      w.all_upper=candidate; w.all_seed=seed; w.all_ready=true;
    }
    if (candidate.objective > result.objective) result = std::move(candidate);
  };
  if ((w.proposal_mask & 2) && w.history_count) {
    w.original_weights();
    for (std::size_t k = 0; k < w.by_level[0].size(); ++k) {
      double sum = 0;
      for (std::size_t slot = 0; slot < w.history_count; ++slot) sum += w.history[slot][k];
      w.weights[w.by_level[0][k]] = sum / w.history_count;
    }
    consider("mean");
  }
  if (w.proposal_mask & 4) {
    w.guide(&dual_selected);
    consider("selected_upper");
    if (w.cache_static && w.all_ready) {
      if (std::find(seen.begin(),seen.end(),w.all_seed)==seen.end()) {
        seen.push_back(w.all_seed);
        if (w.all_upper.objective > result.objective) result=w.all_upper;
      }
    } else {
      w.guide(nullptr);
      consider("all_upper",true);
    }
    // Keep the unpenalized proposals too: opportunity estimates can be too
    // pessimistic, especially for signed contacts or overlapping upper pairs.
    w.guide(&dual_selected, true);
    consider("upper_opportunity");
  }
  result.proposals = proposals;
  return result;
}
