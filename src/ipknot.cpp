/*
 * $Id$
 * 
 * Copyright (C) 2010 Kengo Sato
 *
 * This file is part of IPknot.
 *
 * IPknot is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * IPknot is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with IPknot.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <cmath>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <iostream>
#include <iterator>
#include <list>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <cstdint>

#include "ipknot.h"
#include "ip.h"
#include "dd_constrained.h"
#include "dd_constraint_builder.h"
#include "bpseq.h"
#include "nmr_pair_index.h"
#include "nmr_sequence_index.h"
#include "nmr_candidate_filter.h"

#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
//#include "spdlog/sinks/basic_file_sink.h"

// Utility functions for base pair type normalization
char normalize_base(char c) {
  c = std::toupper(c);
  return (c == 'T') ? 'U' : c;
}

std::string normalize_base_pair_type(char a, char b) {
  // Normalize individual bases
  a = normalize_base(a);
  b = normalize_base(b);

  // Use the conventional names for canonical/wobble pairs.
  if ((a == 'G' && b == 'C') || (a == 'C' && b == 'G')) return "GC";
  if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A')) return "AU";
  if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G')) return "GU";

  // Give all remaining (non-canonical) types a direction-independent name.
  if (a > b) std::swap(a, b);

  // Return as string
  return std::string{a, b};
}

std::string normalize_base_pair_type(const std::string& bp_type) {
  if (bp_type.size() != 2) {
    throw std::invalid_argument("Base pair type must be exactly 2 characters: " + bp_type);
  }
  return normalize_base_pair_type(bp_type[0], bp_type[1]);
}

bool is_canonical_base_pair(const std::string& bp_type) {
  return bp_type == "GC" || bp_type == "AU" || bp_type == "GU";
}

bool BPConstraints::has_noncanonical_constraints() const {
  for (const auto& [bp_type, count] : constraints) {
    if (count >= 0 && !is_canonical_base_pair(bp_type)) {
      return true;
    }
  }
  return false;
}

using Pair = std::pair<int, int>;
struct NMRPairHash {
  size_t operator()(const Pair& pair) const {
    const auto key = (uint64_t(uint32_t(pair.first)) << 32) | uint32_t(pair.second);
    return std::hash<uint64_t>{}(key);
  }
};

static bool pair_encloses(const Pair& outer, const Pair& inner) {
  return outer.first < inner.first && inner.second < outer.second;
}

static bool pair_crosses(const Pair& a, const Pair& b) {
  return (a.first < b.first && b.first < a.second && a.second < b.second) ||
         (b.first < a.first && a.first < b.second && b.second < a.second);
}

static std::vector<Pair>
collect_candidate_pairs(const VVSVI& v_l) {
  std::vector<Pair> pairs;
  for (const auto& level : v_l) {
    for (int i = 0; i < static_cast<int>(level.size()); ++i) {
      for (const auto& [j, var] : level[i]) {
        if (i < static_cast<int>(j)) pairs.emplace_back(i, static_cast<int>(j));
      }
    }
  }
  std::sort(pairs.begin(), pairs.end());
  pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
  return pairs;
}

// Enumerate flush coaxial-stacking candidates.  Only two-base-pair NMR
// constraints are eligible: longer patterns continue to mean a conventional
// stack/bulge.  Supporting third-helix pairs are taken from the ordinary BPP
// candidates so the NMR observation does not invent an arbitrary branch.
static void
find_coaxial_instances(const std::string& seq, const VVSVI& v_l,
                       StackConstraints& stack_constraints) {
  std::set<std::string> observed_types;
  for (const auto& constraint : stack_constraints.constraints) {
    if (constraint.size() == 2) {
      observed_types.insert(constraint.bp_types.begin(), constraint.bp_types.end());
    }
  }
  if (observed_types.empty()) {
    spdlog::info("Found 0 flush coaxial-stacking instances");
    return;
  }
  const auto ordinary_pairs = collect_candidate_pairs(v_l);
  const NMRPairRangeIndex support_index(ordinary_pairs);
  auto find_supports = [&](int left_low, int left_high, int right_low, int right_high) {
    std::vector<size_t> indices;
    support_index.append(left_low, left_high, right_low, right_high, indices);
    std::vector<Pair> supports;
    supports.reserve(indices.size());
    for (size_t index : indices) supports.push_back(ordinary_pairs[index]);
    return supports;
  };

  const NMRSequenceIndex sequence_index(seq, normalize_base,
      [](char a, char b) { return normalize_base_pair_type(a, b); });
  struct ObservedPairs {
    std::vector<Pair> pairs;
    std::map<int, std::vector<Pair>> by_left, by_right;
  };
  std::map<std::string, ObservedPairs> observed;
  for (const auto& type : observed_types) {
    auto& index = observed[type];
    sequence_index.for_each_pair({NMRSequenceIndex::type_code(type)}, 4,
        [&](size_t i, size_t j) {
          const Pair pair{static_cast<int>(i), static_cast<int>(j)};
          index.pairs.push_back(pair);
          index.by_left[pair.first].push_back(pair);
          index.by_right[pair.second].push_back(pair);
        });
  }

  using InstanceKey = std::tuple<int, int, int, int, int, int>;
  std::set<InstanceKey> seen;
  auto add_instance = [&](int constraint_id, CoaxialKind kind,
                          const Pair& pair1, HelixFace face1,
                          const Pair& pair2, HelixFace face2,
                          std::vector<Pair> support_pairs) {
    if (support_pairs.empty()) return;
    const auto key = std::make_tuple(constraint_id, static_cast<int>(kind),
                                     pair1.first, pair1.second,
                                     pair2.first, pair2.second);
    if (!seen.insert(key).second) return;

    CoaxialInstance instance{pair1, pair2, face1, face2, kind,
                             constraint_id, std::move(support_pairs)};
    stack_constraints.coaxial_instances.push_back(std::move(instance));
  };

  for (int constraint_id = 0;
       constraint_id < static_cast<int>(stack_constraints.constraints.size());
       ++constraint_id) {
    const auto& constraint = stack_constraints.constraints[constraint_id];
    if (constraint.size() != 2) continue;

    for (int direction = 0; direction < 2; ++direction) {
      const auto& type1 = constraint.bp_types[direction == 0 ? 0 : 1];
      const auto& type2 = constraint.bp_types[direction == 0 ? 1 : 0];
      if (direction == 1 && type1 == type2) continue;
      const auto& pairs1 = observed.at(type1).pairs;
      const auto& pairs2_by_left = observed.at(type2).by_left;
      const auto& pairs2_by_right = observed.at(type2).by_right;

      // Closing helix coaxially stacked with its first direct child.
      for (const auto& closing : pairs1) {
        const auto children = pairs2_by_left.find(closing.first + 1);
        if (children == pairs2_by_left.end()) continue;
        for (const auto& child : children->second) {
          if (!pair_encloses(closing, child)) continue;
          auto supports = find_supports(child.second, closing.second,
                                       NMRPairRangeIndex::low, closing.second);
          add_instance(constraint_id, CoaxialKind::CLOSING_FIRST_CHILD,
                       closing, HelixFace::INNER, child, HelixFace::OUTER,
                       std::move(supports));
        }
      }

      // Two adjacent direct child helices with a common closing pair.
      for (const auto& child1 : pairs1) {
        const auto children = pairs2_by_left.find(child1.second + 1);
        if (children == pairs2_by_left.end()) continue;
        for (const auto& child2 : children->second) {
          auto supports = find_supports(NMRPairRangeIndex::low, child1.first,
                                       child2.second, NMRPairRangeIndex::high);
          add_instance(constraint_id, CoaxialKind::ADJACENT_CHILDREN,
                       child1, HelixFace::OUTER, child2, HelixFace::OUTER,
                       std::move(supports));
        }
      }

      // Last direct child coaxially stacked with the closing helix.
      for (const auto& child : pairs1) {
        const auto closings = pairs2_by_right.find(child.second + 1);
        if (closings == pairs2_by_right.end()) continue;
        for (const auto& closing : closings->second) {
          if (!pair_encloses(closing, child)) continue;
          auto supports = find_supports(closing.first, child.first,
                                       NMRPairRangeIndex::low, child.first);
          add_instance(constraint_id, CoaxialKind::LAST_CHILD_CLOSING,
                       child, HelixFace::OUTER, closing, HelixFace::INNER,
                       std::move(supports));
        }
      }
    }
  }
  spdlog::info("Found {} flush coaxial-stacking instances",
               stack_constraints.coaxial_instances.size());
}
//#include "spdlog/stopwatch.h"

IPknot::IPknot(uint pk_level, const float* alpha,
         bool levelwise, bool stacking_constraints, int n_th,
         bool require_canonical_neighbor,
         bool allow_coaxial_stacking,
         NMRConstraintOptions nmr_options,
         DDOptions dd_options)
    : pk_level_(pk_level),
      alpha_(alpha, alpha+pk_level_),
      levelwise_(levelwise),
      stacking_constraints_(stacking_constraints),
      n_th_(n_th),
      require_canonical_neighbor_(require_canonical_neighbor),
      allow_coaxial_stacking_(allow_coaxial_stacking),
      nmr_options_(nmr_options),
      dd_options_(dd_options)
{
    dd_options_.validate();
    if (dd_options_.enabled && !levelwise_)
      throw std::invalid_argument("DD requires levelwise prediction; remove --no-levelwise or use --decoder ilp");
}

#if 0
void IPknot::solve(const std::string& seq, const VF& bp, const VI& offset,
             const VF& th, VI& bpseq, VI& plevel, bool constraint,
             const BPConstraints& bp_constraints) const
{
    uint L = seq.size();
    IP ip(IP::MAX, n_th_);
    VVSVI v_l(pk_level_, VSVI(L));
    VVSVI v_r(pk_level_, VSVI(L));
    VI c_l(L, 0), c_r(L, 0);
    uint n=0;

    // make objective variables with their weights
    for (auto j=1; j!=L; ++j)
    {
      for (auto i=j-1; i!=-1u; --i)
      {
        const float& p=bp[offset[i+1]+(j+1)];
        for (auto lv=0; lv!=pk_level_; ++lv)
          if (p>th[lv])
          {
            const auto v_ij = ip.make_variable((p-th[lv])*alpha_[lv]);
            v_l[lv][i].emplace_back(j, v_ij);
            v_r[lv][j].emplace_back(i, v_ij);
            c_l[i]++; c_r[j]++;
            n++;
          }
      }
    }
    ip.update();

    if (n>0)
      solve(seq, ip, v_l, v_r, c_l, c_r, th, bpseq, plevel, constraint, bp_constraints);
    else
    {
      bpseq.resize(L);
      std::fill(std::begin(bpseq), std::end(bpseq), -1);
      plevel.resize(L);
      std::fill(std::begin(plevel), std::end(plevel), -1);
    }
  }
#endif

void IPknot::solve(const std::string& seq, const VSVF& bp,
             const VF& th, VI& bpseq, VI& plevel, bool constraint,
             const BPConstraints& bp_constraints,
             const StackConstraints& stack_constraints) const
{
    solve_with_penalty(seq, bp, th, bpseq, plevel, constraint,
                       bp_constraints, stack_constraints);
}

double IPknot::solve_with_penalty(const std::string& seq, const VSVF& bp,
             const VF& th, VI& bpseq, VI& plevel, bool constraint,
             const BPConstraints& bp_constraints,
             const StackConstraints& stack_constraints) const
{
    if (dd_options_.enabled && !constraint && !bp_constraints.has_constraints() && !stack_constraints.has_constraints()) {
      if (th.size() != pk_level_ || bp.size() != seq.size() + 1)
        throw std::invalid_argument("Invalid sparse posterior or threshold dimensions for DD");
      std::vector<DDPair> pairs;
      for (unsigned i = 1; i < bp.size(); ++i)
        for (const auto& [j, p] : bp[i]) {
          if (j == 0 || j >= bp.size() || j == i || !std::isfinite(p) || p < 0)
            throw std::invalid_argument("Invalid sparse posterior pair for DD");
          if (i >= j) continue;
          for (unsigned level = 0; level < pk_level_; ++level) {
            if (p > th[level]) pairs.push_back({static_cast<int>(i - 1), static_cast<int>(j - 1),
                static_cast<int>(level), (p - th[level]) * alpha_[level]});
          }
        }
      const auto result = solve_dual_decomposition(seq.size(), pairs, pk_level_,
          stacking_constraints_, dd_options_);
      bpseq = result.bpseq; plevel = result.levels;
      spdlog::info("DD: iterations={}, pairs={}, support_rows={}, contacts={}, "
                   "crossing_beam_drops={}, witness_drops={}, objective={:.12g}, "
                   "graph_upper_bound={:.12g}, stop={}",
                   result.iterations, result.pairs, result.support_rows, result.contacts,
                   result.crossing_beam_drops, result.witness_drops,
                   result.objective, result.upper_bound, result.stop_reason);
      return 0.0;
    }
    if (dd_options_.enabled && dd_options_.linear_constraints)
      return decode_linear_constraints(seq, bp, th, alpha_, pk_level_, stacking_constraints_,
          require_canonical_neighbor_, allow_coaxial_stacking_, nmr_options_,
          dd_options_, bpseq, plevel, constraint, bp_constraints, stack_constraints);
    uint L = seq.size();
    StackConstraints effective_stack_constraints = stack_constraints;
    effective_stack_constraints.coaxial_instances.clear();
    IPModel dd_model;
    auto ip_owner = dd_options_.enabled ? std::make_unique<IP>(dd_model)
                                       : std::make_unique<IP>(IP::MAX, n_th_);
    IP& ip = *ip_owner;
    if (dd_options_.enabled) {
      if (th.size() != pk_level_ || bp.size() != seq.size() + 1 ||
          (constraint && bpseq.size() != seq.size()))
        throw std::invalid_argument("Invalid constrained DD input dimensions");
      for (unsigned i = 1; i < bp.size(); ++i) for (const auto& [j, p] : bp[i])
        if (j == 0 || j >= bp.size() || j == i || !std::isfinite(p) || p < 0)
          throw std::invalid_argument("Invalid sparse posterior pair for DD");
      if (constraint) for (int i = 0; i < static_cast<int>(seq.size()); ++i) {
        if (bpseq[i] >= static_cast<int>(seq.size()) || bpseq[i] == i || bpseq[i] < BPSEQ::LR ||
            (bpseq[i] >= 0 && bpseq[bpseq[i]] != i))
          throw std::invalid_argument("Invalid fixed structure for DD");
      }
    }
    VVSVI v_l(pk_level_, VSVI(L));
    VVSVI v_r(pk_level_, VSVI(L));
    VI c_l(L, 0), c_r(L, 0);
    uint n=0;

    // Collect required non-canonical base pair types from constraints
    std::set<std::string> required_noncanonical_bp_types;

    // From base pair constraints
    for (const auto& [bp_type, count] : bp_constraints.constraints) {
      if (count >= 0 && !is_canonical_base_pair(bp_type)) {
        required_noncanonical_bp_types.insert(bp_type);
      }
    }

    // From stack constraints
    if (stack_constraints.has_constraints()) {
      for (const auto& sc : stack_constraints.constraints) {
        for (const auto& bp_type : sc.bp_types) {
          if (!is_canonical_base_pair(bp_type)) {
            required_noncanonical_bp_types.insert(bp_type);
          }
        }
      }
    }

    // make objective variables with their weights (canonical base pairs first)
    // Allocate a sparse membership index only for the NMR path. No dense L*L
    // table is introduced into LinearPartition's ordinary sparse path.
    const bool need_pair_membership = !required_noncanonical_bp_types.empty() ||
                                     effective_stack_constraints.has_constraints();
    std::vector<std::unordered_set<uint>> pair_membership(need_pair_membership ? L : 0);
    for (auto i=1; i<=L; ++i)
    {
      bool found_constraint_j = false;
      for (const auto [j, p]: bp[i])
        if (i<j)
        {
          for (auto lv=0; lv!=pk_level_; ++lv)
            if (p > th[lv]
                || (constraint && bpseq[i-1]==j-1))
            {
              const auto v_ij = ip.make_variable((p-th[lv])*alpha_[lv]);
              v_l[lv][i-1].emplace_back(j-1, v_ij);
              v_r[lv][j-1].emplace_back(i-1, v_ij);
              c_l[i-1]++; c_r[j-1]++;
              n++;
              if (need_pair_membership && lv == 0) pair_membership[i-1].insert(j-1);
            }
          if (constraint && bpseq[i-1]==j-1) found_constraint_j = true;
        }

      if (constraint && !found_constraint_j && bpseq[i-1]>=0 && i-1<bpseq[i-1])
      {
        const auto j = bpseq[i-1]+1;
        const auto p = 0.0;
        for (auto lv=0; lv!=pk_level_; ++lv)
        {
          const auto v_ij = ip.make_variable((p-th[lv])*alpha_[lv]);
          v_l[lv][i-1].emplace_back(j-1, v_ij);
          v_r[lv][j-1].emplace_back(i-1, v_ij);
          c_l[i-1]++; c_r[j-1]++;
          n++;
          if (need_pair_membership && lv == 0) pair_membership[i-1].insert(j-1);
        }
      }
    }

    if (allow_coaxial_stacking_ && effective_stack_constraints.has_constraints()) {
      find_coaxial_instances(seq, v_l, effective_stack_constraints);
    }

    // Helper function to check if a variable exists for position (i, j) at level 0
    auto has_variable_at = [&pair_membership](uint i_0idx, uint j_0idx) -> bool {
      return pair_membership[i_0idx].count(j_0idx) != 0;
    };

    // Add non-canonical base pairs that are required by constraints
    // This is done after canonical base pairs so we can check for canonical neighbors
    if (!required_noncanonical_bp_types.empty()) {
      const NMRSequenceIndex sequence_index(seq, normalize_base,
          [](char a, char b) { return normalize_base_pair_type(a, b); });
      std::vector<uint16_t> required_types;
      for (const auto& type : required_noncanonical_bp_types)
        required_types.push_back(NMRSequenceIndex::type_code(type));
      // Preserve i/j insertion order, including mutable neighbor checks, but
      // never visit the unrequested pair types or allocate a dense type table.
      sequence_index.for_each_pair(required_types, 4, [&](size_t left, size_t right) {
          const auto i = left + 1, j = right + 1;
          const std::string bp_type = normalize_base_pair_type(seq[left], seq[right]);
            // Check if this non-canonical pair is already added
            if (has_variable_at(i-1, j-1)) {
              return;
            }

            // If require_canonical_neighbor_ is set, check if canonical neighbor variable exists
            if (require_canonical_neighbor_) {
              bool has_canonical_neighbor = false;
              // Check above: (i-1, j+1) - positions (i-2, j) in 0-indexed
              if (i >= 2 && j < L) {
                if (has_variable_at(i-2, j)) {
                  has_canonical_neighbor = true;
                }
              }
              // Check below: (i+1, j-1) - positions (i, j-2) in 0-indexed
              if (!has_canonical_neighbor && i < L && j >= 2) {
                if (has_variable_at(i, j-2)) {
                  has_canonical_neighbor = true;
                }
              }
              if (!has_canonical_neighbor) {
                spdlog::debug("Skipping non-canonical base pair ({},{}) type {}: no canonical neighbor variable", i, j, bp_type);
                return;
              }
            }

            // Add with zero weight (will only be selected if constraint requires it)
            const auto p = 0.0;
            for (auto lv=0; lv!=pk_level_; ++lv) {
              const auto v_ij = ip.make_variable((p-th[lv])*alpha_[lv]);
              v_l[lv][i-1].emplace_back(j-1, v_ij);
              v_r[lv][j-1].emplace_back(i-1, v_ij);
              c_l[i-1]++; c_r[j-1]++;
              n++;
              if (lv == 0) pair_membership[i-1].insert(j-1);
            }
            spdlog::debug("Added non-canonical base pair ({},{}) type {} for constraint", i, j, bp_type);
      });
    }

    // Every pair belonging to a concrete stack/bulge instance must be an IP
    // candidate even when its posterior probability is below the threshold.
    // Otherwise an explicit NMR constraint could be silently skipped when no
    // ordinary candidate variable was created.
    if (effective_stack_constraints.has_constraints()) {
      std::set<std::pair<int, int>> required_instance_pairs;
      for (const auto& instance : effective_stack_constraints.instances) {
        required_instance_pairs.insert(instance.pairs.begin(), instance.pairs.end());
      }
      for (const auto& instance : effective_stack_constraints.coaxial_instances) {
        required_instance_pairs.insert(instance.pair1);
        required_instance_pairs.insert(instance.pair2);
      }
      for (const auto& [i, j] : required_instance_pairs) {
        if (i < 0 || j < 0 || static_cast<uint>(j) >= L || j < i + 4 ||
            has_variable_at(i, j)) {
          continue;
        }
        for (auto lv=0; lv!=pk_level_; ++lv) {
          const auto v_ij = ip.make_variable((0.0-th[lv])*alpha_[lv]);
          v_l[lv][i].emplace_back(j, v_ij);
          v_r[lv][j].emplace_back(i, v_ij);
          c_l[i]++; c_r[j]++;
          n++;
          if (lv == 0) pair_membership[i].insert(j);
        }
        spdlog::debug("Added base pair ({},{}) required by stack/bulge/coaxial constraint", i+1, j+1);
      }
    }
    if (dd_options_.enabled) {
      dd_model.optimize = [&]() {
        std::vector<DDPair> pairs;
        std::vector<int> columns;
        for (unsigned level = 0; level < pk_level_; ++level)
          for (int i = 0; i < static_cast<int>(L); ++i)
            for (const auto& [j, col] : v_l[level][i]) {
              pairs.push_back({i, static_cast<int>(j), static_cast<int>(level), dd_model.variables[col].coefficient});
              columns.push_back(col);
            }
        const auto result = solve_constrained_dd(L, pairs, columns, pk_level_, dd_model, dd_options_);
        if (dd_options_.noe_ilp)
          spdlog::info("DD NOE ILP: variables={}, rows={}, calls={}, cache_hits={}, time={:.6f}s, "
                       "primal_calls={}, primal_feasible={}, primal_time={:.6f}s",
                       result.noe_ilp_variables, result.noe_ilp_rows, result.noe_ilp_calls,
                       result.noe_ilp_cache_hits, result.noe_ilp_seconds,
                       result.noe_primal_calls, result.noe_primal_feasible, result.noe_primal_seconds);
        spdlog::info("DD: iterations={}, pairs={}, constraint_rows={}, objective={:.12g}, "
                     "model_upper_bound={:.12g}, repair_states={}, repair_budget_exhausted={}, stop={}",
                     result.iterations, pairs.size(), dd_model.rows.size(), result.objective,
                     result.upper_bound, result.repair_states, result.repair_budget_exhausted, result.stop_reason);
        return result.objective;
      };
    }
    ip.update();

    // Constraints must still be processed when the posterior threshold leaves
    // no ordinary candidate pairs.  Soft constraints then pay a violation
    // penalty; hard constraints correctly report infeasibility/no instances.
    if (n>0 || constraint || bp_constraints.has_constraints() ||
        effective_stack_constraints.has_constraints())
    {
      const double penalty = solve(seq, ip, v_l, v_r, c_l, c_r, th,
                                   bpseq, plevel, constraint, bp_constraints,
                                   effective_stack_constraints);
      return penalty;
    }
    else
    {
      bpseq.resize(L);
      std::fill(std::begin(bpseq), std::end(bpseq), -1);
      plevel.resize(L);
      std::fill(std::begin(plevel), std::end(plevel), -1);
      return 0.0;
    }
  }

double IPknot::solve(const std::string& seq, IP& ip, const VVSVI& v_l, const VVSVI& v_r, const VI& c_l, const VI& c_r,
             const VF& th, VI& bpseq, VI& plevel, bool constraint,
             const BPConstraints& bp_constraints,
             const StackConstraints& stack_constraints) const
{
    uint L = seq.size();
    struct CountViolationVars {
      std::string bp_type;
      int expected;
      std::vector<int> bp_vars;
      int excess_var;
      int missing_var;
    };
    std::vector<CountViolationVars> count_violation_vars;
    std::vector<std::pair<size_t, int>> stack_violation_vars;

    if (!constraint)
    {
      bpseq.resize(L);
      std::fill(bpseq.begin(), bpseq.end(), -2);
    }
    
    // constraint 1: each s_i is paired with at most one base
    for (auto i=0; i!=L; ++i)
    {
      auto row_l = -1, row_r = -1;
      switch (bpseq[i])
      {
        default:
        case BPSEQ::DOT: // no constraints
          row_l = row_r = ip.make_constraint(IP::UP, 0, 1);
          break;
        case BPSEQ::U: // unpaired
          row_l = row_r = ip.make_constraint(IP::UP, 0, 0);
          break;
        case BPSEQ::LR: // paired with left or right
          if (c_l[i]+c_r[i]>0)
            row_l = row_r = ip.make_constraint(IP::FX, 1, 1);
          break;
        case BPSEQ::L: // paired with right j
          if (c_l[i]>0) 
          {
            row_l = ip.make_constraint(IP::FX, 1, 1);
            row_r = ip.make_constraint(IP::UP, 0, 0);
          }
          break;
        case BPSEQ::R: // paired with left j
          if (c_r[i]>0) 
          {
            row_l = ip.make_constraint(IP::UP, 0, 0);
            row_r = ip.make_constraint(IP::FX, 1, 1);
          }
          break;
      }
      if (row_l<0 || row_r<0)
      {
        if (dd_options_.enabled)
          throw DDInfeasible("Fixed paired-base constraint has no candidate partner");
        spdlog::warn("invalid constraint for the base {}, ignored.", i+1);
        row_l = row_r = ip.make_constraint(IP::UP, 0, 1); // fallback to no constraint
      }
      
      for (auto lv=0; lv!=pk_level_; ++lv)
      {
        for (const auto [j, v_ij]: v_r[lv][i]) 
          ip.add_constraint(row_r, v_ij, 1);
        for (const auto [j, v_ij]: v_l[lv][i])
          ip.add_constraint(row_l, v_ij, 1);
      }

      if (bpseq[i]>=0 && i<bpseq[i]) // paired with j=bpseq[i]
      {
        const auto j = bpseq[i];
        int c=0;
        std::vector<int> vals(pk_level_, -1);
        for (auto lv=0; lv!=pk_level_; ++lv)
        {
          for (auto [temp, v_ij]: v_l[lv][i])
            if (j==temp) 
            { 
              vals[lv] = v_ij; 
              c++;
              break; 
            }
        }
        if (c>0)
        {
          auto row = ip.make_constraint(IP::FX, 1, 1);
          for (auto lv=0; lv!=pk_level_; ++lv)
            if (vals[lv]>=0)
              ip.add_constraint(row, vals[lv], 1);
        }
        else
          spdlog::warn("invalid constraint for the bases {} and {}, ignored.", i+1, bpseq[i]+1);
      }
    }

    if (levelwise_)
    {
      // constraint 2: disallow pseudoknots in x[lv]
      for (auto lv=0; lv!=pk_level_; ++lv)
        for (auto i=0; i<v_l[lv].size(); ++i)
          for (auto [j, v_ij]: v_l[lv][i])
            for (auto k=i+1; k<j; ++k)
              for (auto [l, v_kl]: v_l[lv][k])
                if (j<l)
                {
                  auto row = ip.make_constraint(IP::UP, 0, 1);
                  ip.add_constraint(row, v_ij, 1);
                  ip.add_constraint(row, v_kl, 1);
                }

      // constraint 3: any x[t]_kl must be pseudoknotted with x[u]_ij for t>u
      for (auto lv=1; lv!=pk_level_; ++lv)
        for (auto k=0; k<v_l[lv].size(); ++k)
          for (auto [l, v_kl]: v_l[lv][k])
            for (auto plv=0; plv!=lv; ++plv)
            {
              int row = ip.make_constraint(IP::LO, 0, 0);
              ip.add_constraint(row, v_kl, -1);
              for (auto i=0; i<k; ++i)
                for (auto [j, v_ij]: v_l[plv][i])
                  if (k<j && j<l)
                    ip.add_constraint(row, v_ij, 1);

              for (auto i=k+1; i<l; ++i)
                for (auto [j, v_ij]: v_l[plv][i])
                  if (l<j)
                    ip.add_constraint(row, v_ij, 1);
            }
    }

    const bool has_bulged_instance = std::any_of(
        stack_constraints.instances.begin(), stack_constraints.instances.end(),
        [](const StackInstance& instance) { return instance.has_bulge; });
    const int max_neighbor_distance = has_bulged_instance ? 2 : 1;
    if (stacking_constraints_)
    {
      // Relax ordinary stacking support only when the effective NMR witness
      // set actually contains a one-nucleotide bulge.  In particular,
      // --nmr-bulge-mode none must retain IPknot's original adjacent-pair
      // rule instead of broadening the feasible structure set merely because
      // an NMR stack constraint is present.
      spdlog::debug("Stacking neighbor distance: {}", max_neighbor_distance);
      for (auto lv=0; lv!=pk_level_; ++lv)
      {
        // upstream
        for (auto i=0; i<L; ++i)
        {
          int row = ip.make_constraint(IP::LO, 0, 0);
          for (auto [j, v_ji]: v_r[lv][i])
            ip.add_constraint(row, v_ji, -1);
          for (int d=1; d<=max_neighbor_distance; ++d) {
            if (i>=static_cast<uint>(d))
              for (auto [j, v_ji]: v_r[lv][i-d])
                ip.add_constraint(row, v_ji, 1);
            if (i+d<L)
              for (auto [j, v_ji]: v_r[lv][i+d])
                ip.add_constraint(row, v_ji, 1);
          }
        }

        // downstream
        for (auto i=0; i<L; ++i)
        {
          auto row = ip.make_constraint(IP::LO, 0, 0);
          for (auto [j, v_ij]: v_l[lv][i])
            ip.add_constraint(row, v_ij, -1);
          for (int d=1; d<=max_neighbor_distance; ++d) {
            if (i>=static_cast<uint>(d))
              for (auto [j, v_ij]: v_l[lv][i-d])
                ip.add_constraint(row, v_ij, 1);
            if (i+d<L)
              for (auto [j, v_ij]: v_l[lv][i+d])
                ip.add_constraint(row, v_ij, 1);
          }
        }
      }
    }

    // Add base pair type constraints if specified
    if (bp_constraints.has_constraints() && !seq.empty())
    {

      // Map to collect variables for each base pair type
      std::map<std::string, std::vector<int>> bp_type_vars;

      // Collect variables for each base pair type
      for (auto lv = 0; lv != pk_level_; ++lv) {
        for (auto i = 0; i < L; ++i) {
          for (const auto [j, v_ij] : v_l[lv][i]) {
            if (i < j && i < seq.size() && j < seq.size()) {
              std::string bp_type = normalize_base_pair_type(seq[i], seq[j]);
              bp_type_vars[bp_type].push_back(v_ij);
            }
          }
        }
      }

      // Add constraints for all base pair types
      for (const auto& [bp_type, count] : bp_constraints.constraints) {
        if (count >= 0) {
          auto it = bp_type_vars.find(bp_type);
          const bool lower_bound =
              nmr_options_.count_mode == NMRCountMode::LOWER_BOUND;
          int row = lower_bound
              ? ip.make_constraint(IP::LO, count, 0)
              : ip.make_constraint(IP::FX, count, count);
          std::vector<int> constrained_vars;
          if (it != bp_type_vars.end()) {
            for (int var : it->second) {
              ip.add_constraint(row, var, 1);
              constrained_vars.push_back(var);
            }
          }
          if (nmr_options_.soft) {
            const int max_deviation = std::max(static_cast<int>(L), count);
            int excess_var = -1;
            if (!lower_bound) {
              // exact: actual - excess + missing = expected
              excess_var = ip.make_variable(
                  -nmr_options_.count_penalty, 0, max_deviation);
              ip.add_constraint(row, excess_var, -1);
            }
            // lower bound: actual + missing >= observed
            const int missing_var = ip.make_variable(
                -nmr_options_.count_penalty, 0, max_deviation);
            ip.add_constraint(row, missing_var, 1);
            count_violation_vars.push_back(
                {bp_type, count, std::move(constrained_vars), excess_var,
                 missing_var});
            spdlog::debug(
                "Added soft base-pair constraint {} {} {} with penalty {}",
                bp_type, lower_bound ? ">=" : "=", count,
                nmr_options_.count_penalty);
          } else {
            // An empty sum is zero, so a positive requested count correctly
            // makes the hard model infeasible.
            spdlog::debug("Added hard base-pair constraint {}={}",
                          bp_type, count);
          }
        }
      }

    }

    // Non-canonical base pairs do not require any special stacking constraints
    // They can be used freely just like canonical base pairs

    // Add stack constraints if specified
    if (stack_constraints.has_constraints())
    {
      // Create level-independent base pair variables for stack constraints
      // bp_pair[i][j] = 1 if positions i and j form a base pair at any level
      // Only create variables for base pairs used in stack constraint instances
      std::unordered_map<Pair, int, NMRPairHash> bp_pair_vars;
      spdlog::stopwatch nmr_build_timer;

      // Index once instead of searching sparse rows again for every witness
      // and every topology blocker. Retain the first variable at each level,
      // as the former sparse-row lookup did.
      std::map<Pair, std::vector<int>> pair_level_vars;
      for (auto lv = 0; lv != pk_level_; ++lv) {
        for (int i = 0; i < static_cast<int>(L); ++i) {
          for (const auto& [j, var] : v_l[lv][i]) {
            auto& vars = pair_level_vars[{i, static_cast<int>(j)}];
            if (vars.empty()) vars.resize(pk_level_, -1);
            if (vars[lv] < 0) vars[lv] = var;
          }
        }
      }
      std::vector<Pair> pair_keys;
      pair_keys.reserve(pair_level_vars.size());
      for (const auto& [pair, vars] : pair_level_vars) pair_keys.push_back(pair);
      const NMRCandidateFilter candidate_filter(
          L, pair_keys, stacking_constraints_ ? max_neighbor_distance : 0);
      std::unique_ptr<NMRPairRangeIndex> topology_index;
      if (!stack_constraints.coaxial_instances.empty())
        topology_index = std::make_unique<NMRPairRangeIndex>(pair_keys);
      std::vector<std::vector<const StackInstance*>> stacks_by_observation(stack_constraints.constraints.size());
      std::vector<std::vector<const CoaxialInstance*>> coaxials_by_observation(stack_constraints.constraints.size());
      for (const auto& instance : stack_constraints.instances)
        stacks_by_observation.at(instance.constraint_id).push_back(&instance);
      for (const auto& instance : stack_constraints.coaxial_instances)
        coaxials_by_observation.at(instance.constraint_id).push_back(&instance);

      // First pass: identify all base pairs needed by stack constraint instances
      std::set<std::pair<int,int>> required_pairs;
      for (const auto& instance : stack_constraints.instances) {
        for (const auto& [i, j] : instance.pairs) {
          required_pairs.insert(std::make_pair(i, j));
        }
      }
      for (const auto& instance : stack_constraints.coaxial_instances) {
        required_pairs.insert(instance.pair1);
        required_pairs.insert(instance.pair2);
        for (const auto& support : instance.support_pairs) {
          required_pairs.insert(support);
        }
      }

      // Second pass: create bp_pair_var only for required pairs
      for (const auto& [i, j] : required_pairs) {
        const auto level_it = pair_level_vars.find({i, j});
        // A textual pattern match is not necessarily a feasible RNA base-pair
        // candidate (for example, it may violate the minimum hairpin length).
        if (level_it == pair_level_vars.end()) {
          continue;
        }

        // Create a new binary variable for this base pair position
        int bp_pair_var = ip.make_variable(0.0, 0, 1);
        bp_pair_vars[std::make_pair(i, j)] = bp_pair_var;

        // Add constraint: bp_pair_var = sum of v_ij across all levels
        // This means bp_pair_var = 1 iff the pair exists at any level
        int row = ip.make_constraint(IP::FX, 0, 0);  // bp_pair_var - sum(v_ij) = 0
        ip.add_constraint(row, bp_pair_var, 1);

        // Add all level-specific variables for this position
        for (int var : level_it->second) {
          if (var >= 0) ip.add_constraint(row, var, -1);
        }
      }

      spdlog::info("Created {} level-independent base pair variables for stack constraints", bp_pair_vars.size());

      // Exact-one per observation lets implications be aggregated: at integer
      // points sum(witnesses using pair) <= selected_pair is equivalent to
      // separate implications plus the no-sharing capacity. Its LP relaxation
      // is stronger, and it needs one row per pair rather than per incidence.
      // With sharing enabled the aggregation must be per observation.
      std::vector<int> all_instance_vars;
      std::vector<int> instance_to_constraint_id;
      std::map<int, std::vector<int>> global_pair_witnesses;
      std::map<std::tuple<int,int,int>, std::vector<int>> coaxial_face_vars;
      size_t unaggregated_pair_rows = 0;
      size_t pair_link_rows = 0;
      size_t unaggregated_topology_rows = 0;
      size_t topology_rows = 0;
      size_t coaxial_context_vars = 0;
      size_t pruned_stack_infeasible = 0;
      size_t pruned_coaxial_infeasible = 0;
      size_t pruned_coaxial_adjacent_bulge = 0;

      auto add_pair_sum = [&](int row, const Pair& pair, double coefficient) {
        for (int var : pair_level_vars.at(pair)) {
          if (var >= 0) ip.add_constraint(row, var, coefficient);
        }
      };
      auto add_pair_links = [&](const std::map<int, std::vector<int>>& groups) {
        for (const auto& [bp_var, witnesses] : groups) {
          int row = ip.make_constraint(IP::UP, 0, 0);
          for (int witness : witnesses) ip.add_constraint(row, witness, 1);
          ip.add_constraint(row, bp_var, -1);
          ++pair_link_rows;
        }
      };

      // For each stack constraint, at least one instance must be selected
      for (size_t constraint_id = 0; constraint_id < stack_constraints.constraints.size(); ++constraint_id)
      {
        std::vector<int> instance_vars;  // Variables representing each instance
        std::map<int, std::vector<int>> observation_pair_witnesses;
        std::vector<std::vector<int>> blocker_witnesses(pair_keys.size());
        auto register_pairs = [&](int witness, const std::vector<int>& bp_vars) {
          auto& groups = nmr_options_.allow_shared_stack_pairs
              ? observation_pair_witnesses : global_pair_witnesses;
          for (int bp_var : bp_vars) {
            groups[bp_var].push_back(witness);
            ++unaggregated_pair_rows;
          }
        };

        // For each instance of this constraint (now level-independent)
        for (const auto* instance_pointer : stacks_by_observation[constraint_id])
        {
          const auto& instance = *instance_pointer;

          // Find the level-independent bp_pair variables for this instance
          std::vector<int> bp_vars;
          bool all_pairs_found = true;
          for (const auto& [i, j] : instance.pairs)
          {
            auto pair_key = std::make_pair(i, j);
            auto it = bp_pair_vars.find(pair_key);
            if (it != bp_pair_vars.end())
            {
              bp_vars.push_back(it->second);
            }
            else
            {
              all_pairs_found = false;
              break;
            }
          }

          // Skip this instance if any base pair variable was not found
          if (!all_pairs_found)
          {
            continue;
          }
          if (!candidate_filter.has_endpoint_support(instance.pairs)) {
            ++pruned_stack_infeasible;
            continue;
          }

          // Create a binary variable for this instance (level-independent)
          // This variable is 1 if all base pairs in the stack are selected
          int instance_var = ip.make_variable(0.0, 0, 1);  // Binary variable with no weight
          ip.mark_noe_variable(instance_var);
          instance_vars.push_back(instance_var);

          all_instance_vars.push_back(instance_var);
          register_pairs(instance_var, bp_vars);
          instance_to_constraint_id.push_back(constraint_id);

          // Debug: log the positions for this instance
          if (spdlog::default_logger()->should_log(spdlog::level::debug)) {
            std::ostringstream pos_str;
            for (const auto& [i, j] : instance.pairs) {
              pos_str << "(" << (i+1) << "," << (j+1) << ") ";
            }
            spdlog::debug("Constraint {} instance: {}", constraint_id + 1, pos_str.str());
          }
        }  // end inst_idx loop

        // Flush coaxial-stacking alternatives for this NMR constraint.
        for (const auto* instance_pointer : coaxials_by_observation[constraint_id])
        {
          const auto& instance = *instance_pointer;
          const auto it1 = bp_pair_vars.find(instance.pair1);
          const auto it2 = bp_pair_vars.find(instance.pair2);
          if (it1 == bp_pair_vars.end() || it2 == bp_pair_vars.end()) continue;

          // This optional prior is local to the two observed helices, not a
          // ban on bulges elsewhere or in the third (supporting) helix.
          auto stem_continuation = [](const Pair& pair, HelixFace loop_face) {
            return loop_face == HelixFace::INNER
                ? Pair{pair.first - 1, pair.second + 1}
                : Pair{pair.first + 1, pair.second - 1};
          };
          const Pair continuation1 = stem_continuation(instance.pair1, instance.face1);
          const Pair continuation2 = stem_continuation(instance.pair2, instance.face2);
          const bool straight_stems = nmr_options_.coaxial_no_adjacent_bulge;
          if (straight_stems && (!candidate_filter.contains(continuation1) ||
                                 !candidate_filter.contains(continuation2))) {
            pruned_coaxial_adjacent_bulge += instance.support_pairs.size();
            continue;
          }

          // Terminal-pair blockers are identical for every support choice.
          // Factor that common implication through one context variable,
          // rather than repeating its coefficient for each third helix.
          auto shares_base = [](const Pair& a, const Pair& b) {
            return a.first == b.first || a.first == b.second ||
                   a.second == b.first || a.second == b.second;
          };
          using Index = NMRPairRangeIndex;
          const auto crossing_left = [](const Pair& pair) {
            return Index::Rectangle{Index::low, pair.first, pair.first, pair.second};
          };
          const auto crossing_right = [](const Pair& pair) {
            return Index::Rectangle{pair.first, pair.second, pair.second, Index::high};
          };
          const auto between = [](const Pair& outer, const Pair& inner) {
            return Index::Rectangle{outer.first, inner.first, inner.second, outer.second};
          };
          std::vector<size_t> common_blockers;
          Index::Rectangle common_between{0, 0, 0, 0};
          if (instance.kind == CoaxialKind::CLOSING_FIRST_CHILD)
            common_between = between(instance.pair1, instance.pair2);
          else if (instance.kind == CoaxialKind::LAST_CHILD_CLOSING)
            common_between = between(instance.pair2, instance.pair1);
          topology_index->append_union({crossing_left(instance.pair1), crossing_right(instance.pair1),
              crossing_left(instance.pair2), crossing_right(instance.pair2), common_between}, common_blockers);
          common_blockers.erase(std::remove_if(common_blockers.begin(), common_blockers.end(),
              [&](size_t index) { return shares_base(pair_keys[index], instance.pair1) ||
                                        shares_base(pair_keys[index], instance.pair2); }), common_blockers.end());
          // Exclude common blockers before support queries, not after emitting
          // and sorting them for every witness. Dense contexts often leave only
          // a small sparse remainder. Keep its indices in the global pair order.
          std::vector<Pair> support_candidates;
          std::vector<size_t> support_candidate_indices;
          support_candidates.reserve(pair_keys.size() - common_blockers.size());
          support_candidate_indices.reserve(support_candidates.capacity());
          size_t common_offset = 0;
          for (size_t index = 0; index < pair_keys.size(); ++index) {
            if (common_offset < common_blockers.size() && common_blockers[common_offset] == index) {
              ++common_offset;
              continue;
            }
            if (shares_base(pair_keys[index], instance.pair1) ||
                shares_base(pair_keys[index], instance.pair2)) continue;
            support_candidates.push_back(pair_keys[index]);
            support_candidate_indices.push_back(index);
          }
          const Index support_index(support_candidates);
          std::vector<int> context_witnesses;
          std::vector<size_t> support_blockers;

          // Make a separate witness for every possible third helix/closing
          // pair.  This lets the topology blockers refer to the exact
          // multibranch-loop context selected by the solver.
          for (const auto& support : instance.support_pairs) {
            const auto support_it = bp_pair_vars.find(support);
            if (support_it == bp_pair_vars.end()) continue;

            using Topology = NMRCandidateFilter::CoaxialTopology;
            const Topology topology = instance.kind == CoaxialKind::CLOSING_FIRST_CHILD
                ? Topology{instance.pair1, instance.pair2, support}
                : instance.kind == CoaxialKind::ADJACENT_CHILDREN
                    ? Topology{support, instance.pair1, instance.pair2}
                    : Topology{instance.pair2, support, instance.pair1};
            const std::vector<Pair> required{instance.pair1, instance.pair2, support};
            if (!candidate_filter.has_endpoint_support(required, &topology)) {
              ++pruned_coaxial_infeasible;
              continue;
            }
            if (straight_stems &&
                (!NMRCandidateFilter::compatible(continuation1, required, &topology) ||
                 !NMRCandidateFilter::compatible(continuation2, required, &topology))) {
              ++pruned_coaxial_adjacent_bulge;
              continue;
            }

            int instance_var = ip.make_variable(0.0, 0, 1);
            ip.mark_noe_variable(instance_var);
            instance_vars.push_back(instance_var);
            all_instance_vars.push_back(instance_var);
            context_witnesses.push_back(instance_var);
            register_pairs(instance_var,
                {it1->second, it2->second, support_it->second});
            instance_to_constraint_id.push_back(constraint_id);

            // If this witness is selected, its three helix-terminal pairs
            // form one planar context.  Reject any intervening pair that
            // would make an observed/supporting child indirect to the
            // selected multibranch-loop closing pair.
            support_blockers.clear();
            Index::Rectangle support_between1{0, 0, 0, 0};
            Index::Rectangle support_between2{0, 0, 0, 0};
            if (instance.kind == CoaxialKind::CLOSING_FIRST_CHILD)
              support_between1 = between(instance.pair1, support);
            else if (instance.kind == CoaxialKind::ADJACENT_CHILDREN) {
              support_between1 = between(support, instance.pair1);
              support_between2 = between(support, instance.pair2);
            } else
              support_between1 = between(instance.pair2, support);
            support_index.append_union({crossing_left(support), crossing_right(support),
                support_between1, support_between2}, support_blockers);
            for (size_t candidate_index : support_blockers) {
              const auto& blocker = support_candidates[candidate_index];
              // Base uniqueness already excludes any pair sharing a base
              // with a required terminal/support pair.
              if (shares_base(blocker, support)) continue;
              blocker_witnesses[support_candidate_indices[candidate_index]].push_back(instance_var);
              ++unaggregated_topology_rows;
            }

            if (nmr_options_.allow_shared_stack_pairs) {
              coaxial_face_vars[{instance.pair1.first, instance.pair1.second,
                                 static_cast<int>(instance.face1)}].push_back(instance_var);
              coaxial_face_vars[{instance.pair2.first, instance.pair2.second,
                                 static_cast<int>(instance.face2)}].push_back(instance_var);
            }
            spdlog::debug("Constraint {} coaxial witness: ({},{}) ({},{}) support ({},{})",
                          constraint_id + 1,
                          instance.pair1.first + 1, instance.pair1.second + 1,
                          instance.pair2.first + 1, instance.pair2.second + 1,
                          support.first + 1, support.second + 1);
          }
          if (!context_witnesses.empty() && (!common_blockers.empty() || straight_stems)) {
            const int context_var = ip.make_variable(0.0, 0, 1);
            ip.mark_noe_variable(context_var);
            int context_row = ip.make_constraint(IP::FX, 0, 0);
            ip.add_constraint(context_row, context_var, -1);
            for (int witness : context_witnesses) {
              ip.add_constraint(context_row, witness, 1);
            }
            ++coaxial_context_vars;
            if (straight_stems) {
              // A presence-only filter is not enough: the selected RNA must
              // actually contain the two unbulged continuation pairs.
              for (const Pair& continuation : {continuation1, continuation2}) {
                const int stem_row = ip.make_constraint(IP::UP, 0, 0);
                ip.add_constraint(stem_row, context_var, 1);
                add_pair_sum(stem_row, continuation, -1);
              }
            }
            for (size_t blocker_index : common_blockers) {
              blocker_witnesses[blocker_index].push_back(context_var);
              unaggregated_topology_rows += context_witnesses.size();
            }
          }
        }

        // Select exactly one explanation for each observation.  Additional
        // witnesses have no structural meaning: a witness only certifies that
        // the observation is satisfied by one concrete instance.  Exact-one
        // removes symmetric witness assignments without restricting the
        // underlying selected base pairs.  In soft mode the violation variable
        // is the alternative explanation when no instance is selected.
        if (instance_vars.empty() && !nmr_options_.soft) {
          spdlog::error("Stack constraint {} has no valid instances - cannot satisfy constraint",
                      constraint_id + 1);
          if (dd_options_.enabled)
            throw DDInfeasible("Stack constraint cannot be satisfied: no valid instances found");
          throw std::runtime_error("Stack constraint cannot be satisfied: no valid instances found");
        }

        int row = ip.make_constraint(IP::FX, 1, 1);
        for (int inst_var : instance_vars) ip.add_constraint(row, inst_var, 1);
        if (nmr_options_.soft) {
          const int violation_var = ip.make_variable(
              -nmr_options_.stack_penalty, 0, 1);
          ip.mark_noe_variable(violation_var);
          ip.add_constraint(row, violation_var, 1);
          stack_violation_vars.emplace_back(constraint_id, violation_var);
        }
        if (nmr_options_.allow_shared_stack_pairs) {
          add_pair_links(observation_pair_witnesses);
        }
        // One observation selects at most one witness, so the individual
        // w + selected_blocker <= 1 rows can be summed over its blocked
        // alternatives without changing any integer-feasible assignment.
        for (size_t blocker_index = 0; blocker_index < blocker_witnesses.size(); ++blocker_index) {
          const auto& witnesses = blocker_witnesses[blocker_index];
          if (witnesses.empty()) continue;
          int blocker_row = ip.make_constraint(IP::UP, 0, 1);
          for (int witness : witnesses) ip.add_constraint(blocker_row, witness, 1);
          add_pair_sum(blocker_row, pair_keys[blocker_index], 1);
          ++topology_rows;
        }
        spdlog::info(
            "Added {} stack constraint {} with {} possible instances "
            "(total instances so far: {})",
            nmr_options_.soft ? "soft" : "hard", constraint_id + 1,
            instance_vars.size(), all_instance_vars.size());
      }

      if (!nmr_options_.allow_shared_stack_pairs) {
        add_pair_links(global_pair_witnesses);
      }

      // One loop-facing end of a helix can have at most one coaxial partner.
      // With sharing disabled, the pair-link rows already imply this rule.
      for (const auto& [face, vars] : coaxial_face_vars) {
        if (vars.size() < 2) continue;
        int row = ip.make_constraint(IP::UP, 0, 1);
        for (int var : vars) ip.add_constraint(row, var, 1);
      }

      spdlog::info(
          "NMR candidate pruning: {} stack/bulge and {} coaxial proven infeasible; "
          "{} coaxial removed by experimental no-adjacent-bulge assumption",
          pruned_stack_infeasible, pruned_coaxial_infeasible,
          pruned_coaxial_adjacent_bulge);
      spdlog::info(
          "NMR formulation: {} witnesses, {} pair-link rows ({} separate), "
          "{} topology rows ({} separate), {} factored coaxial contexts; "
          "built in {:.6f}s",
          all_instance_vars.size(), pair_link_rows, unaggregated_pair_rows,
          topology_rows, unaggregated_topology_rows, coaxial_context_vars,
          nmr_build_timer.elapsed().count());

      // Update IP solver
      ip.update();

      // execute optimization
      spdlog::stopwatch solver_timer;
      const double objective = ip.solve();
      spdlog::debug("{} objective: {:.12g}", dd_options_.enabled ? "DD" : "IP", objective);
      spdlog::info("{} optimization finished in {:.6f}s", dd_options_.enabled ? "DD" : "IP", solver_timer.elapsed().count());

      // Log which stack constraint instances were selected
      for (size_t i = 0; i < all_instance_vars.size(); ++i)
      {
        double val = ip.get_value(all_instance_vars[i]);
        if (val > 0.5)
        {
          spdlog::info("Instance {} (constraint {}) selected with value {}",
                      i, instance_to_constraint_id[i] + 1, val);
        }
      }
    }
    else
    {
      // All backends must see auxiliary columns before optimization.
      ip.update();
      spdlog::stopwatch solver_timer;
      const double objective = ip.solve();
      spdlog::debug("{} objective: {:.12g}", dd_options_.enabled ? "DD" : "IP", objective);
      spdlog::info("{} optimization finished in {:.6f}s", dd_options_.enabled ? "DD" : "IP", solver_timer.elapsed().count());
    }

    double nmr_violation_penalty = 0.0;
    for (const auto& violation : count_violation_vars) {
      double actual = 0.0;
      for (int bp_var : violation.bp_vars) actual += ip.get_value(bp_var);
      const double excess = violation.excess_var >= 0
          ? ip.get_value(violation.excess_var) : 0.0;
      const double missing = ip.get_value(violation.missing_var);
      const double cost = nmr_options_.count_penalty * (excess + missing);
      nmr_violation_penalty += cost;
      if (excess > 0.5 || missing > 0.5) {
        spdlog::info(
            "Soft NMR base-pair constraint {} violated: expected {}, "
            "actual {:.0f}, missing {:.0f}, excess {:.0f}, penalty {}",
            violation.bp_type, violation.expected, actual, missing, excess,
            cost);
      } else {
        spdlog::info("Soft NMR base-pair constraint {}={} satisfied",
                     violation.bp_type, violation.expected);
      }
    }
    for (const auto& [constraint_id, violation_var] : stack_violation_vars) {
      const double violation = ip.get_value(violation_var);
      const double cost = nmr_options_.stack_penalty * violation;
      nmr_violation_penalty += cost;
      if (violation > 0.5) {
        spdlog::info(
            "Soft NMR stack constraint {} violated, penalty {}",
            constraint_id + 1, cost);
      } else {
        spdlog::info("Soft NMR stack constraint {} satisfied",
                     constraint_id + 1);
      }
    }

    // build the result
    bpseq.resize(L);
    std::fill(bpseq.begin(), bpseq.end(), -1);
    plevel.resize(L);
    std::fill(plevel.begin(), plevel.end(), -1);
    for (auto lv=0; lv!=pk_level_; ++lv)
      for (auto i=0; i<L; ++i)
        for (const auto [j, v_ij]: v_l[lv][i])
          if (ip.get_value(v_ij)>0.5)
          {
            bpseq[i]=j; bpseq[j]=i;
            plevel[i]=plevel[j]=lv;
          }

    if (!levelwise_)
      decompose_plevel(bpseq, plevel);

    return nmr_violation_penalty;
  }

auto IPknot::solve(const std::string& seq, const VSVF& bp,
             EnumParam<float>& ep, VI& bpseq, VI& plevel, bool constraint,
             const BPConstraints& bp_constraints,
             const StackConstraints& stack_constraints) const -> std::pair<float,float>
{
    uint L = seq.size();
    std::vector<float> th(ep.size());
    VI bpseq_temp, plevel_temp;
    VI max_bpseq, max_plevel;
    float max_fval=-100.0, max_fval_pk=-100.0;
    double max_nmr_penalty = 0.0;
    double max_selection_score = -std::numeric_limits<double>::infinity();
    spdlog::info("Search for the best thresholds by pseudo expected F-value:");
    spdlog::info("NMR automatic-threshold penalty scale: {}",
                 nmr_options_.threshold_penalty_scale);
    std::vector<DDCrossingEvidence> dd_evidence;
    if (dd_options_.enabled) dd_evidence = dd_crossing_evidence(bp, dd_options_.crossing_beam);
    double dd_sump = 0;
    for (const auto& contact : dd_evidence) dd_sump += contact.product;
    const auto sump = dd_options_.enabled ? static_cast<float>(dd_sump) : compute_sump_pk(bp);
    do {
      ep.get(th);
      uint i;
      for (i=1; i!=th.size(); i++)
        if (th[i-1]<th[i]) break;
      if (i!=th.size()) continue;
      bpseq_temp = bpseq;
      plevel_temp = plevel;
      double nmr_penalty;
      try {
        nmr_penalty = solve_with_penalty(seq, bp, th, bpseq_temp, plevel_temp, constraint,
            bp_constraints, stack_constraints);
      } catch (const DDInfeasible& error) {
        spdlog::info("Skipping infeasible DD thresholds: {}", error.what());
        continue;
      }
      const auto [sen, ppv, mcc, fval] = compute_expected_accuracy(bpseq_temp, bp);
      const auto pk_accuracy = [&]() {
        if (!dd_options_.enabled) return compute_expected_accuracy_pk(bpseq_temp, bp, sump);
        double etp = 0;
        int selected_contacts = 0;
        for (const auto& c : dd_evidence)
          if (bpseq_temp[c.left1] == c.right1 && bpseq_temp[c.left2] == c.right2) {
            etp += c.product; ++selected_contacts;
          }
        const double total = static_cast<double>(bpseq_temp.size()) * (bpseq_temp.size() - 1) / 2;
        return compute_expected_accuracy(etp, total - selected_contacts - sump + etp,
                                         selected_contacts - etp, sump - etp);
      };
      const auto [sen_pk, ppv_pk, mcc_pk, fval_pk] = pk_accuracy();
      if (spdlog::get_level() <= spdlog::level::info)
      {
        std::ostringstream th_ss;
        std::copy(th.begin(), th.end(), std::ostream_iterator<float>(th_ss, ","));
        spdlog::info("th={} pF={}, pF_pk={}, NMR penalty={}",
                     th_ss.str(), fval, fval_pk, nmr_penalty);
      }
      const double selection_score = fval + fval_pk -
          nmr_options_.threshold_penalty_scale * nmr_penalty;
      if (selection_score > max_selection_score)
      {
        max_selection_score = selection_score;
        max_fval = fval;
        max_fval_pk = fval_pk;
        max_nmr_penalty = nmr_penalty;
        max_bpseq = bpseq_temp;
        max_plevel = plevel_temp;
      }
    } while (!ep.succ());
    if (dd_options_.enabled && !std::isfinite(max_selection_score))
      throw DDInfeasible("No feasible structure for any DD threshold combination");
    bpseq = max_bpseq;
    plevel = max_plevel;
    spdlog::info("max pF={}, pF_pk={}, NMR penalty={}, selection score={}",
                 max_fval, max_fval_pk, max_nmr_penalty,
                 max_selection_score);

    return {max_fval, max_fval_pk};
  }

template < class T >
IPknot::EnumParam<T>::EnumParam(const std::vector<std::vector<T> >& p)
  : p_(p), m_(p.size()), v_(p.size(), 0)
{
  for (uint i=0; i!=p.size(); ++i)
    m_[i] = p[i].size();
}

template < class T >
uint IPknot::EnumParam<T>::size() const { return m_.size(); }

template < class T >
void IPknot::EnumParam<T>::get(std::vector<T>& q) const
{
  for (uint i=0; i!=v_.size(); ++i)
    q[i] = p_[i][v_[i]];
}

template < class T >
bool IPknot::EnumParam<T>::succ()
{
  return succ(m_.size(), &m_[0], &v_[0]);
}

template < class T >
bool IPknot::EnumParam<T>::succ(int n, const int* m, int* v)
{
  if (n==0) return true;
  if (++(*v)==*m)
  {
    *v=0;
    return succ(n-1, ++m, ++v);
  }
  return false;
}

int IPknot::decompose_plevel(const std::vector<int>& bpseq, std::vector<int>& plevel)
{
    // resolve the symbol of parenthsis by the graph coloring problem
    uint L=bpseq.size();
    
    // make an adjacent graph, in which pseudoknotted base-pairs are connected.
    std::vector< std::vector<int> > g(L);
    for (uint i=0; i!=L; ++i)
    {
      if (bpseq[i]<0 || bpseq[i]<=(int)i) continue;
      uint j=bpseq[i];
      for (uint k=i+1; k!=L; ++k)
      {
        uint l=bpseq[k];
        if (bpseq[k]<0 || bpseq[k]<=(int)k) continue;
        if (k<j && j<l)
        {
          g[i].push_back(k);
          g[k].push_back(i);
        }
      }
    }
    // vertices are indexed by the position of the left base
    std::vector<int> v;
    for (uint i=0; i!=bpseq.size(); ++i)
      if (bpseq[i]>=0 && (int)i<bpseq[i]) 
        v.push_back(i);
    // sort vertices by degree
    std::sort(v.begin(), v.end(), [&](int x, int y) { return g[y].size() < g[x].size(); });

    // determine colors
    std::vector<int> c(L, -1);
    int max_color=0;
    for (uint i=0; i!=v.size(); ++i)
    {
      // find the smallest color that is unused
      std::vector<int> used;
      for (uint j=0; j!=g[v[i]].size(); ++j)
        if (c[g[v[i]][j]]>=0) used.push_back(c[g[v[i]][j]]);
      std::sort(used.begin(), used.end());
      used.erase(std::unique(used.begin(), used.end()), used.end());
      int j=0;
      for (j=0; j!=(int)used.size(); ++j)
        if (used[j]!=j) break;
      c[v[i]]=j;
      max_color=std::max(max_color, j);
    }

    // renumber colors in decentant order by the number of base-pairs for each color
    std::vector<int> count(max_color+1, 0);
    for (uint i=0; i!=c.size(); ++i)
      if (c[i]>=0) count[c[i]]++;
    std::vector<int> idx(count.size());
    for (uint i=0; i!=idx.size(); ++i) idx[i]=i;
    sort(idx.begin(), idx.end(), [&](int x, int y) { return count[y] < count[x]; });
    std::vector<int> rev(idx.size());
    for (uint i=0; i!=rev.size(); ++i) rev[idx[i]]=i;
    plevel.resize(L);
    for (uint i=0; i!=c.size(); ++i)
      plevel[i]= c[i]>=0 ? rev[c[i]] : -1;

    return max_color+1;
  }

auto IPknot::compute_expected_accuracy(float etp, float etn, float efp, float efn) -> std::tuple<float,float,float,float>
{
    float sen, ppv, mcc, f;
    sen = ppv = mcc = f = 0;
    if (etp+efn!=0) sen = etp / (etp + efn);
    if (etp+efp!=0) ppv = etp / (etp + efp);
    if (etp+efp!=0 && etp+efn!=0 && etn+efp!=0 && etn+efn!=0)
      mcc = (etp*etn-efp*efn) / std::sqrt((etp+efp)*(etp+efn)*(etn+efp)*(etn+efn));
    if (sen+ppv!=0) f = 2*sen*ppv/(sen+ppv);

    return {sen, ppv, mcc, f};
  }

#if 0
auto IPknot::compute_expected_accuracy(const VI& bpseq, const VF& bp, const VI& offset) -> std::tuple<float,float,float,float>
{
    int L  = bpseq.size();
    int L2 = L*(L-1)/2;
    int N = 0;

    float sump = 0.0;
    float etp  = 0.0;

    for (uint i=0; i!=bp.size(); ++i) sump += bp[i];

    for (uint i=0; i!=bpseq.size(); ++i)
    {
      if (bpseq[i]!=-1 && bpseq[i]>(int)i)
      {
        etp += bp[offset[i+1]+bpseq[i]+1];
        N++;
      }
    }

    float etn = L2 - N - sump + etp;
    float efp = N - etp;
    float efn = sump - etp;

    return compute_expected_accuracy(etp, etn, efp, efn);
  }
#endif

auto IPknot::compute_expected_accuracy(const VI& bpseq, const VSVF& bp) -> std::tuple<float,float,float,float>
{
    int L  = bpseq.size();
    int L2 = L*(L-1)/2;
    int N = 0;

    float sump = 0.0;
    float etp  = 0.0;

    for (uint i=1; i!=bp.size(); ++i) 
      for (const auto [j, p]: bp[i])
        if (i<j)
        {
            sump += p; 
            if (bpseq[i-1]==j-1)
            {
              etp += p;
              N++;
            }
        }

    float etn = L2 - N - sump + etp;
    float efp = N - etp;
    float efn = sump - etp;

    return compute_expected_accuracy(etp, etn, efp, efn);
  }

auto IPknot::compute_expected_accuracy_pk(const VI& bpseq, const VSVF& bp) -> std::tuple<float,float,float,float>
{
    int L  = bpseq.size();
    int L2 = L*(L-1)/2;
    int N = 0;

    float sump = 0.0;
    float etp  = 0.0;

    for (uint i=1; i!=bp.size(); ++i) 
      for (const auto [j, p]: bp[i])
        if (i<j)
          for (uint k=i+1; k!=j; ++k)
            for (const auto [l, q]: bp[k])
              if (/*k<l &&*/ j<l)
              {
                sump += p*q;
                if (bpseq[i-1]==j-1 && bpseq[k-1]==l-1)
                {
                  etp += p*q;
                  N++;
                }
              } 

    float etn = L2 - N - sump + etp;
    float efp = N - etp;
    float efn = sump - etp;

    return compute_expected_accuracy(etp, etn, efp, efn);
  }

auto IPknot::compute_sump_pk(const VSVF& bp) -> float
{
    float sump = 0.0;

    for (uint i=1; i!=bp.size(); ++i) 
      for (const auto [j, p]: bp[i])
        if (i<j /*&& p>=DEFAULT_THRESHOLD*/)
          for (uint k=i+1; k!=j; ++k)
            for (const auto [l, q]: bp[k])
              if (/*k<l &&*/ j<l /*&& q>=DEFAULT_THRESHOLD*/)
                sump += p*q;

    return sump;
  }

auto IPknot::compute_expected_accuracy_pk(const VI& bpseq, const VSVF& bp, float sump) -> std::tuple<float,float,float,float>
{
    int L  = bpseq.size();
    int L2 = L*(L-1)/2;
    int N = 0;

    float etp  = 0.0;
    for (auto i=0; i!=bpseq.size(); ++i) 
    {
      const auto j=bpseq[i];
      if (j>=0 && i<j)
      {
        const auto r = std::lower_bound(std::begin(bp[i+1]), std::end(bp[i+1]), std::make_pair<uint,float>(j+1, 0.0));
        const auto p = r!=std::end(bp[i+1]) && r->first==j+1 ? r->second : 0.;
        for (auto k=i+1; k!=j; ++k)
        {
          const auto l=bpseq[k];
          if (k>=0 && k<l && j<l)
          {
            const auto r = std::lower_bound(std::begin(bp[k+1]), std::end(bp[k+1]), std::make_pair<uint,float>(l+1, 0.0));
            const auto q = r!=std::end(bp[k+1]) && r->first==l+1 ? r->second : 0.;
            etp += p*q;
            N++;
          }
        }
      }
    }

    float etn = L2 - N - sump + etp;
    float efp = N - etp;
    float efn = sump - etp;

    return compute_expected_accuracy(etp, etn, efp, efn);
  }

auto IPknot::check_pseudoknots(const VI& bpseq) -> VI
{
    std::vector<std::pair<int,int>> st;
    VI bpseq_pk(bpseq.size(), -1);
    for (auto i=0; i!=bpseq.size(); i++)
    {
      auto j = bpseq[i];
      if (i<j) 
      {
        st.emplace_back(i, j);
      }
      else if (j>=0)
      {
        for (auto it=std::rbegin(st); it!=std::rend(st); ++it) 
        {
          const auto [k, l] = *it;
          if (k==j && l==j)
          {
            st.erase((++it).base());
            break;
          }
          bpseq_pk[i]=j; bpseq_pk[j]=i;
          bpseq_pk[k]=l; bpseq_pk[l]=k;
        }
      }
    }
    return bpseq_pk;
  }

uint IPknot::length(const std::string& seq) { return seq.size(); }
uint IPknot::length(const std::list<std::string>& aln) { return aln.front().size(); }

// Explicit template instantiation
template class IPknot::EnumParam<float>;
