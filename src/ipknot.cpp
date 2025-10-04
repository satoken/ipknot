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
#include <sstream>

#include "ipknot.h"
#include "ip.h"
#include "bpseq.h" 

#include "spdlog/spdlog.h"
//#include "spdlog/sinks/basic_file_sink.h"
//#include "spdlog/stopwatch.h"

IPknot::IPknot(uint pk_level, const float* alpha,
         bool levelwise, bool stacking_constraints, int n_th)
    : pk_level_(pk_level),
      alpha_(alpha, alpha+pk_level_),
      levelwise_(levelwise),
      stacking_constraints_(stacking_constraints),
      n_th_(n_th)
{
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
    uint L = seq.size();
    IP ip(IP::MAX, n_th_);
    VVSVI v_l(pk_level_, VSVI(L));
    VVSVI v_r(pk_level_, VSVI(L));
    VI c_l(L, 0), c_r(L, 0);
    uint n=0;

    // Helper to get base pair type
    auto get_bp_type = [](char a, char b) -> std::string {
      a = std::toupper(a);
      b = std::toupper(b);
      if (a > b) std::swap(a, b);
      if ((a == 'G' && b == 'C') || (a == 'C' && b == 'G')) return "GC";
      if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A')) return "AU";
      if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G')) return "GU";
      if (a == 'U' && b == 'U') return "UU";
      if ((a == 'C' && b == 'U') || (a == 'U' && b == 'C')) return "CU";
      std::string result;
      result += a;
      result += b;
      return result;
    };

    // Helper to check if a base pair type is canonical
    auto is_canonical = [](const std::string& bp_type) -> bool {
      return bp_type == "GC" || bp_type == "AU" || bp_type == "GU";
    };

    // Check if non-canonical base pairs are needed
    bool has_noncanonical_bp_constraint = bp_constraints.has_noncanonical_constraints();
    bool has_noncanonical_stack_constraint = false;
    if (stack_constraints.has_constraints()) {
      for (const auto& sc : stack_constraints.constraints) {
        for (const auto& bp_type : sc.bp_types) {
          if (!is_canonical(bp_type)) {
            has_noncanonical_stack_constraint = true;
            break;
          }
        }
        if (has_noncanonical_stack_constraint) break;
      }
    }

    bool add_noncanonical_pairs = has_noncanonical_bp_constraint || has_noncanonical_stack_constraint;

    // make objective variables with their weights
    for (auto i=1; i<=L; ++i)
    {
      bool found_constraint_j = false;
      for (const auto [j, p]: bp[i])
        if (i<j)
        {
          for (auto lv=0; lv!=pk_level_; ++lv)
            if (p>th[lv] || (constraint && bpseq[i-1]==j-1))
            {
              const auto v_ij = ip.make_variable((p-th[lv])*alpha_[lv]);
              v_l[lv][i-1].emplace_back(j-1, v_ij);
              v_r[lv][j-1].emplace_back(i-1, v_ij);
              c_l[i-1]++; c_r[j-1]++;
              n++;
            }
          if (constraint && bpseq[i-1]==j-1) found_constraint_j = true;
        }

      if (constraint && !found_constraint_j && bpseq[i-1]>=0)
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
        }
      }

      // Add non-canonical base pairs if needed
      if (add_noncanonical_pairs) {
        for (auto j=i+1; j<=L; ++j) {
          std::string bp_type = get_bp_type(seq[i-1], seq[j-1]);
          if (!is_canonical(bp_type)) {
            // Check if this non-canonical pair is already added
            bool already_added = false;
            for (const auto [jj, v_ij]: v_l[0][i-1]) {
              if (jj == j-1) {
                already_added = true;
                break;
              }
            }

            if (!already_added) {
              // Add with zero weight (will only be selected if constraint requires it)
              const auto p = 0.0;
              for (auto lv=0; lv!=pk_level_; ++lv) {
                const auto v_ij = ip.make_variable((p-th[lv])*alpha_[lv]);
                v_l[lv][i-1].emplace_back(j-1, v_ij);
                v_r[lv][j-1].emplace_back(i-1, v_ij);
                c_l[i-1]++; c_r[j-1]++;
                n++;
              }
              spdlog::debug("Added non-canonical base pair ({},{}) type {} for constraint", i, j, bp_type);
            }
          }
        }
      }
    }
    ip.update();

    if (n>0)
      solve(seq, ip, v_l, v_r, c_l, c_r, th, bpseq, plevel, constraint, bp_constraints, stack_constraints);
    else
    {
      bpseq.resize(L);
      std::fill(std::begin(bpseq), std::end(bpseq), -1);
      plevel.resize(L);
      std::fill(std::begin(plevel), std::end(plevel), -1);
    }
  }

void IPknot::solve(const std::string& seq, IP& ip, const VVSVI& v_l, const VVSVI& v_r, const VI& c_l, const VI& c_r,
             const VF& th, VI& bpseq, VI& plevel, bool constraint,
             const BPConstraints& bp_constraints,
             const StackConstraints& stack_constraints) const
{
    uint L = seq.size();
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

    if (stacking_constraints_)
    {
      for (auto lv=0; lv!=pk_level_; ++lv)
      {
        // upstream
        for (auto i=0; i<L; ++i)
        {
          int row = ip.make_constraint(IP::LO, 0, 0);
          for (auto [j, v_ji]: v_r[lv][i])
            ip.add_constraint(row, v_ji, -1);
          if (i>0)
            for (auto [j, v_ji]: v_r[lv][i-1])
              ip.add_constraint(row, v_ji, 1);
          if (i+1<L)
            for (auto [j, v_ji]: v_r[lv][i+1])
              ip.add_constraint(row, v_ji, 1);
        }

        // downstream
        for (auto i=0; i<L; ++i)
        {
          auto row = ip.make_constraint(IP::LO, 0, 0);
          for (auto [j, v_ij]: v_l[lv][i])
            ip.add_constraint(row, v_ij, -1);
          if (i>0)
            for (auto [j, v_ij]: v_l[lv][i-1])
              ip.add_constraint(row, v_ij, 1);
          if (i+1<L)
            for (auto [j, v_ij]: v_l[lv][i+1])
              ip.add_constraint(row, v_ij, 1);
        }
      }
    }

    // Helper function to get base pair type (used by both bp and stacking constraints)
    auto get_bp_type = [](char a, char b) -> std::string {
      a = std::toupper(a);
      b = std::toupper(b);
      if (a > b) std::swap(a, b);
      if ((a == 'G' && b == 'C') || (a == 'C' && b == 'G')) return "GC";
      if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A')) return "AU";
      if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G')) return "GU";
      if (a == 'U' && b == 'U') return "UU";
      if ((a == 'C' && b == 'U') || (a == 'U' && b == 'C')) return "CU";
      std::string result;
      result += a;
      result += b;
      return result;
    };

    // Helper to check if a base pair type is canonical
    auto is_canonical = [](const std::string& bp_type) -> bool {
      return bp_type == "GC" || bp_type == "AU" || bp_type == "GU";
    };

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
              std::string bp_type = get_bp_type(seq[i], seq[j]);
              bp_type_vars[bp_type].push_back(v_ij);
            }
          }
        }
      }

      // Add constraints for GC, AU, GU, UU
      if (bp_constraints.GC >= 0) {
        auto it = bp_type_vars.find("GC");
        if (it != bp_type_vars.end() && !it->second.empty()) {
          int row = ip.make_constraint(IP::FX, bp_constraints.GC, bp_constraints.GC);
          for (int var : it->second) {
            ip.add_constraint(row, var, 1);
          }
        }
      }

      if (bp_constraints.AU >= 0) {
        auto it = bp_type_vars.find("AU");
        if (it != bp_type_vars.end() && !it->second.empty()) {
          int row = ip.make_constraint(IP::FX, bp_constraints.AU, bp_constraints.AU);
          for (int var : it->second) {
            ip.add_constraint(row, var, 1);
          }
        }
      }

      if (bp_constraints.GU >= 0) {
        auto it = bp_type_vars.find("GU");
        if (it != bp_type_vars.end() && !it->second.empty()) {
          int row = ip.make_constraint(IP::FX, bp_constraints.GU, bp_constraints.GU);
          for (int var : it->second) {
            ip.add_constraint(row, var, 1);
          }
        }
      }

      if (bp_constraints.UU >= 0) {
        auto it = bp_type_vars.find("UU");
        if (it != bp_type_vars.end() && !it->second.empty()) {
          int row = ip.make_constraint(IP::FX, bp_constraints.UU, bp_constraints.UU);
          for (int var : it->second) {
            ip.add_constraint(row, var, 1);
          }
        }
      }

      // Add constraints for other non-canonical base pair types
      for (const auto& [bp_type, count] : bp_constraints.other_bp_types) {
        if (count >= 0) {
          auto it = bp_type_vars.find(bp_type);
          if (it != bp_type_vars.end() && !it->second.empty()) {
            int row = ip.make_constraint(IP::FX, count, count);
            for (int var : it->second) {
              ip.add_constraint(row, var, 1);
            }
            spdlog::info("Added constraint for base pair type {}: {} pairs", bp_type, count);
          }
        }
      }

    }

    // Non-canonical base pairs do not require any special stacking constraints
    // They can be used freely just like canonical base pairs

    // Add stack constraints if specified
    if (stack_constraints.has_constraints())
    {
      // For each stack constraint, at least one instance must be selected
      for (size_t constraint_id = 0; constraint_id < stack_constraints.constraints.size(); ++constraint_id)
      {
        std::vector<int> instance_vars;  // Variables representing each instance

        // For each instance of this constraint
        for (size_t inst_idx = 0; inst_idx < stack_constraints.instances.size(); ++inst_idx)
        {
          const auto& instance = stack_constraints.instances[inst_idx];
          if (instance.constraint_id != (int)constraint_id) continue;

          // Find the IP variables corresponding to the base pairs in this instance
          std::vector<int> bp_vars;
          bool all_pairs_found = true;
          for (const auto& [i, j] : instance.pairs)
          {
            // Find the variable v_ij in v_l
            int found_var = -1;
            for (auto lv = 0; lv != pk_level_; ++lv)
            {
              for (const auto& [jj, v_ij] : v_l[lv][i])
              {
                if (jj == j)
                {
                  found_var = v_ij;
                  break;
                }
              }
              if (found_var >= 0) break;
            }

            if (found_var >= 0)
            {
              bp_vars.push_back(found_var);
            }
            else
            {
              //spdlog::warn("Stack constraint: base pair ({},{}) not found in variables, skipping this instance", i+1, j+1);
              all_pairs_found = false;
              break;
            }
          }

          // Skip this instance if any base pair variable was not found
          if (!all_pairs_found)
          {
            continue;
          }

          // Create a binary variable for this instance
          // This variable is 1 if all base pairs in the stack are selected
          int instance_var = ip.make_variable(0.0);  // No weight, just a helper variable
          instance_vars.push_back(instance_var);

          // Add constraints: instance_var = 1 if and only if all bp_vars = 1
          // This is implemented as:
          // 1. instance_var <= bp_var_k for all k (if any bp is 0, instance must be 0)
          // 2. instance_var >= sum(bp_vars) - n + 1 (if all bp are 1, instance must be 1)
          const int n_pairs = bp_vars.size();
          if (n_pairs > 0)
          {
            // Constraint 1: instance_var <= each bp_var
            for (int bp_var : bp_vars)
            {
              int row = ip.make_constraint(IP::LO, 0, 0);  // instance_var - bp_var <= 0
              ip.add_constraint(row, instance_var, 1);
              ip.add_constraint(row, bp_var, -1);
            }

            // Constraint 2: instance_var >= sum(bp_vars) - n + 1
            // Equivalent to: sum(bp_vars) - instance_var <= n - 1
            int row = ip.make_constraint(IP::UP, 0, n_pairs - 1);
            for (int bp_var : bp_vars)
            {
              ip.add_constraint(row, bp_var, 1);
            }
            ip.add_constraint(row, instance_var, -1);
          }
        }

        // At least one instance of this constraint must be selected
        if (!instance_vars.empty())
        {
          int row = ip.make_constraint(IP::LO, 1, 0);  // sum(instance_vars) >= 1
          for (int inst_var : instance_vars)
          {
            ip.add_constraint(row, inst_var, 1);
          }
          spdlog::info("Added stack constraint {} with {} possible instances",
                       constraint_id + 1, instance_vars.size());
        }
        else
        {
          // No valid instances found for this constraint - cannot satisfy
          spdlog::error("Stack constraint {} has no valid instances - cannot satisfy constraint",
                      constraint_id + 1);
          throw std::runtime_error("Stack constraint cannot be satisfied: no valid instances found");
        }
      }
    }

    // execute optimization
    ip.solve();

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
    spdlog::info("Search for the best thresholds by pseudo expected F-value:");
    const auto sump = compute_sump_pk(bp);
    do {
      ep.get(th);
      uint i;
      for (i=1; i!=th.size(); i++)
        if (th[i-1]<th[i]) break;
      if (i!=th.size()) continue;
      bpseq_temp = bpseq;
      plevel_temp = plevel;
      solve(seq, bp, th, bpseq_temp, plevel_temp, constraint, bp_constraints, stack_constraints);
      const auto [sen, ppv, mcc, fval] = compute_expected_accuracy(bpseq_temp, bp);
      const auto [sen_pk, ppv_pk, mcc_pk, fval_pk] = compute_expected_accuracy_pk(bpseq_temp, bp, sump);
      if (spdlog::get_level() <= spdlog::level::info)
      {
        std::ostringstream th_ss;
        std::copy(th.begin(), th.end(), std::ostream_iterator<float>(th_ss, ","));
        spdlog::info("th={} pF={}, pF_pk={}", th_ss.str(), fval, fval_pk);
      }
      if (fval+fval_pk>max_fval+max_fval_pk)
      {
        max_fval = fval;
        max_fval_pk = fval_pk;
        max_bpseq = bpseq_temp;
        max_plevel = plevel_temp;
      }
    } while (!ep.succ());
    bpseq = max_bpseq;
    plevel = max_plevel;
    spdlog::info("max pF={}, pF_pk={}", max_fval, max_fval_pk);

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
