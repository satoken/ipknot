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
#include "ipknot/decoder.h"
#include "ilp_model.h"
#include "decoder_internal.h"
#include "ip.h"
#include "dd_constrained.h"
#include "dd_constraint_builder.h"
#include "bpseq.h"

#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
//#include "spdlog/sinks/basic_file_sink.h"

IPknot::IPknot(uint pk_level, const float* alpha,
         bool levelwise, bool stacking_constraints, int n_th,
         DDOptions dd_options)
    : pk_level_(pk_level),
      alpha_(alpha, alpha+pk_level_),
      levelwise_(levelwise),
      stacking_constraints_(stacking_constraints),
      n_th_(n_th),
      dd_options_(dd_options)
{
    dd_options_.validate();
    if (dd_options_.enabled && !levelwise_)
      throw std::invalid_argument("DD requires levelwise prediction; remove --no-levelwise or use --decoder ilp");
}

void IPknot::solve(const std::string& seq, const VSVF& bp,
             const VF& th, VI& bpseq, VI& plevel, bool constraint) const
{
    if (levelwise_ && !constraint) {
      ipknot::DecoderOptions options;
      options.backend = dd_options_.enabled ? ipknot::Backend::DD : ipknot::Backend::ILP;
      options.thresholds = th;
      options.weights = alpha_;
      options.no_lonely_pairs = stacking_constraints_;
      options.threads = n_th_;
      options.dd = dd_options_;
      if (bp.size() != seq.size() + 1)
        throw std::invalid_argument("Invalid sparse posterior dimensions");
      // Refinement accumulates posterior contributions and can exceed one.
      spdlog::stopwatch decoder_timer;
      const auto decoded = ipknot::detail::decode_posterior_scores(bp, std::move(options));
      bpseq = decoded.bpseq; plevel = decoded.levels;
      if (!decoded.dd) {
        spdlog::debug("IP objective: {:.12g}", decoded.objective);
        spdlog::info("IP optimization finished in {:.6f}s", decoder_timer.elapsed().count());
        return;
      }
      const auto& result = *decoded.dd;
      spdlog::info("DD: iterations={}, pairs={}, support_rows={}, contacts={}, "
                   "crossing_beam_drops={}, witness_drops={}, objective={:.12g}, "
                   "graph_upper_bound={:.12g}, stop={}",
                   result.iterations, result.pairs, result.support_rows, result.contacts,
                   result.crossing_beam_drops, result.witness_drops,
                   result.objective, result.upper_bound, result.stop_reason);
      return;
    }
    if (dd_options_.enabled && dd_options_.linear_constraints) {
      decode_linear_constraints(seq, bp, th, alpha_, pk_level_, stacking_constraints_,
                                dd_options_, bpseq, plevel, constraint);
      return;
    }
    uint L = seq.size();
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

    // Make objective variables, including explicitly fixed pairs.
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
        }
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
        spdlog::info("DD: iterations={}, pairs={}, constraint_rows={}, objective={:.12g}, "
                     "model_upper_bound={:.12g}, repair_states={}, repair_budget_exhausted={}, stop={}",
                     result.iterations, pairs.size(), dd_model.rows.size(), result.objective,
                     result.upper_bound, result.repair_states, result.repair_budget_exhausted, result.stop_reason);
        return result.objective;
      };
    }
    ip.update();

    // Fixed constraints must also be checked when no posterior pair remains.
    if (n > 0 || constraint)
      solve(seq, ip, v_l, v_r, c_l, c_r, th, bpseq, plevel, constraint);
    else
    {
      bpseq.resize(L);
      std::fill(std::begin(bpseq), std::end(bpseq), -1);
      plevel.resize(L);
      std::fill(std::begin(plevel), std::end(plevel), -1);
      return;
    }
  }

void IPknot::solve(const std::string& seq, IP& ip, const VVSVI& v_l, const VVSVI& v_r, const VI& c_l, const VI& c_r,
             const VF& th, VI& bpseq, VI& plevel, bool constraint) const
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

    if (levelwise_) ipknot::detail::add_level_constraints(ip, v_l);

    if (stacking_constraints_)
      ipknot::detail::add_stacking_constraints(ip, v_l, v_r);

    ip.update();
    spdlog::stopwatch solver_timer;
    const double objective = ip.solve();
    spdlog::debug("{} objective: {:.12g}", dd_options_.enabled ? "DD" : "IP", objective);
    spdlog::info("{} optimization finished in {:.6f}s", dd_options_.enabled ? "DD" : "IP", solver_timer.elapsed().count());

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
             EnumParam<float>& ep, VI& bpseq, VI& plevel, bool constraint) const -> std::pair<float,float>
{
    uint L = seq.size();
    std::vector<float> th(ep.size());
    VI bpseq_temp, plevel_temp;
    VI max_bpseq, max_plevel;
    float max_fval=-100.0, max_fval_pk=-100.0;
    double max_selection_score = -std::numeric_limits<double>::infinity();
    spdlog::info("Search for the best thresholds by pseudo expected F-value:");
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
      try {
        solve(seq, bp, th, bpseq_temp, plevel_temp, constraint);
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
        spdlog::info("th={} pF={}, pF_pk={}", th_ss.str(), fval, fval_pk);
      }
      const double selection_score = fval + fval_pk;
      if (selection_score > max_selection_score)
      {
        max_selection_score = selection_score;
        max_fval = fval;
        max_fval_pk = fval_pk;
        max_bpseq = bpseq_temp;
        max_plevel = plevel_temp;
      }
    } while (!ep.succ());
    if (dd_options_.enabled && !std::isfinite(max_selection_score))
      throw DDInfeasible("No feasible structure for any DD threshold combination");
    bpseq = max_bpseq;
    plevel = max_plevel;
    spdlog::info("max pF={}, pF_pk={}, selection score={}",
                 max_fval, max_fval_pk, max_selection_score);

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
