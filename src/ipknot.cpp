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
#include <unistd.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
//#include <sys/resource.h>
#include <strings.h>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <memory>

#include "ip.h"
#include "fa.h"
#include "aln.h"
#include "fold.h"
#include "nupack/nupack.h"
#include "bpseq.h"
#include "cxxopts.hpp"

using uint = unsigned int;
using VF = std::vector<float>;
using VVF = std::vector<VF>;
using VI = std::vector<int>;
using VVI = std::vector<VI>;
using VVVI = std::vector<VVI>;
using SVI = std::vector<std::pair<uint,int>>;
using VSVI = std::vector<SVI>; 
using VVSVI = std::vector<VSVI>;

class IPknot
{
public:
  template < class T > class EnumParam;

public:
  IPknot(uint pk_level, const float* alpha,
         bool levelwise, bool stacking_constraints, int n_th)
    : pk_level_(pk_level),
      alpha_(alpha, alpha+pk_level_),
      levelwise_(levelwise),
      stacking_constraints_(stacking_constraints),
      n_th_(n_th)
  {
  }

public:
  void solve(uint L, const std::vector<float>& bp, const std::vector<int>& offset,
             const std::vector<float>& th, std::vector<int>& bpseq, std::vector<int>& plevel, bool constraint) const
  {
    IP ip(IP::MAX, n_th_);
    VVSVI v_l(pk_level_, VSVI(L));
    VVSVI v_r(pk_level_, VSVI(L));
    VI c_l(L, 0), c_r(L, 0);
    uint n=0;

    // make objective variables with their weights
    for (uint j=1; j!=L; ++j)
    {
      for (uint i=j-1; i!=-1u; --i)
      {
        const float& p=bp[offset[i+1]+(j+1)];
        for (uint lv=0; lv!=pk_level_; ++lv)
          if (p>th[lv])
          {
            const auto v_ij = ip.make_variable(p*alpha_[lv]);
            v_l[lv][i].emplace_back(j, v_ij);
            v_r[lv][j].emplace_back(i, v_ij);
            c_l[i]++; c_r[j]++;
            n++;
          }
      }
    }
    ip.update();

    if (n>0)
      solve(L, ip, v_l, v_r, c_l, c_r, th, bpseq, plevel, constraint);
    else
    {
      bpseq.resize(L);
      std::fill(std::begin(bpseq), std::end(bpseq), -1);
      plevel.resize(L);
      std::fill(std::begin(plevel), std::end(plevel), -1);
    }
  }

  void solve(uint L, const std::vector<std::vector<std::pair<uint, float>>>& bp,
             const std::vector<float>& th, std::vector<int>& bpseq, std::vector<int>& plevel, bool constraint) const
  {
    IP ip(IP::MAX, n_th_);
    VVSVI v_l(pk_level_, VSVI(L));
    VVSVI v_r(pk_level_, VSVI(L));
    VI c_l(L, 0), c_r(L, 0);
    uint n=0;

    // make objective variables with their weights
    for (uint i=1; i<=L; ++i)
    {
      for (const auto [j, p]: bp[i])
        if (i<j)
          for (uint lv=0; lv!=pk_level_; ++lv)
            if (p>th[lv])
            {
              const auto v_ij = ip.make_variable(p*alpha_[lv]);
              v_l[lv][i-1].emplace_back(j-1, v_ij);
              v_r[lv][j-1].emplace_back(i-1, v_ij);
              c_l[i-1]++; c_r[j-1]++;
              n++;
            }
    }
    ip.update();

    if (n>0)
      solve(L, ip, v_l, v_r, c_l, c_r, th, bpseq, plevel, constraint);
    else
    {
      bpseq.resize(L);
      std::fill(std::begin(bpseq), std::end(bpseq), -1);
      plevel.resize(L);
      std::fill(std::begin(plevel), std::end(plevel), -1);
    }
  }

private:
  void solve(uint L, IP& ip, const VVSVI& v_l, const VVSVI& v_r, const VI& c_l, const VI& c_r,
             const std::vector<float>& th, std::vector<int>& bpseq, std::vector<int>& plevel, bool constraint) const
  {
    if (!constraint)
    {
      bpseq.resize(L);
      std::fill(bpseq.begin(), bpseq.end(), -2);
    }
    
    // constraint 1: each s_i is paired with at most one base
    for (uint i=0; i!=L; ++i)
    {
      int row_l = -1, row_r = -1;
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
        std::cerr << "invalid constraint for the base " << i+1 << ", ignored." << std::endl;
        row_l = row_r = ip.make_constraint(IP::UP, 0, 1); // fallback to no constraint
      }
      
      for (uint lv=0; lv!=pk_level_; ++lv)
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
        for (uint lv=0; lv!=pk_level_; ++lv)
        {
          for (const auto [temp, v_ij]: v_l[lv][i])
            if (j==temp) 
            { 
              vals[lv] = v_ij; 
              c++;
              break; 
            }
        }
        if (c>0)
        {
          int row = ip.make_constraint(IP::FX, 1, 1);
          for (uint lv=0; lv!=pk_level_; ++lv)
            if (vals[lv]>=0)
              ip.add_constraint(row, vals[lv], 1);
        }
        else
          std::cerr << "invalid constraint for the bases " << i+1 << " and " << bpseq[i]+1 << ", ignored." << std::endl;
      }
    }

    if (levelwise_)
    {
      // constraint 2: disallow pseudoknots in x[lv]
      for (uint lv=0; lv!=pk_level_; ++lv)
        for (uint i=0; i<v_l[lv].size(); ++i)
          for (const auto [j, v_ij]: v_l[lv][i])
            for (uint k=i+1; k<j; ++k)
              for (const auto [l, v_kl]: v_l[lv][k])
                if (j<l)
                {
                  int row = ip.make_constraint(IP::UP, 0, 1);
                  ip.add_constraint(row, v_ij, 1);
                  ip.add_constraint(row, v_kl, 1);
                }

      // constraint 3: any x[t]_kl must be pseudoknotted with x[u]_ij for t>u
      for (uint lv=1; lv!=pk_level_; ++lv)
        for (uint k=0; k<v_l[lv].size(); ++k)
          for (const auto [l, v_kl]: v_l[lv][k])
            for (uint plv=0; plv!=lv; ++plv)
            {
              int row = ip.make_constraint(IP::LO, 0, 0);
              ip.add_constraint(row, v_kl, -1);
              for (uint i=0; i<k; ++i)
                for (const auto [j, v_ij]: v_l[plv][i])
                  if (k<j && j<l)
                    ip.add_constraint(row, v_ij, 1);

              for (uint i=k+1; i<l; ++i)
                for (const auto [j, v_ij]: v_l[plv][i])
                  if (l<j)
                    ip.add_constraint(row, v_ij, 1);
            }
    }

    if (stacking_constraints_)
    {
      for (uint lv=0; lv!=pk_level_; ++lv)
      {
        // upstream
        for (uint i=0; i<L; ++i)
        {
          int row = ip.make_constraint(IP::LO, 0, 0);
          for (const auto [j, v_ji]: v_r[lv][i])
            ip.add_constraint(row, v_ji, -1);
          if (i>0)
            for (const auto [j, v_ji]: v_r[lv][i-1])
              ip.add_constraint(row, v_ji, 1);
          if (i+1<L)
            for (const auto [j, v_ji]: v_r[lv][i+1])
              ip.add_constraint(row, v_ji, 1);
        }

        // downstream
        for (uint i=0; i<L; ++i)
        {
          int row = ip.make_constraint(IP::LO, 0, 0);
          for (const auto [j, v_ij]: v_l[lv][i])
            ip.add_constraint(row, v_ij, -1);
          if (i>0)
            for (const auto [j, v_ij]: v_l[lv][i-1])
              ip.add_constraint(row, v_ij, 1);
          if (i+1<L)
            for (const auto [j, v_ij]: v_l[lv][i+1])
              ip.add_constraint(row, v_ij, 1);
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
    for (uint lv=0; lv!=pk_level_; ++lv)
      for (uint i=0; i<L; ++i)
        for (const auto [j, v_ij]: v_l[lv][i])
          if (ip.get_value(v_ij)>0.5)
          {
            bpseq[i]=j; bpseq[j]=i;
            plevel[i]=plevel[j]=lv;
          }

    if (!levelwise_)
      decompose_plevel(bpseq, plevel);
  }

public:
  void solve(uint L, const std::vector<float>& bp, const std::vector<int>& offset,
             EnumParam<float>& ep, std::vector<int>& bpseq, std::vector<int>& plevel, bool constraint) const
  {
    std::vector<float> th(ep.size());
    std::vector<int> bpseq_temp;
    std::vector<int> plevel_temp;
    float max_fval=-100.0;
    std::cerr << "Search for the best thresholds by pseudo expected F-value:" << std::endl;
    do {
      ep.get(th);
      std::cerr << "th=";
      std::copy(th.begin(), th.end(), std::ostream_iterator<float>(std::cerr, ","));
      solve(L, bp, offset, th, bpseq_temp, plevel_temp, constraint);
      const auto [sen, ppv, mcc, fval] = compute_expected_accuracy(bpseq_temp, bp, offset);
      std::cerr << " pF=" << fval << std::endl;
      if (fval>max_fval)
      {
        max_fval = fval;
        bpseq = bpseq_temp;
        plevel = plevel_temp;
      }
    } while (!ep.succ());
    std::cerr << "max pF=" << max_fval << std::endl << std::endl;
  }

  void solve(uint L, const std::vector<std::vector<std::pair<uint, float>>>& bp,
             EnumParam<float>& ep, std::vector<int>& bpseq, std::vector<int>& plevel, bool constraint) const
  {
    std::vector<float> th(ep.size());
    std::vector<int> bpseq_temp;
    std::vector<int> plevel_temp;
    float max_fval=-100.0;
    std::cerr << "Search for the best thresholds by pseudo expected F-value:" << std::endl;
    do {
      ep.get(th);
      std::cerr << "th=";
      std::copy(th.begin(), th.end(), std::ostream_iterator<float>(std::cerr, ","));
      solve(L, bp, th, bpseq_temp, plevel_temp, constraint);
      const auto [sen, ppv, mcc, fval] = compute_expected_accuracy(bpseq_temp, bp);
      std::cerr << " pF=" << mcc << std::endl;
      if (fval>max_fval)
      {
        max_fval = fval;
        bpseq = bpseq_temp;
        plevel = plevel_temp;
      }
    } while (!ep.succ());
    std::cerr << "max pF=" << max_fval << std::endl << std::endl;
  }

public:
  template < class T >
  class EnumParam
  {
  public:
    EnumParam(const std::vector<std::vector<T> >& p)
      : p_(p), m_(p.size()), v_(p.size(), 0)
    {
      for (uint i=0; i!=p.size(); ++i)
        m_[i] = p[i].size();
    }

    uint size() const { return m_.size(); }
    
    void get(std::vector<T>& q) const
    {
      for (uint i=0; i!=v_.size(); ++i)
        q[i] = p_[i][v_[i]];
    }
  
    bool succ()
    {
      return succ(m_.size(), &m_[0], &v_[0]);
    }

  private:
    static bool succ(int n, const int* m, int* v)
    {
      if (n==0) return true;
      if (++(*v)==*m)
      {
        *v=0;
        return succ(n-1, ++m, ++v);
      }
      return false;
    }

  private:
    const std::vector<std::vector<T> >& p_;
    std::vector<int> m_;
    std::vector<int> v_;
  };

private:
  struct cmp_by_degree : public std::less<int>
  {
    cmp_by_degree(const std::vector< std::vector<int> >& g) : g_(g) {}
    bool operator()(int x, int y) const { return g_[y].size()<g_[x].size(); }
    const std::vector< std::vector<int> >& g_;
  };

  struct cmp_by_count : public std::less<int>
  {
    cmp_by_count(const std::vector<int>& count) : count_(count) { }
    bool operator()(int x, int y) const { return count_[y]<count_[x]; }
    const std::vector<int>& count_;
  };

public:
  static int
  decompose_plevel(const std::vector<int>& bpseq, std::vector<int>& plevel)
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
    std::sort(v.begin(), v.end(), cmp_by_degree(g));

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
    sort(idx.begin(), idx.end(), cmp_by_count(count));
    std::vector<int> rev(idx.size());
    for (uint i=0; i!=rev.size(); ++i) rev[idx[i]]=i;
    plevel.resize(L);
    for (uint i=0; i!=c.size(); ++i)
      plevel[i]= c[i]>=0 ? rev[c[i]] : -1;

    return max_color+1;
  }

private:
  static auto
  compute_expected_accuracy(const std::vector<int>& bpseq,
                            const std::vector<float>& bp, const std::vector<int>& offset) -> std::tuple<float,float,float,float>
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

    float sen, ppv, mcc, f;
    sen = ppv = mcc = f = 0;
    if (etp+efn!=0) sen = etp / (etp + efn);
    if (etp+efp!=0) ppv = etp / (etp + efp);
    if (etp+efp!=0 && etp+efn!=0 && etn+efp!=0 && etn+efn!=0)
      mcc = (etp*etn-efp*efn) / std::sqrt((etp+efp)*(etp+efn)*(etn+efp)*(etn+efn));
    if (sen+ppv!=0) f = 2*sen*ppv/(sen+ppv);

    return {sen, ppv, mcc, f};
  }

  static auto
  compute_expected_accuracy(const std::vector<int>& bpseq,
                            const std::vector<std::vector<std::pair<uint, float>>>& bp) -> std::tuple<float,float,float,float>
  {
    int L  = bpseq.size();
    int L2 = L*(L-1)/2;
    int N = 0;

    float sump = 0.0;
    float etp  = 0.0;

    for (uint i=1; i!=bp.size(); ++i) 
      for (const auto [j, p]: bp[i])
        if (i<j) sump += p;

    for (uint i=0; i!=bpseq.size(); ++i)
    {
      if (bpseq[i]!=-1 && bpseq[i]>(int)i)
      {
        const auto j=bpseq[i];
        auto res = std::find_if(std::cbegin(bp[i+1]), std::cend(bp[i+1]),
                          [&](const auto& x) { return x.first==j+1; });
        if (res!=std::cend(bp[i+1])) {
          etp += res->second;
          N++;
        }
      }
    }

    float etn = L2 - N - sump + etp;
    float efp = N - etp;
    float efn = sump - etp;

    float sen, ppv, mcc, f;
    sen = ppv = mcc = f = 0;
    if (etp+efn!=0) sen = etp / (etp + efn);
    if (etp+efp!=0) ppv = etp / (etp + efp);
    if (etp+efp!=0 && etp+efn!=0 && etn+efp!=0 && etn+efn!=0)
      mcc = (etp*etn-efp*efn) / std::sqrt((etp+efp)*(etp+efn)*(etn+efp)*(etn+efn));
    if (sen+ppv!=0) f = 2*sen*ppv/(sen+ppv);

    return {sen, ppv, mcc, f};
  }

  static uint length(const std::string& seq) { return seq.size(); }
  static uint length(const std::list<std::string>& aln) { return aln.front().size(); }

private:
  // options
  uint pk_level_;
  std::vector<float> alpha_;
  bool levelwise_;
  bool stacking_constraints_;
  int n_th_;
};

std::string
make_parenthsis(const std::vector<int>& bpseq, const std::vector<int>& plevel)
{
  const int n_support_parens=4;
  const char* left_paren="([{<";
  const char* right_paren=")]}>";

  std::string r(bpseq.size(), '.');
  for (int i=0; i!=(int)bpseq.size(); ++i)
  {
    if (bpseq[i]>=0 && i<bpseq[i])
    {
      int j=bpseq[i];
      if (plevel[i]<n_support_parens)
      {
        r[i]=left_paren[plevel[i]];
        r[j]=right_paren[plevel[i]];
      }
      else if (plevel[i]<n_support_parens+'Z'-'A'+1)
      {
        r[i]='A'+plevel[i]-n_support_parens;
        r[j]='a'+plevel[i]-n_support_parens;
      }
    }
  }
  return r;
}

#if 0
double
timing()
{
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return ru.ru_utime.tv_sec+ru.ru_utime.tv_usec*1e-6;
}
#endif

template < class SEQ, class EN >
void
update_bpm(uint pk_level, const SEQ& seq, EN& en,
           const std::vector<int>& bpseq, const std::vector<int>& plevel,
           std::vector<float>& bp, std::vector<int>& offset)
{
  // update the base-pairing probability matrix by the previous result
  uint L=bpseq.size();
  bp.resize((L+1)*(L+2)/2, 0.0);
  std::fill(bp.begin(), bp.end(), 0.0);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
  
  std::vector<float> bpl;
  std::vector<int> offsetl;
  for (uint l=0; l!=pk_level; ++l)
  {
    // make the constraint string
    std::string str(L, '?');
    for (uint i=0; i!=bpseq.size(); ++i)
    {
      switch (bpseq[i])
      {
        case BPSEQ::U: str[i] = '.'; break;
        case BPSEQ::L: str[i] = '<'; break;
        case BPSEQ::R: str[i] = '>'; break;
        case BPSEQ::LR: str[i] = '|'; break;
        case BPSEQ::DOT: str[i] = '?'; break;
        default:
          if (bpseq[i]>=0 && (int)i<bpseq[i])
          {
            if ((int)l==plevel[i])
            {
              str[i]='('; str[bpseq[i]]=')';
            }
            else
            {
              str[i]=str[bpseq[i]]='.';
            }
          }
          break;
      }
    }

    // re-folding the seq with the constraint
    std::fill(bpl.begin(), bpl.end(), 0.0);
    en.calculate_posterior(seq, str, bpl, offsetl);
    assert(bp.size()==bpl.size());
    // update the base-pairing probability matrix
#if 0 // original behaivior
    for (uint k=0; k!=bp.size(); ++k) bp[k]+=bpl[k];
#else
    for (uint j=1; j!=L; ++j)
      for (uint i=j-1; i!=-1u; --i)
        if (bpseq[i]>=0)
          bp[offset[i+1]+(j+1)] += bpl[offset[i+1]+(j+1)];
        else
          bp[offset[i+1]+(j+1)] += bpl[offset[i+1]+(j+1)] / pk_level;
#endif
  }
#ifndef NDEBUG
  for (uint k=0; k!=bp.size(); ++k) assert(bp[k]<=1.0);
#endif
}

template < class SEQ, class EN>
void
update_bpm(uint pk_level, const SEQ& seq, EN& en,
           const std::vector<int>& bpseq, const std::vector<int>& plevel,
           std::vector<std::vector<std::pair<uint, float>>>& sbp)
{
  // update the base-pairing probability matrix by the previous result
  uint L=bpseq.size();
  sbp.resize(L+1);
  
  for (uint l=0; l!=pk_level; ++l)
  {
    // make the constraint string
    std::string str(L, '?');
    for (uint i=0; i!=bpseq.size(); ++i)
    {
      switch (bpseq[i])
      {
        case BPSEQ::U: str[i] = '.'; break;
        case BPSEQ::L: str[i] = '<'; break;
        case BPSEQ::R: str[i] = '>'; break;
        case BPSEQ::LR: str[i] = '|'; break;
        case BPSEQ::DOT: str[i] = '?'; break;
        default:
          if (bpseq[i]>=0 && (int)i<bpseq[i])
          {
            if ((int)l==plevel[i])
            {
              str[i]='('; str[bpseq[i]]=')';
            }
            else
            {
              str[i]=str[bpseq[i]]='.';
            }
          }
          break;
      }
    }

    // re-folding the seq with the constraint
    auto sbpl = en.calculate_posterior(seq, str);
    assert(sbp.size()==sbpl.size());

    // update the base-pairing probability matrix
    for (uint i=1; i!=sbpl.size(); i++) 
    {
      for (const auto [jl, vl]: sbpl[i]) 
      {
        auto v = bpseq[i-1]>=0 ? vl : vl / pk_level;
        auto re = std::find_if(std::begin(sbp[i]), std::end(sbp[i]),
                            [&](const auto& x) { return x.first == jl; });
        if (re != std::end(sbp[i]))
          re->second += v;
        else
          sbp[i].emplace_back(jl, v);
      }
    }
  }
}

static
void
read_constraints(const char* filename, VI& bpseq)
{
  std::ifstream is(filename);
  std::string s;
  int i;
  while (is >> i >> s)
  {
    if (i<=0 && i>bpseq.size())
      std::cerr << "invalid format base number i=" << i << ", ignored." << std::endl;
    else switch (s[0]) {
      default: 
        if (std::isdigit(s[0]))
        {
          int j = std::atoi(s.c_str());
          if (j>0 && j<=bpseq.size()) {
            bpseq[i-1] = j-1; bpseq[j-1] = i-1;
          } 
          else
            std::cerr << "invalid format base number j=" << j << ", ignored." << std::endl;
        }
        break;
      case '.': bpseq[i-1] = BPSEQ::DOT; break;
      case 'x': case 'X': bpseq[i-1] = BPSEQ::U; break;
      case '|': bpseq[i-1] = BPSEQ::LR; break;
      case '<': bpseq[i-1] = BPSEQ::L; break;
      case '>': bpseq[i-1] = BPSEQ::R; break;
    }
    
  }
}

static
void
output_fa(std::ostream& os,
          const std::string& desc, const std::string& seq,
          const std::vector<int>& bpseq, const std::vector<int>& plevel, 
          bool output_energy)
{
  if (output_energy)
  {
    Nupack<float> nupack;
    nupack.load_default_parameters();
    nupack.load_sequence(seq);
    nupack.load_constraints(bpseq);
    long double e = nupack.calculate_minimum_free_energy();
    os << ">" << desc << " (e=" << e << ")" << std::endl;
  }
  else
  {
    os << ">" << desc << std::endl; 
  }
  os << seq << std::endl
     << make_parenthsis(bpseq, plevel) << std::endl;
}

static
void
output_mfa(std::ostream& os, const Aln& aln,
           const std::vector<int>& bpseq, const std::vector<int>& plevel)
{
  os << ">SS_cons" << std::endl
     << make_parenthsis(bpseq, plevel) << std::endl;
  std::list<std::string>::const_iterator name=aln.name().begin();
  std::list<std::string>::const_iterator seq=aln.seq().begin();
  while (name!=aln.name().end() && seq!=aln.seq().end())
  {
    os << ">" << *name << std::endl
       << *seq << std::endl;
    ++seq; ++name;
  }
}

static
void
output_bpseq(std::ostream& os,
             const std::string& desc, const std::string& seq,
             const std::vector<int>& bpseq, const std::vector<int>& plevel)
{
  os << "# " << desc << std::endl;
  for (uint i=0; i!=bpseq.size(); ++i)
    os << i+1 << " " << seq[i] << " " << bpseq[i]+1 << std::endl;
}

template <class T>
std::vector<T>
parse_csv_line(const char* l, char delim=',')
{
  std::string s;
  std::vector<T> r;
  std::istringstream ss(l);
  while (std::getline(ss, s, delim))
    r.push_back(atof(s.c_str()));
  return r;
}

auto
build_engine_seq(const char* model, const char* param, uint beam_size=100)
{
  std::unique_ptr<BPEngineSeq> en;
  if (model==nullptr || strcasecmp(model, "McCaskill")==0)
    en = std::make_unique<RNAfoldModel>(param);
  else if (strcasecmp(model, "CONTRAfold")==0)
    en = std::make_unique<CONTRAfoldModel>();
#if 0
  else if (strcasecmp(model, "nupack")==0)
    if (param)
      en = std::make_unique<NupackModel>(param);
    else
      en = std::make_unique<NupackModel>(2);
  else if (strcasecmp(model, "nupack03")==0)
    en = std::make_unique<NupackModel>(0);
  else if (strcasecmp(model, "nupack09")==0)
    en = std::make_unique<NupackModel>(1);
#else
  else if (strcasecmp(model, "nupack")==0)
    en = std::make_unique<NupackModel>(param);
#endif
  else if (strcasecmp(model, "LinearPartition")==0 || strcasecmp(model, "lpc")==0)
    en = std::make_unique<LinearPartitionModel>(false, beam_size);
  else if (strcasecmp(model, "LinearPartitionV")==0 || strcasecmp(model, "lpv")==0)
    en = std::make_unique<LinearPartitionModel>(true, beam_size);
  return en;
}

auto
build_engine_aln(const std::vector<std::string>& model, const char* param, uint beam_size=100)
{
  std::unique_ptr<BPEngineAln> mix_en;
  std::vector<std::unique_ptr<BPEngineAln>> en_a;
  if (model.empty())
  {
    auto e = std::make_unique<RNAfoldModel>(param);
    en_a.push_back(std::make_unique<AveragedModel>(std::move(e)));
    en_a.push_back(std::make_unique<AlifoldModel>(param));
    mix_en = std::make_unique<MixtureModel>(std::move(en_a));
  }
  else
  {
    for (const auto mo : model) 
    {
      auto m = mo.c_str();
      if (strcasecmp(m, "McCaskill")==0)
      {
        auto e = std::make_unique<RNAfoldModel>(param);
        en_a.push_back(std::make_unique<AveragedModel>(std::move(e)));
      }
      else if (strcasecmp(m, "CONTRAfold")==0)
      {
        auto e = std::make_unique<CONTRAfoldModel>();
        en_a.push_back(std::make_unique<AveragedModel>(std::move(e)));
      }
      else if (strcasecmp(m, "Alifold")==0)
      {
        en_a.push_back(std::make_unique<AlifoldModel>(param));
      }
      else if (strcasecmp(m, "LinearPartition")==0 || strcasecmp(m, "lpc")==0)
      {
        auto e = std::make_unique<LinearPartitionModel>(false, beam_size);
        en_a.push_back(std::make_unique<AveragedModel>(std::move(e)));
      }
      else if (strcasecmp(m, "LinearPartitionV")==0 || strcasecmp(m, "lpv")==0)
      {
        auto e = std::make_unique<LinearPartitionModel>(true, beam_size);
        en_a.push_back(std::make_unique<AveragedModel>(std::move(e)));
      }
      else
        return std::unique_ptr<BPEngineAln>();
    }
    if (en_a.size()>1)
      mix_en = std::make_unique<MixtureModel>(std::move(en_a));
  }
  return std::move(mix_en ? mix_en : en_a[0]);
}

template <typename ... Args>
std::string format(const std::string& fmt, Args ... args )
{
    size_t len = std::snprintf( nullptr, 0, fmt.c_str(), args ... );
    std::vector<char> buf(len + 1);
    std::snprintf(&buf[0], len + 1, fmt.c_str(), args ... );
    return std::string(&buf[0], &buf[0] + len);
}

int
main(int argc, char* argv[])
{
  char* progname=argv[0];
  // parse options
  uint pk_level=0;
  std::vector< std::vector<float> > th;
  std::vector<float> alpha;
  bool isolated_bp=false;
  std::vector<std::string> model;
  int n_th=1;
  int n_refinement=0;
  std::string param;
  bool aux=false;
  bool levelwise=true;
  bool max_pfval=false;
  bool output_energy=false;
  std::ostream *os_bpseq=nullptr;
  std::ostream *os_mfa=nullptr;
  std::string constraint;
  uint beam_size;
  std::string input;

  cxxopts::Options options{progname, format("IPknot version %s", PACKAGE_VERSION)};
  options.add_options()
    ("input", "FASTA-formatted file or ALN-formatted file",
      cxxopts::value<std::string>(), "FASTA_OR_ALN")
    ("e,model", "Probabilistic model",
      cxxopts::value<std::vector<std::string>>()->default_value("McCaskill"), "MODEL")
    ("r,refinement", "The number of the iterative refinement",
      cxxopts::value<int>()->default_value("0"), "N")
    ("a,alpha", "The weight for each level",
      cxxopts::value<std::vector<float>>(), "ALPHA")
    // ("m,mcc", "Select thresholds that maximize pseudo MCC",
    //   cxxopts::value<bool>()->default_value("false"))
    ("pseudo-fval", "Select thresholds that maximize pseudo F-value",
      cxxopts::value<bool>()->default_value("false"))
    ("t,threshold", "The threshold of base-pairing probabilities for each level",
      cxxopts::value<std::vector<std::string>>(), "TH")
    ("g,gamma", "The weight for true base-pairs equivalent to '-t 1/(gamma+1)'",
      cxxopts::value<std::vector<std::string>>(), "G")
    ("i,allow-isolated", "Allow isolated base-pairs",
      cxxopts::value<bool>()->default_value("false"))
    ("b,bpseq", "Output the prediction by BPSEQ format",
      cxxopts::value<bool>()->default_value("false"))
    ("B,bpseq-file", "Output file for BPSEQ format",
      cxxopts::value<std::string>(), "FILE")
#ifndef WITH_GLPK
    ("n,threads", "The number of threads for the available solvers",
      cxxopts::value<uint>()->default_value("1"), "N")
#endif
    ("P,param", "Read the energy parameter file for Vienna RNA package",
      cxxopts::value<std::string>(), "FILE")
    ("x,aux", "Import an auxiliary file for base-pairing probabilities",
      cxxopts::value<bool>()->default_value("false"))
    ("u,levelwise", "Perform the levelwise prediction",
      cxxopts::value<bool>()->default_value("false"))
    ("l,mfa", "Output the prediction with the given mulple alignment",
      cxxopts::value<bool>()->default_value("false"))
    ("L,mfa-file", "Output file for the multiple alignment",
      cxxopts::value<std::string>(), "FILE")
    ("E,energy", "Output with the free energy",
      cxxopts::value<bool>()->default_value("false"))
    ("c,constraint", "Specify the structure constraint by a BPSEQ formatted file",
      cxxopts::value<std::string>(), "FILE")
    ("beam-size", "Beam size for LinearPartition algorithm",
      cxxopts::value<uint>()->default_value("100"), "N")
    ("h,help", "Print usage"); 
  options.parse_positional({"input"});
  options
    .positional_help("FASTA_OR_ALN")
    .show_positional_help();

  auto res = options.parse(argc, argv);
  if (res.count("help") || res.count("input")==0)
  {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  model = res["model"].as<std::vector<std::string>>();
  n_refinement = res["refinement"].as<int>();
  max_pfval = res["pseudo-fval"].as<bool>();
  isolated_bp = res["allow-isolated"].as<bool>();
  if (res.count("param")) param = res["param"].as<std::string>();
  aux = res["aux"].as<bool>();
  levelwise = res["levelwise"].as<bool>();
#ifndef WITH_GLPK
  n_th = res["threads"].as<uint>();
#endif
  output_energy = res["energy"].as<bool>();
  if (res.count("constraint")) constraint = res["constraint"].as<std::string>();
  beam_size = res["beam-size"].as<uint>();
  input = res["input"].as<std::string>();

  if (res.count("bpseq"))
    os_bpseq = &std::cout;
  if (res.count("bpseq-file"))
  {
    auto f = res["bpseq-file"].as<std::string>().c_str();
    os_bpseq = new std::ofstream(f);
    if (!dynamic_cast<std::ofstream*>(os_bpseq)->is_open())
    {
      perror(f);
      return 1;
    }
  }

  if (res.count("mfa"))
    os_mfa = &std::cout;
  if (res.count("mfa-file"))
  {
    auto f = res["mfa-file"].as<std::string>().c_str();
    os_mfa = new std::ofstream(f);
    if (!dynamic_cast<std::ofstream*>(os_mfa)->is_open())
    {
      perror(f);
      return 1;
    }
  }

  if (res.count("threshold"))
  {
    const auto& arg_th = res["threshold"].as<std::vector<std::string>>();
    th.resize(arg_th.size());
    for (uint i=0; i!=th.size(); ++i)
    {
      th[i] = parse_csv_line<float>(arg_th[i].c_str(), '_');
    }
  }
  else if (res.count("gamma"))
  {
    const auto& arg_gamma = res["gamma"].as<std::vector<std::string>>();
    th.resize(arg_gamma.size());
    for (uint i=0; i!=th.size(); ++i)
    {
      auto temp = parse_csv_line<float>(arg_gamma[i].c_str(), '_');
      th[i].resize(temp.size());
      std::transform(std::cbegin(temp), std::cend(temp), std::begin(th[i]),
        [&](auto v) { return 1./(v+1.); });
      // for (uint j=0; j!=temp.size(); ++j)
      //   temp[j] = 1./(temp[i]+1.);
    }
  }
  else // default
  {
    th.resize(2);
    if (n_refinement==0)
    {
      th[0].resize(1, 1/(2.0+1)); // -g 2
      th[1].resize(1, 1/(4.0+1)); // -g 4
    }
    else
    {
      th[0].resize(1, 1/(1.0+1)); // -g 1
      th[1].resize(1, 1/(1.0+1)); // -g 1
    }
  }

  if (res.count("alpha"))
  {
    alpha = res["alpha"].as<std::vector<float>>();
  }
  else
  {
    alpha.resize(th.size());
    std::fill(std::begin(alpha), std::end(alpha), 1./alpha.size());
  }
  pk_level = alpha.size();

  try
  {
    IPknot ipknot(pk_level, &alpha[0], levelwise, !isolated_bp, n_th);
    std::vector<int> bpseq;
    std::vector<int> plevel;

    IPknot::EnumParam<float> ep(th);
    std::vector<float> t(th.size());
    ep.get(t);

    std::list<Fasta> f;
    std::list<Aln> a;

    if (aux)
    {
      std::vector<float> bp;
      std::vector<int> offset;
      AuxModel aux;
      std::string seq;
      aux.calculate_posterior(input.c_str(), seq, bp, offset);
      if (max_pfval)
        ipknot.solve(seq.size(), bp, offset, ep, bpseq, plevel, false);
      else
        ipknot.solve(seq.size(), bp, offset, t, bpseq, plevel, false);
      if (os_bpseq)
        output_bpseq(*os_bpseq, input.c_str(), seq, bpseq, plevel);
      if (os_bpseq!=&std::cout)
        output_fa(std::cout, input.c_str(), seq, bpseq, plevel, output_energy);
    }
    else if (Fasta::load(f, input.c_str())>0)
    {
      auto en = build_engine_seq(model.empty() ? nullptr : model[0].c_str(), 
                                  param.empty() ? nullptr : param.c_str(), beam_size);
      if (!en) 
      {
        std::cout << options.help() << std::endl;
        return 1;
      }

      while (!f.empty())
      {
        std::list<Fasta>::iterator fa = f.begin();

        std::vector<std::vector<std::pair<uint, float>>> sbp;
        if (constraint.empty())
          sbp = en->calculate_posterior(fa->seq());
        else
        { // constraint folding
          bpseq.resize(fa->size(), BPSEQ::DOT);
          read_constraints(constraint.c_str(), bpseq);
          int pl = IPknot::decompose_plevel(bpseq, plevel);
          update_bpm(pl, fa->seq(), *en, bpseq, plevel, sbp);
        }

        if (max_pfval)
          ipknot.solve(fa->size(), sbp, ep, bpseq, plevel, !constraint.empty());
        else
          ipknot.solve(fa->size(), sbp, t, bpseq, plevel, !constraint.empty());

        for (int i=0; i!=n_refinement; ++i) // iterative refinement
        {
          update_bpm(pk_level, fa->seq(), *en, bpseq, plevel, sbp);
          if (max_pfval)
            ipknot.solve(fa->size(), sbp, ep, bpseq, plevel, !constraint.empty());
          else
            ipknot.solve(fa->size(), sbp, t, bpseq, plevel, !constraint.empty());
        }

        if (os_bpseq)
          output_bpseq(*os_bpseq, fa->name(), fa->seq(), bpseq, plevel);
        if (os_bpseq!=&std::cout)
          output_fa(std::cout, fa->name(), fa->seq(), bpseq, plevel, output_energy);
        f.erase(fa);
      }
    }
    else if (Aln::load(a, input.c_str())>0)
    {
      auto en = build_engine_aln(model, param.empty() ? nullptr : param.c_str(), beam_size);
      if (!en) 
      {
        std::cout << options.help() << std::endl;
        return 1;
      }

      while (!a.empty())
      {
        std::list<Aln>::iterator aln = a.begin();

        std::vector<std::vector<std::pair<uint, float>>> sbp;
        if (constraint.empty()) // default behaiviro
          sbp = en->calculate_posterior(aln->seq());
        else
        { // constraint folding
          bpseq.resize(aln->size(), BPSEQ::DOT);
          read_constraints(constraint.c_str(), bpseq);
          int pl = IPknot::decompose_plevel(bpseq, plevel);
          update_bpm(pl, aln->seq(), *en, bpseq, plevel, sbp);
        }
        
        if (max_pfval)
          ipknot.solve(aln->size(), sbp, ep, bpseq, plevel, !constraint.empty());
        else
          ipknot.solve(aln->size(), sbp, t, bpseq, plevel, !constraint.empty());

        for (int i=0; i!=n_refinement; ++i)
        {
          update_bpm(pk_level, aln->seq(), *en, bpseq, plevel, sbp);
          if (max_pfval)
            ipknot.solve(aln->size(), sbp, ep, bpseq, plevel, !constraint.empty());
          else
            ipknot.solve(aln->size(), sbp, t, bpseq, plevel, !constraint.empty());
        }

        if (os_bpseq)
          output_bpseq(*os_bpseq, aln->name().front(), aln->consensus(), bpseq, plevel);
        if (os_mfa)
          output_mfa(*os_mfa, *aln, bpseq, plevel);
        if (os_bpseq!=&std::cout && os_mfa!=&std::cout)
          output_fa(std::cout, aln->name().front(), aln->consensus(), bpseq, plevel, output_energy);
        a.erase(aln);
      }
    }
    else
    {
      throw (input+": Format error").c_str();
    }
  }
  catch (const char* msg)
  {
    std::cout << msg << std::endl;
  }
  catch (std::logic_error err)
  {
    std::cout << err.what() << std::endl;
  }

  if (os_bpseq!=&std::cout) delete os_bpseq;
  if (os_mfa!=&std::cout) delete os_mfa;

  return 0;
}
