/*
 * $Id:$
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

#include "config.h"
#include <unistd.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/multi_array.hpp>
#include <boost/algorithm/string.hpp>

#include "ip.h"
#include "fa.h"
#include "aln.h"
extern "C" {
#include "new_param.h"
};

#include "contrafold/SStruct.hpp"
#include "contrafold/InferenceEngine.hpp"
#include "contrafold/Defaults.ipp"

namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/aln_util.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/PS_dot.h>
  extern void read_parameter_file(const char fname[]);
};
};

typedef unsigned int uint;

const uint n_support_parens=4;
const char* left_paren="([{<";
const char* right_paren=")]}>";

double
timing()
{
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return ru.ru_utime.tv_sec+ru.ru_utime.tv_usec*1e-6;
}

// The base class for calculating base-pairing probabilities of an indivisual sequence
class BPEngineSeq
{
public:
  BPEngineSeq() { }
  virtual ~BPEngineSeq() {}

  virtual void calculate_posterior(const std::string& seq,
                                   std::vector<float>& bp, std::vector<int>& offset) const = 0;

  virtual void calculate_posterior(const std::string& seq, const std::string& paren,
                                   std::vector<float>& bp, std::vector<int>& offset) const = 0;
};

// The base class for calculating base-pairing probabilities of aligned sequences
class BPEngineAln
{
public:
  BPEngineAln() { }
  virtual ~BPEngineAln() {}

  virtual void calculate_posterior(const std::list<std::string>& aln,
                                   std::vector<float>& bp, std::vector<int>& offset) const = 0;
  virtual void calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                                   std::vector<float>& bp, std::vector<int>& offset) const = 0;
};

class CONTRAfoldModel : public BPEngineSeq
{
public:
  CONTRAfoldModel()
    : BPEngineSeq()
  {
  }

  void calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const
  {
    SStruct ss("unknown", seq);
    ParameterManager<float> pm;
    InferenceEngine<float> en(false);
    std::vector<float> w = GetDefaultComplementaryValues<float>();
    bp.resize((seq.size()+1)*(seq.size()+2)/2, 0.0);
    en.RegisterParameters(pm);
    en.LoadValues(w);
    en.LoadSequence(ss);
    en.ComputeInside();
    en.ComputeOutside();
    en.ComputePosterior();
    en.GetPosterior(0, bp, offset);
  }

  void calculate_posterior(const std::string& seq, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const
  {
    SStruct ss("unknown", seq, paren);
    ParameterManager<float> pm;
    InferenceEngine<float> en(false);
    std::vector<float> w = GetDefaultComplementaryValues<float>();
    bp.resize((seq.size()+1)*(seq.size()+2)/2, 0.0);
    en.RegisterParameters(pm);
    en.LoadValues(w);
    en.LoadSequence(ss);
    en.UseConstraints(ss.GetMapping());
    en.ComputeInside();
    en.ComputeOutside();
    en.ComputePosterior();
    en.GetPosterior(0, bp, offset);
  }
};

class RNAfoldModel : public BPEngineSeq
{
public:
  RNAfoldModel(const char* param)
  {
    copy_new_parameters();
    if (param) Vienna::read_parameter_file(param);
  }
  
  void calculate_posterior(const std::string& seq, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const
  {
    std::string p(paren);
    std::replace(p.begin(), p.end(), '.', 'x');
    std::replace(p.begin(), p.end(), '?', '.');

    int bk = Vienna::fold_constrained;
    Vienna::fold_constrained = 1;

    uint L=seq.size();
    bp.resize((L+1)*(L+2)/2);
    offset.resize(L+1);
    for (uint i=0; i<=L; ++i)
      offset[i] = i*((L+1)+(L+1)-i-1)/2;
#if 0
    std::string str(seq.size()+1, '.');
    float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &str[0]);
    float sfact = 1.07;
    float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
    Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
#else
    Vienna::pf_scale = -1;
#endif
    Vienna::init_pf_fold(L);
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), &p[0]);
    for (uint i=0; i!=L-1; ++i)
      for (uint j=i+1; j!=L; ++j)
        bp[offset[i+1]+(j+1)] = Vienna::pr[Vienna::iindx[i+1]-(j+1)];
    Vienna::free_pf_arrays();
    Vienna::fold_constrained = bk;
  }
    
  void calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const
  {
    uint L=seq.size();
    bp.resize((L+1)*(L+2)/2);
    offset.resize(L+1);
    for (uint i=0; i<=L; ++i)
      offset[i] = i*((L+1)+(L+1)-i-1)/2;
#if 0
    std::string str(seq.size()+1, '.');
    float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &str[0]);
    float sfact = 1.07;
    float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
    Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
#else
    Vienna::pf_scale = -1;
#endif
    Vienna::init_pf_fold(L);
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
    for (uint i=0; i!=L-1; ++i)
      for (uint j=i+1; j!=L; ++j)
        bp[offset[i+1]+(j+1)] = Vienna::pr[Vienna::iindx[i+1]-(j+1)];
    Vienna::free_pf_arrays();
  }
};

class AlifoldModel : public BPEngineAln
{
public:
  AlifoldModel(const char* param)
  {
    if (param) Vienna::read_parameter_file(param);
  }

  void calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const
  {
    std::string p(paren);
    std::replace(p.begin(), p.end(), '.', 'x');
    std::replace(p.begin(), p.end(), '?', '.');

    int bk = Vienna::fold_constrained;
    Vienna::fold_constrained = 1;

    //uint N=aln.size();
    uint L=aln.front().size();
    bp.resize((L+1)*(L+2)/2, 0.0);
    offset.resize(L+1);
    for (uint i=0; i<=L; ++i)
      offset[i] = i*((L+1)+(L+1)-i-1)/2;

    char** seqs=alloc_aln(aln);
    std::string res(p);
    // scaling parameters to avoid overflow
    double min_en = Vienna::alifold(seqs, &res[0]);
    double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
    Vienna::pf_scale = exp(-(1.07*min_en)/kT/L);
    Vienna::free_alifold_arrays();

#ifdef HAVE_VIENNA18
    Vienna::plist* pi;
#else
    Vienna::pair_info* pi;
#endif
    Vienna::alipf_fold(seqs, &p[0], &pi);
    for (uint k=0; pi[k].i!=0; ++k)
      bp[offset[pi[k].i]+pi[k].j]=pi[k].p;
    free(pi);

    Vienna::free_alipf_arrays();
    free_aln(seqs);
    Vienna::fold_constrained = bk;
  }

  void calculate_posterior(const std::list<std::string>& aln,
                           std::vector<float>& bp, std::vector<int>& offset) const
  {

    //uint N=aln.size();
    uint L=aln.front().size();
    bp.resize((L+1)*(L+2)/2, 0.0);
    offset.resize(L+1);
    for (uint i=0; i<=L; ++i)
      offset[i] = i*((L+1)+(L+1)-i-1)/2;

    char** seqs=alloc_aln(aln);
    std::string res(L+1, ' ');
    // scaling parameters to avoid overflow
    double min_en = Vienna::alifold(seqs, &res[0]);
    double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
    Vienna::pf_scale = exp(-(1.07*min_en)/kT/L);
    Vienna::free_alifold_arrays();

#ifdef HAVE_VIENNA18
    Vienna::plist* pi;
#else
    Vienna::pair_info* pi;
#endif
    Vienna::alipf_fold(seqs, NULL, &pi);
    for (uint k=0; pi[k].i!=0; ++k)
      bp[offset[pi[k].i]+pi[k].j]=pi[k].p;
    free(pi);

    Vienna::free_alipf_arrays();
    free_aln(seqs);
  }

private:
  static char** alloc_aln(const std::list<std::string>& aln)
  {
    uint L=aln.front().size();
    char** seqs = new char*[aln.size()+1];
    seqs[aln.size()]=NULL;
    std::list<std::string>::const_iterator x;
    uint i=0;
    for (x=aln.begin(); x!=aln.end(); ++x)
    {
      seqs[i] = new char[L+1];
      strcpy(seqs[i++], x->c_str());
    }
    return seqs;
  }

  static void free_aln(char** seqs)
  {
    for (uint i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];
    delete[] seqs;
  }
};

class AveragedModel : public BPEngineAln
{
public:
  AveragedModel(BPEngineSeq* en) : en_(en) { }
  ~AveragedModel() { }

  void calculate_posterior(const std::list<std::string>& aln,
                           std::vector<float>& bp, std::vector<int>& offset) const
  {
    uint N=aln.size();
    uint L=aln.front().size();
    bp.resize((L+1)*(L+2)/2, 0.0);
    offset.resize(L+1);
    for (uint i=0; i<=L; ++i)
      offset[i] = i*((L+1)+(L+1)-i-1)/2;
    for (std::list<std::string>::const_iterator s=aln.begin(); s!=aln.end(); ++s)
    {
      std::vector<float> lbp;
      std::vector<int> loffset;
      std::string seq;
      std::vector<int> idx;
      for (uint i=0; i!=s->size(); ++i)
      {
        if ((*s)[i]!='-')
        {
          seq.push_back((*s)[i]);
          idx.push_back(i);
        }
      }
      en_->calculate_posterior(seq, lbp, loffset);
      for (uint i=0; i!=seq.size()-1; ++i)
        for (uint j=i+1; j!=seq.size(); ++j)
          bp[offset[idx[i]+1]+(idx[j]+1)] += lbp[loffset[i+1]+(j+1)]/N;
    }
  }

  void calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const
  {
    uint N=aln.size();
    uint L=aln.front().size();
    bp.resize((L+1)*(L+2)/2, 0.0);
    offset.resize(L+1);
    for (uint i=0; i<=L; ++i)
      offset[i] = i*((L+1)+(L+1)-i-1)/2;
    std::vector<int> p = bpseq(paren);
    for (std::list<std::string>::const_iterator s=aln.begin(); s!=aln.end(); ++s)
    {
      std::vector<float> lbp;
      std::vector<int> loffset;
      std::string seq;
      std::vector<int> idx;
      std::vector<int> rev(s->size(), -1);
      for (uint i=0; i!=s->size(); ++i)
      {
        if ((*s)[i]!='-')
        {
          seq.push_back((*s)[i]);
          idx.push_back(i);
          rev[i]=seq.size()-1;
        }
      }
      std::string lparen(seq.size(), '.');
      for (uint i=0; i!=p.size(); ++i)
      {
        if (rev[i]>=0)
        {
          if (p[i]<0 || rev[p[i]]>=0)
            lparen[rev[i]] = paren[i];
          else
            lparen[rev[i]] = '.';
        }
      }
      
      en_->calculate_posterior(seq, lparen, lbp, loffset);
      for (uint i=0; i!=seq.size()-1; ++i)
        for (uint j=i+1; j!=seq.size(); ++j)
          bp[offset[idx[i]+1]+(idx[j]+1)] += lbp[loffset[i+1]+(j+1)]/N;
    }
  }

private:
  std::vector<int> bpseq(const std::string& paren) const
  {
    std::vector<int> ret(paren.size(), -1);
    std::stack<int> st;
    for (uint i=0; i!=paren.size(); ++i)
    {
      if (paren[i]=='(')
      {
        st.push(i);
      }
      else if (paren[i]==')')
      {
        ret[st.top()]=i;
        ret[i]=st.top();
        st.pop();
      }
    }
    return ret;
  }

private:
  BPEngineSeq* en_;
};

class AuxModel
{
public:
  AuxModel() { }

  bool calculate_posterior(const char* filename, std::string& seq,
                           std::vector<float>& bp, std::vector<int>& offset) const
  {
    std::string l;
    uint L=0;
    std::ifstream in(filename);
    if (!in) return false;
    while (std::getline(in, l)) ++L;
    bp.resize((L+1)*(L+2)/2, 0.0);
    offset.resize(L+1);
    for (uint i=0; i<=L; ++i)
      offset[i] = i*((L+1)+(L+1)-i-1)/2;
    seq.resize(L);
    in.clear();
    in.seekg(0, std::ios::beg);
    while (std::getline(in, l))
    {
      std::vector<std::string> v;
      boost::algorithm::split(v, l, boost::is_space(), boost::algorithm::token_compress_on);
      uint up = atoi(v[0].c_str());
      seq[up-1] = v[1][0];
      for (uint i=2; i!=v.size(); ++i)
      {
        uint down;
        float p;
        if (sscanf(v[i].c_str(), "%u:%f", &down, &p)==2)
          bp[offset[up]+down]=p;
      }
    }

    return true;
  }
};

class IPknot
{
public:
  IPknot(uint pk_level, const float* th, const float* alpha,
         bool stacking_constraints, int n_th)
    : pk_level_(pk_level),
      th_(th, th+pk_level),
      alpha_(alpha, alpha+pk_level),
      stacking_constraints_(stacking_constraints),
      n_th_(n_th)
  {
  }

  void set_parameters(uint pk_level, const float* th, const float* alpha,
                      bool stacking_constraints)
  {
    pk_level_ = pk_level;
    th_.resize(pk_level_);
    std::copy(th, th+pk_level_, th_.begin());
    alpha_.resize(pk_level_);
    std::copy(alpha, alpha+pk_level_, alpha_.begin());
    stacking_constraints_ = stacking_constraints;
  }

  template < class SEQ, class EN >
  void solve(const SEQ& seq, EN& en, std::string& r, std::vector<int>& bpseq, int fill=0) const
  {
    uint L=length(seq);
    IP* ip = NULL;
    boost::multi_array<int, 3> v(boost::extents[pk_level_][L][L]);
    boost::multi_array<std::vector<int>, 2> w(boost::extents[pk_level_][L]);
    std::fill(v.data(), v.data()+v.num_elements(), -1);

    std::vector<float> bp;
    std::vector<int> offset;
    en.calculate_posterior(seq, bp, offset);
    ip = new IP(IP::MAX, n_th_);
    solve(length(seq), bp, offset, *ip, v, w);

    while (--fill>=0)
    {
      // re-calculate the base-pairing probability matrix
      std::vector<float> bpl;
      std::vector<int> offsetl;
      std::fill(bp.begin(), bp.end(), 0.0);

      for (uint l=0; l!=pk_level_; ++l)
      {
        std::string str(L, '?');
        for (uint lv=0; lv!=pk_level_; ++lv)
          for (uint i=0; i<L; ++i)
            for (uint j=i+1; j<L; ++j)
              if (v[lv][i][j]>=0 && ip->get_value(v[lv][i][j])>0.5)
              {
                str[i]= l==lv ? '(' : '.';
                str[j]= l==lv ? ')' : '.';
              }

        std::fill(bpl.begin(), bpl.end(), 0.0);
        en.calculate_posterior(seq, str, bpl, offsetl);
        assert(bp.size()==bpl.size());
        for (uint k=0; k!=bp.size(); ++k) bp[k]+=bpl[k];
      }
#ifndef NDEBUG
      for (uint k=0; k!=bp.size(); ++k) assert(bp[k]<=1.0);
#endif
      delete ip;
      ip = new IP(IP::MAX, n_th_);

      // solve IP under the updated base-pairing probability matrix
      std::fill(v.data(), v.data()+v.num_elements(), -1);
      std::fill(w.data(), w.data()+w.num_elements(), std::vector<int>());
      solve(L, bp, offset, *ip, v, w);
    }

    // build the resultant structure
    r.resize(L);
    std::fill(r.begin(), r.end(), '.');
    bpseq.resize(L);
    std::fill(bpseq.begin(), bpseq.end(), -1);
    for (uint lv=0; lv!=pk_level_; ++lv)
      for (uint i=0; i<L; ++i)
        for (uint j=i+1; j<L; ++j)
          if (v[lv][i][j]>=0 && ip->get_value(v[lv][i][j])>0.5)
          {
            assert(r[i]=='.'); assert(r[j]=='.');
            bpseq[i]=j; bpseq[j]=i;
            if (lv<n_support_parens)
            {
              r[i]=left_paren[lv]; r[j]=right_paren[lv];
            }
          }
    delete ip;
  }

  void solve(uint L, const std::vector<float>& bp, const std::vector<int>& offset,
             std::string& r, std::vector<int>& bpseq) const
  {
    IP ip(IP::MAX, n_th_);
    boost::multi_array<int, 3> v(boost::extents[pk_level_][L][L]);
    boost::multi_array<std::vector<int>, 2> w(boost::extents[pk_level_][L]);
    std::fill(v.data(), v.data()+v.num_elements(), -1);

    solve(L, bp, offset, ip, v, w);

    // build the resultant structure
    r.resize(L);
    std::fill(r.begin(), r.end(), '.');
    bpseq.resize(L);
    std::fill(bpseq.begin(), bpseq.end(), -1);
    for (uint lv=0; lv!=pk_level_; ++lv)
      for (uint i=0; i<L; ++i)
        for (uint j=i+1; j<L; ++j)
          if (v[lv][i][j]>=0 && ip.get_value(v[lv][i][j])>0.5)
          {
            assert(r[i]=='.'); assert(r[j]=='.');
            bpseq[i]=j; bpseq[j]=i;
            if (lv<n_support_parens)
            {
              r[i]=left_paren[lv]; r[j]=right_paren[lv];
            }
          }
  }

  void solve(uint L, const std::vector<float>& bp, const std::vector<int>& offset,
             IP& ip, boost::multi_array<int, 3>& v,
             boost::multi_array<std::vector<int>, 2>& w) const
  {
    // make objective variables with their weights
    for (uint j=1; j!=L; ++j)
    {
      for (uint i=j-1; i!=-1u; --i)
      {
        const float& p=bp[offset[i+1]+(j+1)];
        for (uint lv=0; lv!=pk_level_; ++lv)
          if (p>th_[lv])
          {
            v[lv][i][j] = ip.make_variable(p*alpha_[lv]);
            w[lv][i].push_back(j);
          }
      }
    }
    ip.update();

    // constraint 1: each s_i is paired with at most one base
    for (uint i=0; i!=L; ++i)
    {
      int row = ip.make_constraint(IP::UP, 0, 1);
      for (uint lv=0; lv!=pk_level_; ++lv)
      {
        for (uint j=0; j<i; ++j)
          if (v[lv][j][i]>=0)
            ip.add_constraint(row, v[lv][j][i], 1);
        for (uint j=i+1; j<L; ++j)
          if (v[lv][i][j]>=0)
            ip.add_constraint(row, v[lv][i][j], 1);
      }
    }

    // constraint 2: disallow pseudoknots in x[lv]
    for (uint lv=0; lv!=pk_level_; ++lv)
      for (uint i=0; i<w[lv].size(); ++i)
        for (uint p=0; p<w[lv][i].size(); ++p)
        {
          uint j=w[lv][i][p];
          for (uint k=i+1; k<j; ++k)
            for (uint q=0; q<w[lv][k].size(); ++q)
            {
              uint l=w[lv][k][q];
              if (j<l)
              {
                int row = ip.make_constraint(IP::UP, 0, 1);
                ip.add_constraint(row, v[lv][i][j], 1);
                ip.add_constraint(row, v[lv][k][l], 1);
              }
            }
        }

    // constraint 3: any x[t]_kl must be pseudoknotted with x[u]_ij for t>u
    for (uint lv=1; lv!=pk_level_; ++lv)
      for (uint k=0; k<w[lv].size(); ++k)
        for (uint q=0; q<w[lv][k].size(); ++q)
        {
          uint l=w[lv][k][q];
          for (uint plv=0; plv!=lv; ++plv)
          {
            int row = ip.make_constraint(IP::LO, 0, 0);
            ip.add_constraint(row, v[lv][k][l], -1);
            for (uint i=0; i<k; ++i)
              for (uint p=0; p<w[plv][i].size(); ++p)
              {
                uint j=w[plv][i][p];
                if (k<j && j<l)
                  ip.add_constraint(row, v[plv][i][j], 1);
              }
            for (uint i=k+1; i<l; ++i)
              for (uint p=0; p<w[plv][i].size(); ++p)
              {
                uint j=w[plv][i][p];
                if (l<j)
                  ip.add_constraint(row, v[plv][i][j], 1);
              }
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
          for (uint j=0; j<i; ++j)
            if (v[lv][j][i]>=0)
              ip.add_constraint(row, v[lv][j][i], -1);
          if (i>0)
            for (uint j=0; j<i-1; ++j)
              if (v[lv][j][i-1]>=0)
                ip.add_constraint(row, v[lv][j][i-1], 1);
          if (i+1<L)
            for (uint j=0; j<i+1; ++j)
              if (v[lv][j][i+1]>=0)
                ip.add_constraint(row, v[lv][j][i+1], 1);
        }

        // downstream
        for (uint i=0; i<L; ++i)
        {
          int row = ip.make_constraint(IP::LO, 0, 0);
          for (uint j=i+1; j<L; ++j)
            if (v[lv][i][j]>=0)
              ip.add_constraint(row, v[lv][i][j], -1);
          if (i>0)
            for (uint j=i; j<L; ++j)
              if (v[lv][i-1][j]>=0)
                ip.add_constraint(row, v[lv][i-1][j], 1);
          if (i+1<L)
            for (uint j=i+2; j<L; ++j)
              if (v[lv][i+1][j]>=0)
                ip.add_constraint(row, v[lv][i+1][j], 1);
        }
      }
    }

    // execute optimization
    ip.solve();
  }

private:
  uint length(const std::string& seq) const { return seq.size(); }
  uint length(const std::list<std::string>& aln) const { return aln.front().size(); }

private:
  // options
  uint pk_level_;
  std::vector<float> th_;
  std::vector<float> alpha_;
  bool stacking_constraints_;
  int n_th_;
};

void
compute_expected_accuracy(const std::vector<int>& bpseq,
                          const std::vector<float>& bp, const std::vector<int>& offset,
                          double& sen, double& ppv, double& mcc)
{
  int L  = bpseq.size();
  int L2 = L*(L-1)/2;
  int N;

  double sump = 0.0;
  double etp  = 0.0;

  for (uint i=0; i!=bp.size(); ++i) sump += bp[i];

  for (uint i=0; i!=bpseq.size(); ++i)
  {
    if (bpseq[i]!=-1 && bpseq[i]>(int)i)
    {
      etp += bp[offset[i+1]+bpseq[i]+1];
      N++;
    }
  }

  double etn = L2 - N - sump + etp;
  double efp = N - etp;
  double efn = sump - etp;

  sen = ppv = mcc = 0;
  if (etp+efn!=0) sen = etp / (etp + efn);
  if (etp+efp!=0) ppv = etp / (etp + efp);
  if (etp+efp!=0 && etp+efn!=0 && etn+efp!=0 && etn+efn!=0)
    mcc = (etp*etn-efp*efn) / std::sqrt((etp+efp)*(etp+efn)*(etn+efp)*(etn+efn));
}

void
usage(const char* progname)
{
  std::cout << progname << ": [options] fasta" << std::endl
            << " -h:       show this message" << std::endl
    //      << " -a alpha: weight for each level" << std::endl
            << " -t th:    threshold of base-pairing probabilities for each level" << std::endl
            << " -g gamma: weight for true base-pairs equivalent to -t 1/(gamma+1)" << std::endl
            << "           (default: -g 4 -g 8)" << std::endl
            << " -m model: probaility distribution model" << std::endl
            << " -i:       allow isolated base-pairs" << std::endl
            << " -b:       output the prediction via BPSEQ format" << std::endl
            << " -P param: read the energy parameter file for the Vienna RNA package" << std::endl
#ifndef WITH_GLPK
            << " -n n_th:  specify the number of threads (default: 1)" << std::endl
#endif
    ;
}

int
main(int argc, char* argv[])
{
  char* progname=argv[0];
  // parse options
  char ch;
  std::vector<float> th;
  std::vector<float> alpha;
  bool isolated_bp=false;
  char* model=NULL;
  bool use_bpseq=false;
  int n_th=1;
  int n_fill=0;
  const char* param=NULL;
  bool aux=false;
  std::vector<float> temp;
  while ((ch=getopt(argc, argv, "a:t:g:m:f:ibn:P:xh"))!=-1)
  {
    switch (ch)
    {
      case 'm':
        model=optarg;
        break;
      case 'f':
        n_fill=atoi(optarg);
        break;
      case 'a':
        alpha.push_back(atof(optarg));
        break;
      case 't':
        th.push_back(atof(optarg));
        break;
      case 'g':
        th.push_back(1/(atof(optarg)+1));
        break;
      case 'i':
        isolated_bp=true;
        break;
      case 'b':
        use_bpseq=true;
        break;
      case 'n':
        n_th=atoi(optarg);
        break;
      case 'P':
        param=optarg;
        break;
      case 'x':
        aux=true;
        break;
      case 'h': case '?': default:
        usage(progname);
        return 1;
        break;
    }
  }
  argc -= optind;
  argv += optind;

  if (argc!=1) { usage(progname); return 1; }

  // set default parameters
  if (th.empty())
  {
    th.resize(2);
    th[0]=1/(4.0+1);
    th[1]=1/(8.0+1);
  }
  if (alpha.empty())
  {
    alpha.resize(2);
    for (uint i=0; i!=alpha.size(); ++i)
      alpha[i]=1.0/alpha.size();
  }
  IPknot ipknot(th.size(), &th[0], &alpha[0], !isolated_bp, n_th);
  std::string r;
  std::vector<int> bpseq;
  
  std::list<Fasta> f;
  std::list<Aln> a;
  if (aux)
  {
    std::vector<float> bp;
    std::vector<int> offset;
    AuxModel aux;
    std::string seq;
    aux.calculate_posterior(argv[0], seq, bp, offset);
    ipknot.solve(seq.size(), bp, offset, r, bpseq);
    if (!use_bpseq && th.size()<n_support_parens)
    {
      std::cout << ">" << argv[0] << std::endl
                << seq << std::endl << r << std::endl;
    }
    else
    {
      std::cout << "# " << argv[0] << std::endl;
      for (uint i=0; i!=bpseq.size(); ++i)
        std::cout << i+1 << " " << seq[i] << " " << bpseq[i]+1 << std::endl;
    }
  }    
  else if (Fasta::load(f, argv[0])>0)
  {
    BPEngineSeq* en=NULL;
    if (model==NULL || strcasecmp(model, "McCaskill")==0)
      en = new RNAfoldModel(param);
    else if (strcasecmp(model, "CONTRAfold")==0)
      en = new CONTRAfoldModel();
    else
    {
      usage(progname);
      return 1;
    }
    
    while (!f.empty())
    {
      std::list<Fasta>::iterator fa = f.begin();
      ipknot.solve(fa->seq(), *en, r, bpseq, n_fill);
      if (!use_bpseq && th.size()<n_support_parens)
      {
        std::cout << ">" << fa->name() << std::endl
                  << fa->seq() << std::endl << r << std::endl;
      }
      else
      {
        //double sen, ppv, mcc;
        //compute_expected_accuracy(bpseq, bp, offset, sen, ppv, mcc);
        std::cout << "# " << fa->name() << std::endl;
        for (uint i=0; i!=bpseq.size(); ++i)
          std::cout << i+1 << " " << fa->seq()[i] << " " << bpseq[i]+1 << std::endl;
      }
      f.erase(fa);
    }

    delete en;
  }
  else if (Aln::load(a, argv[0])>0)
  {
    BPEngineAln* en=NULL;
    BPEngineSeq* en_s=NULL;
    if (model==NULL || strcasecmp(model, "McCaskill")==0)
    {
      en_s = new RNAfoldModel(param);
      en = new AveragedModel(en_s);
    }
    else if (strcasecmp(model, "CONTRAfold")==0)
    {
      en_s = new CONTRAfoldModel();
      en = new AveragedModel(en_s);
    }
    else if (strcasecmp(model, "Alifold")==0)
      en = new AlifoldModel(param);
    else
    {
      usage(progname);
      return 1;
    }

    while (!a.empty())
    {
      std::list<Aln>::iterator aln = a.begin();
      std::string consensus(aln->consensus());
      ipknot.solve(aln->seq(), *en, r, bpseq, n_fill);
      if (!use_bpseq && th.size()<n_support_parens)
      {
        std::cout << ">" << aln->name().front() << std::endl
                  << consensus << std::endl << r << std::endl;
      }
      else
      {
        std::cout << "# " << aln->name().front() << std::endl;
        for (uint i=0; i!=bpseq.size(); ++i)
          std::cout << i+1 << " " << consensus[i] << " " << bpseq[i]+1 << std::endl;
      }
      a.erase(aln);
    }

    if (en) delete en;
    if (en_s) delete en_s;
  }

  return 0;
}
