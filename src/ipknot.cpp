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
#include <iostream>
#include <string>
#include <vector>
#include <boost/multi_array.hpp>

#include "ip.h"
#include "fa.h"
#include "contrafold/SStruct.hpp"
#include "contrafold/InferenceEngine.hpp"
#include "contrafold/DuplexEngine.hpp"
#include "contrafold/Defaults.ipp"

namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
};
};

typedef unsigned int uint;

const int n_support_parens=4;
const char* left_paren="([{<";
const char* right_paren=")]}>";

class IPknot
{
public:
  IPknot(uint pk_level, const float* th, const float* alpha,
         bool use_contrafold, bool stacking_constraints, int n_th)
    : ip_(IP::MAX, n_th),
      pk_level_(pk_level),
      th_(th, th+pk_level),
      alpha_(alpha, alpha+pk_level),
      use_contrafold_(use_contrafold),
      stacking_constraints_(stacking_constraints)
  {
  }

  void solve(const std::string& seq, std::string& res, std::vector<int>& ct);

private:
  void contrafold(const std::string& seq);
  void rnafold(const std::string& seq);

private:
  IP ip_;
  
  // options
  uint pk_level_;
  std::vector<float> th_;
  std::vector<float> alpha_;
  bool use_contrafold_;        // use CONTRAfold model or not
  bool stacking_constraints_;

  // index
  boost::multi_array<int, 3> v_;
};

void
IPknot::
contrafold(const std::string& seq)
{
  float th_min=1.0;
  for (uint i=0; i!=th_.size(); ++i) th_min=std::min(th_min, th_[i]);
  SStruct ss("unknown", seq);
  ParameterManager<float> pm;
  InferenceEngine<float> en(false);
  std::vector<float> w = GetDefaultComplementaryValues<float>();
  std::vector<float> bp((seq.size()+1)*(seq.size()+2)/2, 0.0);
  en.RegisterParameters(pm);
  en.LoadValues(w);
  en.LoadSequence(ss);
  en.ComputeInside();
  en.ComputeOutside();
  en.ComputePosterior();
  en.GetPosterior(th_min, bp);

  for (uint j=1; j!=seq.size(); ++j)
  {
    for (uint i=j-1; i!=-1u; --i)
    {
      float p=bp[en.GetOffset(i+1)+(j+1)];
      for (uint lv=0; lv!=pk_level_; ++lv)
        if (p>th_[lv])
          v_[lv][i][j] = ip_.make_variable(p*alpha_[lv]);
    }
  }
}

void
IPknot::
rnafold(const std::string& seq)
{
  Vienna::pf_scale = -1;
  Vienna::init_pf_fold(seq.size());
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  for (uint i=0; i!=seq.size()-1; ++i)
  {
    for (uint j=i+1; j!=seq.size(); ++j)
    {
      float p=Vienna::pr[Vienna::iindx[i+1]-(j+1)];
      for (uint lv=0; lv!=pk_level_; ++lv)
        if (p>th_[lv])
          v_[lv][i][j] = ip_.make_variable(p*alpha_[lv]);
    }
  }
  Vienna::free_pf_arrays();
}

void
IPknot::
solve(const std::string& s, std::string& r, std::vector<int>& bpseq)
{
  v_.resize(boost::extents[pk_level_][s.size()][s.size()]);
  std::fill(v_.data(), v_.data()+v_.num_elements(), -1);

  // make objective variavles with their weights
  if (use_contrafold_) contrafold(s);
  else rnafold(s);
  ip_.update();

  // constraint 1: each s_i is paired with at most one base
  for (uint i=0; i!=s.size(); ++i)
  {
    int row = ip_.make_constraint(IP::UP, 0, 1);
    for (uint lv=0; lv!=pk_level_; ++lv)
    {
      for (uint j=0; j<i; ++j)
        if (v_[lv][j][i]>=0)
          ip_.add_constraint(row, v_[lv][j][i], 1);
      for (uint j=i+1; j<s.size(); ++j)
        if (v_[lv][i][j]>=0)
          ip_.add_constraint(row, v_[lv][i][j], 1);
    }
  }

  // constraint 2: disallow pseudoknots in x[lv]
  for (uint lv=0; lv!=pk_level_; ++lv)
    for (uint i=0; i<s.size(); ++i)
      for (uint k=i+1; k<s.size(); ++k)
        for (uint j=k+1; j<s.size(); ++j)
          if (v_[lv][i][j]>=0)
            for (uint l=j+1; l<s.size(); ++l)
              if (v_[lv][k][l]>=0)
              {
                int row = ip_.make_constraint(IP::UP, 0, 1);
                ip_.add_constraint(row, v_[lv][i][j], 1);
                ip_.add_constraint(row, v_[lv][k][l], 1);
              }

  // constraint 3: any x[t]_kl must be pseudoknotted with x[u]_ij for t>u
  for (uint lv=1; lv!=pk_level_; ++lv)
    for (uint k=0; k<s.size(); ++k)
      for (uint l=k+1; l<s.size(); ++l)
        if (v_[lv][k][l]>=0)
          for (uint plv=0; plv!=lv; ++plv)
          {
            int row = ip_.make_constraint(IP::LO, 0, 0);
            ip_.add_constraint(row, v_[lv][k][l], -1);
            for (uint i=0; i<k; ++i)
              for (uint j=k+1; j<l; ++j)
                if (v_[plv][i][j]>=0)
                  ip_.add_constraint(row, v_[plv][i][j], 1);
            for (uint i=k+1; i<l; ++i)
              for (uint j=l+1; j<s.size(); ++j)
                if (v_[plv][i][j]>=0)
                  ip_.add_constraint(row, v_[plv][i][j], 1);
          }

  if (stacking_constraints_)
  {
    for (uint lv=0; lv!=pk_level_; ++lv)
    {
      // upstream
      for (uint i=0; i<s.size(); ++i)
      {
        int row = ip_.make_constraint(IP::LO, 0, 0);
        for (uint j=0; j<i; ++j)
          if (v_[lv][j][i]>=0)
            ip_.add_constraint(row, v_[lv][j][i], -1);
        if (i>0)
          for (uint j=0; j<i-1; ++j)
            if (v_[lv][j][i-1]>=0)
              ip_.add_constraint(row, v_[lv][j][i-1], 1);
        if (i+1<s.size())
          for (uint j=0; j<i+1; ++j)
            if (v_[lv][j][i+1]>=0)
              ip_.add_constraint(row, v_[lv][j][i+1], 1);
      }

      // downstream
      for (uint i=0; i<s.size(); ++i)
      {
        int row = ip_.make_constraint(IP::LO, 0, 0);
        for (uint j=i+1; j<s.size(); ++j)
          if (v_[lv][i][j]>=0)
            ip_.add_constraint(row, v_[lv][i][j], -1);
        if (i>0)
          for (uint j=i; j<s.size(); ++j)
            if (v_[lv][i-1][j]>=0)
              ip_.add_constraint(row, v_[lv][i-1][j], 1);
        if (i+1<s.size())
          for (uint j=i+2; j<s.size(); ++j)
            if (v_[lv][i+1][j]>=0)
              ip_.add_constraint(row, v_[lv][i+1][j], 1);
      }
    }
  }

  // execute optimization
  ip_.solve();

  // build the resultant structure
  r.resize(s.size());
  std::fill(r.begin(), r.end(), '.');
  bpseq.resize(s.size());
  std::fill(bpseq.begin(), bpseq.end(), -1);
  for (uint lv=0; lv!=pk_level_; ++lv)
    for (uint i=0; i<s.size(); ++i)
      for (uint j=i+1; j<s.size(); ++j)
        if (v_[lv][i][j]>=0 && ip_.get_value(v_[lv][i][j])>0.5)
        {
          assert(r[i]=='.'); assert(r[j]=='.');
          r[i]=left_paren[lv]; r[j]=right_paren[lv];
          bpseq[i]=j; bpseq[j]=i;
        }
}

void
usage(const char* progname)
{
  std::cout << progname << ": [options] fasta" << std::endl
            << " -h:       show this message" << std::endl
            << " -a alpha: weight for hybridation probabilities (default: 0.5)" << std::endl
            << " -t th:    threshold of base-pairing probabilities (default: 0.5)" << std::endl
            << " -m:       use McCaskill model (default: CONTRAfold model)" << std::endl
            << " -i:       allow isolated base-pairs" << std::endl
            << " -b:       output the prediction via BPSEQ format" << std::endl
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
  bool use_contrafold=true;
  bool use_bpseq=false;
  int n_th=1;
  while ((ch=getopt(argc, argv, "a:t:mibn:r:h"))!=-1)
  {
    switch (ch)
    {
      case 'm':
        use_contrafold=false;
        break;
      case 'a':
        alpha.push_back(atof(optarg));
        break;
      case 't':
        th.push_back(atof(optarg));
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
      case 'h': case '?': default:
        usage(progname);
        return 1;
        break;
    }
  }
  argc -= optind;
  argv += optind;

  if (argc!=1) { usage(progname); return 1; }
  std::list<Fasta> f;
  Fasta::load(f, argv[0]);

  IPknot ipknot(th.size(), &th[0], &alpha[0], use_contrafold, !isolated_bp, n_th);
  while (!f.empty())
  {
    std::list<Fasta>::iterator fa = f.begin();
    std::string r;
    std::vector<int> bpseq;
    ipknot.solve(fa->seq(), r, bpseq);
    if (!use_bpseq)
    {
      std::cout << ">" << fa->name() << std::endl
                << fa->seq() << std::endl << r << std::endl;
    }
    else
    {
      std::cout << "# " << fa->name() << std::endl;
      for (uint i=0; i!=bpseq.size(); ++i)
        std::cout << i+1 << " " << fa->seq()[i] << " " << bpseq[i]+1 << std::endl;
    }
    f.erase(fa);
  }

  return 0;
}
