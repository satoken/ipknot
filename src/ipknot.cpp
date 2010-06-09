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
#include <ctime>
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
    : pk_level_(pk_level),
      th_(th, th+pk_level),
      alpha_(alpha, alpha+pk_level),
      use_contrafold_(use_contrafold),
      stacking_constraints_(stacking_constraints),
      n_th_(n_th)
  {
  }

  void solve(const std::string& seq, std::string& res, std::vector<int>& ct) const;

  void solve(const std::string& seq, const std::vector<float>& bp, const std::vector<int>& offset,
             std::string& res, std::vector<int>& ct) const;

  void calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const
  {
    use_contrafold_ ? contrafold(seq, bp, offset) : rnafold(seq, bp, offset);
  }

private:
  void contrafold(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const;
  void rnafold(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const;

private:
  // options
  uint pk_level_;
  std::vector<float> th_;
  std::vector<float> alpha_;
  bool use_contrafold_;        // use CONTRAfold model or not
  bool stacking_constraints_;
  int n_th_;
};

void
IPknot::
contrafold(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const
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

void
IPknot::
rnafold(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const
{
  bp.resize((seq.size()+1)*(seq.size()+2)/2);
  offset.resize(seq.size()+1);
  for (uint i=0; i<=offset.size(); ++i)
    offset[i] = i*(offset.size()+offset.size()-i-1)/2;
  
  Vienna::pf_scale = -1;
  Vienna::init_pf_fold(seq.size());
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  for (uint i=0; i!=seq.size()-1; ++i)
    for (uint j=i+1; j!=seq.size(); ++j)
      bp[offset[i+1]+(j+1)] = Vienna::pr[Vienna::iindx[i+1]-(j+1)];
  Vienna::free_pf_arrays();
}

void
IPknot::
solve(const std::string& s, std::string& r, std::vector<int>& bpseq) const
{
  std::vector<float> bp;
  std::vector<int> offset;

  clock_t t1 = clock();
  std::cerr << "Calculating base-pairing probabilities ...";
  calculate_posterior(s, bp, offset);
  clock_t t2 = clock();
  std::cerr << " done ("
            << static_cast<float>(t2-t1)/CLOCKS_PER_SEC
            << "s)." << std::endl;
  
  solve(s, bp, offset, r, bpseq);
}

void
IPknot::
solve(const std::string& s, const std::vector<float>& bp, const std::vector<int>& offset,
      std::string& r, std::vector<int>& bpseq) const
{
  IP ip(IP::MAX, n_th_);
  
  boost::multi_array<int, 3> v(boost::extents[pk_level_][s.size()][s.size()]);
  boost::multi_array<std::vector<int>, 2> w(boost::extents[pk_level_][s.size()]);
  std::fill(v.data(), v.data()+v.num_elements(), -1);

  // make objective variavles with their weights
  clock_t t1 = clock();
  std::cerr << "Making variables ...";
  for (uint j=1; j!=s.size(); ++j)
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
  clock_t t2 = clock();
  std::cerr << " done ("
            << static_cast<float>(t2-t1)/CLOCKS_PER_SEC
            << "s)." << std::endl;

  // constraint 1: each s_i is paired with at most one base
  t1 = clock();
  std::cerr << "Making constraints 1 ...";
  for (uint i=0; i!=s.size(); ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint lv=0; lv!=pk_level_; ++lv)
    {
      for (uint j=0; j<i; ++j)
        if (v[lv][j][i]>=0)
          ip.add_constraint(row, v[lv][j][i], 1);
      for (uint j=i+1; j<s.size(); ++j)
        if (v[lv][i][j]>=0)
          ip.add_constraint(row, v[lv][i][j], 1);
    }
  }
  t2 = clock();
  std::cerr << " done (" 
            << static_cast<float>(t2-t1)/CLOCKS_PER_SEC
            << "s)." << std::endl;

  // constraint 2: disallow pseudoknots in x[lv]
  t1 = clock();
  std::cerr << "Making constraints 2 ...";
#if 0
  for (uint lv=0; lv!=pk_level_; ++lv)
    for (uint i=0; i<s.size(); ++i)
      for (uint k=i+1; k<s.size(); ++k)
        for (uint j=k+1; j<s.size(); ++j)
          if (v[lv][i][j]>=0)
            for (uint l=j+1; l<s.size(); ++l)
              if (v[lv][k][l]>=0)
              {
                int row = ip.make_constraint(IP::UP, 0, 1);
                ip.add_constraint(row, v[lv][i][j], 1);
                ip.add_constraint(row, v[lv][k][l], 1);
              }
#else
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
#endif
  t2 = clock();
  std::cerr << " done (" 
            << static_cast<float>(t2-t1)/CLOCKS_PER_SEC
            << "s)." << std::endl;

  // constraint 3: any x[t]_kl must be pseudoknotted with x[u]_ij for t>u
  t1 = clock();
  std::cerr << "Making constraints 3 ...";
#if 0
  for (uint lv=1; lv!=pk_level_; ++lv)
    for (uint k=0; k<s.size(); ++k)
      for (uint l=k+1; l<s.size(); ++l)
        if (v[lv][k][l]>=0)
          for (uint plv=0; plv!=lv; ++plv)
          {
            int row = ip.make_constraint(IP::LO, 0, 0);
            ip.add_constraint(row, v[lv][k][l], -1);
            for (uint i=0; i<k; ++i)
              for (uint j=k+1; j<l; ++j)
                if (v[plv][i][j]>=0)
                  ip.add_constraint(row, v[plv][i][j], 1);
            for (uint i=k+1; i<l; ++i)
              for (uint j=l+1; j<s.size(); ++j)
                if (v[plv][i][j]>=0)
                  ip.add_constraint(row, v[plv][i][j], 1);
          }
#else
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
#endif
  t2 = clock();
  std::cerr << " done (" 
            << static_cast<float>(t2-t1)/CLOCKS_PER_SEC
            << "s)." << std::endl;

  if (stacking_constraints_)
  {
    t1 = clock();
    std::cerr << "Making stacking constraints ...";
    for (uint lv=0; lv!=pk_level_; ++lv)
    {
      // upstream
      for (uint i=0; i<s.size(); ++i)
      {
        int row = ip.make_constraint(IP::LO, 0, 0);
        for (uint j=0; j<i; ++j)
          if (v[lv][j][i]>=0)
            ip.add_constraint(row, v[lv][j][i], -1);
        if (i>0)
          for (uint j=0; j<i-1; ++j)
            if (v[lv][j][i-1]>=0)
              ip.add_constraint(row, v[lv][j][i-1], 1);
        if (i+1<s.size())
          for (uint j=0; j<i+1; ++j)
            if (v[lv][j][i+1]>=0)
              ip.add_constraint(row, v[lv][j][i+1], 1);
      }

      // downstream
      for (uint i=0; i<s.size(); ++i)
      {
        int row = ip.make_constraint(IP::LO, 0, 0);
        for (uint j=i+1; j<s.size(); ++j)
          if (v[lv][i][j]>=0)
            ip.add_constraint(row, v[lv][i][j], -1);
        if (i>0)
          for (uint j=i; j<s.size(); ++j)
            if (v[lv][i-1][j]>=0)
              ip.add_constraint(row, v[lv][i-1][j], 1);
        if (i+1<s.size())
          for (uint j=i+2; j<s.size(); ++j)
            if (v[lv][i+1][j]>=0)
              ip.add_constraint(row, v[lv][i+1][j], 1);
      }
    }
    t2 = clock();
    std::cerr << " done (" 
              << static_cast<float>(t2-t1)/CLOCKS_PER_SEC
              << "s)." << std::endl;
  }

  // execute optimization
  t1 = clock();
  std::cerr << "Solving IP problem...";
  ip.solve();
  t2 = clock();
  std::cerr << " done (" 
            << static_cast<float>(t2-t1)/CLOCKS_PER_SEC
            << "s)." << std::endl;

  // build the resultant structure
  r.resize(s.size());
  std::fill(r.begin(), r.end(), '.');
  bpseq.resize(s.size());
  std::fill(bpseq.begin(), bpseq.end(), -1);
  for (uint lv=0; lv!=pk_level_; ++lv)
    for (uint i=0; i<s.size(); ++i)
      for (uint j=i+1; j<s.size(); ++j)
        if (v[lv][i][j]>=0 && ip.get_value(v[lv][i][j])>0.5)
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
            << " -a alpha: weight for hybridation probabilities" << std::endl
            << " -t th:    threshold of base-pairing probabilities" << std::endl
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
  while ((ch=getopt(argc, argv, "a:t:g:mibn:r:h"))!=-1)
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
