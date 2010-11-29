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

#include "config.h"
#include <unistd.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <strings.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/multi_array.hpp>

#include "ip.h"
#include "fa.h"
#include "aln.h"
#include "fold.h"

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

class IPknot
{
public:
  IPknot(uint pk_level, const float* th, const float* alpha,
         bool levelwise, bool stacking_constraints, int n_th)
    : pk_level_(pk_level),
      th_(th, th+pk_level_),
      alpha_(alpha, alpha+pk_level_),
      levelwise_(levelwise),
      stacking_constraints_(stacking_constraints),
      n_th_(n_th)
  {
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
    make_result(L, *ip, v, w, r, bpseq);
    delete ip;
  }

  void solve(uint L, const std::vector<float>& bp, const std::vector<int>& offset,
             std::string& r, std::vector<int>& bpseq) const
  {
    IP ip(IP::MAX, n_th_);
    boost::multi_array<int, 3> v(boost::extents[pk_level_][L][L]);
    boost::multi_array<std::vector<int>, 2> w(boost::extents[pk_level_][L]);
    std::fill(v.data(), v.data()+v.num_elements(), -1);

    // solve the IP problem
    solve(L, bp, offset, ip, v, w);

    // build the resultant structure
    make_result(L, ip, v, w, r, bpseq);
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

    if (levelwise_)
    {
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

  void make_result(uint L, const IP& ip, boost::multi_array<int, 3>& v,
                   boost::multi_array<std::vector<int>, 2>& w,
                   std::string& r, std::vector<int>& bpseq) const
  {
    r.resize(L);
    std::fill(r.begin(), r.end(), '.');
    bpseq.resize(L);
    std::fill(bpseq.begin(), bpseq.end(), -1);

    if (levelwise_)
    {
      for (uint lv=0; lv!=pk_level_; ++lv)
      {
        for (uint i=0; i<L; ++i)
          for (uint j=i+1; j<L; ++j)
            if (v[lv][i][j]>=0 && ip.get_value(v[lv][i][j])>0.5)
            {
              assert(r[i]=='.'); assert(r[j]=='.');
              bpseq[i]=j; bpseq[j]=i;
              if (lv<n_support_parens)
              {
                r[i]=left_paren[lv];
                r[j]=right_paren[lv];
              }
              else if (lv<n_support_parens+26)
              {
                r[i]='A'+lv-n_support_parens;
                r[j]='a'+lv-n_support_parens;
              }
            }
      }
    }
    else
    {
      // make BPSEQ
      for (uint lv=0; lv!=pk_level_; ++lv)
        for (uint i=0; i<L; ++i)
          for (uint j=i+1; j<L; ++j)
            if (v[lv][i][j]>=0 && ip.get_value(v[lv][i][j])>0.5)
            {
              bpseq[i]=j; bpseq[j]=i;
            }

      // resolve the symbol of parenthsis by the graph coloring problem
      // make an adjacent graph, in which pseudoknotted base-pairs are connected.
      std::vector< std::vector<int> > g(bpseq.size());
      for (uint i=0; i!=bpseq.size(); ++i)
      {
        if (bpseq[i]<0 || bpseq[i]<=(int)i) continue;
        uint j=bpseq[i];
        for (uint k=i+1; k!=bpseq.size(); ++k)
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
      // sort vertices by degree
      std::vector<int> v;
      for (uint i=0; i!=bpseq.size(); ++i)
        if (bpseq[i]>=0 && (int)i<bpseq[i]) 
          v.push_back(i);
      std::sort(v.begin(), v.end(), cmp_by_degree(g));

      // determine colors
      std::vector<int> c(g.size(), -1);
      int max_color=0;
      for (uint i=0; i!=v.size(); ++i)
      {
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
      for (uint i=0; i!=c.size(); ++i)
        if (c[i]>=0) c[i]=rev[c[i]];

      // make the parenthsis string
      for (uint i=0; i!=bpseq.size(); ++i)
      {
        if (bpseq[i]<0 || bpseq[i]<(int)i) continue;
        uint j=bpseq[i];
        assert(r[i]=='.'); assert(r[j]=='.');
        if (c[i]<(int)n_support_parens)
        {
          r[i]=left_paren[c[i]];
          r[j]=right_paren[c[i]];
        }
        else if (c[i]<(int)n_support_parens+26)
        {
          r[i]='A'+c[i]-n_support_parens;
          r[j]='a'+c[i]-n_support_parens;
        }
      }
    }
  }

private:
  uint length(const std::string& seq) const { return seq.size(); }
  uint length(const std::list<std::string>& aln) const { return aln.front().size(); }

private:
  // options
  uint pk_level_;
  std::vector<float> th_;
  std::vector<float> alpha_;
  bool levelwise_;
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
            << " -m model: probabilistic model (default: McCaskill)" << std::endl
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
  std::vector<const char*> model;
  bool output_bpseq=false;
  int n_th=1;
  int n_fill=0;
  const char* param=NULL;
  bool aux=false;
  bool levelwise=true;
  while ((ch=getopt(argc, argv, "a:t:g:m:f:ibn:P:xuh"))!=-1)
  {
    switch (ch)
    {
      case 'm':
        model.push_back(optarg);
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
        output_bpseq=true;
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
      case 'u':
        levelwise=false;
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
    th[0]=1/(4.0+1);            // -g 4
    th[1]=1/(8.0+1);            // -g 8
  }
  if (alpha.empty())
  {
    alpha.resize(th.size());
    for (uint i=0; i!=alpha.size(); ++i)
      alpha[i]=1.0/alpha.size();
  }

  IPknot ipknot(th.size(), &th[0], &alpha[0], levelwise, !isolated_bp, n_th);
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
    if (!output_bpseq /*&& th.size()<n_support_parens*/)
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
    if (model.empty() || strcasecmp(model[0], "McCaskill")==0)
      en = new RNAfoldModel(param);
    else if (strcasecmp(model[0], "CONTRAfold")==0)
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
      if (!output_bpseq /*&& th.size()<n_support_parens*/)
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
    BPEngineAln* mix_en=NULL;
    std::vector<BPEngineSeq*> en_s;
    std::vector<BPEngineAln*> en_a;
    if (model.empty())
    {
      BPEngineSeq* e = new RNAfoldModel(param);
      en_s.push_back(e);
      en_a.push_back(new AveragedModel(e));
      en_a.push_back(new AlifoldModel(param));
      mix_en = new MixtureModel(en_a);
    }
    else
    {
      for (uint i=0; i!=model.size(); ++i)
      {
        if (strcasecmp(model[i], "McCaskill")==0)
        {
          BPEngineSeq* e = new RNAfoldModel(param);
          en_s.push_back(e);
          en_a.push_back(new AveragedModel(e));
        }
        else if (strcasecmp(model[0], "CONTRAfold")==0)
        {
          BPEngineSeq* e = new CONTRAfoldModel();
          en_s.push_back(e);
          en_a.push_back(new AveragedModel(e));
        }
        else if (strcasecmp(model[0], "Alifold")==0)
        {
          en_a.push_back(new AlifoldModel(param));
        }
        else
        {
          usage(progname);
          return 1;
        }
      }
      if (en_a.size()>1)
        mix_en = new MixtureModel(en_a);
    }

    while (!a.empty())
    {
      std::list<Aln>::iterator aln = a.begin();
      std::string consensus(aln->consensus());
      ipknot.solve(aln->seq(), (mix_en ? *mix_en : *en_a[0]), r, bpseq, n_fill);
      if (!output_bpseq /*&& th.size()<n_support_parens*/)
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

    if (mix_en) delete mix_en;
    for (uint i=0; i!=en_s.size(); ++i) delete en_s[i];
    for (uint i=0; i!=en_a.size(); ++i) delete en_a[i];
  }

  return 0;
}
