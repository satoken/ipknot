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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iterator>
#include <algorithm>

#include "ip.h"
#include "fa.h"
#include "aln.h"
#include "fold.h"
#include "nupack/nupack.h"
#include "bpseq.h"

typedef unsigned int uint;
typedef std::vector<float> VF;
typedef std::vector<VF> VVF;
typedef std::vector<int> VI;
typedef std::vector<VI> VVI;
typedef std::vector<VVI> VVVI;

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

  void solve(uint L, const std::vector<float>& bp, const std::vector<int>& offset,
             const std::vector<float>& th, std::vector<int>& bpseq, std::vector<int>& plevel, bool constraint) const
  {
    IP ip(IP::MAX, n_th_);
    VVVI v(pk_level_, VVI(L, VI(L, -1)));
    VVVI w(pk_level_, VVI(L));
    VI c_l(L, 0), c_r(L, 0);

    if (!constraint)
    {
      bpseq.resize(L);
      std::fill(bpseq.begin(), bpseq.end(), -2);
    }
    
    // make objective variables with their weights
    for (uint j=1; j!=L; ++j)
    {
      for (uint i=j-1; i!=-1u; --i)
      {
        const float& p=bp[offset[i+1]+(j+1)];
        for (uint lv=0; lv!=pk_level_; ++lv)
          if (p>th[lv])
          {
            v[lv][i][j] = ip.make_variable(p*alpha_[lv]);
            w[lv][i].push_back(j);
            c_l[i]++; c_r[j]++;
          }
      }
    }
    ip.update();

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
          //std::cout << "L " << i << " " << c_l[i] << std::endl;
          if (c_l[i]>0) 
          {
            row_l = ip.make_constraint(IP::FX, 1, 1);
            row_r = ip.make_constraint(IP::UP, 0, 0);
          }
          break;
        case BPSEQ::R: // paired with left j
          //std::cout << "R " << i << " " << c_l[i] << std::endl;
          if (c_l[i]>0) 
          {
            row_l = ip.make_constraint(IP::UP, 0, 0);
            row_r = ip.make_constraint(IP::FX, 1, 1);
          }
          break;
      }
      if (row_l<0 || row_r<0)
      {
        std::cerr << "invalid constraint for the base " << i+1 << std::endl;
        row_l = row_r = ip.make_constraint(IP::UP, 0, 1); // fallback to no constraint
      }
      
      for (uint lv=0; lv!=pk_level_; ++lv)
      {
        for (uint j=0; j<i; ++j)
          if (v[lv][j][i]>=0)
            ip.add_constraint(row_r, v[lv][j][i], 1);
        for (uint j=i+1; j<L; ++j)
          if (v[lv][i][j]>=0)
            ip.add_constraint(row_l, v[lv][i][j], 1);
      }

      if (bpseq[i]>=0 && i<bpseq[i]) // paired with j=bpseq[i]
      {
        int c=0;
        for (uint lv=0; lv!=pk_level_; ++lv)
          if (v[lv][i][bpseq[i]]>=0) c++;
        if (c>0)
        {
          int row = ip.make_constraint(IP::FX, 1, 1);
          for (uint lv=0; lv!=pk_level_; ++lv)
            if (v[lv][i][bpseq[i]]>=0)  
              ip.add_constraint(row, v[lv][i][bpseq[i]], 1);
        }
        else
          std::cerr << "invalid constraint for the bases " << i+1 << " and " << bpseq[i]+1 << std::endl;
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

    // build the result
    bpseq.resize(L);
    std::fill(bpseq.begin(), bpseq.end(), -1);
    plevel.resize(L);
    std::fill(plevel.begin(), plevel.end(), -1);
    for (uint lv=0; lv!=pk_level_; ++lv)
    {
      for (uint i=0; i<L; ++i)
        for (uint j=i+1; j<L; ++j)
          if (v[lv][i][j]>=0 && ip.get_value(v[lv][i][j])>0.5)
          {
            bpseq[i]=j; bpseq[j]=i;
            plevel[i]=plevel[j]=lv;
          }
    }

    if (!levelwise_)
      decompose_plevel(bpseq, plevel);
  }

  void solve(uint L, const std::vector<float>& bp, const std::vector<int>& offset,
             EnumParam<float>& ep, std::vector<int>& bpseq, std::vector<int>& plevel, bool constraint) const
  {
    std::vector<float> th(ep.size());
    std::vector<int> bpseq_temp;
    std::vector<int> plevel_temp;
    float max_mcc=-100.0;
    std::cerr << "Search for the best thresholds by pseudo expected MCC:" << std::endl;
    do {
      ep.get(th);
      std::cerr << "th=";
      std::copy(th.begin(), th.end(), std::ostream_iterator<float>(std::cerr, ","));
      solve(L, bp, offset, th, bpseq_temp, plevel_temp, constraint);
      float sen, ppv, mcc;
      compute_expected_accuracy(bpseq_temp, bp, offset, sen, ppv, mcc);
      std::cerr << " pMCC=" << mcc << std::endl;
      if (mcc>max_mcc)
      {
        max_mcc = mcc;
        bpseq = bpseq_temp;
        plevel = plevel_temp;
      }
    } while (!ep.succ());
    std::cerr << "max pMCC=" << max_mcc << std::endl << std::endl;
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

  static void
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
  }

  static void
  compute_expected_accuracy(const std::vector<int>& bpseq,
                            const std::vector<float>& bp, const std::vector<int>& offset,
                            float& sen, float& ppv, float& mcc)
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

    sen = ppv = mcc = 0;
    if (etp+efn!=0) sen = etp / (etp + efn);
    if (etp+efp!=0) ppv = etp / (etp + efp);
    if (etp+efp!=0 && etp+efn!=0 && etn+efp!=0 && etn+efn!=0)
      mcc = (etp*etn-efp*efn) / std::sqrt((etp+efp)*(etp+efn)*(etn+efp)*(etn+efn));
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
  //bp.resize(L);
  std::fill(bp.begin(), bp.end(), 0.0);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
  
  std::vector<float> bpl;
  std::vector<int> offsetl;
  for (uint l=0; l!=pk_level; ++l)
  {
    // make the constraint string
    std::string str(L, '?');
    for (uint i=0; i!=bpseq.size(); ++i)
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

    // re-folding the seq with the constraint
    std::fill(bpl.begin(), bpl.end(), 0.0);
    en.calculate_posterior(seq, str, bpl, offsetl);
    assert(bp.size()==bpl.size());
    // update the base-pairing probability matrix
    for (uint k=0; k!=bp.size(); ++k) bp[k]+=bpl[k];
  }
#ifndef NDEBUG
  for (uint k=0; k!=bp.size(); ++k) assert(bp[k]<=1.0);
#endif
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
parse_csv_line(const char* l)
{
  std::string s;
  std::vector<T> r;
  std::istringstream ss(l);
  while (std::getline(ss, s, ','))
    r.push_back(atof(s.c_str()));
  return r;
}

void
usage(const char* progname)
{
  std::cout << "IPknot version " << PACKAGE_VERSION << std::endl
            << "  Available probabilistic models: McCaskill, CONTRAfold, Alifold, NUPACK"
            << std::endl << std::endl;

  std::cout << progname << ": [options] fasta" << std::endl
            << ""
            << " -h:       show this message" << std::endl
    //      << " -a alpha: weight for each level" << std::endl
            << " -t th:    threshold of base-pairing probabilities for each level" << std::endl
            << " -g gamma: weight for true base-pairs equivalent to -t 1/(gamma+1)" << std::endl
            << "           (default: -g 2 -g 4)" << std::endl
    //      << " -m:       select thresholds that maxmize pseudo MCC" << std::endl
            << " -e model: probabilistic model (default: McCaskill)" << std::endl
            << " -r n:     the number of the iterative refinement (default: 0)" << std::endl
            << " -i:       allow isolated base-pairs" << std::endl
            << " -b:       output the prediction by BPSEQ format" << std::endl
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
  uint pk_level=0;
  char ch;
  std::vector< std::vector<float> > th;
  std::vector<float> alpha;
  bool isolated_bp=false;
  std::vector<const char*> model;
  int n_th=1;
  int n_refinement=0;
  const char* param=NULL;
  bool aux=false;
  bool levelwise=true;
  bool max_pmcc=false;
  bool output_energy=false;
  std::ostream *os_bpseq=NULL;
  std::ostream *os_mfa=NULL;
  std::string constraint;
  while ((ch=getopt(argc, argv, "a:t:g:me:f:r:ibB:n:P:xulL:Ec:h"))!=-1)
  {
    switch (ch)
    {
      case 'e':
        model.push_back(optarg);
        break;
      case 'f': case 'r':
        n_refinement=atoi(optarg);
        break;
      case 'a':
        alpha.push_back(atof(optarg));
        break;
      case 'm':
        max_pmcc=true;
        break;
      case 't':
        th.push_back(parse_csv_line<float>(optarg));
        break;
      case 'g':
      {
        std::vector<float> temp = parse_csv_line<float>(optarg);
        for (uint i=0; i!=temp.size(); ++i)
          temp[i] = 1/(temp[i]+1);
        th.push_back(temp);
        break;
      }
      case 'i':
        isolated_bp=true;
        break;
      case 'b':
        os_bpseq=&std::cout;
        break;
      case 'B':
        os_bpseq=new std::ofstream(optarg);
        if (!dynamic_cast<std::ofstream*>(os_bpseq)->is_open())
        {
          perror(optarg);
          return 1;
        }
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
      case 'l':
        os_mfa=&std::cout;
        break;
      case 'L':
        os_mfa=new std::ofstream(optarg);
        if (!dynamic_cast<std::ofstream*>(os_mfa)->is_open())
        {
          perror(optarg);
          return 1;
        }
        break;
      case 'E':
        output_energy = true;
        break;
      case 'c':
        constraint = optarg;
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
  if (alpha.empty())
  {
    alpha.resize(th.size());
    for (uint i=0; i!=alpha.size(); ++i)
      alpha[i]=1.0/alpha.size();
  }
  pk_level=alpha.size();

  try
  {
    IPknot ipknot(pk_level, &alpha[0], levelwise, !isolated_bp, n_th);
    std::vector<float> bp;
    std::vector<int> offset;
    std::vector<int> bpseq;
    std::vector<int> plevel;

    IPknot::EnumParam<float> ep(th);
    std::vector<float> t(th.size());
    ep.get(t);

    std::list<Fasta> f;
    std::list<Aln> a;

    if (!constraint.empty()) {
      BPSEQ c;
      c.load(constraint.c_str());
      bpseq.resize(c.seq().size());
      for (uint i=0; i!=bpseq.size(); i++) {
        bpseq[i] = c.bp()[i+1];
        if (bpseq[i]>=0) bpseq[i]--;
      }
    }

    if (aux)
    {
      AuxModel aux;
      std::string seq;
      aux.calculate_posterior(argv[0], seq, bp, offset);
      if (max_pmcc)
        ipknot.solve(seq.size(), bp, offset, ep, bpseq, plevel, !constraint.empty());
      else
        ipknot.solve(seq.size(), bp, offset, t, bpseq, plevel, !constraint.empty());
      if (os_bpseq)
        output_bpseq(*os_bpseq, argv[0], seq, bpseq, plevel);
      if (os_bpseq!=&std::cout)
        output_fa(std::cout, argv[0], seq, bpseq, plevel, output_energy);
    }
    else if (Fasta::load(f, argv[0])>0)
    {
      BPEngineSeq* en=NULL;
      if (model.empty() || strcasecmp(model[0], "McCaskill")==0)
        en = new RNAfoldModel(param);
      else if (strcasecmp(model[0], "CONTRAfold")==0)
        en = new CONTRAfoldModel();
#if 0
      else if (strcasecmp(model[0], "nupack")==0)
        if (param)
          en = new NupackModel(param);
        else
          en = new NupackModel(2);
      else if (strcasecmp(model[0], "nupack03")==0)
        en = new NupackModel(0);
      else if (strcasecmp(model[0], "nupack09")==0)
        en = new NupackModel(1);
#else
      else if (strcasecmp(model[0], "nupack")==0)
        en = new NupackModel(param);
#endif
      else
      {
        usage(progname);
        return 1;
      }

      while (!f.empty())
      {
        std::list<Fasta>::iterator fa = f.begin();
        en->calculate_posterior(fa->seq(), bp, offset);
        if (max_pmcc)
          ipknot.solve(fa->size(), bp, offset, ep, bpseq, plevel, !constraint.empty());
        else
          ipknot.solve(fa->size(), bp, offset, t, bpseq, plevel, !constraint.empty());
        for (int i=0; i!=n_refinement; ++i)
        {
          update_bpm(pk_level, fa->seq(), *en, bpseq, plevel, bp, offset);
          if (max_pmcc)
            ipknot.solve(fa->size(), bp, offset, ep, bpseq, plevel, !constraint.empty());
          else
            ipknot.solve(fa->size(), bp, offset, t, bpseq, plevel, !constraint.empty());
        }
        if (os_bpseq)
          output_bpseq(*os_bpseq, fa->name(), fa->seq(), bpseq, plevel);
        if (os_bpseq!=&std::cout)
          output_fa(std::cout, fa->name(), fa->seq(), bpseq, plevel, output_energy);
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
      BPEngineAln* en= mix_en ? mix_en : en_a[0];

      while (!a.empty())
      {
        std::list<Aln>::iterator aln = a.begin();
        en->calculate_posterior(aln->seq(), bp, offset);
        if (max_pmcc)
          ipknot.solve(aln->size(), bp, offset, ep, bpseq, plevel, !constraint.empty());
        else
          ipknot.solve(aln->size(), bp, offset, t, bpseq, plevel, !constraint.empty());
        for (int i=0; i!=n_refinement; ++i)
        {
          update_bpm(pk_level, aln->seq(), *en, bpseq, plevel, bp, offset);
          if (max_pmcc)
            ipknot.solve(aln->size(), bp, offset, ep, bpseq, plevel, !constraint.empty());
          else
            ipknot.solve(aln->size(), bp, offset, t, bpseq, plevel, !constraint.empty());
        }
        if (os_bpseq)
          output_bpseq(*os_bpseq, aln->name().front(), aln->consensus(), bpseq, plevel);
        if (os_mfa)
          output_mfa(*os_mfa, *aln, bpseq, plevel);
        if (os_bpseq!=&std::cout && os_mfa!=&std::cout)
          output_fa(std::cout, aln->name().front(), aln->consensus(), bpseq, plevel, output_energy);
        a.erase(aln);
      }

      if (mix_en) delete mix_en;
      for (uint i=0; i!=en_s.size(); ++i) delete en_s[i];
      for (uint i=0; i!=en_a.size(); ++i) delete en_a[i];
    }
    else
    {
      throw (std::string(argv[0])+": Format error").c_str();
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
