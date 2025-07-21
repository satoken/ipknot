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
#include <strings.h>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <memory>

#include "ipknot.h"
#include "ip.h"
#include "fa.h"
#include "aln.h"
#include "fold.h"
#include "nupack/nupack.h"
#include "bpseq.h"
#include "cxxopts.hpp"

std::string
make_parenthsis(const VI& bpseq, const VI& plevel)
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

template < class SEQ, class EN >
void
update_bpm(uint pk_level, const SEQ& seq, EN& en,
           const VI& bpseq, const VI& plevel, VF& bp, VI& offset)
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
update_bpm(uint pk_level, const SEQ& seq, EN& en, const VI& bpseq, const VI& plevel, VSVF& sbp)
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
                            [&, &jl=jl](const auto& x) { return x.first == jl; });
        if (re != std::end(sbp[i]))
          re->second += v;
        else
          sbp[i].emplace_back(jl, v);
      }
    }
  }

  VSVF sbp_temp(L+1);
  for (uint i=1; i!=sbp.size(); i++)
    for (const auto [j, v]: sbp[i]) 
      if (v>=DEFAULT_THRESHOLD) 
        sbp_temp[i].emplace_back(j, v);
  std:swap(sbp, sbp_temp);
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
          const VI& bpseq, const VI& plevel, bool output_energy)
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
output_mfa(std::ostream& os, const Aln& aln, const VI& bpseq, const VI& plevel)
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
             const VI& bpseq, const VI& plevel, 
             bool max_pfval, float fval, float fval_pk)
{
  os << "# " << desc; 
  if (max_pfval)
    os << " (max pF=" << fval << "," << fval_pk << ")";
  os << std::endl;
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
  if (model==nullptr || strcasecmp(model, "McCaskill")==0 || strcasecmp(model, "Boltzmann")==0)
    en = std::make_unique<RNAfoldModel>(param);
  else if (strcasecmp(model, "ViennaRNA")==0)
    en = std::make_unique<RNAfoldModel>("default");
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
  else if (strcasecmp(model, "LinearPartition-C")==0 || strcasecmp(model, "lpc")==0)
    en = std::make_unique<LinearPartitionModel>(false, beam_size);
  else if (strcasecmp(model, "LinearPartition-V")==0 || strcasecmp(model, "lpv")==0)
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
      if (strcasecmp(m, "McCaskill")==0 || strcasecmp(m, "Boltzmann")==0)
      {
        auto e = std::make_unique<RNAfoldModel>(param);
        en_a.push_back(std::make_unique<AveragedModel>(std::move(e)));
      }
      else if (strcasecmp(m, "ViennaRNA")==0)
      {
        auto e = std::make_unique<RNAfoldModel>("default");
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
      else if (strcasecmp(m, "LinearPartition-C")==0 || strcasecmp(m, "lpc")==0)
      {
        auto e = std::make_unique<LinearPartitionModel>(false, beam_size);
        en_a.push_back(std::make_unique<AveragedModel>(std::move(e)));
      }
      else if (strcasecmp(m, "LinearPartition-V")==0 || strcasecmp(m, "lpv")==0)
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
  bool verbose = false;

  cxxopts::Options options{progname, format("IPknot version %s", PACKAGE_VERSION)};
  options.add_options()
    ("input", "FASTA-formatted file or ALN-formatted file",
      cxxopts::value<std::string>(), "FASTA_OR_ALN")
    ("e,model", "Probabilistic model",
      cxxopts::value<std::vector<std::string>>()->default_value("LinearPartition-C"), "MODEL")
    ("r,refinement", "The number of the iterative refinement",
      cxxopts::value<int>()->default_value("1"), "N")
#if 0
    ("a,alpha", "The weight for each level",
      cxxopts::value<std::vector<float>>(), "ALPHA")
#endif
    ("t,threshold", "The threshold of base-pairing probabilities for each level",
      cxxopts::value<std::vector<std::string>>()->default_value("auto,auto"), "TH")
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
    ("u,no-levelwise", "Do not perform the levelwise prediction",
      cxxopts::value<bool>()->default_value("false"))
    ("l,mfa", "Output the prediction with the given mulple alignment",
      cxxopts::value<bool>()->default_value("false"))
    ("L,mfa-file", "Output file for the multiple alignment",
      cxxopts::value<std::string>(), "FILE")
    ("E,energy", "Output with the free energy",
      cxxopts::value<bool>()->default_value("false"))
    ("c,constraint", "Specify the structure constraint by a BPSEQ formatted file",
      cxxopts::value<std::string>(), "FILE")
    ("V,verbose", "Verbose output")
    ("beam-size", "Beam size for LinearPartition algorithm",
      cxxopts::value<uint>()->default_value("100"), "N")
    ("version", "Print version")
    ("h,help", "Print usage"); 
  options.parse_positional({"input"});
  options
    .positional_help("FASTA_OR_ALN")
    .show_positional_help();

  auto res = options.parse(argc, argv);
  if (res.count("version")) 
  {
    std::cout << format("IPknot version %s", PACKAGE_VERSION) << std::endl;
    exit(0);
  }
  if (res.count("help") || res.count("input")==0)
  {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  model = res["model"].as<std::vector<std::string>>();
  n_refinement = res["refinement"].as<int>();
  isolated_bp = res["allow-isolated"].as<bool>();
  if (res.count("param")) param = res["param"].as<std::string>();
  aux = res["aux"].as<bool>();
  levelwise = !res["no-levelwise"].as<bool>();
#ifndef WITH_GLPK
  n_th = res["threads"].as<uint>();
#endif
  output_energy = res["energy"].as<bool>();
  if (res.count("constraint")) constraint = res["constraint"].as<std::string>();
  beam_size = res["beam-size"].as<uint>();
  verbose = res["verbose"].as<bool>();
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

  if (res.count("gamma"))
  {
    const auto& arg_gamma = res["gamma"].as<std::vector<std::string>>();
    th.resize(arg_gamma.size());
    for (uint i=0; i!=th.size(); ++i)
    {
      if (arg_gamma[i]=="auto")
      {
        th[i] = VF{0.5, 0.25, 0.125, 0.0625};
        max_pfval = true;
      }
      else
      {
        auto temp = parse_csv_line<float>(arg_gamma[i].c_str(), '_');
        th[i].resize(temp.size());
        std::transform(std::cbegin(temp), std::cend(temp), std::begin(th[i]),
          [&](auto v) { return 1./(v+1.); });
        if (th[i].size()>1) max_pfval = true;
      }
    }
  }
  else
  {
    const auto& arg_th = res["threshold"].as<std::vector<std::string>>();
    th.resize(arg_th.size());
    for (uint i=0; i!=th.size(); ++i)
    {
      if (arg_th[i]=="auto") 
      {
        th[i] = VF{0.5, 0.25, 0.125, 0.0625};
        max_pfval = true;
      }
      else
      {
        th[i] = parse_csv_line<float>(arg_th[i].c_str(), '_');
        if (th[i].size()>1) max_pfval = true;
      }
    }
  }

#if 0
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
#endif

#if 0
  if (res.count("alpha"))
  {
    alpha = res["alpha"].as<std::vector<float>>();
  }
  else
  {
#endif
    alpha.resize(th.size());
    std::fill(std::begin(alpha), std::end(alpha), 1./alpha.size());
#if 0
  }
#endif
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
      AuxModel aux;
      std::string seq;
      float fval, fval_pk;
      auto sbp = aux.calculate_posterior(input.c_str(), seq);
      if (max_pfval)
        std::tie(fval, fval_pk) = ipknot.solve(seq.size(), sbp, ep, bpseq, plevel, false, verbose);
      else
        ipknot.solve(seq.size(), sbp, t, bpseq, plevel, false);
      if (os_bpseq)
        output_bpseq(*os_bpseq, input.c_str(), seq, bpseq, plevel, max_pfval, fval, fval_pk);
      if (os_bpseq!=&std::cout)
        output_fa(std::cout, input.c_str(), seq, bpseq, plevel, output_energy);
    }
    else if (Fasta::load(f, input.c_str())>0)
    {
      float fval, fval_pk;
      auto en = build_engine_seq(model.empty() ? nullptr : model[0].c_str(), 
                                  param.empty() ? nullptr : param.c_str(), beam_size);
      if (!en) 
      {
        std::cout << options.help() << std::endl;
        return 1;
      }

      if (verbose)
      {
        std::cerr << "Model: ";
        std::copy(std::begin(model), std::end(model), std::ostream_iterator<std::string>(std::cerr, ", "));
        std::cerr << std::endl;
        std::cerr << "Thresholds: ";
        for (auto t: th) 
        {
          std::cerr << "(";
          std::copy(std::begin(t), std::end(t), std::ostream_iterator<float>(std::cerr, ", "));
          std::cerr << "), ";
        }
        std::cerr << std::endl << std::endl;
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
          std::tie(fval, fval_pk) = ipknot.solve(fa->size(), sbp, ep, bpseq, plevel, !constraint.empty(), verbose);
        else
          ipknot.solve(fa->size(), sbp, t, bpseq, plevel, !constraint.empty());

        for (int i=0; i!=n_refinement; ++i) // iterative refinement
        {
          update_bpm(pk_level, fa->seq(), *en, bpseq, plevel, sbp);
          if (max_pfval)
            std::tie(fval, fval_pk) = ipknot.solve(fa->size(), sbp, ep, bpseq, plevel, !constraint.empty(), verbose);
          else
            ipknot.solve(fa->size(), sbp, t, bpseq, plevel, !constraint.empty());
        }

        if (os_bpseq)
          output_bpseq(*os_bpseq, fa->name(), fa->seq(), bpseq, plevel, max_pfval, fval, fval_pk);
        if (os_bpseq!=&std::cout)
          output_fa(std::cout, fa->name(), fa->seq(), bpseq, plevel, output_energy);
        f.erase(fa);
      }
    }
    else if (Aln::load(a, input.c_str())>0)
    {
      float fval, fval_pk;
      auto en = build_engine_aln(model, param.empty() ? nullptr : param.c_str(), beam_size);
      if (!en) 
      {
        std::cout << options.help() << std::endl;
        return 1;
      }

      if (verbose)
      {
        std::cerr << "Model: ";
        std::copy(std::begin(model), std::end(model), std::ostream_iterator<std::string>(std::cerr, ", "));
        std::cerr << std::endl;
        std::cerr << "Thresholds: ";
        for (auto t: th) 
        {
          std::cerr << "(";
          std::copy(std::begin(t), std::end(t), std::ostream_iterator<float>(std::cerr, ", "));
          std::cerr << "), ";
        }
        std::cerr << std::endl << std::endl;
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
          std::tie(fval, fval_pk) = ipknot.solve(aln->size(), sbp, ep, bpseq, plevel, !constraint.empty(), verbose);
        else
          ipknot.solve(aln->size(), sbp, t, bpseq, plevel, !constraint.empty());

        for (int i=0; i!=n_refinement; ++i)
        {
          update_bpm(pk_level, aln->seq(), *en, bpseq, plevel, sbp);
          if (max_pfval)
            std::tie(fval, fval_pk) = ipknot.solve(aln->size(), sbp, ep, bpseq, plevel, !constraint.empty(), verbose);
          else
            ipknot.solve(aln->size(), sbp, t, bpseq, plevel, !constraint.empty());
        }

        if (os_bpseq)
          output_bpseq(*os_bpseq, aln->name().front(), aln->consensus(), bpseq, plevel, max_pfval, fval, fval_pk);
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