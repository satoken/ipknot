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
#include <cctype>

#include "ipknot.h"
#include "ip.h"
#include "fa.h"
#include "aln.h"
#include "fold.h"
#include "nupack/nupack.h"
#include "bpseq.h"
#include "cxxopts.hpp"

// Function to normalize base pair type to canonical form
std::string normalize_base_pair_type(const std::string& bp_type) {
  std::string normalized = bp_type;
  
  // Convert to uppercase
  std::transform(normalized.begin(), normalized.end(), normalized.begin(), ::toupper);
  
  // Normalize to canonical order
  if (normalized == "GC" || normalized == "CG") {
    return "GC";
  } else if (normalized == "AU" || normalized == "AT" || normalized == "UA" || normalized == "TA") {
    return "AU";
  } else if (normalized == "GU" || normalized == "GT" || normalized == "UG" || normalized == "TG") {
    return "GU";
  } else if (normalized == "UU" || normalized == "TT") {
    return "UU";
  } else {
    throw std::invalid_argument("Unknown base pair type: " + bp_type);
  }
}

// Function to parse base pair constraints from string
BPConstraints parse_base_pair_constraints(const std::string& str) {
  BPConstraints constraints;
  
  std::istringstream ss(str);
  std::string pair;
  
  while (std::getline(ss, pair, ',')) {
    // Remove whitespace
    pair.erase(std::remove_if(pair.begin(), pair.end(), ::isspace), pair.end());
    
    size_t eq_pos = pair.find('=');
    if (eq_pos == std::string::npos) {
      throw std::invalid_argument("Invalid base pair constraint format: " + pair);
    }
    
    std::string bp_type = pair.substr(0, eq_pos);
    int count = std::stoi(pair.substr(eq_pos + 1));
    
    // Normalize base pair type
    std::string normalized_bp_type = normalize_base_pair_type(bp_type);
    
    if (normalized_bp_type == "GC") {
      constraints.GC = count;
    } else if (normalized_bp_type == "AU") {
      constraints.AU = count;
    } else if (normalized_bp_type == "GU") {
      constraints.GU = count;
    } else if (normalized_bp_type == "UU") {
      constraints.UU = count;
    }
  }
  
  return constraints;
}

// Function to determine base pair type
std::string get_base_pair_type(char a, char b) {
  // Normalize to uppercase and canonical order
  if (a > b) std::swap(a, b);
  a = std::toupper(a);
  b = std::toupper(b);
  
  if ((a == 'G' && b == 'C') || (a == 'C' && b == 'G')) {
    return "GC";
  } else if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A')) {
    return "AU";
  } else if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G')) {
    return "GU";
  } else if (a == 'U' && b == 'U') {
    return "UU";
  }
  return "UNKNOWN";
}

// Function to count base pairs in a structure
BPConstraints count_base_pairs(const std::string& seq, const VI& bpseq) {
  BPConstraints counts;
  counts.GC = counts.AU = counts.GU = counts.UU = 0;
  
  for (uint i = 0; i < bpseq.size(); ++i) {
    if (bpseq[i] >= 0 && (int)i < bpseq[i]) {  // Only count each pair once
      int j = bpseq[i];
      std::string bp_type = get_base_pair_type(seq[i], seq[j]);
      
      if (bp_type == "GC") {
        counts.GC++;
      } else if (bp_type == "AU") {
        counts.AU++;
      } else if (bp_type == "GU") {
        counts.GU++;
      } else if (bp_type == "UU") {
        counts.UU++;
      }
    }
  }
  
  return counts;
}

// Function to check if base pair counts satisfy constraints
bool satisfies_constraints(const BPConstraints& counts, const BPConstraints& constraints) {
  if (constraints.GC >= 0 && counts.GC != constraints.GC) return false;
  if (constraints.AU >= 0 && counts.AU != constraints.AU) return false;
  if (constraints.GU >= 0 && counts.GU != constraints.GU) return false;
  if (constraints.UU >= 0 && counts.UU != constraints.UU) return false;
  return true;
}

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
  BPConstraints bp_constraints;

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
    ("base-pairs", "Specify base pair count constraints (e.g., GC=1,AU=3,GU=1,UU=1)",
      cxxopts::value<std::string>(), "CONSTRAINTS")
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
  
  // Parse base pair constraints
  if (res.count("base-pairs")) {
    try {
      bp_constraints = parse_base_pair_constraints(res["base-pairs"].as<std::string>());
      if (verbose) {
        std::cerr << "Base pair constraints: ";
        if (bp_constraints.GC >= 0) std::cerr << "GC=" << bp_constraints.GC << " ";
        if (bp_constraints.AU >= 0) std::cerr << "AU=" << bp_constraints.AU << " ";
        if (bp_constraints.GU >= 0) std::cerr << "GU=" << bp_constraints.GU << " ";
        if (bp_constraints.UU >= 0) std::cerr << "UU=" << bp_constraints.UU << " ";
        std::cerr << std::endl;
      }
    } catch (const std::exception& e) {
      std::cerr << "Error parsing base pair constraints: " << e.what() << std::endl;
      return 1;
    }
  }

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
        std::tie(fval, fval_pk) = ipknot.solve(seq, sbp, ep, bpseq, plevel, false, verbose, bp_constraints);
      else
        ipknot.solve(seq, sbp, t, bpseq, plevel, false, bp_constraints);
      if (os_bpseq)
        output_bpseq(*os_bpseq, input.c_str(), seq, bpseq, plevel, max_pfval, fval, fval_pk);
      if (os_bpseq!=&std::cout)
        output_fa(std::cout, input.c_str(), seq, bpseq, plevel, output_energy);
    }
    else if (Fasta::load(f, input.c_str())>0)
    {
      float fval, fval_pk;
      auto en = BPEngineSeq::build(model.empty() ? nullptr : model[0].c_str(), 
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
          en->update_bpm(pl, fa->seq(), bpseq, plevel, sbp);
        }

        if (max_pfval)
          std::tie(fval, fval_pk) = ipknot.solve(fa->seq(), sbp, ep, bpseq, plevel, !constraint.empty(), verbose, bp_constraints);
        else
          ipknot.solve(fa->seq(), sbp, t, bpseq, plevel, !constraint.empty(), bp_constraints);

        for (int i=0; i!=n_refinement; ++i) // iterative refinement
        {
          en->update_bpm(pk_level, fa->seq(), bpseq, plevel, sbp);
          if (max_pfval)
            std::tie(fval, fval_pk) = ipknot.solve(fa->seq(), sbp, ep, bpseq, plevel, !constraint.empty(), verbose, bp_constraints);
          else
            ipknot.solve(fa->seq(), sbp, t, bpseq, plevel, !constraint.empty(), bp_constraints);
        }

        // Count and display base pairs if constraints are specified or verbose mode
        if (bp_constraints.has_constraints() || verbose) {
          BPConstraints actual_counts = count_base_pairs(fa->seq(), bpseq);
          if (verbose || bp_constraints.has_constraints()) {
            std::cerr << "Base pair counts in predicted structure: ";
            std::cerr << "GC=" << actual_counts.GC << " ";
            std::cerr << "AU=" << actual_counts.AU << " ";
            std::cerr << "GU=" << actual_counts.GU << " ";
            std::cerr << "UU=" << actual_counts.UU << std::endl;
          }
          
          if (bp_constraints.has_constraints()) {
            if (satisfies_constraints(actual_counts, bp_constraints)) {
              std::cerr << "✓ Base pair constraints are satisfied." << std::endl;
            } else {
              std::cerr << "✗ Base pair constraints are NOT satisfied." << std::endl;
            }
          }
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
      auto en = BPEngineAln::build(model, param.empty() ? nullptr : param.c_str(), beam_size);
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
          en->update_bpm(pl, aln->seq(), bpseq, plevel, sbp);
        }
        
        // For alignments, ignore bp_constraints (use empty constraints)
        BPConstraints empty_constraints;
        if (bp_constraints.has_constraints() && verbose) {
          std::cerr << "Warning: Base pair constraints are ignored for alignment input." << std::endl;
        }
        if (max_pfval)
          std::tie(fval, fval_pk) = ipknot.solve(aln->consensus(), sbp, ep, bpseq, plevel, !constraint.empty(), verbose, empty_constraints);
        else
          ipknot.solve(aln->consensus(), sbp, t, bpseq, plevel, !constraint.empty(), empty_constraints);

        for (int i=0; i!=n_refinement; ++i)
        {
          en->update_bpm(pk_level, aln->seq(), bpseq, plevel, sbp);
          if (max_pfval)
            std::tie(fval, fval_pk) = ipknot.solve(aln->consensus(), sbp, ep, bpseq, plevel, !constraint.empty(), verbose, empty_constraints);
          else
            ipknot.solve(aln->consensus(), sbp, t, bpseq, plevel, !constraint.empty(), empty_constraints);
        }

        // Count and display base pairs if constraints are specified or verbose mode
        if (bp_constraints.has_constraints() || verbose) {
          BPConstraints actual_counts = count_base_pairs(aln->consensus(), bpseq);
          if (verbose || bp_constraints.has_constraints()) {
            std::cerr << "Base pair counts in predicted structure: ";
            std::cerr << "GC=" << actual_counts.GC << " ";
            std::cerr << "AU=" << actual_counts.AU << " ";
            std::cerr << "GU=" << actual_counts.GU << " ";
            std::cerr << "UU=" << actual_counts.UU << std::endl;
          }
          
          if (bp_constraints.has_constraints()) {
            if (satisfies_constraints(actual_counts, bp_constraints)) {
              std::cerr << "✓ Base pair constraints are satisfied." << std::endl;
            } else {
              std::cerr << "✗ Base pair constraints are NOT satisfied." << std::endl;
            }
          }
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