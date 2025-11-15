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
#ifdef WITH_MXFOLD2
#include "mxfold2.h"
#endif
#include "bpseq.h"

#include "cxxopts.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/stopwatch.h"

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

    // Use the new set_constraint method
    constraints.set_constraint(normalized_bp_type, count);
  }

  return constraints;
}


// Function to count base pairs in a structure
BPConstraints count_base_pairs(const std::string& seq, const VI& bpseq) {
  BPConstraints counts;

  for (uint i = 0; i < bpseq.size(); ++i) {
    if (bpseq[i] >= 0 && (int)i < bpseq[i]) {  // Only count each pair once
      int j = bpseq[i];
      std::string bp_type = normalize_base_pair_type(seq[i], seq[j]);

      // Increment count for this base pair type
      int current_count = counts.get_constraint(bp_type);
      counts.set_constraint(bp_type, current_count < 0 ? 1 : current_count + 1);
    }
  }

  return counts;
}

// Function to check if base pair counts satisfy constraints
bool satisfies_constraints(const BPConstraints& counts, const BPConstraints& bp_constraints) {
  // Check all constraint requirements
  for (const auto& [bp_type, expected_count] : bp_constraints.constraints) {
    if (expected_count >= 0) {
      int actual_count = counts.get_constraint(bp_type);
      if (actual_count < 0) actual_count = 0;  // Not found means 0
      if (actual_count != expected_count) return false;
    }
  }
  return true;
}

// Function to check if predicted structure satisfies stack constraints
bool satisfies_stack_constraints(const std::string& seq, const VI& bpseq,
                                  const StackConstraints& stack_constraints) {
  if (!stack_constraints.has_constraints()) return true;

  // For each constraint, check if at least one instance is fully satisfied
  for (size_t constraint_id = 0; constraint_id < stack_constraints.constraints.size(); ++constraint_id) {
    bool constraint_satisfied = false;

    // Check each instance of this constraint
    for (const auto& instance : stack_constraints.instances) {
      if (instance.constraint_id != (int)constraint_id) continue;

      // Check if all base pairs in this instance are present in bpseq
      bool instance_satisfied = true;
      for (const auto& [i, j] : instance.pairs) {
        if (bpseq[i] != j || bpseq[j] != i) {
          instance_satisfied = false;
          break;
        }
      }

      if (instance_satisfied) {
        constraint_satisfied = true;
        std::ostringstream pos_ss;
        for (const auto& [l, r] : instance.pairs) {
          pos_ss << "(" << l+1 << "," << r+1 << ") ";
        }
        spdlog::info("Stack constraint {} satisfied by instance: {}", constraint_id + 1, pos_ss.str());
        break;  // At least one instance is satisfied, move to next constraint
      }
    }

    if (!constraint_satisfied) {
      spdlog::warn("Stack constraint {} is NOT satisfied", constraint_id + 1);
      return false;
    }
  }

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
      spdlog::warn("invalid format base number i={}, ignored.", i);
    else switch (s[0]) {
      default:
        if (std::isdigit(s[0]))
        {
          int j = std::atoi(s.c_str());
          if (j>0 && j<=bpseq.size()) {
            bpseq[i-1] = j-1; bpseq[j-1] = i-1;
          }
          else
            spdlog::warn("invalid format base number j={}, ignored.", j);
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

// Find all instances of stack patterns in the sequence
// A stack pattern can be matched in two directions (forward or reverse)
static
void
find_stack_instances(const std::string& seq, StackConstraints& stack_constraints)
{
  stack_constraints.clear_instances();

  for (size_t constraint_id = 0; constraint_id < stack_constraints.constraints.size(); ++constraint_id) {
    const auto& constraint = stack_constraints.constraints[constraint_id];
    const size_t n = constraint.size();
    const size_t L = seq.size();

    // Try all possible positions in the sequence
    // For each position i, try to match the stack pattern starting from (i, j) pairs
    for (size_t i = 0; i + n - 1 < L; ++i) {
      for (size_t j = i + n; j < L; ++j) {
        // Try forward direction: (i, j), (i+1, j-1), (i+2, j-2), ...
        bool forward_match = true;
        StackInstance forward_instance(constraint_id);

        for (size_t k = 0; k < n; ++k) {
          size_t left = i + k;
          size_t right = j - k;

          if (left >= right) {
            forward_match = false;
            break;
          }

          std::string actual_bp = normalize_base_pair_type(seq[left], seq[right]);
          if (actual_bp != constraint.bp_types[k]) {
            forward_match = false;
            break;
          }

          forward_instance.add_pair(left, right);
        }

        if (forward_match) {
          stack_constraints.add_instance(forward_instance);
          std::ostringstream pos_ss;
          for (const auto& [l, r] : forward_instance.pairs) {
            pos_ss << "(" << l+1 << "," << r+1 << ") ";
          }
          spdlog::debug("Found stack instance (forward): {}", pos_ss.str());
        }

        // Try reverse direction: match the pattern in reverse order
        bool reverse_match = true;
        StackInstance reverse_instance(constraint_id);

        for (size_t k = 0; k < n; ++k) {
          size_t left = i + k;
          size_t right = j - k;

          if (left >= right) {
            reverse_match = false;
            break;
          }

          std::string actual_bp = normalize_base_pair_type(seq[left], seq[right]);
          // Match in reverse order
          if (actual_bp != constraint.bp_types[n - 1 - k]) {
            reverse_match = false;
            break;
          }

          reverse_instance.add_pair(left, right);
        }

        // Only add reverse instance if it's different from forward
        if (reverse_match && !forward_match) {
          stack_constraints.add_instance(reverse_instance);
          std::ostringstream pos_ss;
          for (const auto& [l, r] : reverse_instance.pairs) {
            pos_ss << "(" << l+1 << "," << r+1 << ") ";
          }
          spdlog::debug("Found stack instance (reverse): {}", pos_ss.str());
        }
      }
    }
  }

  spdlog::info("Found {} stack instances in sequence", stack_constraints.instances.size());
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

static
void
output_bpp(std::ostream& os,
            const std::string& desc, const std::string& seq,
            const VSVF& sbp)
{
  os << "# " << desc << std::endl; 
  for (uint i=1; i!=sbp.size(); ++i)
  {
    os << i << " " << seq[i-1];
    for (const auto [j, v]: sbp[i])
      if (i<j)
        os << " " << j << ":" << v;
    os << std::endl;
  }
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
  std::ostream *os_bpp=nullptr;
  std::ostream *os_mfa=nullptr;
  std::string constraint;
  std::vector<std::string> stack_constraint_args;
  uint beam_size;
  std::string input;
  bool verbose = false;
  BPConstraints bp_constraints;
  StackConstraints stack_constraints;

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
    ("bpp", "Output base-pairing probabilities",
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
    ("stack-constraint", "Specify stack constraint as space-separated base pairs (e.g., 'GC AU GU'). Can be specified multiple times for multiple constraints.",
      cxxopts::value<std::vector<std::string>>(), "\"BP1 BP2 ...\"")
    ("V,verbose", "Verbose output")
    ("loglevel", "Set the logging level (trace, debug, info, warn, error, critical)",
      cxxopts::value<std::string>()->default_value("warn"), "LEVEL")
    ("beam-size", "Beam size for LinearPartition algorithm",
      cxxopts::value<uint>()->default_value("100"), "N")
    ("base-pairs", "Specify base pair count constraints (e.g., GC=1,AU=3,GU=1,UU=1)",
      cxxopts::value<std::string>(), "CONSTRAINTS")
#ifdef WITH_MXFOLD2
    ("mxfold2-config", "config file for MXfold2 model",
      cxxopts::value<std::string>()->default_value(""), "FILE")
    ("mxfold2-gpu", "Use GPU for MXfold2 model (default: -1 for CPU)",
      cxxopts::value<int>()->default_value("-1"), "GPUID")
#endif
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
  if (res.count("stack-constraint")) stack_constraint_args = res["stack-constraint"].as<std::vector<std::string>>();
  beam_size = res["beam-size"].as<uint>();
  verbose = res["verbose"].as<bool>();
  spdlog::set_level(spdlog::level::warn); // Default log level
  if (verbose) 
    spdlog::set_level(spdlog::level::info);
  if (res.count("loglevel")) {
    const auto& loglevel = res["loglevel"].as<std::string>();
    if (loglevel == "trace") {
      spdlog::set_level(spdlog::level::trace);
    } else if (loglevel == "debug") {
      spdlog::set_level(spdlog::level::debug);
    } else if (loglevel == "info") {
      spdlog::set_level(spdlog::level::info);
    } else if (loglevel == "warn") {
      spdlog::set_level(spdlog::level::warn);
    } else if (loglevel == "error") {
      spdlog::set_level(spdlog::level::err);
    } else if (loglevel == "critical") {
      spdlog::set_level(spdlog::level::critical);
    }
  }
  input = res["input"].as<std::string>();
#ifdef WITH_MXFOLD2
  auto mxfold2_config = res["mxfold2-config"].as<std::string>();
  auto mxfold2_gpu = res["mxfold2-gpu"].as<int>();
#else
  std::string mxfold2_config;
  int mxfold2_gpu = -1;
#endif
  
  // Parse base pair constraints
  if (res.count("base-pairs")) {
    try {
      bp_constraints = parse_base_pair_constraints(res["base-pairs"].as<std::string>());
      std::ostringstream bp_ss;
      bool first = true;
      for (const auto& [bp_type, count] : bp_constraints.constraints) {
        if (!first) bp_ss << ", ";
        bp_ss << bp_type << "=" << count;
        first = false;
      }
      spdlog::info("Base pair constraints: {}", bp_ss.str());
    } catch (const std::exception& e) {
      spdlog::error("Error parsing base pair constraints: {}", e.what());
      return 1;
    }
  }

  // Parse stack constraints from command-line arguments
  if (res.count("stack-constraint")) {
    try {
      for (const auto& constraint_str : stack_constraint_args) {
        spdlog::debug("Parsing stack constraint string: '{}'", constraint_str);
        StackConstraint constraint;
        std::istringstream ss(constraint_str);
        std::string bp_type;

        // Parse space-separated base pair types
        while (ss >> bp_type) {
          spdlog::debug("Parsed bp_type: '{}'", bp_type);

          if (!bp_type.empty()) {
            try {
              std::string normalized = normalize_base_pair_type(bp_type);
              constraint.add_bp_type(normalized);
              spdlog::debug("Added normalized bp_type: '{}'", normalized);
            } catch (const std::exception& e) {
              spdlog::error("Invalid base pair type '{}': {}", bp_type, e.what());
              return 1;
            }
          }
        }

        spdlog::debug("Constraint size: {}", constraint.size());
        if (constraint.is_valid()) {
          stack_constraints.add_constraint(constraint);
          std::ostringstream bp_ss;
          for (const auto& bp : constraint.bp_types) {
            bp_ss << bp << " ";
          }
          spdlog::info("Added stack constraint: {}", bp_ss.str());
        } else {
          spdlog::error("Stack constraint '{}' must have at least 2 base pairs", constraint_str);
          return 1;
        }
      }

      if (stack_constraints.has_constraints()) {
        spdlog::info("Loaded {} stack constraints", stack_constraints.constraints.size());
      }
    } catch (const std::exception& e) {
      spdlog::error("Error parsing stack constraints: {}", e.what());
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
  if (res.count("bpp"))
  {
    auto f = res["bpp"].as<std::string>().c_str();
    os_bpp = new std::ofstream(f);
    if (!dynamic_cast<std::ofstream*>(os_bpp)->is_open())
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
        std::tie(fval, fval_pk) = ipknot.solve(seq, sbp, ep, bpseq, plevel, false, bp_constraints);
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
                                   param.empty() ? nullptr : param.c_str(), beam_size, n_th, mxfold2_config, mxfold2_gpu);
      if (!en) 
      {
        std::cout << options.help() << std::endl;
        return 1;
      }

      if (spdlog::get_level() <= spdlog::level::info)
      {
        std::ostringstream model_ss;
        std::copy(std::begin(model), std::end(model), std::ostream_iterator<std::string>(model_ss, ", "));
        spdlog::info("Model: {}", model_ss.str());
        
        std::ostringstream th_ss;
        for (auto t: th) 
        {
          th_ss << "(";
          std::copy(std::begin(t), std::end(t), std::ostream_iterator<float>(th_ss, ", "));
          th_ss << "), ";
        }
        spdlog::info("Thresholds: {}", th_ss.str());
      }

      while (!f.empty())
      {
        std::list<Fasta>::iterator fa = f.begin();

        // Find stack constraint instances in the sequence
        if (stack_constraints.has_constraints()) {
          find_stack_instances(fa->seq(), stack_constraints);
        }

        std::vector<std::vector<std::pair<uint, float>>> sbp;
        if (constraint.empty())
          sbp = en->calculate_posterior(fa->seq());
        else
        { // constraint folding
          bpseq.resize(fa->size(), BPSEQ::DOT);
          read_constraints(constraint.c_str(), bpseq);
          // int pl = IPknot::decompose_plevel(bpseq, plevel);
          // en->update_bpm(pl, fa->seq(), bpseq, plevel, sbp);
          sbp = en->calculate_posterior(fa->seq());
        }

        if (max_pfval)
          std::tie(fval, fval_pk) = ipknot.solve(fa->seq(), sbp, ep, bpseq, plevel, !constraint.empty(), bp_constraints, stack_constraints);
        else
          ipknot.solve(fa->seq(), sbp, t, bpseq, plevel, !constraint.empty(), bp_constraints, stack_constraints);

        for (int i=0; i!=n_refinement; ++i) // iterative refinement
        {
          en->update_bpm(pk_level, fa->seq(), bpseq, plevel, sbp);
          if (max_pfval)
            std::tie(fval, fval_pk) = ipknot.solve(fa->seq(), sbp, ep, bpseq, plevel, !constraint.empty(), bp_constraints, stack_constraints);
          else
            ipknot.solve(fa->seq(), sbp, t, bpseq, plevel, !constraint.empty(), bp_constraints, stack_constraints);
        }

        // Count and display base pairs if constraints are specified or verbose mode
        if (bp_constraints.has_constraints() && spdlog::get_level() <= spdlog::level::info) {
          BPConstraints actual_counts = count_base_pairs(fa->seq(), bpseq);
          std::ostringstream counts_ss;
          bool first = true;
          for (const auto& [bp_type, count] : actual_counts.constraints) {
            if (!first) counts_ss << " ";
            counts_ss << bp_type << "=" << count;
            first = false;
          }
          spdlog::info("Base pair counts in predicted structure: {}", counts_ss.str());
          if (satisfies_constraints(actual_counts, bp_constraints)) {
            spdlog::info("Base pair constraints are satisfied.");
          } else {
            spdlog::warn("Base pair constraints are NOT satisfied.");
          }
        }

        // Check stack constraints
        if (stack_constraints.has_constraints() && spdlog::get_level() <= spdlog::level::info) {
          if (satisfies_stack_constraints(fa->seq(), bpseq, stack_constraints)) {
            spdlog::info("All stack constraints are satisfied.");
          } else {
            spdlog::warn("Some stack constraints are NOT satisfied.");
          }
        }
        
        if (os_bpseq)
          output_bpseq(*os_bpseq, fa->name(), fa->seq(), bpseq, plevel, max_pfval, fval, fval_pk);
        if (os_bpseq!=&std::cout)
          output_fa(std::cout, fa->name(), fa->seq(), bpseq, plevel, output_energy);
        if (os_bpp)
          output_bpp(*os_bpp, fa->name(), fa->seq(), sbp);
        f.erase(fa);
      }
    }
    else if (Aln::load(a, input.c_str())>0)
    {
      float fval, fval_pk;
      auto en = BPEngineAln::build(model, param.empty() ? nullptr : param.c_str(), beam_size, n_th, mxfold2_config, mxfold2_gpu);
      if (!en) 
      {
        std::cout << options.help() << std::endl;
        return 1;
      }

      if (spdlog::get_level() <= spdlog::level::info)
      {
        std::ostringstream model_ss;
        std::copy(std::begin(model), std::end(model), std::ostream_iterator<std::string>(model_ss, ", "));
        spdlog::info("Model: {}", model_ss.str());
        
        std::ostringstream th_ss;
        for (auto t: th) 
        {
          th_ss << "(";
          std::copy(std::begin(t), std::end(t), std::ostream_iterator<float>(th_ss, ", "));
          th_ss << "), ";
        }
        spdlog::info("Thresholds: {}", th_ss.str());
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
          // int pl = IPknot::decompose_plevel(bpseq, plevel);
          // en->update_bpm(pl, aln->seq(), bpseq, plevel, sbp);
          sbp = en->calculate_posterior(aln->seq());
        }
        
        // For alignments, ignore bp_constraints (use empty constraints)
        BPConstraints empty_constraints;
        if (bp_constraints.has_constraints()) {
          spdlog::warn("Base pair constraints are ignored for alignment input.");
        }
        if (max_pfval)
          std::tie(fval, fval_pk) = ipknot.solve(aln->consensus(), sbp, ep, bpseq, plevel, !constraint.empty(), empty_constraints);
        else
          ipknot.solve(aln->consensus(), sbp, t, bpseq, plevel, !constraint.empty(), empty_constraints);

        for (int i=0; i!=n_refinement; ++i)
        {
          en->update_bpm(pk_level, aln->seq(), bpseq, plevel, sbp);
          if (max_pfval)
            std::tie(fval, fval_pk) = ipknot.solve(aln->consensus(), sbp, ep, bpseq, plevel, !constraint.empty(), empty_constraints);
          else
            ipknot.solve(aln->consensus(), sbp, t, bpseq, plevel, !constraint.empty(), empty_constraints);
        }

        // Count and display base pairs if constraints are specified or verbose mode
#if 0
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
#endif
        
        if (os_bpseq)
          output_bpseq(*os_bpseq, aln->name().front(), aln->consensus(), bpseq, plevel, max_pfval, fval, fval_pk);
        if (os_mfa)
          output_mfa(*os_mfa, *aln, bpseq, plevel);
        if (os_bpseq!=&std::cout && os_mfa!=&std::cout)
          output_fa(std::cout, aln->name().front(), aln->consensus(), bpseq, plevel, output_energy);
        if (os_bpp)
          output_bpp(*os_bpp, aln->name().front(), aln->consensus(), sbp);
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
  catch (std::runtime_error err)
  {
    std::cout << err.what() << std::endl;
  }

  if (os_bpseq!=&std::cout) delete os_bpseq;
  if (os_bpp) delete os_bpp;
  if (os_mfa!=&std::cout) delete os_mfa;

  return 0;
}
