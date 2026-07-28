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
#include <cmath>
#include <functional>
#include <set>

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
    const std::string count_text = pair.substr(eq_pos + 1);
    if (bp_type.empty() || count_text.empty()) {
      throw std::invalid_argument("Invalid base pair constraint format: " + pair);
    }
    size_t parsed_chars = 0;
    int count = std::stoi(count_text, &parsed_chars);
    if (parsed_chars != count_text.size()) {
      throw std::invalid_argument("Invalid base pair count: " + count_text);
    }
    if (count < 0) {
      throw std::invalid_argument("Base pair count must be non-negative: " + pair);
    }

    // Normalize base pair type
    std::string normalized_bp_type = normalize_base_pair_type(bp_type);

    if (constraints.get_constraint(normalized_bp_type) >= 0) {
      throw std::invalid_argument("Duplicate base pair constraint: " + bp_type);
    }

    constraints.set_constraint(normalized_bp_type, count);
  }

  return constraints;
}

// Parse either the regular stacking notation ("GU GC AU" / "GU,GC,AU")
// or the compact NMR notation ("(GU)-G-U").  In the NMR notation the
// parenthesized token is an explicitly observed base pair, while each bare
// G or U denotes its canonical pair (G-C or U-A, respectively).
StackConstraint parse_stack_constraint(const std::string& str) {
  StackConstraint constraint;
  const bool nmr_notation = str.find('(') != std::string::npos ||
                            str.find(')') != std::string::npos ||
                            str.find('-') != std::string::npos;

  if (!nmr_notation) {
    std::string normalized_str = str;
    std::replace(normalized_str.begin(), normalized_str.end(), ',', ' ');
    std::istringstream ss(normalized_str);
    std::string bp_type;
    while (ss >> bp_type) {
      constraint.add_bp_type(normalize_base_pair_type(bp_type));
    }
    return constraint;
  }

  std::istringstream ss(str);
  std::string token;
  bool saw_parenthesized_pair = false;
  while (std::getline(ss, token, '-')) {
    token.erase(std::remove_if(token.begin(), token.end(), ::isspace), token.end());
    if (token.empty()) {
      throw std::invalid_argument("Empty token in NMR stack constraint: " + str);
    }

    if (token.front() == '(' || token.back() == ')') {
      if (token.size() != 4 || token.front() != '(' || token.back() != ')') {
        throw std::invalid_argument("Invalid parenthesized base pair in NMR stack constraint: " + token);
      }
      constraint.add_bp_type(normalize_base_pair_type(token.substr(1, 2)));
      saw_parenthesized_pair = true;
      continue;
    }

    if (token.size() != 1) {
      throw std::invalid_argument("Invalid base in NMR stack constraint: " + token);
    }
    const char base = normalize_base(token[0]);
    if (base == 'G') {
      constraint.add_bp_type("GC");
    } else if (base == 'U') {
      constraint.add_bp_type("AU");
    } else {
      throw std::invalid_argument("NMR stack constraint supports bare G or U only: " + token);
    }
  }

  if (!saw_parenthesized_pair) {
    throw std::invalid_argument("NMR stack constraint requires a parenthesized base pair: " + str);
  }
  return constraint;
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
bool satisfies_flush_coaxial_constraint(const std::string& seq, const VI& bpseq,
                                         const StackConstraint& constraint) {
  if (constraint.size() != 2) return false;
  std::vector<std::pair<int,int>> selected_pairs;
  for (int i = 0; i < static_cast<int>(bpseq.size()); ++i)
    if (i < bpseq[i]) selected_pairs.emplace_back(i, bpseq[i]);

  auto encloses = [](const auto& outer, const auto& inner) {
    return outer.first < inner.first && inner.second < outer.second;
  };
  auto crosses = [](const auto& a, const auto& b) {
    return (a.first < b.first && b.first < a.second && a.second < b.second) ||
           (b.first < a.first && a.first < b.second && b.second < a.second);
  };
  auto types_match = [&](const auto& a, const auto& b) {
    const auto ta = normalize_base_pair_type(seq[a.first], seq[a.second]);
    const auto tb = normalize_base_pair_type(seq[b.first], seq[b.second]);
    return (ta == constraint.bp_types[0] && tb == constraint.bp_types[1]) ||
           (ta == constraint.bp_types[1] && tb == constraint.bp_types[0]);
  };

  for (const auto& closing : selected_pairs) {
    bool planar = true;
    for (const auto& p : selected_pairs)
      if (crosses(closing, p)) { planar = false; break; }
    if (!planar) continue;

    std::vector<std::pair<int,int>> children;
    for (const auto& p : selected_pairs) {
      if (!encloses(closing, p)) continue;
      bool direct = true;
      for (const auto& q : selected_pairs) {
        if (q != p && encloses(closing, q) && encloses(q, p)) {
          direct = false;
          break;
        }
      }
      if (direct) children.push_back(p);
    }
    std::sort(children.begin(), children.end());
    if (children.size() < 2) continue;

    if (children.front().first == closing.first + 1 &&
        types_match(closing, children.front())) return true;
    for (size_t i = 1; i < children.size(); ++i) {
      if (children[i].first == children[i-1].second + 1 &&
          types_match(children[i-1], children[i])) return true;
    }
    if (closing.second == children.back().second + 1 &&
        types_match(children.back(), closing)) return true;
  }
  return false;
}

bool satisfies_stack_constraints(const std::string& seq, const VI& bpseq,
                                  const StackConstraints& stack_constraints,
                                  bool allow_coaxial_stacking) {
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
      if (allow_coaxial_stacking &&
          satisfies_flush_coaxial_constraint(seq, bpseq,
                                              stack_constraints.constraints[constraint_id])) {
        spdlog::info("Stack constraint {} satisfied by flush coaxial stacking",
                     constraint_id + 1);
        constraint_satisfied = true;
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

// Find all instances of stack patterns in the sequence.  NMR does not reveal
// a one-nucleotide bulge reliably, so consecutive constrained base pairs may
// be directly stacked or separated by one bulged base on either strand.
enum class NMRBulgeMode {
  NONE,
  FALLBACK,
  ALL
};

static
void
find_stack_instances(const std::string& seq, StackConstraints& stack_constraints,
                     NMRBulgeMode bulge_mode)
{
  stack_constraints.clear_instances();

  for (size_t constraint_id = 0; constraint_id < stack_constraints.constraints.size(); ++constraint_id) {
    const auto& constraint = stack_constraints.constraints[constraint_id];
    const size_t L = seq.size();
    std::set<std::vector<std::pair<int, int>>> unique_instances;
    const size_t first_instance = stack_constraints.instances.size();

    auto enumerate_direction = [&](const std::vector<std::string>& bp_types) {
      std::function<void(size_t, size_t, size_t, bool,
                         std::vector<std::pair<int, int>>&)> extend;
      extend = [&](size_t type_index, size_t left, size_t right,
                   bool has_bulge,
                   std::vector<std::pair<int, int>>& pairs) {
        // IPknot does not admit sharp hairpins with fewer than three enclosed
        // nucleotides, so such textual matches can never become witnesses.
        if (left >= right || right < left + 4 ||
            normalize_base_pair_type(seq[left], seq[right]) != bp_types[type_index]) {
          return;
        }

        pairs.emplace_back(left, right);
        if (type_index + 1 == bp_types.size()) {
          if (unique_instances.insert(pairs).second) {
            StackInstance instance(constraint_id, has_bulge);
            for (const auto& [l, r] : pairs) instance.add_pair(l, r);
            stack_constraints.add_instance(instance);

            std::ostringstream pos_ss;
            for (const auto& [l, r] : pairs) pos_ss << "(" << l+1 << "," << r+1 << ") ";
            spdlog::debug("Found stack/bulge instance: {}", pos_ss.str());
          }
        } else {
          // Direct stack, one-base bulge on the left strand, or one-base
          // bulge on the right strand.
          const std::pair<size_t, size_t> steps_with_bulge[] = {
              {1, 1}, {2, 1}, {1, 2}};
          const std::pair<size_t, size_t> direct_step[] = {{1, 1}};
          const bool enumerate_bulges = bulge_mode != NMRBulgeMode::NONE;
          const auto* steps = enumerate_bulges ? steps_with_bulge : direct_step;
          const size_t step_count = enumerate_bulges ? 3 : 1;
          for (size_t step_index = 0; step_index < step_count; ++step_index) {
            const auto [left_step, right_step] = steps[step_index];
            if (left + left_step < L && right >= right_step) {
              const size_t next_left = left + left_step;
              const size_t next_right = right - right_step;
              if (next_left < next_right) {
                extend(type_index + 1, next_left, next_right,
                       has_bulge || left_step != 1 || right_step != 1, pairs);
              }
            }
          }
        }
        pairs.pop_back();
      };

      for (size_t i = 0; i < L; ++i) {
        for (size_t j = i + 1; j < L; ++j) {
          std::vector<std::pair<int, int>> pairs;
          extend(0, i, j, false, pairs);
        }
      }
    };

    enumerate_direction(constraint.bp_types);
    std::vector<std::string> reverse_types(constraint.bp_types.rbegin(), constraint.bp_types.rend());
    enumerate_direction(reverse_types);

    if (bulge_mode == NMRBulgeMode::FALLBACK) {
      const auto begin = stack_constraints.instances.begin() + first_instance;
      const bool has_direct_instance = std::any_of(
          begin, stack_constraints.instances.end(),
          [](const StackInstance& instance) { return !instance.has_bulge; });
      if (has_direct_instance) {
        const size_t before = stack_constraints.instances.size();
        stack_constraints.instances.erase(
            std::remove_if(begin, stack_constraints.instances.end(),
                           [](const StackInstance& instance) {
                             return instance.has_bulge;
                           }),
            stack_constraints.instances.end());
        const size_t direct_count = stack_constraints.instances.size() - first_instance;
        spdlog::info(
            "NMR bulge fallback: constraint {} uses {} direct instance(s); "
            "discarded {} bulged instance(s)",
            constraint_id + 1, direct_count,
            before - stack_constraints.instances.size());
      } else {
        spdlog::info(
            "NMR bulge fallback: constraint {} has no direct instance; "
            "retained {} bulged instance(s)",
            constraint_id + 1,
            stack_constraints.instances.size() - first_instance);
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
  // cxxopts uses commas as the delimiter for vector-valued options.  For a
  // stack pattern, however, commas separate base-pair types within one
  // pattern.  Normalize only stack-constraint arguments before parsing so
  // both "GC AU GU" and "GC,AU,GU" retain the same grouping.
  std::vector<std::string> normalized_args;
  normalized_args.reserve(argc);
  for (int i = 0; i < argc; ++i) {
    std::string arg = argv[i];
    if (i > 0 && std::string(argv[i - 1]) == "--stack-constraint") {
      std::replace(arg.begin(), arg.end(), ',', ' ');
    } else if (arg.rfind("--stack-constraint=", 0) == 0) {
      std::replace(arg.begin() + arg.find('=') + 1, arg.end(), ',', ' ');
    }
    normalized_args.push_back(std::move(arg));
  }
  std::vector<const char*> normalized_argv;
  normalized_argv.reserve(normalized_args.size());
  for (const auto& arg : normalized_args) normalized_argv.push_back(arg.c_str());

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
  bool require_canonical_neighbor = false;
  bool allow_coaxial_stacking = false;
  NMRBulgeMode nmr_bulge_mode = NMRBulgeMode::ALL;
  NMRConstraintOptions nmr_options;
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
    ("stack-constraint", "Specify stack constraint using base-pair types (e.g., 'GU GC AU') or compact NMR notation (e.g., '(GU)-G-U'). Can be specified multiple times.",
      cxxopts::value<std::vector<std::string>>(), "\"BP1 BP2 ...\"")
    ("V,verbose", "Verbose output")
    ("loglevel", "Set the logging level (trace, debug, info, warn, error, critical)",
      cxxopts::value<std::string>()->default_value("warn"), "LEVEL")
    ("beam-size", "Beam size for LinearPartition algorithm",
      cxxopts::value<uint>()->default_value("100"), "N")
    ("base-pairs", "Specify base pair count constraints (e.g., GC=1,AU=3,GU=1,UU=1)",
      cxxopts::value<std::string>(), "CONSTRAINTS")
    ("without-canonical-neighbor", "Add non-canonical base pairs even without a canonical neighbor above or below",
      cxxopts::value<bool>()->default_value("false"))
    ("coaxial-stacking", "Allow NMR stacking constraints to match flush coaxial stacking in multibranch loops",
      cxxopts::value<bool>()->default_value("false"))
    ("without-nmr-bulge", "Require directly adjacent base pairs for NMR stack constraints (disable one-nucleotide bulges)",
      cxxopts::value<bool>()->default_value("false"))
    ("nmr-bulge-mode", "NMR stack matching mode: none, fallback (use bulges only when no direct instance exists), or all",
      cxxopts::value<std::string>()->default_value("all"), "MODE")
    ("nmr-soft", "Treat NMR base-pair counts and stacking observations as soft constraints",
      cxxopts::value<bool>()->default_value("false"))
    ("nmr-count-penalty", "Penalty per missing or excess base pair in soft NMR mode",
      cxxopts::value<double>()->default_value("1.0"), "WEIGHT")
    ("nmr-stack-penalty", "Penalty per unsatisfied stacking observation in soft NMR mode",
      cxxopts::value<double>()->default_value("1.0"), "WEIGHT")
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

  auto res = options.parse(argc, normalized_argv.data());
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
  require_canonical_neighbor = !res["without-canonical-neighbor"].as<bool>();
  allow_coaxial_stacking = res["coaxial-stacking"].as<bool>();
  const auto nmr_bulge_mode_arg = res["nmr-bulge-mode"].as<std::string>();
  if (nmr_bulge_mode_arg == "none") {
    nmr_bulge_mode = NMRBulgeMode::NONE;
  } else if (nmr_bulge_mode_arg == "fallback") {
    nmr_bulge_mode = NMRBulgeMode::FALLBACK;
  } else if (nmr_bulge_mode_arg == "all") {
    nmr_bulge_mode = NMRBulgeMode::ALL;
  } else {
    spdlog::error("NMR bulge mode must be 'none', 'fallback', or 'all' (got '{}')",
                  nmr_bulge_mode_arg);
    return 1;
  }
  // Retain the original flag as a backward-compatible alias.
  if (res["without-nmr-bulge"].as<bool>()) {
    nmr_bulge_mode = NMRBulgeMode::NONE;
  }
  nmr_options.soft = res["nmr-soft"].as<bool>();
  nmr_options.count_penalty = res["nmr-count-penalty"].as<double>();
  nmr_options.stack_penalty = res["nmr-stack-penalty"].as<double>();
  if (!std::isfinite(nmr_options.count_penalty) ||
      nmr_options.count_penalty <= 0.0 ||
      !std::isfinite(nmr_options.stack_penalty) ||
      nmr_options.stack_penalty <= 0.0) {
    spdlog::error("NMR soft-constraint penalties must be finite and positive");
    return 1;
  }
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
  if (nmr_options.soft) {
    spdlog::info(
        "NMR soft-constraint mode: count penalty={}, stack penalty={}",
        nmr_options.count_penalty, nmr_options.stack_penalty);
  } else {
    spdlog::info("NMR hard-constraint mode");
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
        StackConstraint constraint = parse_stack_constraint(constraint_str);

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

  int exit_code = 0;
  try
  {
    IPknot ipknot(pk_level, &alpha[0], levelwise, !isolated_bp, n_th,
                  require_canonical_neighbor, allow_coaxial_stacking,
                  nmr_options);
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
          find_stack_instances(fa->seq(), stack_constraints, nmr_bulge_mode);
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
            spdlog::warn("Base pair constraints are NOT satisfied{}.",
                         nmr_options.soft ? " (allowed in soft NMR mode)" : "");
          }
        }

        // Check stack constraints
        if (stack_constraints.has_constraints() && spdlog::get_level() <= spdlog::level::info) {
          if (satisfies_stack_constraints(fa->seq(), bpseq, stack_constraints,
                                          allow_coaxial_stacking)) {
            spdlog::info("All stack constraints are satisfied.");
          } else {
            spdlog::warn("Some stack constraints are NOT satisfied{}.",
                         nmr_options.soft ? " (allowed in soft NMR mode)" : "");
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
      throw std::runtime_error(input + ": Format error");
    }
  }
  catch (const char* msg)
  {
    std::cerr << msg << std::endl;
    exit_code = 1;
  }
  catch (const std::logic_error& err)
  {
    std::cerr << err.what() << std::endl;
    exit_code = 1;
  }
  catch (const std::runtime_error& err)
  {
    std::cerr << err.what() << std::endl;
    exit_code = 1;
  }

  if (os_bpseq!=&std::cout) delete os_bpseq;
  if (os_bpp) delete os_bpp;
  if (os_mfa!=&std::cout) delete os_mfa;

  return exit_code;
}
