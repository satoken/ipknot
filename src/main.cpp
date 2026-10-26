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
#include <cstdint>
#include <cmath>
#include <functional>
#include <set>
#include <limits>

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
    if (i<=0 || i>static_cast<int>(bpseq.size()))
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
  const auto previous_precision = os.precision();
  const auto previous_flags = os.flags();
  // Cached probabilities must round-trip to the same float used by the
  // decoder. Six significant digits can change threshold/objective ties.
  os.precision(std::numeric_limits<float>::max_digits10);
  os.setf(std::ios::fmtflags(0), std::ios::floatfield);
  os << "# " << desc << std::endl; 
  for (uint i=1; i!=sbp.size(); ++i)
  {
    os << i << " " << seq[i-1];
    for (const auto [j, v]: sbp[i])
      if (i<j)
        os << " " << j << ":" << v;
    os << std::endl;
  }
  os.precision(previous_precision);
  os.flags(previous_flags);
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
  uint beam_size;
  std::string input;
  bool verbose = false;

  cxxopts::Options options{progname, format("IPknot version %s", PACKAGE_VERSION)};
  options.add_options()
    ("input", "FASTA-formatted file or ALN-formatted file",
      cxxopts::value<std::string>(), "FASTA_OR_ALN")
    ("e,model", "Probabilistic model",
      cxxopts::value<std::vector<std::string>>()->default_value("LinearPartition-C"), "MODEL")
    ("decoder", "Decoder: dd (default), ilp, or auto (ILP when linked, DD otherwise)",
      cxxopts::value<std::string>()->default_value("dd"), "NAME")
    ("dd-constraints", "Constraint model: linear (bounded candidates/work), full (diagnostic)",
      cxxopts::value<std::string>()->default_value("linear"))
    ("dd-constraint-passes", "Maximum propagation sweeps per repair state in linear mode",
      cxxopts::value<int>()->default_value("8"))
    ("dd-constraint-states", "Total visited repair states per threshold/refinement solve (0: unlimited, full mode)",
      cxxopts::value<int>()->default_value("2048"))
    ("dd-constraint-recovery-every", "Resume constrained primal recovery every N iterations (0: final-only); shares total state budget",
      cxxopts::value<int>()->default_value("0"), "N")
    ("dd-max-iter", "Maximum DD iterations per threshold/refinement",
      cxxopts::value<int>()->default_value("50"), "N")
    ("dd-dp", "Per-level DP: beam, improved-beam (suffix dominance), nussinov (sparse exact), or wavefront (exact interval DP)",
      cxxopts::value<std::string>()->default_value("improved-beam"), "NAME")
    ("dd-beam", "Beam DP width (0: exact); ignored by nussinov/wavefront",
      cxxopts::value<int>()->default_value("100"), "N")
    ("dd-threads", "Threads for independent DD level oracles (0: automatic; requires OpenMP for >1)",
      cxxopts::value<int>()->default_value("1"), "N")
    ("dd-dp-threads", "Within-DP threads for beam/improved-beam/wavefront (0: automatic); overrides --dd-threads",
      cxxopts::value<int>()->default_value("1"), "N")
    ("dd-crossing-beam", "Active crossing candidates per level (0: unbounded)",
      cxxopts::value<int>()->default_value("100"), "N")
    ("dd-witnesses", "Crossing witnesses per pair/lower-level row (0: all)",
      cxxopts::value<int>()->default_value("16"), "N")
    ("dd-step", "Polyak relaxation (0 < value < 2)",
      cxxopts::value<double>()->default_value("1.5"), "VALUE")
    ("dd-patience", "Stop after this many iterations without primal improvement (0: disabled)",
      cxxopts::value<int>()->default_value("0"), "N")
    ("dd-projected-norm", "Exclude outward boundary gradients from the DD step denominator",
      cxxopts::value<bool>()->default_value("true"), "BOOL")
    ("dd-unpruned-bound", "Use an exact oracle bound when no DP state was pruned",
      cxxopts::value<bool>()->default_value("false"))
    ("dd-bound-block", "Exact internal-block certificate width (0: off; maximum 64)",
      cxxopts::value<int>()->default_value("0"), "N")
    ("dd-bound-every", "Evaluate block bounds every N iterations, also at first and last",
      cxxopts::value<int>()->default_value("10"), "N")
    ("dd-global-bound", "Add a global base-capacity certificate, computed once per solve",
      cxxopts::value<bool>()->default_value("false"))
    ("dd-bound-shift", "Take the minimum of two block partitions shifted by half a window",
      cxxopts::value<bool>()->default_value("false"))
    ("dd-bound-strict-stack", "Borrow only available external stack candidates in block bounds",
      cxxopts::value<bool>()->default_value("false"))
    ("dd-joint-bound", "Joint integer window certificate width (0: off; maximum 12)",
      cxxopts::value<int>()->default_value("0"), "N")
    ("dd-joint-clusters", "Group distant PK/support endpoints instead of consecutive sequence windows",
      cxxopts::value<bool>()->default_value("false"))
    ("dd-joint-shift", "Also evaluate joint windows shifted by half their width",
      cxxopts::value<bool>()->default_value("false"))
    ("dd-joint-matching", "Use matching-only joint certificates for comparison",
      cxxopts::value<bool>()->default_value("false"))
    ("dd-joint-states", "Search states per certificate window (0: unlimited diagnostic)",
      cxxopts::value<int>()->default_value("20000"), "N")
    ("dd-exchange", "Feasible joint local exchange width (0: off; maximum 12)",
      cxxopts::value<int>()->default_value("0"), "N")
    ("dd-exchange-passes", "Exchange passes with alternating partitions (1..4)",
      cxxopts::value<int>()->default_value("1"), "N")
    ("dd-exchange-every", "Additional exchange interval (0: first and before stopping only)",
      cxxopts::value<int>()->default_value("0"), "N")
    ("dd-exchange-states", "Search states per exchange window (0: unlimited diagnostic)",
      cxxopts::value<int>()->default_value("20000"), "N")
    ("dd-recovery-cache", "Cache the static all-upper recovery proposal",
      cxxopts::value<bool>()->default_value("true"))
    ("dd-recovery-share", "Reuse dual DP buffers for extra recovery",
      cxxopts::value<bool>()->default_value("true"))
    ("dd-recovery-every", "Extra feasible recovery every N iterations and before stopping (0: off)",
      cxxopts::value<int>()->default_value("0"), "N")
    ("dd-recovery-target", "Polyak target for extra recovery: baseline or best (experimental)",
      cxxopts::value<std::string>()->default_value("baseline"))
    ("dd-recovery-mode", "Primal proposals: original, mean, lookahead or all",
      cxxopts::value<std::string>()->default_value("all"), "MODE")
    ("dd-schedule", "DD step schedule: polyak or diminishing (exponent 0.75)",
      cxxopts::value<std::string>()->default_value("polyak"), "SCHEDULE")
    ("dd-trace", "Write per-solve DD convergence diagnostics as JSONL (overwrites FILE)",
      cxxopts::value<std::string>(), "FILE")
    ("dd-trace-state", "Include graph and oracle states for offline exact diagnostics",
      cxxopts::value<bool>()->default_value("false"))
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
    ("n,threads", "The number of threads for the available solvers",
      cxxopts::value<uint>()->default_value("1"), "N")
    ("P,param", "Read an energy parameter file for the selected model",
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
    ("loglevel", "Set the logging level (trace, debug, info, warn, error, critical)",
      cxxopts::value<std::string>()->default_value("warn"), "LEVEL")
    ("beam-size", "Beam size for LinearPartition / LinearNUPACK (0 disables pruning)",
      cxxopts::value<uint>()->default_value("100"), "N")
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
  n_th = res["threads"].as<uint>();
  output_energy = res["energy"].as<bool>();
  if (res.count("constraint")) constraint = res["constraint"].as<std::string>();
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
    DDOptions dd_options;
    const auto decoder = res["decoder"].as<std::string>();
    if (decoder != "auto" && decoder != "ilp" && decoder != "dd")
      throw std::invalid_argument("Decoder must be auto, ilp, or dd");
    dd_options.enabled = decoder == "dd";
    if (!IP::available()) {
      if (decoder == "ilp") throw std::invalid_argument("No ILP solver is linked; use --decoder dd");
      dd_options.enabled = true;
    }
    const auto constraint_model = res["dd-constraints"].as<std::string>();
    if (constraint_model != "linear" && constraint_model != "full")
      throw std::invalid_argument("DD constraints must be linear or full");
    dd_options.linear_constraints = constraint_model == "linear";
    dd_options.constraint_passes = res["dd-constraint-passes"].as<int>();
    const int constraint_states = res["dd-constraint-states"].as<int>();
    if (constraint_states < 0) throw std::invalid_argument("DD constraint state budget must be nonnegative");
    dd_options.constraint_states = constraint_states;
    dd_options.constraint_recovery_every = res["dd-constraint-recovery-every"].as<int>();
    dd_options.max_iterations = res["dd-max-iter"].as<int>();
    dd_options.beam = res["dd-beam"].as<int>();
    dd_options.threads = res["dd-threads"].as<int>();
    dd_options.dp_threads = res["dd-dp-threads"].as<int>();
    const auto dd_dp = res["dd-dp"].as<std::string>();
    if (dd_dp != "beam" && dd_dp != "improved-beam" && dd_dp != "nussinov" && dd_dp != "wavefront")
      throw std::invalid_argument("DD DP must be beam, improved-beam, nussinov or wavefront");
    dd_options.wavefront_dp = dd_dp == "wavefront";
    dd_options.nussinov_dp = dd_dp == "nussinov";
    dd_options.improved_beam = dd_dp == "improved-beam";
    dd_options.crossing_beam = res["dd-crossing-beam"].as<int>();
    dd_options.witnesses = res["dd-witnesses"].as<int>();
    dd_options.step = res["dd-step"].as<double>();
    dd_options.patience = res["dd-patience"].as<int>();
    dd_options.projected_norm = res["dd-projected-norm"].as<bool>();
    dd_options.unpruned_bound = res["dd-unpruned-bound"].as<bool>();
    dd_options.bound_block = res["dd-bound-block"].as<int>();
    dd_options.bound_every = res["dd-bound-every"].as<int>();
    dd_options.global_bound = res["dd-global-bound"].as<bool>();
    dd_options.bound_shift = res["dd-bound-shift"].as<bool>();
    dd_options.bound_strict_stack = res["dd-bound-strict-stack"].as<bool>();
    dd_options.joint_bound_width = res["dd-joint-bound"].as<int>();
    dd_options.joint_bound_clusters = res["dd-joint-clusters"].as<bool>();
    dd_options.joint_bound_shift = res["dd-joint-shift"].as<bool>();
    dd_options.joint_bound_matching = res["dd-joint-matching"].as<bool>();
    dd_options.exchange_width = res["dd-exchange"].as<int>();
    dd_options.exchange_passes = res["dd-exchange-passes"].as<int>();
    dd_options.exchange_every = res["dd-exchange-every"].as<int>();
    const int joint_states = res["dd-joint-states"].as<int>();
    const int exchange_states = res["dd-exchange-states"].as<int>();
    if (joint_states < 0 || exchange_states < 0)
      throw std::invalid_argument("DD window state budgets must be nonnegative");
    dd_options.joint_bound_states = joint_states;
    dd_options.exchange_states = exchange_states;
    dd_options.recovery_cache = res["dd-recovery-cache"].as<bool>();
    dd_options.recovery_share = res["dd-recovery-share"].as<bool>();
    dd_options.recovery_every = res["dd-recovery-every"].as<int>();
    const auto recovery_target = res["dd-recovery-target"].as<std::string>();
    if (recovery_target != "baseline" && recovery_target != "best")
      throw std::runtime_error("Invalid --dd-recovery-target: " + recovery_target);
    dd_options.recovery_target_best = recovery_target == "best";
    const auto recovery_mode = res["dd-recovery-mode"].as<std::string>();
    if (recovery_mode == "original") dd_options.recovery_mode = 1;
    else if (recovery_mode == "mean") dd_options.recovery_mode = 2;
    else if (recovery_mode == "lookahead") dd_options.recovery_mode = 4;
    else if (recovery_mode == "all") dd_options.recovery_mode = 7;
    else throw std::invalid_argument("DD recovery mode must be original, mean, lookahead or all");
    const auto dd_schedule = res["dd-schedule"].as<std::string>();
    if (dd_schedule != "polyak" && dd_schedule != "diminishing")
      throw std::invalid_argument("DD schedule must be polyak or diminishing");
    dd_options.diminishing_step = dd_schedule == "diminishing";
    if (res.count("dd-trace")) dd_options.trace_file = res["dd-trace"].as<std::string>();
    dd_options.trace_state = res["dd-trace-state"].as<bool>();
    dd_options.validate();
    if (!dd_options.trace_file.empty()) {
      if (!dd_options.enabled) throw std::invalid_argument("DD tracing requires --decoder dd");
      std::ofstream trace(dd_options.trace_file);
      if (!trace) throw std::runtime_error("Cannot create DD trace: " + dd_options.trace_file);
    }
    spdlog::info("Decoder: {}", dd_options.enabled ? "dd" : "ilp");
    if (dd_options.enabled)
      spdlog::info("DD DP: {}, beam={}, max iterations={}", res["dd-dp"].as<std::string>(),
                   dd_options.dp_beam(), dd_options.max_iterations);
    IPknot ipknot(pk_level, &alpha[0], levelwise, !isolated_bp, n_th, dd_options);
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
      if (!constraint.empty()) {
        bpseq.assign(seq.size(), BPSEQ::DOT);
        read_constraints(constraint.c_str(), bpseq);
      }
      if (max_pfval)
        std::tie(fval, fval_pk) = ipknot.solve(seq, sbp, ep, bpseq, plevel, !constraint.empty());
      else
        ipknot.solve(seq, sbp, t, bpseq, plevel, !constraint.empty());
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

        const auto input_constraints = constraint.empty() ? VI{} : bpseq;
        if (max_pfval)
          std::tie(fval, fval_pk) = ipknot.solve(fa->seq(), sbp, ep, bpseq, plevel, !constraint.empty());
        else
          ipknot.solve(fa->seq(), sbp, t, bpseq, plevel, !constraint.empty());

        for (int i=0; i!=n_refinement; ++i) // iterative refinement
        {
          en->update_bpm(pk_level, fa->seq(), bpseq, plevel, sbp);
          if (!constraint.empty()) bpseq = input_constraints;
          if (max_pfval)
            std::tie(fval, fval_pk) = ipknot.solve(fa->seq(), sbp, ep, bpseq, plevel, !constraint.empty());
          else
            ipknot.solve(fa->seq(), sbp, t, bpseq, plevel, !constraint.empty());
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
        
        const auto input_constraints = constraint.empty() ? VI{} : bpseq;
        if (max_pfval)
          std::tie(fval, fval_pk) = ipknot.solve(aln->consensus(), sbp, ep, bpseq, plevel, !constraint.empty());
        else
          ipknot.solve(aln->consensus(), sbp, t, bpseq, plevel, !constraint.empty());

        for (int i=0; i!=n_refinement; ++i)
        {
          en->update_bpm(pk_level, aln->seq(), bpseq, plevel, sbp);
          if (!constraint.empty()) bpseq = input_constraints;
          if (max_pfval)
            std::tie(fval, fval_pk) = ipknot.solve(aln->consensus(), sbp, ep, bpseq, plevel, !constraint.empty());
          else
            ipknot.solve(aln->consensus(), sbp, t, bpseq, plevel, !constraint.empty());
        }

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
