#include "dd_constrained.h"
#include "dd_noe.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <deque>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <unordered_map>

namespace {
constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double tolerance = 1e-8;
struct Range { double lower, upper; };
struct Row {
  double lower, upper;
  std::vector<std::pair<int, double>> terms;
};
struct DualRow {
  double bound, multiplier = 0, gradient = 0;
  bool equality;
  std::vector<std::pair<int, double>> terms;
};

// Interval propagation preserves all feasible assignments. It is also used
// on auxiliary counts, witness links, topology blockers and PK score columns.
class Repair {
public:
  const IPModel& model;
  std::vector<Row> rows;
  std::vector<std::vector<int>> incident;
  std::vector<Range> root;
  std::vector<double> best;
  double objective = -infinity;
  std::size_t states = 0, calls = 0;
  bool exhausted = false;
  struct Decision { int parent, column; Range range; };
  struct State { std::vector<Range> domain; int changed; int path = -1; };
  std::vector<Decision> decisions;
  std::vector<State> pending;
  bool search_started = false;
  bool search_complete() const { return search_started && pending.empty(); }
  void restart_search() {
    // A second branching order shares the total node budget and incumbent.
    // The completed-search case is excluded by the caller.
    pending.clear(); decisions.clear(); search_started = false;
  }
  const std::vector<DDPair>& pairs;
  const std::vector<int>& columns;
  int length, levels;
  const DDOptions& options;
  DDNoEOracle* noe = nullptr;
  bool noe_leaf = true;
  mutable std::size_t propagation_work = 0, structural_work = 0;
  std::size_t nonzeros = 0;

  Repair(const IPModel& m, int n, const std::vector<DDPair>& p,
         const std::vector<int>& c, int k, const DDOptions& o)
      : model(m), incident(m.variables.size()), pairs(p), columns(c), length(n), levels(k), options(o) {
    for (const auto& v : m.variables) {
      if (!std::isfinite(v.coefficient) || !std::isfinite(v.lower) ||
          !std::isfinite(v.upper) || v.lower > v.upper ||
          (v.integer && (v.lower != std::floor(v.lower) || v.upper != std::floor(v.upper))))
        throw std::invalid_argument("Invalid constrained DD variable");
      root.push_back({v.lower, v.upper});
    }
    for (const auto& source : m.rows) {
      Row row{-infinity, infinity, source.terms};
      if (source.bound == IP::LO || source.bound == IP::DB || source.bound == IP::FX)
        row.lower = source.lower;
      if (source.bound == IP::UP || source.bound == IP::DB) row.upper = source.upper;
      if (source.bound == IP::FX) row.upper = source.lower;
      if (std::isnan(row.lower) || std::isnan(row.upper) || row.lower > row.upper)
        throw std::invalid_argument("Invalid constrained DD row bounds");
      // Preserve first-incidence order; hashing avoids sorting long count,
      // endpoint and NMR blocker rows in the linear path.
      std::unordered_map<int, std::size_t> slot;
      std::vector<std::pair<int, double>> merged;
      for (const auto& term : row.terms) {
        if (term.first < 0 || term.first >= static_cast<int>(m.variables.size()) ||
            !std::isfinite(term.second)) throw std::invalid_argument("Invalid constrained DD row term");
        auto found = slot.emplace(term.first, merged.size());
        if (found.second) merged.push_back(term);
        else merged[found.first->second].second += term.second;
      }
      merged.erase(std::remove_if(merged.begin(), merged.end(),
          [](const auto& term) { return term.second == 0; }), merged.end());
      row.terms = std::move(merged);
      nonzeros += row.terms.size();
      int continuous = 0;
      for (const auto& [col, coefficient] : row.terms) {
        incident[col].push_back(rows.size());
        continuous += !m.variables[col].integer;
      }
      // The crossing/blocks PK formulation has at most one continuous score
      // per row. Once integers are fixed its exact feasible interval is known.
      if (continuous > 1)
        throw std::invalid_argument("Constrained DD requires independent continuous score columns");
      rows.push_back(std::move(row));
    }
    if (!propagate(root, -1)) throw DDInfeasible(o.linear_constraints ? "Retained constrained DD graph is infeasible (constraint propagation); increase candidate/witness beams or use --dd-constraints full" : "Constrained DD model is infeasible (constraint propagation)");
  }

  // A candidate crosses a forced planar matching exactly when its free
  // endpoints have different innermost forced parents. Build these labels by
  // a linear sweep, avoiding a matrix of pair-pair conflicts.
  bool planar(std::vector<Range>& domain, bool prune, bool* tightened = nullptr) const {
    std::vector<int> starts(length * levels, -1), ends(length * levels, -1), parents(length * levels, -1);
    structural_work += 3 * static_cast<std::size_t>(length) * levels + pairs.size();
    for (std::size_t id = 0; id < pairs.size(); ++id) if (domain[columns[id]].lower > .5) {
      const auto& p = pairs[id]; const int left = p.level * length + p.left, right = p.level * length + p.right;
      if (starts[left] >= 0 || ends[left] >= 0 || starts[right] >= 0 || ends[right] >= 0) return false;
      starts[left] = id; ends[right] = id;
    }
    for (int level = 0; level < levels; ++level) {
      std::vector<int> stack;
      for (int position = 0; position < length; ++position) {
        const int pos = level * length + position;
        if (ends[pos] >= 0) {
          if (stack.empty() || stack.back() != ends[pos]) return false;
          stack.pop_back();
        }
        parents[pos] = stack.empty() ? -1 : stack.back();
        if (starts[pos] >= 0) stack.push_back(starts[pos]);
      }
    }
    if (prune) for (std::size_t id = 0; id < pairs.size(); ++id) {
      auto& d = domain[columns[id]]; const auto& p = pairs[id];
      if (d.lower < .5 && d.upper > .5 && parents[p.level * length + p.left] != parents[p.level * length + p.right]) {
        d.upper = 0; if (tightened) *tightened = true;
      }
    }
    return true;
  }

  bool propagate_linear(std::vector<Range>& domain) const {
    for (int pass = 0; pass < options.constraint_passes; ++pass) {
      bool tightened = false;
      for (const auto& row : rows) {
        double minimum = 0, maximum = 0;
        propagation_work += 2 * row.terms.size();
        for (const auto& [col, a] : row.terms) {
          minimum += a * (a > 0 ? domain[col].lower : domain[col].upper);
          maximum += a * (a > 0 ? domain[col].upper : domain[col].lower);
        }
        if (minimum > row.upper + tolerance || maximum < row.lower - tolerance) return false;
        // Sums are a snapshot. Terms are unique, so each column still has its
        // snapshot bounds when visited; previous tightenings only make the
        // other-column interval smaller. All updates remain conservative.
        for (const auto& [col, a] : row.terms) {
          auto& d = domain[col];
          const double other_min = minimum - a * (a > 0 ? d.lower : d.upper);
          const double other_max = maximum - a * (a > 0 ? d.upper : d.lower);
          double lower = a > 0 ? (row.lower - other_max) / a : (row.upper - other_min) / a;
          double upper = a > 0 ? (row.upper - other_min) / a : (row.lower - other_max) / a;
          if (model.variables[col].integer) { lower = std::ceil(lower - tolerance); upper = std::floor(upper + tolerance); }
          lower = std::max(lower, d.lower); upper = std::min(upper, d.upper);
          if (lower > upper + tolerance) return false;
          if (lower > upper) lower = upper;
          if (lower > d.lower + 1e-10 || upper < d.upper - 1e-10) { d = {lower, upper}; tightened = true; }
        }
      }
      if (!planar(domain, true, &tightened)) return false;
      if (!tightened) break;
    }
    // Stopping before a fixed point cannot remove feasible assignments.
    // Every exposed incumbent still undergoes complete row/structure checks.
    return true;
  }

  bool propagate(std::vector<Range>& domain, int changed) const {
    if (options.linear_constraints) return propagate_linear(domain);
    std::deque<int> queue;
    std::vector<unsigned char> queued(rows.size());
    const auto enqueue = [&](int r) {
      if (!queued[r]) { queued[r] = 1; queue.push_back(r); }
    };
    if (changed < 0) for (int r = 0; r < static_cast<int>(rows.size()); ++r) enqueue(r);
    else for (int r : incident[changed]) enqueue(r);
    while (!queue.empty()) {
      const int id = queue.front(); queue.pop_front(); queued[id] = 0;
      const auto& row = rows[id];
      double minimum = 0, maximum = 0;
      for (const auto& [col, a] : row.terms) {
        minimum += a * (a > 0 ? domain[col].lower : domain[col].upper);
        maximum += a * (a > 0 ? domain[col].upper : domain[col].lower);
      }
      if (minimum > row.upper + tolerance || maximum < row.lower - tolerance) return false;
      for (const auto& [col, a] : row.terms) {
        auto& d = domain[col];
        const double other_min = minimum - a * (a > 0 ? d.lower : d.upper);
        const double other_max = maximum - a * (a > 0 ? d.upper : d.lower);
        double lower = a > 0 ? (row.lower - other_max) / a : (row.upper - other_min) / a;
        double upper = a > 0 ? (row.upper - other_min) / a : (row.lower - other_max) / a;
        if (model.variables[col].integer) {
          lower = std::ceil(lower - tolerance); upper = std::floor(upper + tolerance);
        }
        lower = std::max(d.lower, lower); upper = std::min(d.upper, upper);
        if (lower > upper + tolerance) return false;
        if (lower > upper) lower = upper; // continuous roundoff only
        if (lower > d.lower + 1e-10 || upper < d.upper - 1e-10) {
          d = {lower, upper};
          // Recompute the row sums before using the tightened domain again.
          for (int r : incident[col]) enqueue(r);
          break;
        }
      }
    }
    return true;
  }

  bool feasible(const std::vector<double>& value) const {
    if (value.size() != root.size()) return false;
    for (std::size_t col = 0; col < value.size(); ++col)
      if (!std::isfinite(value[col]) || value[col] < root[col].lower - tolerance ||
          value[col] > root[col].upper + tolerance ||
          (model.variables[col].integer && std::abs(value[col] - std::round(value[col])) > tolerance)) return false;
    std::vector<Range> fixed; fixed.reserve(value.size());
    for (double v : value) fixed.push_back({v,v});
    if (!planar(fixed, false)) return false;
    for (const auto& row : rows) {
      double sum = 0;
      for (const auto& [col, a] : row.terms) sum += a * value[col];
      if (sum < row.lower - tolerance || sum > row.upper + tolerance) return false;
    }
    return true;
  }

  double consider(const std::vector<double>& value) {
    if (!feasible(value)) return -infinity;
    double score = 0;
    for (std::size_t col = 0; col < value.size(); ++col) score += model.variables[col].coefficient * value[col];
    if (score > objective) { objective = score; best = value; }
    return score;
  }

  // Bounded primal repair; an optional NOE ILP handles explanation choices.
  // Other integer bounds branch by bisection;
  // continuous score columns are optimized only after integers are fixed.
  // This search does not supply the dual certificate or claim optimality
  // when its budget is exhausted.
  void search(const std::vector<double>& preference, std::size_t budget, bool first_feasible = false) {
    exhausted = false;
    // Linear mode stores only branch decisions. Reconstructing a domain for
    // each state keeps memory O(model size + state budget), rather than
    // retaining a complete model-sized snapshot at every DFS depth.
    // Keep the frontier and branch decisions between calls. A later call
    // may use new DD preferences for fresh branches. A scheduled second
    // branching order can restart the frontier, but never the incumbent.
    if (!search_started) { pending.push_back({root, -1}); search_started = true; }
    if (pending.empty()) return;
    ++calls;
    std::size_t visited = 0;
    while (!pending.empty()) {
      if (budget && visited >= budget) { exhausted = true; break; }
      State state = std::move(pending.back()); pending.pop_back(); ++visited; ++states;
      if (options.linear_constraints && state.path >= 0) {
        state.domain = root;
        for (int path = state.path; path >= 0; path = decisions[path].parent) {
          const auto& cut = decisions[path]; auto& d = state.domain[cut.column];
          d.lower = std::max(d.lower, cut.range.lower); d.upper = std::min(d.upper, cut.range.upper);
        }
      }
      if (!propagate(state.domain, state.changed)) continue;
      // A partially forced assignment can already have a feasible all-lower
      // completion. Checking it avoids branching on thousands of irrelevant
      // binary zeros just to materialize that assignment.
      if (options.linear_constraints) {
        std::vector<double> completion; completion.reserve(root.size());
        for (const auto& d : state.domain) completion.push_back(d.lower);
        consider(completion);
        if (noe && noe_leaf) {
          double candidate_bound = 0;
          for (std::size_t col = 0; col < completion.size(); ++col) {
            const auto c = model.variables[col].coefficient;
            candidate_bound += c * (noe->contains(col)
                ? (c > 0 ? state.domain[col].upper : state.domain[col].lower) : completion[col]);
          }
          if (candidate_bound > objective + 1e-10 && noe->recover(model, completion)) consider(completion);
        }
      }
      double upper = 0;
      int branch = -1;
      double priority = -1;
      for (int col = 0; col < static_cast<int>(root.size()); ++col) {
        const auto& v = model.variables[col]; const auto& d = state.domain[col];
        upper += v.coefficient * (v.coefficient > 0 ? d.upper : d.lower);
        if (v.integer && d.lower < d.upper && !(noe && noe_leaf && noe->contains(col))) {
          const double score = (1 + incident[col].size()) / (d.upper - d.lower);
          if (score > priority) { priority = score; branch = col; }
        }
      }
      if (upper <= objective + 1e-10) {
        if (first_feasible && std::isfinite(objective)) return;
        continue;
      }
      if (branch < 0) {
        std::vector<double> value(root.size());
        for (std::size_t col = 0; col < value.size(); ++col)
          value[col] = model.variables[col].coefficient > 0 ? state.domain[col].upper : state.domain[col].lower;
        if (!(noe && noe_leaf) || noe->recover(model, value)) consider(value);
        if (first_feasible && std::isfinite(objective)) return;
        continue;
      }
      const auto d = state.domain[branch];
      const double middle = std::floor(d.lower + (d.upper - d.lower) / 2);
      auto other = state.domain;
      const bool high_first = preference[branch] > 0;
      if (high_first) { state.domain[branch].lower = middle + 1; other[branch].upper = middle; }
      else { state.domain[branch].upper = middle; other[branch].lower = middle + 1; }
      if (options.linear_constraints) {
        const int low = decisions.size(); decisions.push_back({state.path, branch, other[branch]});
        const int high = decisions.size(); decisions.push_back({state.path, branch, state.domain[branch]});
        pending.push_back({{}, branch, low}); pending.push_back({{}, branch, high});
      } else {
        pending.push_back({std::move(other), branch});
        pending.push_back({std::move(state.domain), branch});
      }
      // Finish this node before pausing, so even a feasible partial
      // completion does not lose the current node's remaining subtree.
      if (first_feasible && std::isfinite(objective)) return;
    }
  }
};
}

DDConstrainedResult solve_constrained_dd(int length,
    const std::vector<DDPair>& pairs, const std::vector<int>& columns,
    int levels, IPModel& model, const DDOptions& options) {
  options.validate();
  if (options.linear_constraints && !options.constraint_states)
    throw std::invalid_argument("Linear constrained DD needs a positive repair budget; use --dd-constraints full for unlimited diagnostics");
  if (options.bound_block || options.global_bound || options.joint_bound_width ||
      options.exchange_width || options.recovery_every)
    throw std::invalid_argument("Constrained DD uses --dd-constraint-states for repair; extra unconstrained bound/recovery options are unsupported");
  if (length < 0 || levels < 1 || pairs.size() != columns.size())
    throw std::invalid_argument("Invalid constrained DD dimensions");
  const int count = model.variables.size();
  std::vector<int> pair_column(count, -1);
  for (int id = 0; id < static_cast<int>(pairs.size()); ++id) {
    const int col = columns[id];
    if (col < 0 || col >= count || pair_column[col] >= 0 || pairs[id].left < 0 || pairs[id].right >= length ||
        pairs[id].left >= pairs[id].right || pairs[id].level < 0 || pairs[id].level >= levels ||
        model.variables[col].lower != 0 || model.variables[col].upper != 1)
      throw std::invalid_argument("Invalid constrained DD pair column");
    pair_column[col] = id;
  }
  Repair repair(model, length, pairs, columns, levels, options);
  std::unique_ptr<DDNoEOracle> noe;
  if (options.noe_ilp && !model.noe_columns.empty()) {
    std::vector<double> lower, upper;
    for (const auto& domain : repair.root) {
      lower.push_back(domain.lower); upper.push_back(domain.upper);
    }
    noe = std::make_unique<DDNoEOracle>(dd_noe_factor(model, pairs, columns), lower, upper);
    repair.noe = noe.get();
    // Real NOE links contain only integer RNA/NOE columns. If a generic
    // recorded model mixes NOE and continuous scores, preserve exhaustive
    // branching so the continuous interval is optimized at fixed integers.
    for (const auto& row : repair.rows) {
      bool has_noe = false, continuous = false;
      for (const auto& [col, a] : row.terms) {
        has_noe = has_noe || noe->contains(col);
        continuous = continuous || !model.variables[col].integer;
      }
      if (has_noe && continuous) repair.noe_leaf = false;
    }
  }
  std::vector<DualRow> dual;
  for (int r = 0; r < static_cast<int>(repair.rows.size()); ++r) {
    if (noe && noe->retains(r)) continue;
    const auto& row = repair.rows[r];
    if (row.lower == row.upper) dual.push_back({row.lower, 0, 0, true, row.terms});
    else {
      if (std::isfinite(row.upper)) dual.push_back({row.upper, 0, 0, false, row.terms});
      if (std::isfinite(row.lower)) {
        DualRow lower{-row.lower, 0, 0, false, row.terms};
        for (auto& term : lower.terms) term.second = -term.second;
        dual.push_back(std::move(lower));
      }
    }
  }
  std::vector<std::unique_ptr<DDNussinov>> decoders;
  for (int level = 0; level < levels; ++level)
    decoders.push_back(std::make_unique<DDNussinov>(length, pairs, level, options.dp_beam(), false, options.improved_beam));
  std::vector<double> weights(pairs.size()), adjusted(count), selected(count), preference(count);
  std::vector<unsigned char> allowed(pairs.size());
  for (std::size_t id = 0; id < pairs.size(); ++id) allowed[id] = repair.root[columns[id]].upper > 0;
  for (int col = 0; col < count; ++col) preference[col] = model.variables[col].coefficient;
  const std::size_t first_budget = options.constraint_states ? std::max<std::size_t>(1, options.constraint_states / 2) : 0;
  repair.consider(model.solution); // optional, fully validated initial incumbent
  if (!std::isfinite(repair.objective)) repair.search(preference, first_budget, true);
  if (!std::isfinite(repair.objective) && !repair.exhausted)
    throw DDInfeasible(options.linear_constraints ? "Retained constrained DD graph is infeasible (exhaustive primal repair); increase candidate/witness beams or use --dd-constraints full" : "Constrained DD model is infeasible (exhaustive primal repair)");
  std::ofstream trace;
  static std::atomic<unsigned long long> next_id{0};
  const auto solve_id = ++next_id;
  if (!options.trace_file.empty()) {
    trace.open(options.trace_file, std::ios::app);
    if (!trace) throw std::runtime_error("Cannot open DD trace: " + options.trace_file);
    trace.exceptions(std::ios::badbit | std::ios::failbit);
    trace << std::setprecision(17);
    trace << "{\"event\":\"problem\",\"solve_id\":" << solve_id
          << ",\"constrained\":true,\"dp\":\"" << (options.nussinov_dp ? "nussinov" : options.improved_beam ? "improved-beam" : "beam")
          << "\",\"beam\":" << options.dp_beam()
          << ",\"repair_state_budget\":" << options.constraint_states
          << ",\"constraint_recovery_every\":" << options.constraint_recovery_every
          << ",\"recovery_target\":\"" << (options.recovery_target_best ? "best" : "baseline") << "\""
          << ",\"variables\":" << count << ",\"row_count\":" << dual.size()
          << ",\"noe_solver\":\"" << (options.noe_ilp ? "ilp" : "relaxed") << "\""
          << ",\"noe_ilp_variables\":" << (noe ? noe->variables() : 0)
          << ",\"noe_ilp_rows\":" << (noe ? noe->rows() : 0);
    if (options.trace_state) {
      trace << ",\"variable_domains\":[";
      for (int col = 0; col < count; ++col) {
        if (col) trace << ',';
        const auto& v = model.variables[col];
        trace << '[' << v.coefficient << ',' << repair.root[col].lower << ','
              << repair.root[col].upper << ',' << int(v.integer) << ']';
      }
      trace << "],\"pairs\":[";
      for (std::size_t id = 0; id < pairs.size(); ++id) {
        if (id) trace << ',';
        const auto& p = pairs[id];
        trace << '[' << p.left << ',' << p.right << ',' << p.level << ',' << columns[id] << ']';
      }
      trace << "],\"rows\":[";
      for (std::size_t id = 0; id < dual.size(); ++id) {
        if (id) trace << ',';
        const auto& row = dual[id];
        trace << '[' << row.bound << ',' << int(row.equality) << ",[";
        for (std::size_t k = 0; k < row.terms.size(); ++k) {
          if (k) trace << ',';
          trace << '[' << row.terms[k].first << ',' << row.terms[k].second << ']';
        }
        trace << "]]";
      }
      trace << ']';
      if (noe) {
        trace << ",\"noe_columns\":[";
        const auto& f = noe->factor();
        for (std::size_t k = 0; k < f.columns.size(); ++k) { if (k) trace << ','; trace << f.columns[k]; }
        trace << "],\"noe_rows\":[";
        for (std::size_t r = 0; r < f.model.rows.size(); ++r) {
          if (r) trace << ',';
          const auto& row = f.model.rows[r];
          trace << '[' << int(row.bound) << ',' << row.lower << ',' << row.upper << ",[";
          for (std::size_t k = 0; k < row.terms.size(); ++k) {
            if (k) trace << ',';
            trace << '[' << f.columns[row.terms[k].first] << ',' << row.terms[k].second << ']';
          }
          trace << "]]";
        }
        trace << ']';
      }
    }
    trace << "}\n";
  }
  DDConstrainedResult result;
  // Match the existing recovery-target control: preserve the original DD
  // step trajectory by default, while still retaining every recovered
  // incumbent for certificates, stopping and final output. "best" opts in
  // to using checkpoint improvements in the Polyak target as well.
  double baseline_objective = repair.objective;
  // Four early resumable slices use about one fifth of the finite
  // budget. The fifth checkpoint spends the remainder with a fresh branch
  // order based on later DD preferences and the same incumbent. The schedule
  // is independent of max_iterations: at the default interval, 50/100/500
  // iterations share the recovery prefix through 50. Shorter solves spend
  // the remainder at their final recovery. Zero is unlimited full-mode work.
  const auto available = options.constraint_states - std::min(options.constraint_states, repair.states);
  const auto checkpoint_budget = available ? std::max<std::size_t>(1, available / 20) : 0;
  const auto recover = [&](int iteration, bool periodic) {
    const auto remaining = options.constraint_states - std::min(options.constraint_states, repair.states);
    if (repair.search_complete() || (options.constraint_states && !remaining)) return;
    const bool finish_budget = !periodic || iteration / options.constraint_recovery_every >= 5;
    const bool restart = finish_budget && options.constraint_states && result.periodic_repair_calls > 0 && remaining > checkpoint_budget;
    if (restart) repair.restart_search();
    const auto before = repair.states;
    repair.search(preference, !finish_budget && options.constraint_states ? std::min(checkpoint_budget, remaining) : remaining);
    if (periodic) ++result.periodic_repair_calls;
    if (trace.is_open()) {
      trace << "{\"event\":\"recovery\",\"solve_id\":" << solve_id
            << ",\"constrained\":true,\"iteration\":" << iteration
            << ",\"phase\":\"" << (periodic ? "periodic" : "final")
            << "\",\"states_visited\":" << repair.states - before
            << ",\"repair_states\":" << repair.states
            << ",\"restart\":" << (restart ? "true" : "false") << ",\"lower_bound\":";
      if (std::isfinite(repair.objective)) trace << repair.objective; else trace << "null";
      trace << ",\"search_complete\":" << (repair.search_complete() ? "true" : "false") << "}\n";
    }
  };
  double best_bound = infinity, scale = 1e-3;
  for (const auto& v : model.variables) scale = std::max(scale, std::abs(v.coefficient));
  int stagnant = 0;
  result.stop_reason = "max_iterations";
  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    for (int col = 0; col < count; ++col) adjusted[col] = model.variables[col].coefficient;
    double constant = 0;
    for (const auto& row : dual) {
      constant += row.multiplier * row.bound;
      for (const auto& [col, a] : row.terms) adjusted[col] -= row.multiplier * a;
    }
    std::fill(selected.begin(), selected.end(), 0);
    double value = constant, bound = constant;
    for (int col = 0; col < count; ++col) if (pair_column[col] < 0 && !(noe && noe->contains(col))) {
      selected[col] = adjusted[col] > 0 ? repair.root[col].upper : repair.root[col].lower;
      value += adjusted[col] * selected[col]; bound += adjusted[col] * selected[col];
    }
    if (noe) {
      const auto oracle = noe->solve(adjusted, selected);
      value += oracle.value; bound += oracle.upper_bound;
    }
    for (std::size_t id = 0; id < pairs.size(); ++id) weights[id] = adjusted[columns[id]];
    for (int level = 0; level < levels; ++level) {
      std::vector<int> decoded;
      const double score = decoders[level]->decode(weights, allowed, decoded);
      value += score;
      result.dp_pruned_states += decoders[level]->pruned_states();
      if (!options.dp_beam() || !decoders[level]->pruned_states()) bound += score;
      else {
        std::vector<double> maxima(length);
        for (std::size_t id = 0; id < pairs.size(); ++id) if (allowed[id] && pairs[id].level == level) {
          maxima[pairs[id].left] = std::max(maxima[pairs[id].left], weights[id]);
          maxima[pairs[id].right] = std::max(maxima[pairs[id].right], weights[id]);
        }
        bound += std::accumulate(maxima.begin(), maxima.end(), 0.0) / 2;
      }
      for (int id : decoded) selected[columns[id]] = 1;
    }
    if (!std::isfinite(value) || !std::isfinite(bound))
      throw std::overflow_error("Constrained DD objective overflow");
    best_bound = std::min(best_bound, bound);
    for (int col = 0; col < count; ++col)
      preference[col] += noe && noe->contains(col) ? scale * (2 * selected[col] - 1) : adjusted[col];
    const double previous = repair.objective;
    const double raw_objective = repair.consider(selected);
    baseline_objective = std::max(baseline_objective, raw_objective);
    if (noe && !std::isfinite(raw_objective)) {
      auto proposal = selected;
      if (noe->recover(model, proposal))
        baseline_objective = std::max(baseline_objective, repair.consider(proposal));
    }
    if (options.constraint_recovery_every &&
        (iteration + 1) % options.constraint_recovery_every == 0 &&
        (!std::isfinite(repair.objective) ||
         best_bound - repair.objective > tolerance * std::max(1.0, std::abs(best_bound))))
      recover(iteration + 1, true);
    stagnant = repair.objective > previous ? 0 : stagnant + 1;
    double norm = 0;
    for (auto& row : dual) {
      row.gradient = row.bound;
      for (const auto& [col, a] : row.terms) row.gradient -= a * selected[col];
      if (!options.projected_norm || row.equality || row.multiplier > 0 || row.gradient < 0)
        norm += row.gradient * row.gradient;
    }
    result.iterations = iteration + 1;
    const double gap_tolerance = tolerance * std::max(1.0, std::abs(best_bound));
    std::string stop;
    if (std::isfinite(repair.objective) && best_bound - repair.objective <= gap_tolerance) stop = "bound_gap";
    else if (norm == 0) stop = options.dp_beam() ? "beam_stationary" : "exact_stationary";
    else if (options.patience && stagnant >= options.patience && std::isfinite(repair.objective)) stop = "patience";
    else if (iteration + 1 == options.max_iterations) stop = "max_iterations";
    const double target_objective = options.recovery_target_best ? repair.objective : baseline_objective;
    const double target = std::isfinite(target_objective) ? target_objective : value - scale;
    const double eta = !stop.empty() ? 0 : !options.diminishing_step && value > target + tolerance
        ? options.step * (value - target) / norm
        : options.step * scale / (std::pow(iteration + 1.0, .75) * std::sqrt(norm));
    if (trace.is_open()) {
      trace << "{\"event\":\"iteration\",\"solve_id\":" << solve_id << ",\"constrained\":true,\"iteration\":"
            << iteration + 1 << ",\"oracle_value\":" << value << ",\"upper_bound\":" << best_bound
            << ",\"lower_bound\":";
      if (std::isfinite(repair.objective)) trace << repair.objective; else trace << "null";
      trace << ",\"baseline_lower_bound\":";
      if (std::isfinite(baseline_objective)) trace << baseline_objective; else trace << "null";
      trace << ",\"eta\":" << eta << ",\"stop\":\"" << stop << "\"";
      if (options.trace_state) {
        trace << ",\"selected\":[";
        for (int col = 0; col < count; ++col) { if (col) trace << ','; trace << selected[col]; }
        trace << "],\"recovered\":[";
        for (std::size_t col = 0; col < repair.best.size(); ++col) { if (col) trace << ','; trace << repair.best[col]; }
        trace << ']';
      }
      trace << "}\n";
    }
    if (!stop.empty()) { result.stop_reason = stop; break; }
    for (auto& row : dual) {
      row.multiplier -= eta * row.gradient;
      if (!row.equality) row.multiplier = std::max(0.0, row.multiplier);
    }
  }
  if (!std::isfinite(repair.objective) || best_bound - repair.objective > tolerance * std::max(1.0, std::abs(best_bound))) {
    recover(result.iterations, false);
  }
  if (!std::isfinite(repair.objective)) {
    if (repair.exhausted)
      throw std::runtime_error("Constrained DD found no feasible structure within the repair budget; increase --dd-constraint-states, or use --dd-constraints full --dd-constraint-states 0 for unlimited diagnostics");
    throw DDInfeasible(options.linear_constraints ? "Retained constrained DD graph is infeasible (exhaustive primal repair); increase candidate/witness beams or use --dd-constraints full" : "Constrained DD model is infeasible (exhaustive primal repair)");
  }
  if (!repair.feasible(repair.best)) throw std::logic_error("Constrained DD primal validation failed");
  model.solution = repair.best;
  result.objective = repair.objective; result.upper_bound = best_bound;
  result.repair_states = repair.states;
  result.nonzeros = repair.nonzeros; result.propagation_work = repair.propagation_work;
  result.structural_work = repair.structural_work;
  result.repair_calls = repair.calls;
  if (noe) {
    result.noe_ilp_variables = noe->variables(); result.noe_ilp_rows = noe->rows();
    result.noe_ilp_calls = noe->calls(); result.noe_ilp_cache_hits = noe->cache_hits();
    result.noe_ilp_seconds = noe->seconds();
    result.noe_primal_calls = noe->primal_calls(); result.noe_primal_feasible = noe->primal_feasible();
    result.noe_primal_seconds = noe->primal_seconds();
    result.noe_primal_cache_hits = noe->primal_cache_hits();
  }
  result.repair_budget_exhausted = options.constraint_states && repair.states >= options.constraint_states && !repair.search_complete();
  if (best_bound - repair.objective <= tolerance * std::max(1.0, std::abs(best_bound))) result.stop_reason = "bound_gap";
  if (trace.is_open() && options.trace_state) {
    trace << "{\"event\":\"repair\",\"solve_id\":" << solve_id << ",\"constrained\":true,\"objective\":"
          << result.objective << ",\"selected\":[";
    for (int col = 0; col < count; ++col) { if (col) trace << ','; trace << model.solution[col]; }
    trace << "]}\n";
  }
  if (trace.is_open())
    trace << "{\"event\":\"summary\",\"solve_id\":" << solve_id << ",\"constrained\":true,\"iterations\":"
          << result.iterations << ",\"lower_bound\":" << result.objective << ",\"upper_bound\":" << best_bound
          << ",\"repair_states\":" << result.repair_states << ",\"repair_calls\":" << result.repair_calls
          << ",\"periodic_repair_calls\":" << result.periodic_repair_calls
          << ",\"noe_ilp_variables\":" << result.noe_ilp_variables
          << ",\"noe_ilp_rows\":" << result.noe_ilp_rows
          << ",\"noe_ilp_calls\":" << result.noe_ilp_calls
          << ",\"noe_ilp_cache_hits\":" << result.noe_ilp_cache_hits
          << ",\"noe_ilp_seconds\":" << result.noe_ilp_seconds
          << ",\"noe_primal_calls\":" << result.noe_primal_calls
          << ",\"noe_primal_feasible\":" << result.noe_primal_feasible
          << ",\"noe_primal_seconds\":" << result.noe_primal_seconds
          << ",\"noe_primal_cache_hits\":" << result.noe_primal_cache_hits
          << ",\"repair_budget_exhausted\":" << (result.repair_budget_exhausted ? "true" : "false")
          << ",\"dp_pruned_states\":" << result.dp_pruned_states
          << ",\"exact_dp\":" << (!result.dp_pruned_states ? "true" : "false")
          << ",\"nonzeros\":" << result.nonzeros
          << ",\"propagation_work\":" << result.propagation_work << ",\"structural_work\":" << result.structural_work << ",\"stop\":\"" << result.stop_reason << "\"}\n";
  return result;
}
