#include "dd_noe.h"
#include "dd_constrained.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <tuple>

namespace {
using Terms = std::vector<std::pair<int, double>>;
Terms merge(const Terms& source) {
  std::map<int, double> sum;
  for (const auto& [col, a] : source) sum[col] += a;
  Terms terms;
  for (const auto& [col, a] : sum) if (a != 0) terms.emplace_back(col, a);
  return terms;
}
}

DDNoEFactor dd_noe_factor(const IPModel& source,
    const std::vector<DDPair>& pairs, const std::vector<int>& columns) {
  DDNoEFactor factor;
  const int n = source.variables.size();
  std::vector<int> local(n, -1);
  std::map<int, std::pair<int, int>> physical;
  std::map<int, std::vector<int>> expression;
  std::map<std::pair<int, int>, std::set<int>> pair_columns;
  for (std::size_t id = 0; id < pairs.size(); ++id) {
    physical[columns.at(id)] = {pairs[id].left, pairs[id].right};
    expression[columns[id]] = {columns[id]};
    pair_columns[physical[columns[id]]].insert(columns[id]);
  }
  for (int col : source.noe_columns) {
    if (col < 0 || col >= n || local[col] >= 0 || physical.count(col))
      throw std::invalid_argument("Invalid NOE factor column");
    const auto& v = source.variables[col];
    if (!v.integer || v.lower < 0 || v.upper > 1)
      throw std::invalid_argument("NOE factor requires binary explanation variables");
    local[col] = factor.columns.size();
    factor.columns.push_back(col); factor.model.variables.push_back(v);
  }
  if (factor.columns.empty()) return factor;
  // Full diagnostic formulations use a binary proxy equal to the physical
  // pair's sum across levels. Canonicalize that proxy to the same pair key.
  for (const auto& row : source.rows) if (row.bound == IP::FX && row.lower == 0) {
    int proxy = -1; std::pair<int, int> pair{-1, -1}; bool valid = true;
    std::vector<int> expansion;
    for (const auto& [col, a] : merge(row.terms)) {
      if (local.at(col) >= 0) { valid = false; break; }
      const auto found = physical.find(col);
      if (a == 1 && found == physical.end() && proxy < 0) proxy = col;
      else if (a == -1 && found != physical.end()) {
        const auto& e = expression.at(col); expansion.insert(expansion.end(), e.begin(), e.end());
        if (pair.first < 0) pair = found->second;
        else if (pair != found->second) valid = false;
      } else valid = false;
    }
    if (valid && proxy >= 0 && pair.first >= 0) {
      std::sort(expansion.begin(), expansion.end());
      if (std::adjacent_find(expansion.begin(), expansion.end()) == expansion.end()) {
        physical[proxy] = pair; expression[proxy] = std::move(expansion);
      }
    }
  }
  // Certify physical-pair capacity from the recorded model, rather than
  // assuming it for arbitrary models. A nonnegative <=1 row covering every
  // level choice implies that any subset of those choices has sum <=1.
  std::set<std::pair<int, int>> capacity;
  for (const auto& [p, cols] : pair_columns)
    if (cols.size() == 1 && source.variables[*cols.begin()].upper <= 1) capacity.insert(p);
  for (const auto& row : source.rows) {
    if (row.bound != IP::UP && row.bound != IP::DB && row.bound != IP::FX) continue;
    const double upper = row.bound == IP::FX ? row.lower : row.upper;
    if (upper > 1) continue;
    std::map<std::pair<int, int>, std::set<int>> covered;
    bool valid = true;
    for (const auto& [col, a] : merge(row.terms)) {
      if (a < 0 || source.variables[col].lower < 0) { valid = false; break; }
      if (a >= 1 && physical.count(col) && expression[col].size() == 1)
        covered[physical[col]].insert(expression[col].front());
    }
    if (valid) for (const auto& [p, cols] : covered)
      if (cols == pair_columns[p]) capacity.insert(p);
  }
  using Key = std::tuple<int, double, double, Terms>;
  std::set<Key> seen;
  auto add = [&](IP::BoundType bound, double lower, double upper, Terms terms) {
    terms = merge(terms);
    if (!seen.emplace(int(bound), lower, upper, terms).second) return false;
    factor.model.rows.push_back({bound, lower, upper, std::move(terms)});
    return true;
  };
  using Pair = std::pair<int, int>;
  std::map<std::vector<int>, std::vector<Terms>> requirements;
  std::vector<std::pair<std::vector<int>, Terms>> blockers;
  for (int r = 0; r < static_cast<int>(source.rows.size()); ++r) {
    const auto& row = source.rows[r];
    Terms noe, other;
    for (const auto& [col, a] : merge(row.terms))
      (local.at(col) >= 0 ? noe : other).emplace_back(local[col] >= 0 ? local[col] : col, a);
    if (noe.empty()) continue;
    if (other.empty()) {
      add(row.bound, row.lower, row.upper, noe); factor.retained_rows.push_back(r);
      continue;
    }
    if (row.bound != IP::UP ||
        !std::all_of(noe.begin(), noe.end(), [](auto t) { return t.second == 1; })) continue;
    Pair pair{-1, -1}; bool valid = true;
    std::vector<int> expr;
    const double sign = other.front().second;
    if (sign != -1 && sign != 1) continue;
    for (const auto& [col, a] : other) {
      const auto found = physical.find(col);
      if (a != sign || found == physical.end()) { valid = false; break; }
      if (pair.first < 0) pair = found->second;
      else if (pair != found->second) valid = false;
      const auto& e = expression.at(col); expr.insert(expr.end(), e.begin(), e.end());
    }
    if (!valid) continue;
    std::sort(expr.begin(), expr.end());
    if (std::adjacent_find(expr.begin(), expr.end()) != expr.end()) continue;
    if (sign == -1 && row.upper == 0) {
      // sum witnesses <= physical pair <= 1, including multiple levels.
      if (capacity.count(pair) && noe.size() > 1 && add(IP::UP, 0, 1, noe)) ++factor.sharing_rows;
      requirements[expr].push_back(noe);
    } else if (sign == 1 && row.upper == 1) blockers.emplace_back(std::move(expr), noe);
  }
  for (const auto& [pair, blocked] : blockers)
    for (const auto& required : requirements[pair]) {
      auto conflict = blocked;
      conflict.insert(conflict.end(), required.begin(), required.end());
      // b + x <= 1 and r <= x imply b + r <= 1. No RNA column remains.
      if (add(IP::UP, 0, 1, std::move(conflict))) ++factor.blocker_rows;
    }
  return factor;
}

struct DDNoEOracle::Impl {
  DDNoEFactor factor;
  std::vector<unsigned char> member, retained;
  std::vector<double> costs, values;
  std::size_t calls = 0, cache_hits = 0;
  std::size_t primal_calls = 0, primal_feasible = 0;
  std::size_t primal_cache_hits = 0;
  const IPModel* primal_source = nullptr;
  std::vector<double> primal_fixed, primal_values;
  bool primal_valid = false;
  double seconds = 0, primal_seconds = 0;
  Result last{0, 0};
  std::unique_ptr<IPModelSolver> solver;
  Impl(DDNoEFactor f, const std::vector<double>& lower, const std::vector<double>& upper)
      : factor(std::move(f)), member(lower.size()), costs(factor.columns.size()), values(costs.size()) {
    if (lower.size() != upper.size()) throw std::invalid_argument("Invalid NOE oracle domains");
    for (std::size_t k = 0; k < factor.columns.size(); ++k) {
      int col = factor.columns[k]; member.at(col) = 1;
      factor.model.variables[k].lower = lower.at(col);
      factor.model.variables[k].upper = upper.at(col);
    }
    for (int r : factor.retained_rows) {
      if (r >= static_cast<int>(retained.size())) retained.resize(r + 1);
      retained[r] = 1;
    }
    if (costs.empty()) return;
    if (!IP::available()) throw std::invalid_argument("NOE ILP needs a linked ILP solver");
    solver = std::make_unique<IPModelSolver>(factor.model, IP::MAX, 1);
  }
  Result solve(const std::vector<double>& adjusted, std::vector<double>& selected) {
    const auto started = std::chrono::steady_clock::now();
    std::vector<double> next;
    for (int col : factor.columns) {
      const double c = adjusted.at(col);
      if (!std::isfinite(c)) throw std::overflow_error("NOE ILP coefficient overflow");
      next.push_back(c);
    }
    if (calls && next == costs) ++cache_hits;
    else if (!next.empty()) {
      costs = std::move(next); ++calls;
      IPModelSolver::Result result;
      try { result = solver->solve(costs); }
      catch (const IPInfeasible&) { throw DDInfeasible("Retained NOE ILP is infeasible"); }
      last = {result.objective, result.bound};
      for (std::size_t k = 0; k < costs.size(); ++k) values[k] = solver->get_value(k);

      if (values.size() != costs.size()) throw std::runtime_error("Invalid NOE ILP solution size");
      double score = 0;
      for (std::size_t k = 0; k < costs.size(); ++k) {
        const auto& v = factor.model.variables[k];
        if (!std::isfinite(values[k]) || values[k] < v.lower - 1e-6 || values[k] > v.upper + 1e-6 ||
            std::abs(values[k] - std::round(values[k])) > 1e-6)
          throw std::runtime_error("Invalid NOE ILP integer assignment");
        values[k] = std::round(values[k]); score += costs[k] * values[k];
      }
      for (const auto& row : factor.model.rows) {
        double activity = 0;
        for (const auto& [col, a] : row.terms) activity += a * values[col];
        if (((row.bound == IP::LO || row.bound == IP::DB || row.bound == IP::FX) && activity < row.lower - 1e-6) ||
            ((row.bound == IP::UP || row.bound == IP::DB) && activity > row.upper + 1e-6) ||
            (row.bound == IP::FX && activity > row.lower + 1e-6))
          throw std::runtime_error("NOE ILP violates an internal row");
      }
      if (!std::isfinite(last.value) || !std::isfinite(last.upper_bound) ||
          std::abs(last.value - score) > 1e-6 * std::max(1.0, std::abs(score)))
        throw std::runtime_error("Invalid NOE ILP objective or upper bound");
      last.value = score; last.upper_bound = std::max(score, last.upper_bound);
    }
    for (std::size_t k = 0; k < factor.columns.size(); ++k) selected.at(factor.columns[k]) = values[k];
    seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    return last;
  }
  bool recover(const IPModel& source, std::vector<double>& selected) {
    const auto started = std::chrono::steady_clock::now();
    std::vector<double> fixed = selected;
    for (int col : factor.columns) fixed[col] = 0;
    const auto finish = [&](bool feasible) {
      primal_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
      if (feasible) ++primal_feasible;
      return feasible;
    };
    if (primal_source == &source && fixed == primal_fixed) {
      ++primal_cache_hits;
      if (primal_valid) for (std::size_t k=0;k<factor.columns.size();++k)
        selected[factor.columns[k]]=primal_values[k];
      return finish(primal_valid);
    }
    primal_source = &source; primal_fixed = std::move(fixed); primal_valid = false;
    std::vector<int> local(member.size(), -1);
    for (std::size_t k=0;k<factor.columns.size();++k) local[factor.columns[k]]=k;
    for (std::size_t col=0;col<member.size();++col) if (!member[col]) {
      const auto& v=source.variables[col]; const double x=selected.at(col);
      if (!std::isfinite(x) || x<v.lower-1e-8 || x>v.upper+1e-8 ||
          (v.integer && std::abs(x-std::round(x))>1e-8)) return finish(false);
    }
    std::vector<IPModel::Row> rows;
    for (const auto& row:source.rows) {
      double fixed=0; Terms terms;
      for (const auto& [col,a]:row.terms) {
        if (member[col]) terms.emplace_back(local[col],a);
        else fixed+=a*selected[col];
      }
      if (terms.empty()) {
        if (((row.bound==IP::LO || row.bound==IP::DB || row.bound==IP::FX) && fixed<row.lower-1e-8) ||
            ((row.bound==IP::UP || row.bound==IP::DB) && fixed>row.upper+1e-8) ||
            (row.bound==IP::FX && fixed>row.lower+1e-8)) return finish(false);
      } else rows.push_back({row.bound,row.lower-fixed,row.upper-fixed,merge(terms)});
    }
    ++primal_calls;
    IPModel conditional;
    conditional.variables = factor.model.variables;
    conditional.rows = std::move(rows);
    IPModelSolver solver(conditional, IP::MAX, 1);
    std::vector<double> objective;
    for (const auto& v : conditional.variables) objective.push_back(v.coefficient);
    try { solver.solve(objective); }
    catch(const IPInfeasible&) { return finish(false); }
    primal_values.clear();
    for (std::size_t k=0;k<factor.columns.size();++k) {
      const double x=solver.get_value(k);
      if (!std::isfinite(x) || std::abs(x-std::round(x))>1e-6)
        throw std::runtime_error("Invalid conditional NOE integer assignment");
      primal_values.push_back(std::round(x));
      selected[factor.columns[k]]=primal_values.back();
    }
    primal_valid = true;
    return finish(true);
  }
};

DDNoEOracle::DDNoEOracle(DDNoEFactor f, const std::vector<double>& lo, const std::vector<double>& hi)
    : impl_(std::make_unique<Impl>(std::move(f), lo, hi)) {}
DDNoEOracle::~DDNoEOracle() = default;
bool DDNoEOracle::contains(int c) const { return impl_->member.at(c); }
bool DDNoEOracle::retains(int r) const { return r < int(impl_->retained.size()) && impl_->retained[r]; }
std::size_t DDNoEOracle::variables() const { return impl_->factor.columns.size(); }
std::size_t DDNoEOracle::rows() const { return impl_->factor.model.rows.size(); }
std::size_t DDNoEOracle::calls() const { return impl_->calls; }
std::size_t DDNoEOracle::cache_hits() const { return impl_->cache_hits; }
double DDNoEOracle::seconds() const { return impl_->seconds; }
std::size_t DDNoEOracle::primal_calls() const { return impl_->primal_calls; }
std::size_t DDNoEOracle::primal_feasible() const { return impl_->primal_feasible; }
std::size_t DDNoEOracle::primal_cache_hits() const { return impl_->primal_cache_hits; }
double DDNoEOracle::primal_seconds() const { return impl_->primal_seconds; }
const DDNoEFactor& DDNoEOracle::factor() const { return impl_->factor; }
DDNoEOracle::Result DDNoEOracle::solve(const std::vector<double>& c, std::vector<double>& x) { return impl_->solve(c, x); }
bool DDNoEOracle::recover(const IPModel& m, std::vector<double>& x) { return impl_->recover(m, x); }
