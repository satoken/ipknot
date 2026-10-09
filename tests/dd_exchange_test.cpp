#include "dd_exchange.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}
bool crossing(const DDPair& a, const DDPair& b) {
  return (a.left < b.left && b.left < a.right && a.right < b.right) ||
         (b.left < a.left && a.left < b.right && b.right < a.right);
}
bool adjacent(const DDPair& a, const DDPair& b) {
  return a.level == b.level &&
    ((a.left + 1 == b.left && a.right - 1 == b.right) ||
     (b.left + 1 == a.left && b.right - 1 == a.right));
}
struct Problem {
  int n, levels; bool no_lonely = false;
  std::vector<DDPair> pairs;
  std::vector<unsigned char> allowed;
  std::vector<DDRecoveryRow> rows;
};
bool feasible(const Problem& g, const std::vector<unsigned char>& bits) {
  if (bits.size() != g.pairs.size()) return false;
  std::vector<int> used(g.n, 0);
  for (int p = 0; p < static_cast<int>(bits.size()); ++p) if (bits[p]) {
    const auto& pair = g.pairs[p];
    if (bits[p] != 1 || !g.allowed[p] || used[pair.left]++ || used[pair.right]++) return false;
    bool stack = !g.no_lonely;
    for (int q = 0; q < static_cast<int>(bits.size()); ++q) if (bits[q] && q != p) {
      if (pair.level == g.pairs[q].level && crossing(pair, g.pairs[q])) return false;
      if (adjacent(pair, g.pairs[q])) stack = true;
    }
    if (!stack) return false;
  }
  for (const auto& row : g.rows) if (bits[row.upper]) {
    bool support = false;
    for (const auto& contact : row.contacts) support |= bits[contact.first] != 0;
    if (!support) return false;
  }
  return true;
}
double score(const Problem& g, const std::vector<unsigned char>& bits) {
  double value = 0;
  for (int p = 0; p < static_cast<int>(bits.size()); ++p)
    if (bits[p]) value += g.pairs[p].weight;
  for (const auto& row : g.rows) if (bits[row.upper])
    for (const auto& c : row.contacts) if (bits[c.first]) value += c.second;
  return value;
}
struct Optimum { double value = 0; std::vector<unsigned char> bits; };
Optimum oracle(const Problem& g) {
  std::vector<unsigned char> bits(g.pairs.size(), 0);
  std::vector<int> used(g.n, 0);
  std::vector<std::vector<int>> choices(g.n);
  for (int p = 0; p < static_cast<int>(bits.size()); ++p)
    if (g.allowed[p]) choices[g.pairs[p].left].push_back(p);
  Optimum best{0, bits};
  std::function<void(int)> visit = [&](int base) {
    while (base < g.n && used[base]) ++base;
    if (base == g.n) {
      if (feasible(g, bits)) {
        double value = score(g, bits);
        if (value > best.value) best = {value, bits};
      }
      return;
    }
    visit(base + 1);
    for (int p : choices[base]) if (!used[g.pairs[p].right]) {
      bits[p] = 1; used[base] = used[g.pairs[p].right] = 1;
      visit(base + 1);
      used[base] = used[g.pairs[p].right] = 0; bits[p] = 0;
    }
  };
  visit(0);
  return best;
}
void check(const Problem& g, const std::vector<unsigned char>& incumbent,
           const DDExchangeResult& result) {
  require(feasible(g, result.selected), "exchange returned infeasible structure");
  require(std::abs(score(g, result.selected) - result.objective) < 1e-9,
          "exchange objective differs from original signed objective");
  require(result.objective + 1e-9 >= score(g, incumbent), "exchange decreased incumbent");
}
DDExchangeResult run(const Problem& g, const std::vector<unsigned char>& incumbent,
                     int width = 12, int passes = 1, std::size_t budget = 0) {
  DDLocalExchange solver(g.n, g.pairs, g.levels, g.no_lonely, g.allowed, g.rows);
  auto result = solver.improve(incumbent, width, passes, budget);
  check(g, incumbent, result);
  return result;
}
void boundary_tests() {
  Problem support{7, 2, false, {{1,4,0,-2}, {3,6,1,10}}, {1,1}, {{1,{{0,0}}}}};
  auto r = run(support, {1,1}, 5);
  require(r.objective == 8 && !r.globally_exact, "removed frozen upper's sole witness");
  Problem stack{6, 1, true, {{0,5,0,10}, {1,4,0,-1}}, {1,1}, {}};
  r = run(stack, {1,1}, 5);
  require(r.objective == 9, "removed frozen outer's only stacking neighbor");
  Problem cross{8, 2, false,
    {{2,7,0,1}, {0,4,0,10}, {0,4,1,2}}, {1,1,1}, {{2,{{0,0}}}}};
  r = run(cross, {1,0,0}, 5);
  require(r.objective == 3 && !r.selected[1] && r.selected[2],
          "wrong frozen crossing or level support handling");
  cross.rows[0].contacts[0].second = -5;
  r = run(cross, {1,0,0}, 5);
  require(r.objective == 1, "ignored negative frozen product coefficient");
  Problem nested{12,1,false,{{0,11,0,1},{4,7,0,2}}, {1,1}, {}};
  r = run(nested, {1,0}, 4, 2);
  require(r.objective == 3, "remote enclosing frozen arc blocked nested proposal");

  // A frozen three-level upper requires its witness in each lower level.
  Problem multi{10,3,false,{{1,5,0,-2},{2,6,1,-1},{4,9,2,10}}, {1,1,1},
    {{1,{{0,0}}}, {2,{{0,1}}}, {2,{{1,-1}}}}};
  r = run(multi, {1,1,1}, 7);
  require(r.objective == 7, "failed frozen support in every lower level");
}
void coupled_tests() {
  Problem g{12,3,true,
    {{0,7,0,1},{1,6,0,-1},{2,9,1,2},{3,8,1,-1},{4,11,2,2},{5,10,2,-1}},
    {1,1,1,1,1,1}, {}};
  for (int upper = 0; upper < static_cast<int>(g.pairs.size()); ++upper)
    for (int level = 0; level < g.pairs[upper].level; ++level) {
      DDRecoveryRow row{upper,{}};
      for (int lower = 0; lower < static_cast<int>(g.pairs.size()); ++lower)
        if (g.pairs[lower].level == level && crossing(g.pairs[upper],g.pairs[lower]))
          row.contacts.push_back({lower, (upper + lower) % 2 ? -0.25 : 1.25});
      g.rows.push_back(row);
    }
  auto optimum = oracle(g);
  auto r = run(g, std::vector<unsigned char>(g.pairs.size(),0));
  require(r.globally_exact && r.objective == optimum.value && r.objective > 0,
          "coupled signed three-level exact optimum mismatch");
  auto limited = run(g, optimum.bits, 12, 1, 1);
  require(!limited.globally_exact && limited.budget_windows == 1 &&
          limited.states_visited == 1 && limited.objective == optimum.value,
          "budget must retain feasible incumbent without a proof");
}
Problem random_problem(std::mt19937& rng, int levels, bool no_lonely) {
  Problem g{8,levels,no_lonely,{},{},{}};
  for (int level = 0; level < levels; ++level)
    for (int left = 0; left < g.n; ++left)
      for (int right = left + 1; right < g.n; ++right)
        if (rng() % 100 < 38) {
          g.pairs.push_back({left,right,level,(int(rng()%21)-8)*0.25});
          g.allowed.push_back(1);
        }
  for (int upper = 0; upper < static_cast<int>(g.pairs.size()); ++upper)
    for (int level = 0; level < g.pairs[upper].level; ++level) {
      DDRecoveryRow row{upper,{}};
      for (int lower = 0; lower < static_cast<int>(g.pairs.size()); ++lower)
        if (g.pairs[lower].level == level && crossing(g.pairs[upper],g.pairs[lower]))
          row.contacts.push_back({lower,(int(rng()%17)-8)*0.25});
      if (row.contacts.empty()) g.allowed[upper] = 0;
      g.rows.push_back(row);
    }
  bool changed = true;
  while (changed) {
    changed = false;
    for (const auto& row : g.rows) if (g.allowed[row.upper]) {
      bool found = false;
      for (const auto& c : row.contacts) found |= g.allowed[c.first] != 0;
      if (!found) { g.allowed[row.upper] = 0; changed = true; }
    }
    if (no_lonely)
      for (int p = 0; p < static_cast<int>(g.pairs.size()); ++p) if (g.allowed[p]) {
        bool found = false;
        for (int q = 0; q < static_cast<int>(g.pairs.size()); ++q)
          if (g.allowed[q] && adjacent(g.pairs[p],g.pairs[q])) found = true;
        if (!found) { g.allowed[p] = 0; changed = true; }
      }
  }
  return g;
}
void randomized_tests() {
  std::mt19937 rng(718291);
  for (int sample = 0; sample < 120; ++sample) {
    Problem g = random_problem(rng, 1 + sample % 3, sample % 2);
    auto optimum = oracle(g);
    std::vector<unsigned char> empty(g.pairs.size(),0);
    auto r = run(g, empty);
    require(r.globally_exact && std::abs(r.objective - optimum.value) < 1e-9,
            "independent matching oracle disagrees at sample " + std::to_string(sample));
    auto capped = run(g, empty, 12, 1, 20000);
    require(capped.objective <= optimum.value + 1e-9, "budgeted proposal exceeds optimum");
    auto local = run(g, r.selected, 3, 2, 200);
    require(local.objective == optimum.value && !local.globally_exact,
            "partial local exchange changed exact global incumbent");
    // A feasible but intentionally incomplete incumbent exercises frozen
    // witnesses/stacking while the smaller windows jointly add new pairs.
    auto partial = run(g, empty, 3, 2, 200);
    auto enlarged = run(g, partial.selected, 5, 2, 200);
    require(enlarged.objective <= optimum.value + 1e-9, "local proposal exceeds oracle");
  }
}
void linear_size_test() {
  Problem g{20000,1,true,{},{},{}};
  for (int left = 0; left < g.n; left += 8) {
    g.pairs.push_back({left,left+5,0,2});
    g.pairs.push_back({left+1,left+4,0,-1});
    g.allowed.push_back(1); g.allowed.push_back(1);
  }
  DDLocalExchange solver(g.n,g.pairs,g.levels,true,g.allowed,g.rows);
  auto result = solver.improve(std::vector<unsigned char>(g.pairs.size(),0),8,1,20000);
  require(result.objective == 2500 && result.windows == 2500 &&
          result.improved_windows == 2500 && result.budget_windows == 0 &&
          result.states_visited < 100000 && !result.globally_exact,
          "fixed-width large sparse graph did not scale with windows");
  require(std::find(result.selected.begin(),result.selected.end(),0) == result.selected.end(),
          "large graph dropped a required negative stacking partner");
}
void invalid_tests() {
  Problem g{4,1,false,{{0,3,0,1},{1,2,0,1}}, {1,1}, {}};
  bool threw = false;
  try { run(g,{1,1},13); } catch (const std::invalid_argument&) { threw = true; }
  require(threw,"invalid exchange width accepted");
  Problem missing{6,2,false,{{0,3,0,1},{1,4,1,1}}, {1,1}, {}};
  threw = false;
  try { run(missing,{0,0}); } catch (const std::invalid_argument&) { threw = true; }
  require(threw,"missing lower-level support row accepted");
  g.no_lonely = true;
  threw = false;
  try { run(g,{1,0}); } catch (const std::invalid_argument&) { threw = true; }
  require(threw,"lonely incumbent accepted");
  Problem reused{4,1,false,{{0,2,0,1},{0,3,0,1}}, {1,1}, {}};
  threw = false;
  try { run(reused,{1,1}); } catch (const std::invalid_argument&) { threw = true; }
  require(threw,"incumbent reusing a base accepted");
}
}

int main() {
  try {
    boundary_tests(); coupled_tests(); randomized_tests(); linear_size_test(); invalid_tests();
    std::cout << "dd_exchange_test: independent integer oracle, frozen constraints, signed scores, and linear-size checks passed\n";
  } catch (const std::exception& error) {
    std::cerr << "dd_exchange_test: " << error.what() << '\n';
    return 1;
  }
}
