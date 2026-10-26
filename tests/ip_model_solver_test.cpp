#include "ip.h"
#include <algorithm>
#include <cmath>
#include <limits>

void check(bool condition) {
  if (!condition)
    throw std::runtime_error("IP model solver test failed");
}
template <class Exception, class F> void fails(F operation) {
  try {
    operation();
  } catch (const Exception &) {
    return;
  }
  throw std::runtime_error("Expected IP solver error was not raised");
}
bool feasible(const IPModel &model, const std::vector<double> &values) {
  for (const auto &row : model.rows) {
    double activity = 0;
    for (const auto &[col, a] : row.terms)
      activity += a * values[col];
    if ((row.bound == IP::LO || row.bound == IP::DB || row.bound == IP::FX) &&
        activity < row.lower - 1e-7)
      return false;
    if ((row.bound == IP::UP || row.bound == IP::DB) && activity > row.upper + 1e-7)
      return false;
    if (row.bound == IP::FX && activity > row.lower + 1e-7)
      return false;
  }
  return true;
}

int main() {
  IPModel model;
  IP record(model);
  record.make_variable(0);
  record.make_variable(0);
  record.make_variable(0, 1, 1);
  int cap = record.make_constraint(IP::UP, 0, 1);
  // Recorded rows may contain repeated terms for the same column.
  record.add_constraint(cap, 0, 0.25);
  record.add_constraint(cap, 0, 0.75);
  record.add_constraint(cap, 1, 1);
  int lower = record.make_constraint(IP::LO, 0, 0);
  record.add_constraint(lower, 0, 1);
  record.add_constraint(lower, 1, 1);
  int interval = record.make_constraint(IP::DB, 1, 2);
  for (int col = 0; col < 3; ++col)
    record.add_constraint(interval, col, 1);
  int fixed = record.make_constraint(IP::FX, 1, -100);
  record.add_constraint(fixed, 2, 1); // FX uses lower, ignoring upper.
  int free = record.make_constraint(IP::FR, 100, -100);
  record.add_constraint(free, 0, 1); // FR ignores both bounds.
  const auto original = model.rows[cap].terms;

  for (auto direction : {IP::MAX, IP::MIN}) {
    IPModelSolver solver(model, direction);
    fails<std::logic_error>([&] { solver.get_value(0); });
    for (const auto &cost :
         std::vector<std::vector<double>>{{2, 1, 5}, {1, 3, -2}, {-1, -2, 0}, {2, 1, 5}}) {
      double best = direction == IP::MAX ? -std::numeric_limits<double>::infinity()
                                         : std::numeric_limits<double>::infinity();
      for (int bits = 0; bits < 4; ++bits) {
        const std::vector<double> values{double(bits & 1), double((bits >> 1) & 1), 1};
        if (!feasible(model, values))
          continue;
        double score = 0;
        for (int col = 0; col < 3; ++col)
          score += cost[col] * values[col];
        best = direction == IP::MAX ? std::max(best, score) : std::min(best, score);
      }
      const auto result = solver.solve(cost);
      check(std::abs(result.objective - best) < 1e-7);
      check(direction == IP::MAX ? result.bound >= best - 1e-7 : result.bound <= best + 1e-7);
      std::vector<double> values;
      for (int col = 0; col < 3; ++col)
        values.push_back(solver.get_value(col));
      check(values[2] == 1 && feasible(model, values));
    }
    fails<std::invalid_argument>([&] { solver.solve({1}); });
    fails<std::invalid_argument>(
        [&] { solver.solve({0, 0, std::numeric_limits<double>::infinity()}); });
  }
  check(model.rows[cap].terms == original); // Solver loading leaves the source intact.

  IPModel continuous;
  IP lp(continuous);
  lp.make_continuous_variable(0, 0, 1);
  int limit = lp.make_constraint(IP::UP, 0, 0.25);
  lp.add_constraint(limit, 0, 1);
  IPModelSolver relaxation(continuous, IP::MAX);
  auto result = relaxation.solve({4});
  check(std::abs(result.objective - 1) < 1e-7 && std::abs(result.bound - 1) < 1e-7);
  check(std::abs(relaxation.get_value(0) - 0.25) < 1e-7);
  check(relaxation.solve({-1}).objective == 0);

  IPModel impossible;
  IP hard(impossible);
  hard.make_variable(0);
  int contradiction = hard.make_constraint(IP::LO, 2, 2);
  hard.add_constraint(contradiction, 0, 1);
  IPModelSolver infeasible(impossible, IP::MAX);
  fails<IPInfeasible>([&] { infeasible.solve({1}); });
  fails<std::logic_error>([&] { infeasible.get_value(0); });

  IPModel empty;
  IPModelSolver zero(empty, IP::MAX);
  check(zero.solve({}).objective == 0);
  empty.rows.push_back({IP::FX, 1, 1, {}});
  IPModelSolver constant(empty, IP::MAX);
  fails<IPInfeasible>([&] { constant.solve({}); });
}
