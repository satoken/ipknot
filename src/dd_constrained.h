#ifndef IPKNOT_DD_CONSTRAINED_H
#define IPKNOT_DD_CONSTRAINED_H

#include "dual_decomposition.h"
#include "ip.h"
#include <stdexcept>

class DDInfeasible : public std::runtime_error {
public:
  explicit DDInfeasible(const std::string& reason) : std::runtime_error(reason) {}
};

struct DDConstrainedResult {
  double objective = 0, upper_bound = 0;
  int iterations = 0;
  std::size_t repair_states = 0, nonzeros = 0, propagation_work = 0, structural_work = 0;
  std::size_t repair_calls = 0, periodic_repair_calls = 0, dp_pruned_states = 0;
  std::size_t noe_ilp_variables = 0, noe_ilp_rows = 0, noe_ilp_calls = 0, noe_ilp_cache_hits = 0;
  double noe_ilp_seconds = 0;
  std::size_t noe_primal_calls = 0, noe_primal_feasible = 0, noe_primal_cache_hits = 0;
  double noe_primal_seconds = 0;
  bool repair_budget_exhausted = false;
  std::string stop_reason;
};

// Each pair ID maps to a column in the recorded shared formulation.
// Mixed links are relaxed; optional NOE ILP retains internal NOE relations.
// Per-level RNA matchings remain in the DP.
// Primal repair checks every recorded row before exposing a solution.
DDConstrainedResult solve_constrained_dd(int length,
    const std::vector<DDPair>& pairs, const std::vector<int>& columns,
    int levels, IPModel& model, const DDOptions& options);

#endif
