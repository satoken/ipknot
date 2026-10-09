#ifndef IPKNOT_DD_EXCHANGE_H
#define IPKNOT_DD_EXCHANGE_H

#include "dd_recovery.h"
#include <cstddef>
#include <memory>
#include <vector>

struct DDExchangeResult {
  std::vector<unsigned char> selected;
  double objective = 0;
  std::size_t windows = 0, improved_windows = 0;
  std::size_t states_visited = 0, budget_windows = 0;
  // True only when the whole sequence was searched without a state limit
  // interruption. Completing a proper local window is not a global proof.
  bool globally_exact = false;
};

// Exact, bounded-width integer reoptimization of a feasible incumbent. All
// selected pairs not fully inside a window remain frozen, including crossing
// witnesses and stacking neighbors outside the window. Original signed pair
// and product coefficients are used throughout; these are primal proposals.
// The pair vector must outlive this object; masks and witness rows are copied.
class DDLocalExchange {
public:
  DDLocalExchange(int length, const std::vector<DDPair>& pairs, int levels,
                  bool no_lonely_pairs,
                  const std::vector<unsigned char>& allowed,
                  const std::vector<DDRecoveryRow>& rows);
  ~DDLocalExchange();
  DDLocalExchange(const DDLocalExchange&) = delete;
  DDLocalExchange& operator=(const DDLocalExchange&) = delete;

  // Width 1..12, one to four fixed passes, alternating tiled and half-shifted
  // windows. A positive budget caps DFS nodes per window; zero is unlimited
  // for small diagnostics. A truncated search still keeps its best feasible
  // structure, including the original incumbent, and supplies no certificate.
  DDExchangeResult improve(const std::vector<unsigned char>& incumbent,
                           int width = 8, int passes = 1,
                           std::size_t state_budget = 20000);

private:
  struct Workspace;
  std::unique_ptr<Workspace> workspace_;
};

#endif
