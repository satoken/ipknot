#ifndef IPKNOT_DD_NOE_H
#define IPKNOT_DD_NOE_H

#include "dual_decomposition.h"
#include "ip.h"
#include <memory>

// RNA columns stay outside this factor. Retained rows are original NOE-only
// rows; projected rows are valid consequences of mixed links and pair capacity.
struct DDNoEFactor {
  IPModel model;
  std::vector<int> columns, retained_rows;
  std::size_t sharing_rows = 0, blocker_rows = 0;
};
DDNoEFactor dd_noe_factor(const IPModel& model,
    const std::vector<DDPair>& pairs, const std::vector<int>& columns);

class DDNoEOracle {
public:
  DDNoEOracle(DDNoEFactor factor, const std::vector<double>& lower,
              const std::vector<double>& upper);
  ~DDNoEOracle();
  bool contains(int column) const;
  bool retains(int row) const;
  std::size_t variables() const;
  std::size_t rows() const;
  std::size_t calls() const;
  std::size_t cache_hits() const;
  double seconds() const;
  std::size_t primal_calls() const;
  std::size_t primal_feasible() const;
  std::size_t primal_cache_hits() const;
  double primal_seconds() const;
  const DDNoEFactor& factor() const;
  struct Result { double value, upper_bound; };
  // Adjusted global coefficients in; only NOE columns in selected are changed.
  Result solve(const std::vector<double>& adjusted, std::vector<double>& selected);
  // Keep the DP's RNA choices fixed and optimize the original NOE penalties.
  // The recorded model remains immutable for this oracle's lifetime.
  // Returns false if the fixed non-NOE assignment is infeasible; no dual bound.
  bool recover(const IPModel& model, std::vector<double>& selected);
private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

#endif
