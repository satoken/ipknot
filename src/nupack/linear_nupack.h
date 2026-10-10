// Sparse, left-to-right beam approximation with NUPACK energy parameters.
// See README.md for the restricted pseudoknot grammar and complexity bounds.
#ifndef IPKNOT_LINEAR_NUPACK_H
#define IPKNOT_LINEAR_NUPACK_H

#include "nupack.h"
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

class LinearNupack
{
public:
  struct Statistics
  {
    std::size_t retained_states=0, retained_edges=0, peak_beam=0;
  };

  explicit LinearNupack(unsigned beam=100, int max_loop=30, bool pseudoknots=true);
  bool load_parameters(const char* path) { return energy_.load_parameters(path); }
  // A zero beam disables pruning, but retains the grammar restrictions.
  // Disabling pseudoknots is useful for exact short secondary-structure oracles.
  double calculate(const std::string& sequence, const std::string& constraints="");
  double log_partition_function() const { return log_z_; }
  const Statistics& statistics() const { return statistics_; }
  const std::vector<std::tuple<int,int,double>>& posterior() const { return posterior_; }

private:
  Nupack<long double> energy_;
  unsigned beam_;
  int max_loop_;
  bool pseudoknots_;
  double log_z_=0;
  Statistics statistics_;
  std::vector<std::tuple<int,int,double>> posterior_;
};

#endif
