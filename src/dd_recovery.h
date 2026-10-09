#ifndef IPKNOT_DD_RECOVERY_H
#define IPKNOT_DD_RECOVERY_H

#include "dual_decomposition.h"
#include <memory>
#include <string>
#include <utility>
#include <vector>

// One required lower-level witness row. Scores are the original signed
// product coefficients, never the temporary dual or proposal coefficients.
struct DDRecoveryRow {
  int upper;
  std::vector<std::pair<int, double>> contacts;
};

struct DDRecoveryProposal {
  std::vector<unsigned char> selected;
  double objective = 0;
  std::string method = "empty";
  // Number of distinct level-zero seeds completed on this call.
  int proposals = 0;
};

// Additional feasible lower-bound proposals. With a fixed history width and
// fixed beam, storage and each call are O(n + M + E). The caller controls the
// call interval, e.g. one proposal batch per ten subgradient iterations.
// These scores never supply a dual certificate or a subgradient oracle.
// The pair vector must outlive this object; masks and rows are copied.
class DDPrimalRecovery {
public:
  // proposal_mask: original=1, mean=2, upper-aware=4, all=7.
  DDPrimalRecovery(int length, const std::vector<DDPair>& pairs, int levels,
                   int beam, bool no_lonely_pairs,
                   const std::vector<unsigned char>& allowed,
                   const std::vector<DDRecoveryRow>& rows, int history = 8,
                   unsigned proposal_mask = 7, bool cache_static = true,
                   const std::vector<DDNussinov*>& shared_decoders = {},
                   bool improved_beam = false);
  ~DDPrimalRecovery();
  DDPrimalRecovery(const DDPrimalRecovery&) = delete;
  DDPrimalRecovery& operator=(const DDPrimalRecovery&) = delete;

  // Observe the coefficients used by the dual level decoders, before primal
  // recovery overwrites them. Only level-zero coefficients are retained.
  void observe_adjusted(const std::vector<double>& weights);
  // Shared decoders must outlive this object and match the constructor graph.
  // They may be reused after dual selections/gradients have been saved.
  // Original-score seed (cached), fixed-window mean seed, and three upper-aware
  // seeds. The returned objective always evaluates the original signed model.
  DDRecoveryProposal propose(const std::vector<unsigned char>& dual_selected);

private:
  struct Workspace;
  std::unique_ptr<Workspace> workspace_;
};

#endif
