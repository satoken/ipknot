#ifndef IPKNOT_DD_JOINT_BOUND_H
#define IPKNOT_DD_JOINT_BOUND_H

#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

struct DDPair;
struct DDRecoveryRow;

struct DDJointBoundResult {
  double upper_bound = std::numeric_limits<double>::infinity();
  std::size_t exact_windows = 0, fallback_windows = 0, states = 0;
};

// Original-objective certificate using coupled INTEGER physical matchings
// within disjoint windows/clusters. It can be below the original DD LP.
// Internal signed PK products and structural rules are kept in joint mode;
// matching-only mode majorizes all positive products and drops structure.
// External positive contacts are majorized, external pairs credit-covered.
// A fixed state budget switches a busy window to a certified unary-cover
// relaxation; an unfinished enumeration maximum is never used as an upper
// bound. Width is limited to 12; state_budget=0 is unlimited diagnostics.
// For fixed width/levels/budget, time and memory are O(n+M+E).
// contact_clusters groups possibly distant endpoints by signed contacts,
// then unscored support and pair endpoints. Clusters are disjoint, contain
// at most width physical bases, and retain original-coordinate crossing and
// stacking rules. Offset must be0 for this optional partition.
// Pairs must outlive this object and remain unchanged. Masks/rows are copied.
class DDJointBound {
public:
  DDJointBound(int length, const std::vector<DDPair>& pairs, int levels,
               bool no_lonely_pairs, const std::vector<unsigned char>& allowed,
               const std::vector<DDRecoveryRow>& rows, int width,
               int offset = 0, bool matching_only = false,
               std::size_t state_budget = 20000, bool contact_clusters = false);
  ~DDJointBound();
  DDJointBound(const DDJointBound&) = delete;
  DDJointBound& operator=(const DDJointBound&) = delete;
  DDJointBoundResult evaluate();
private:
  struct Workspace;
  std::unique_ptr<Workspace> workspace_;
};

#endif
