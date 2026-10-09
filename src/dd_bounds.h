#ifndef IPKNOT_DD_BOUNDS_H
#define IPKNOT_DD_BOUNDS_H

#include <memory>
#include <vector>

struct DDPair;

// Certified relaxation of one signed-weight, noncrossing matching oracle.
// Consecutive fixed-width blocks are solved exactly, while external arcs are
// covered by nonnegative endpoint credits. Optional no-lonely states relax
// stacking only for internal pairs whose outer neighbor could lie outside
// the block, because a valid stack can straddle a block boundary.
//
// At fixed width, construction and evaluation use O(n + M) storage and
// O(n * width + M * width) time. The pairs vector must outlive this object.
// Offset shifts the first boundary: 0 means blocks [0,width), ...; otherwise
// the first block is [0,offset). Minima of several offsets remain valid.
// strict_stack additionally requires an actually allowed same-level outer
// neighbor before a boundary pair may borrow support; false keeps the
// geometric boundary relaxation for compatibility.
class DDBlockBound {
public:
  DDBlockBound(int length, const std::vector<DDPair>& pairs, int level,
               int width, int offset = 0, bool no_lonely_pairs = false,
               bool strict_stack = false);
  ~DDBlockBound();
  DDBlockBound(const DDBlockBound&) = delete;
  DDBlockBound& operator=(const DDBlockBound&) = delete;
  double evaluate(const std::vector<double>& weights,
                  const std::vector<unsigned char>& allowed);
private:
  struct Workspace;
  std::unique_ptr<Workspace> workspace_;
};

struct DDScoredContact {
  int upper, lower;
  double score;
};

// An original-objective upper bound shared across ALL levels. Positive
// products are majorized by c*x*y <= c*(x+y)/2; negative products are dropped.
// The resulting unary matching is covered using global physical endpoints.
// It ignores level support/noncrossing/no-lonely constraints and therefore
// bounds both the integer DD graph model and its convex-hull relaxation.
// This certificate is independent of DD multipliers, so it need not bound
// the oracle value L(q) at an individual multiplier vector.
double dd_global_bound(int length, const std::vector<DDPair>& pairs,
                       const std::vector<DDScoredContact>& contacts,
                       const std::vector<unsigned char>& allowed);

#endif
