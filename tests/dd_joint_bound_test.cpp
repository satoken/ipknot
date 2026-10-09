#include "dd_joint_bound.h"
#include "dd_recovery.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>

namespace {
void check(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
bool crosses(const DDPair& a, const DDPair& b) {
  return (a.left < b.left && b.left < a.right && a.right < b.right) ||
         (b.left < a.left && a.left < b.right && b.right < a.right);
}
// Enumerate physical matchings independently, then test the original global
// constraints directly; do not reuse the local certificate's projected rows.
double integer_oracle(int n, const std::vector<DDPair>& pairs,
                      const std::vector<unsigned char>& allowed,
                      const std::vector<DDRecoveryRow>& rows, bool no_lonely) {
  std::vector<std::vector<int>> by_left(n);
  for (int id = 0; id < static_cast<int>(pairs.size()); ++id)
    if (allowed[id]) by_left[pairs[id].left].push_back(id);
  std::vector<unsigned char> selected(pairs.size()), used(n);
  std::vector<int> ids;
  double best = 0;
  std::function<void(int)> visit = [&](int position) {
    if (position == n) {
      for (int a : ids) {
        for (int b : ids) if (pairs[a].level == pairs[b].level && crosses(pairs[a], pairs[b])) return;
        if (no_lonely) {
          bool supported = false;
          for (int b : ids) if (pairs[a].level == pairs[b].level &&
              ((pairs[a].left + 1 == pairs[b].left && pairs[a].right - 1 == pairs[b].right) ||
               (pairs[a].left - 1 == pairs[b].left && pairs[a].right + 1 == pairs[b].right))) supported = true;
          if (!supported) return;
        }
      }
      double value = 0;
      for (int id : ids) value += pairs[id].weight;
      for (const auto& row : rows) if (selected[row.upper]) {
        bool supported = false;
        for (const auto& [lower, score] : row.contacts) if (selected[lower]) {
          supported = true;
          value += score;
        }
        if (!supported) return;
      }
      best = std::max(best, value);
      return;
    }
    visit(position + 1);
    if (used[position]) return;
    for (int id : by_left[position]) if (!used[pairs[id].right]) {
      used[position] = used[pairs[id].right] = selected[id] = 1;
      ids.push_back(id);
      visit(position + 1);
      ids.pop_back();
      used[position] = used[pairs[id].right] = selected[id] = 0;
    }
  };
  visit(0);
  return best;
}
std::vector<DDRecoveryRow> rows_for(const std::vector<DDPair>& pairs, std::mt19937& random) {
  std::vector<DDRecoveryRow> rows;
  for (int upper = 0; upper < static_cast<int>(pairs.size()); ++upper)
    for (int lower_level = 0; lower_level < pairs[upper].level; ++lower_level) {
      DDRecoveryRow row{upper, {}};
      for (int lower = 0; lower < static_cast<int>(pairs.size()); ++lower)
        if (pairs[lower].level == lower_level && crosses(pairs[upper], pairs[lower]))
          row.contacts.push_back({lower, (static_cast<int>(random() % 17) - 8) / 5.0});
      rows.push_back(std::move(row));
    }
  return rows;
}
}

int main() {
  std::mt19937 random(1618033);
  for (int trial = 0; trial < 100; ++trial) {
    const int n = 2 + random() % 8, levels = 1 + random() % 3;
    std::vector<DDPair> pairs;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j)
      for (int level = 0; level < levels; ++level) if (random() % 5 == 0)
        pairs.push_back({i, j, level, (static_cast<int>(random() % 19) - 7) / 7.0});
    std::vector<unsigned char> allowed(pairs.size());
    for (auto& value : allowed) value = random() % 5 != 0;
    const auto rows = rows_for(pairs, random);
    for (bool no_lonely : {false, true}) {
      const double optimum = integer_oracle(n, pairs, allowed, rows, no_lonely);
      for (int width : {1, 2, 4, 6, 8, 12}) for (int offset : {0, width / 2}) {
        DDJointBound exact(n, pairs, levels, no_lonely, allowed, rows, width, offset, false, 0);
        const auto joint = exact.evaluate();
        check(joint.upper_bound >= optimum - 1e-11, "Joint window certificate below original integer optimum");
        check(joint.fallback_windows == 0, "Unlimited joint enumeration unexpectedly fell back");
        if (width >= n && offset == 0)
          check(std::abs(joint.upper_bound - optimum) < 1e-10, "Whole-sequence joint window differs from integer optimum");
        const auto again = exact.evaluate();
        check(again.upper_bound == joint.upper_bound && again.states == joint.states, "Cached static certificate changed");
        DDJointBound matching(n, pairs, levels, no_lonely, allowed, rows, width, offset, true, 0);
        const auto match = matching.evaluate();
        check(match.upper_bound >= joint.upper_bound - 1e-10, "Matching relaxation below joint certificate");
        DDJointBound budget(n, pairs, levels, no_lonely, allowed, rows, width, offset, false, 1);
        const auto fallback = budget.evaluate();
        check(fallback.upper_bound >= optimum - 1e-11, "Budget fallback below original integer optimum");
        check(fallback.upper_bound >= joint.upper_bound - 1e-10, "Budget fallback below complete enumeration certificate");
        check(fallback.states <= fallback.exact_windows + fallback.fallback_windows, "Per-window state cap was exceeded");
      }
      for (int width : {1, 2, 4, 6, 8, 12}) {
        DDJointBound clustered(n, pairs, levels, no_lonely, allowed, rows, width, 0, false, 0, true);
        const auto joint = clustered.evaluate();
        check(joint.upper_bound >= optimum - 1e-11, "Endpoint-cluster certificate below original integer optimum");
        check(joint.fallback_windows == 0, "Unlimited endpoint-cluster enumeration fell back");
        DDJointBound matching(n, pairs, levels, no_lonely, allowed, rows, width, 0, true, 0, true);
        check(matching.evaluate().upper_bound >= joint.upper_bound - 1e-10,
              "Endpoint-cluster matching relaxation below signed joint certificate");
        DDJointBound budget(n, pairs, levels, no_lonely, allowed, rows, width, 0, false, 1, true);
        const auto fallback = budget.evaluate();
        check(fallback.upper_bound >= joint.upper_bound - 1e-10 && fallback.upper_bound >= optimum - 1e-11,
              "Endpoint-cluster budget fallback used unfinished lower value");
        check(fallback.states <= fallback.exact_windows + fallback.fallback_windows,
              "Endpoint-cluster state cap was exceeded");
      }
    }
  }
  // Negative PK creates a strict DD integrality gap: u<=l, w_u=1,w_l=0,c=-2.
  // LP witness u=l=1/2,y=0 has score1/2. For u>1/2, y>=2u-1 makes
  // score<=2-3u<=1/2, proving D*=1/2. Integers choose l alone or empty (0).
  std::vector<DDPair> gap{{0, 2, 0, 0}, {1, 3, 1, 1}};
  std::vector<DDRecoveryRow> gap_rows{{1, {{0, -2}}}};
  std::vector<unsigned char> allowed{1, 1};
  DDJointBound coupled_gap(4, gap, 2, false, allowed, gap_rows, 4, 0, false, 0);
  check(coupled_gap.evaluate().upper_bound == 0, "Joint signed PK certificate did not close explicit DD LP gap");
  DDJointBound relaxed_gap(4, gap, 2, false, allowed, gap_rows, 4, 0, true, 0);
  check(relaxed_gap.evaluate().upper_bound >= 1, "Matching-only relaxation unexpectedly used negative PK constraint");
  // Odd-set matching constraint missing from a degree-only endpoint cover.
  std::vector<DDPair> triangle{{0, 1, 0, 1}, {1, 2, 0, 1}, {0, 2, 0, 1}};
  allowed.assign(3, 1);
  DDJointBound triangle_bound(3, triangle, 1, false, allowed, {}, 3, 0, true, 0);
  check(std::abs(triangle_bound.evaluate().upper_bound - 1) < 1e-12, "Integer matching failed to tighten triangle's degree LP value1.5");
  // External stack borrowing must use an actually allowed adjacent parent.
  std::vector<DDPair> stack{{1, 8, 0, 3}, {2, 7, 0, 2}};
  allowed = {1, 1};
  DDJointBound stack_bound(10, stack, 1, true, allowed, {}, 7, 2, false, 0);
  check(stack_bound.evaluate().upper_bound >= 5, "Internal pair lost actual external stack support");
  allowed[0] = 0;
  DDJointBound masked_stack(10, stack, 1, true, allowed, {}, 7, 2, false, 0);
  check(masked_stack.evaluate().upper_bound == 0, "Masked external stack neighbor was borrowed");
  std::vector<DDPair> witness{{1, 5, 0, 1}, {2, 7, 1, 2}};
  std::vector<DDRecoveryRow> witness_rows{{1, {{0, 0}}}};
  allowed = {1, 1};
  DDJointBound witness_bound(10, witness, 2, false, allowed, witness_rows, 7, 2, false, 0);
  check(witness_bound.evaluate().upper_bound >= 3, "Internal upper lost actual external crossing witness");
  allowed[0] = 0;
  DDJointBound masked_witness(10, witness, 2, false, allowed, witness_rows, 7, 2, false, 0);
  check(masked_witness.evaluate().upper_bound == 0, "Masked exterior crossing witness was borrowed");
  std::vector<DDPair> negative_inner{{0, 5, 0, 3}, {1, 4, 0, -1}};
  allowed = {1, 1};
  DDJointBound negative_stack(6, negative_inner, 1, true, allowed, {}, 6, 0, false, 0);
  check(std::abs(negative_stack.evaluate().upper_bound - 2) < 1e-12, "Joint matching discarded negative stacking support");
  // Extreme signed arithmetic can loosen a certificate, but never yieldNaN.
  std::vector<DDPair> huge{{0, 2, 0, 1e308}, {1, 3, 1, 1e308}};
  std::vector<DDRecoveryRow> huge_rows{{1, {{0, -1e308}}}};
  DDJointBound huge_bound(4, huge, 2, false, allowed, huge_rows, 4, 0, false, 0);
  const double huge_value = huge_bound.evaluate().upper_bound;
  check(!std::isnan(huge_value) && huge_value >= 1e308, "Extreme signed certificate underflowed or becameNaN");
  std::vector<DDPair> dense;
  for (int left = 0; left < 12; ++left) for (int right = left + 1; right < 12; ++right)
    for (int level = 0; level < 3; ++level) dense.push_back({left, right, level, 1});
  allowed.assign(dense.size(), 1);
  const auto dense_rows = rows_for(dense, random);
  DDJointBound dense_cap(12, dense, 3, false, allowed, dense_rows, 12, 0, true, 20000);
  const auto capped_dense = dense_cap.evaluate();
  check(capped_dense.fallback_windows == 1 && capped_dense.states == 20000,
        "Busy width12 window did not obey exact enumeration state budget");
  check(capped_dense.upper_bound >= 6, "Busy matching fallback below six-pair integer structure");
  // Four distant stem arms need eight endpoint slots, even though their span
  // is1,000 bases. Negative signed products have the familiar LP witness
  // selecting each complete stem with probability1/2 and product factors0.
  // Contiguous width12 cannot retain any contact; clusters close the IP gap.
  std::vector<DDPair> distant_h{{0, 500, 0, 0}, {1, 499, 0, 0},
                              {250, 999, 1, .5}, {251, 998, 1, .5}};
  std::vector<DDRecoveryRow> distant_rows{{2, {{0, -2}, {1, -2}}},
                                       {3, {{0, -2}, {1, -2}}}};
  allowed.assign(4, 1);
  DDJointBound contiguous_h(1000, distant_h, 2, true, allowed, distant_rows, 12, 0, false, 0);
  DDJointBound clustered_h(1000, distant_h, 2, true, allowed, distant_rows, 8, 0, false, 0, true);
  check(contiguous_h.evaluate().upper_bound >= 1, "Distant contiguous relaxation unexpectedly charged negative products");
  check(clustered_h.evaluate().upper_bound == 0, "Distant endpoint cluster failed signed H integer bound");
  DDJointBound clustered_h_budget(1000, distant_h, 2, true, allowed, distant_rows, 8, 0, false, 1, true);
  check(clustered_h_budget.evaluate().upper_bound >= 0 && clustered_h_budget.evaluate().fallback_windows > 0,
        "Distant endpoint-cluster fallback was not safely exercised");
  // For noncontiguous clusters either adjacent stack neighbor can be outside.
  std::vector<DDPair> split_inner{{0, 8, 0, 2}, {1, 7, 0, 3}};
  allowed.assign(2, 1);
  DDJointBound inner_borrow(9, split_inner, 1, true, allowed, {}, 2, 0, false, 0, true);
  check(inner_borrow.evaluate().upper_bound >= 5, "Endpoint cluster lost allowed external INNER stack support");
  allowed[1] = 0;
  DDJointBound masked_inner(9, split_inner, 1, true, allowed, {}, 2, 0, false, 0, true);
  check(masked_inner.evaluate().upper_bound == 0, "Endpoint cluster borrowed masked external inner stack support");
  // Original coordinates, rather than consecutive local indices, decide
  // crossing. Zero-score support joins these distant crossing endpoints.
  std::vector<DDPair> distant_witness{{0, 10, 0, 1}, {3, 20, 1, 2}};
  std::vector<DDRecoveryRow> distant_support{{1, {{0, 0}}}};
  allowed.assign(2, 1);
  DDJointBound zero_support(21, distant_witness, 2, false, allowed, distant_support, 4, 0, false, 0, true);
  check(std::abs(zero_support.evaluate().upper_bound - 3) < 1e-12,
        "Endpoint cluster misread original crossing coordinates or zero-score support");
  bool invalid_offset = false;
  try { DDJointBound invalid(21, distant_witness, 2, false, allowed, distant_support, 4, 2, false, 0, true); }
  catch (const std::invalid_argument&) { invalid_offset = true; }
  check(invalid_offset, "Endpoint clusters silently accepted an incompatible positional offset");
  // Long sparse input has a bounded number of tiny independent joint windows.
  const int n = 100000;
  std::vector<DDPair> long_pairs;
  std::vector<DDRecoveryRow> long_rows;
  for (int begin = 0; begin + 5 < n; begin += 8) {
    const int lower = long_pairs.size();
    long_pairs.push_back({begin, begin + 3, 0, 1});
    long_pairs.push_back({begin + 1, begin + 4, 1, 1});
    long_rows.push_back({lower + 1, {{lower, 2}}});
  }
  allowed.assign(long_pairs.size(), 1);
  const auto start = std::chrono::steady_clock::now();
  DDJointBound long_bound(n, long_pairs, 2, false, allowed, long_rows, 8);
  const auto long_result = long_bound.evaluate();
  check(std::abs(long_result.upper_bound - 2 * long_pairs.size()) < 1e-6, "Long independent H windows failed exact joint bound");
  check(long_result.fallback_windows == 0, "Sparse tiny windows exhausted default state budget");
  DDJointBound long_cluster(n, long_pairs, 2, false, allowed, long_rows, 8, 0, false, 20000, true);
  const auto long_cluster_result = long_cluster.evaluate();
  check(std::abs(long_cluster_result.upper_bound - long_result.upper_bound) < 1e-6 &&
        long_cluster_result.fallback_windows == 0, "Long sparse endpoint clusters differed from independent H optima");
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  std::cout << "Independent joint/matching integer bounds, DD LP gap, masks, budgets and100k sparse windows passed (" << seconds << " seconds)\n";
}
