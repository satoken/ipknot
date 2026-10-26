#include "dd_bounds.h"
#include "dual_decomposition.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
bool crosses(const DDPair& a, const DDPair& b) {
  return (a.left < b.left && b.left < a.right && a.right < b.right) ||
         (b.left < a.left && a.left < b.right && b.right < a.right);
}
// Enumerate physical endpoint-disjoint matchings, independently of every DP
// and endpoint-cover recurrence used by the implementation.
double exhaustive(int length, const std::vector<DDPair>& pairs,
                  const std::vector<double>& weights,
                  const std::vector<unsigned char>& allowed,
                  int level, bool no_lonely, bool noncrossing,
                  const std::vector<DDScoredContact>& contacts = {}) {
  std::vector<std::vector<int>> by_left(length);
  for (int id = 0; id < static_cast<int>(pairs.size()); ++id)
    if (allowed[id] && (level < 0 || pairs[id].level == level)) by_left[pairs[id].left].push_back(id);
  std::vector<unsigned char> used(length), chosen(pairs.size());
  std::vector<int> selected;
  double best = 0;
  std::function<void(int)> visit = [&](int position) {
    if (position == length) {
      if (no_lonely) for (int id : selected) {
        bool stacked = false;
        for (int other : selected) if (pairs[id].level == pairs[other].level &&
            ((pairs[id].left + 1 == pairs[other].left && pairs[id].right - 1 == pairs[other].right) ||
             (pairs[id].left - 1 == pairs[other].left && pairs[id].right + 1 == pairs[other].right))) stacked = true;
        if (!stacked) return;
      }
      double score = 0;
      for (int id : selected) score += weights[id];
      for (const auto& contact : contacts)
        if (chosen[contact.upper] && chosen[contact.lower]) score += contact.score;
      best = std::max(best, score);
      return;
    }
    visit(position + 1);
    if (used[position]) return;
    for (int id : by_left[position]) if (!used[pairs[id].right]) {
      bool invalid = false;
      if (noncrossing) for (int other : selected)
        if (pairs[id].level == pairs[other].level && crosses(pairs[id], pairs[other])) invalid = true;
      if (invalid) continue;
      used[position] = used[pairs[id].right] = chosen[id] = 1;
      selected.push_back(id);
      visit(position + 1);
      selected.pop_back();
      used[position] = used[pairs[id].right] = chosen[id] = 0;
    }
  };
  visit(0);
  return best;
}
double endpoint_bound(int length, const std::vector<DDPair>& pairs,
                      const std::vector<double>& weights,
                      const std::vector<unsigned char>& allowed, int level) {
  std::vector<double> maxima(length);
  for (int id = 0; id < static_cast<int>(pairs.size()); ++id)
    if (allowed[id] && pairs[id].level == level) {
      maxima[pairs[id].left] = std::max(maxima[pairs[id].left], weights[id]);
      maxima[pairs[id].right] = std::max(maxima[pairs[id].right], weights[id]);
    }
  double result = 0;
  for (double value : maxima) result += value / 2;
  return result;
}
}

int main() {
  std::mt19937 generator(271828);
  const std::vector<int> widths{1, 2, 3, 4, 8, 16, 32, 64};
  for (int trial = 0; trial < 200; ++trial) {
    const int length = 2 + generator() % 8;
    std::vector<DDPair> pairs;
    for (int left = 0; left < length; ++left) for (int right = left + 1; right < length; ++right)
      if (generator() % 3) pairs.push_back({left, right, static_cast<int>(generator() % 2), 0});
    std::vector<double> weights(pairs.size());
    std::vector<unsigned char> allowed(pairs.size());
    for (int id = 0; id < static_cast<int>(pairs.size()); ++id) {
      weights[id] = (static_cast<int>(generator() % 41) - 20) / 7.0;
      allowed[id] = generator() % 5 != 0;
    }
    for (int level : {0, 1}) {
      const double exact_regular = exhaustive(length, pairs, weights, allowed, level, false, true);
      const double exact_no_lonely = exhaustive(length, pairs, weights, allowed, level, true, true);
      const double old_endpoint = endpoint_bound(length, pairs, weights, allowed, level);
      for (int width : widths) for (int offset : {0, width / 2}) {
        DDBlockBound bound(length, pairs, level, width, offset);
        const double actual = bound.evaluate(weights, allowed);
        require(actual >= exact_regular - 1e-12, "Block certificate below exhaustive regular oracle");
        require(actual >= exact_no_lonely - 1e-12, "Block certificate below exhaustive no-lonely oracle");
        require(actual <= old_endpoint + 1e-10, "Block certificate weaker than endpoint certificate");
        if (width >= length && offset == 0)
          require(std::abs(actual - exact_regular) < 1e-10, "Single-block regular oracle is not exact");
        DDBlockBound stacked_bound(length, pairs, level, width, offset, true);
        const double stacked_actual = stacked_bound.evaluate(weights, allowed);
        require(stacked_actual >= exact_no_lonely - 1e-12, "Boundary-aware block certificate below no-lonely oracle");
        require(stacked_actual <= actual + 1e-10, "No-lonely block relaxation weaker than ordinary relaxation");
        if (width >= length && offset == 0)
          require(std::abs(stacked_actual - exact_no_lonely) < 1e-10, "Single-block no-lonely oracle is not exact");
        DDBlockBound strict_bound(length, pairs, level, width, offset, true, true);
        const double strict_actual = strict_bound.evaluate(weights, allowed);
        require(strict_actual >= exact_no_lonely - 1e-12, "Strict boundary stack certificate below global no-lonely oracle");
        require(strict_actual <= stacked_actual + 1e-10, "Strict boundary stack certificate weakened ordinary boundary relaxation");
        const std::vector<unsigned char> empty(pairs.size());
        require(bound.evaluate(weights, empty) == 0, "Reusable block certificate leaked an active pair");
        require(std::abs(bound.evaluate(weights, allowed) - actual) < 1e-12, "Reusable block certificate changed on repetition");
      }
    }
    std::vector<DDScoredContact> contacts;
    for (int count = 0; count < 8 && pairs.size() > 1; ++count) {
      const int a = generator() % pairs.size(), b = generator() % pairs.size();
      if (a != b) contacts.push_back({a, b, (static_cast<int>(generator() % 31) - 15) / 4.0});
    }
    for (int id = 0; id < static_cast<int>(pairs.size()); ++id) pairs[id].weight = weights[id];
    const double unrestricted_optimum = exhaustive(length, pairs, weights, allowed, -1, false, false, contacts);
    const double global = dd_global_bound(length, pairs, contacts, allowed);
    require(global >= unrestricted_optimum - 1e-12, "Global certificate below independent signed-product matching oracle");
    const std::vector<unsigned char> empty(pairs.size());
    require(dd_global_bound(length, pairs, contacts, empty) == 0, "Global certificate retained masked products");
  }
  // A physical star exposes the half-max certificate's loose leaf credits.
  std::vector<DDPair> star;
  for (int right = 1; right < 18; ++right) star.push_back({0, right, 0, 1});
  std::vector<double> weights(star.size(), 1);
  std::vector<unsigned char> allowed(star.size(), 1);
  DDBlockBound star_bound(18, star, 0, 1);
  require(std::abs(star_bound.evaluate(weights, allowed) - 1) < 1e-12, "Oriented credits failed to tighten a star");
  require(std::abs(dd_global_bound(18, star, {}, allowed) - 1) < 1e-12, "Global endpoint cover failed to tighten a star");
  // A block enforces the crossing exclusion ignored by endpoint-only bounds.
  std::vector<DDPair> crossing_pairs{{0, 2, 0, 1}, {1, 3, 0, 1}};
  weights.assign(2, 1); allowed.assign(2, 1);
  DDBlockBound crossing_bound(4, crossing_pairs, 0, 4);
  require(std::abs(crossing_bound.evaluate(weights, allowed) - 1) < 1e-12, "Internal block did not enforce noncrossing");
  // A valid stack split by a boundary cannot be deleted by an interior
  // no-lonely rule. Our ordinary block relaxation remains above its score.
  std::vector<DDPair> stack{{1, 6, 0, 3}, {2, 5, 0, -1}};
  weights = {3, -1}; allowed.assign(2, 1);
  DDBlockBound stack_bound(8, stack, 0, 4, 2);
  require(stack_bound.evaluate(weights, allowed) >= 2, "Boundary stack was incorrectly disallowed");
  DDBlockBound stacked_stack_bound(8, stack, 0, 4, 2, true);
  require(stacked_stack_bound.evaluate(weights, allowed) >= 2, "Boundary-aware no-lonely stack was incorrectly disallowed");
  std::vector<DDPair> split_stack{{1, 8, 0, 2}, {2, 7, 0, 2}};
  weights = {2, 2}; allowed.assign(2, 1);
  DDBlockBound split_bound(10, split_stack, 0, 7, 2, true);
  require(split_bound.evaluate(weights, allowed) >= 4, "Internal boundary pair lost its external stacking support");
  DDBlockBound strict_split_bound(10, split_stack, 0, 7, 2, true, true);
  require(strict_split_bound.evaluate(weights, allowed) >= 4, "Strict internal boundary pair lost its actual external parent");
  const std::vector<unsigned char> masked_parent{0, 1};
  require(strict_split_bound.evaluate(weights, masked_parent) == 0, "Strict block borrowed a masked external parent");
  require(split_bound.evaluate(weights, masked_parent) >= 2, "Default geometric boundary relaxation changed");
  std::vector<DDPair> duplicate_parent{{1, 8, 0, 2}, {1, 8, 0, 2}, {2, 7, 0, 2}};
  weights.assign(3, 2); allowed = {0, 1, 1};
  DDBlockBound duplicate_strict(10, duplicate_parent, 0, 7, 2, true, true);
  require(duplicate_strict.evaluate(weights, allowed) >= 4, "Strict stack group ignored an allowed duplicate physical parent");
  std::vector<DDPair> negative_inner{{0, 5, 0, 3}, {1, 4, 0, -1}};
  weights = {3, -1}; allowed.assign(2, 1);
  DDBlockBound negative_inner_bound(6, negative_inner, 0, 6, 0, true);
  require(std::abs(negative_inner_bound.evaluate(weights, allowed) - 2) < 1e-12, "No-lonely certificate discarded negative inner support");
  std::vector<DDPair> huge_stack{{0, 5, 0, 1e308}, {1, 4, 0, 1e308}, {2, 3, 0, 1e308}};
  weights.assign(3, 1e308); allowed.assign(3, 1);
  DDBlockBound huge_bound(6, huge_stack, 0, 6, 0, true);
  require(std::isinf(huge_bound.evaluate(weights, allowed)), "Overflowed positive stack lost its upper certificate");
  // Tiny signed values stress outward rounding and masks independently.
  const double tiny = std::numeric_limits<double>::denorm_min();
  std::vector<DDPair> tiny_pairs{{0, 2, 0, tiny}, {1, 3, 0, tiny}};
  weights.assign(2, tiny); allowed.assign(2, 1);
  DDBlockBound tiny_bound(4, tiny_pairs, 0, 1);
  require(tiny_bound.evaluate(weights, allowed) >= 2 * tiny, "Subnormal credits undercovered an external matching");
  // Fixed-width long sparse evaluation exercises iterative buffers and masks
  // without introducing a quadratic chart or a sequence-wide sort.
  const int length = 100000;
  std::vector<DDPair> long_pairs;
  for (int left = 0; left + 4 < length; left += 6) long_pairs.push_back({left, left + 4, 0, 1});
  weights.assign(long_pairs.size(), 1); allowed.assign(long_pairs.size(), 1);
  const auto start = std::chrono::steady_clock::now();
  DDBlockBound long_bound(length, long_pairs, 0, 32);
  const double value = long_bound.evaluate(weights, allowed);
  require(value >= long_pairs.size() && value < long_pairs.size() + 1e-6, "Long sparse isolated-pair certificate mismatch");
  require(dd_global_bound(length, long_pairs, {}, allowed) >= long_pairs.size(), "Long global sparse certificate mismatch");
  const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  std::cout << "Independent exhaustive block/global bound tests and 100k sparse evaluation passed (" << elapsed << " seconds)\n";
}
