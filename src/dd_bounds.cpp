#include "dd_bounds.h"
#include "dual_decomposition.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace {
// Outward rounding avoids turning an upper certificate into a lower value
// through ordinary binary64 addition/subtraction. The compiler must preserve
// IEEE operations (in particular, do not compile this file with fast-math).
double add_up(double a, double b) {
  if (a == 0) return b;
  if (b == 0) return a;
  return std::nextafter(a + b, std::numeric_limits<double>::infinity());
}
double subtract_up(double a, double b) {
  if (b == 0) return a;
  if (a == b) return 0;
  return std::nextafter(a - b, std::numeric_limits<double>::infinity());
}
}

struct DDBlockBound::Workspace {
  int length, width;
  bool no_lonely, strict_stack;
  const std::vector<DDPair>& pairs;
  std::vector<int> selected_level, external, block_of, boundaries;
  std::vector<std::vector<int>> internal_by_right, external_incident;
  std::vector<double> credits, incident_max, left_max, right_max, chart, paired;
  std::vector<int> physical_group, outer_group;
  std::vector<unsigned char> physical_allowed;

  Workspace(int n, const std::vector<DDPair>& p, int level, int w, int offset, bool lonely, bool strict)
      : length(n), width(std::min(n, w)), no_lonely(lonely), strict_stack(strict && lonely), pairs(p), block_of(n),
        internal_by_right(n), external_incident(n), credits(n),
        incident_max(n), left_max(n), right_max(n),
        chart(static_cast<std::size_t>(width) * width), paired(chart.size()),
        physical_group(strict_stack ? p.size() : 0, -1), outer_group(physical_group.size(), -1) {
    boundaries.push_back(0);
    if (offset > 0 && offset < n) boundaries.push_back(offset);
    for (int start = boundaries.back(); start < n;) {
      const int end = start + std::min(w, n - start);
      boundaries.push_back(end);
      start = end;
    }
    for (int block = 0; block + 1 < static_cast<int>(boundaries.size()); ++block)
      for (int i = boundaries[block]; i < boundaries[block + 1]; ++i) block_of[i] = block;
    for (int id = 0; id < static_cast<int>(p.size()); ++id) if (p[id].level == level) {
      selected_level.push_back(id);
      if (block_of[p[id].left] == block_of[p[id].right])
        internal_by_right[p[id].right].push_back(id);
      else {
        external.push_back(id);
        external_incident[p[id].left].push_back(id);
        external_incident[p[id].right].push_back(id);
      }
    }
    if (strict_stack) {
      const auto token = [](int left, int right) {
        return (std::uint64_t(static_cast<unsigned>(left)) << 32) | static_cast<unsigned>(right);
      };
      std::unordered_map<std::uint64_t, int> groups;
      for (int id : selected_level) {
        const int next = static_cast<int>(groups.size());
        const auto entry = groups.emplace(token(p[id].left, p[id].right), next);
        physical_group[id] = entry.first->second;
      }
      physical_allowed.resize(groups.size());
      for (int id : selected_level) if (p[id].left > 0 && p[id].right + 1 < n) {
        const auto parent = groups.find(token(p[id].left - 1, p[id].right + 1));
        if (parent != groups.end()) outer_group[id] = parent->second;
      }
    }
  }

  double certificate(const std::vector<double>& weights,
                     const std::vector<unsigned char>& allowed) {
    double result = 0;
    for (double credit : credits) result = add_up(result, credit);
    for (int block = 0; block + 1 < static_cast<int>(boundaries.size()); ++block) {
      const int begin = boundaries[block], end = boundaries[block + 1], size = end - begin;
      std::fill(chart.begin(), chart.end(), 0);
      if (no_lonely) std::fill(paired.begin(), paired.end(), -std::numeric_limits<double>::infinity());
      // F(i,j) = max(F(i,j-1), max_(k,j) F(i,k-1)+w'(k,j)+F(k+1,j-1)).
      // Each internal sparse edge contributes to at most width states.
      for (int j = 0; j < size; ++j) {
        for (int i = 0; i <= j; ++i)
          chart[static_cast<std::size_t>(i) * width + j] =
              i < j ? chart[static_cast<std::size_t>(i) * width + j - 1] : 0;
        for (int id : internal_by_right[begin + j]) if (allowed[id]) {
          const auto& pair = pairs[id];
          const int k = pair.left - begin;
          const double adjusted = subtract_up(subtract_up(weights[id], credits[pair.left]), credits[pair.right]);
          const double inside = k + 1 < j ? chart[static_cast<std::size_t>(k + 1) * width + j - 1] : 0;
          double component = add_up(adjusted, inside);
          if (no_lonely) {
            const double child = k + 1 < j ? paired[static_cast<std::size_t>(k + 1) * width + j - 1]
                                         : -std::numeric_limits<double>::infinity();
            auto& paired_cell = paired[static_cast<std::size_t>(k) * width + j];
            paired_cell = std::max(paired_cell, add_up(adjusted, std::max(inside, child)));
            const bool outer_can_leave_block = pair.left > 0 && pair.right + 1 < length &&
                                               (k == 0 || begin + j + 1 == end) &&
                                               (!strict_stack || (outer_group[id] >= 0 && physical_allowed[outer_group[id]]));
            component = child != -std::numeric_limits<double>::infinity()
                ? add_up(adjusted, child) : -std::numeric_limits<double>::infinity();
            if (outer_can_leave_block) component = std::max(component, add_up(adjusted, inside));
          }
          if (!(component > 0)) continue;
          for (int i = 0; i <= k; ++i) {
            const double prefix = i < k ? chart[static_cast<std::size_t>(i) * width + k - 1] : 0;
            const double candidate = add_up(prefix, component);
            auto& cell = chart[static_cast<std::size_t>(i) * width + j];
            cell = std::max(cell, candidate);
          }
        }
      }
      if (size) result = add_up(result, chart[size - 1]);
    }
    return result;
  }

  void reduce_credits(const std::vector<double>& weights,
                      const std::vector<unsigned char>& allowed, bool reverse) {
    // This coordinate step only decreases credits and preserves every
    // external edge cover. It cannot increase the mathematical certificate:
    // any internal matching consumes a coordinate at most once.
    for (int position = 0; position < length; ++position) {
      const int i = reverse ? length - position - 1 : position;
      double needed = 0;
      for (int id : external_incident[i]) if (allowed[id]) {
        const auto& pair = pairs[id];
        const int other = pair.left == i ? pair.right : pair.left;
        if (weights[id] > credits[other])
          needed = std::max(needed, subtract_up(weights[id], credits[other]));
      }
      credits[i] = std::min(credits[i], needed);
    }
  }
};

DDBlockBound::DDBlockBound(int length, const std::vector<DDPair>& pairs,
                         int level, int width, int offset, bool no_lonely_pairs, bool strict_stack) {
  if (length < 0 || level < 0 || width < 1 || offset < 0 || offset >= width)
    throw std::invalid_argument("Invalid block-bound dimensions");
  for (const auto& pair : pairs)
    if (pair.left < 0 || pair.right >= length || pair.left >= pair.right || pair.level < 0 || !std::isfinite(pair.weight))
      throw std::invalid_argument("Invalid block-bound pair");
  workspace_ = std::make_unique<Workspace>(length, pairs, level, width, offset, no_lonely_pairs, strict_stack);
}
DDBlockBound::~DDBlockBound() = default;

double DDBlockBound::evaluate(const std::vector<double>& weights,
                            const std::vector<unsigned char>& allowed) {
  auto& w = *workspace_;
  if (weights.size() != w.pairs.size() || allowed.size() != weights.size())
    throw std::invalid_argument("Block-bound weight/mask size mismatch");
  for (int id : w.selected_level)
    if (allowed[id] && !std::isfinite(weights[id]))
      throw std::invalid_argument("Block bound requires finite active weights");
  if (w.strict_stack) {
    std::fill(w.physical_allowed.begin(), w.physical_allowed.end(), 0);
    for (int id : w.selected_level) if (allowed[id]) w.physical_allowed[w.physical_group[id]] = 1;
  }
  std::fill(w.incident_max.begin(), w.incident_max.end(), 0);
  std::fill(w.left_max.begin(), w.left_max.end(), 0);
  std::fill(w.right_max.begin(), w.right_max.end(), 0);
  bool has_positive_external = false;
  for (int id : w.external) if (allowed[id] && weights[id] > 0) {
    has_positive_external = true;
    const auto& pair = w.pairs[id];
    const double weight = weights[id];
    w.incident_max[pair.left] = std::max(w.incident_max[pair.left], weight);
    w.incident_max[pair.right] = std::max(w.incident_max[pair.right], weight);
    w.left_max[pair.left] = std::max(w.left_max[pair.left], weight);
    w.right_max[pair.right] = std::max(w.right_max[pair.right], weight);
  }
  if (!has_positive_external) {
    std::fill(w.credits.begin(), w.credits.end(), 0);
    return w.certificate(weights, allowed);
  }
  double best = std::numeric_limits<double>::infinity();
  // Half-max handles undirected external graphs; the two oriented starts
  // recover much tighter covers on common stars and ordered sparse supports.
  // A fixed three starts and two sweeps retain fixed-budget linear complexity.
  for (int start = 0; start < 3; ++start) {
    if (start == 0) {
      for (int i = 0; i < w.length; ++i)
        w.credits[i] = w.incident_max[i] > 0
            ? std::nextafter(w.incident_max[i] / 2, std::numeric_limits<double>::infinity()) : 0;
    } else w.credits = start == 1 ? w.left_max : w.right_max;
    w.reduce_credits(weights, allowed, false);
    w.reduce_credits(weights, allowed, true);
    best = std::min(best, w.certificate(weights, allowed));
  }
  return best;
}

double dd_global_bound(int length, const std::vector<DDPair>& pairs,
                       const std::vector<DDScoredContact>& contacts,
                       const std::vector<unsigned char>& allowed) {
  if (length < 0 || pairs.size() != allowed.size())
    throw std::invalid_argument("Invalid global-bound dimensions");
  std::vector<double> unary(pairs.size());
  std::vector<std::vector<int>> incident(length);
  for (int id = 0; id < static_cast<int>(pairs.size()); ++id) {
    const auto& pair = pairs[id];
    if (pair.left < 0 || pair.right >= length || pair.left >= pair.right || pair.level < 0 || !std::isfinite(pair.weight))
      throw std::invalid_argument("Invalid global-bound pair");
    unary[id] = pair.weight;
    if (allowed[id]) {
      incident[pair.left].push_back(id);
      incident[pair.right].push_back(id);
    }
  }
  for (const auto& contact : contacts) {
    if (contact.upper < 0 || contact.lower < 0 || contact.upper >= static_cast<int>(pairs.size()) ||
        contact.lower >= static_cast<int>(pairs.size()) || contact.upper == contact.lower || !std::isfinite(contact.score))
      throw std::invalid_argument("Invalid global-bound contact");
    if (contact.score > 0 && allowed[contact.upper] && allowed[contact.lower]) {
      const double half = std::nextafter(contact.score / 2, std::numeric_limits<double>::infinity());
      unary[contact.upper] = add_up(unary[contact.upper], half);
      unary[contact.lower] = add_up(unary[contact.lower], half);
      if (!std::isfinite(unary[contact.upper]) || !std::isfinite(unary[contact.lower]))
        return std::numeric_limits<double>::infinity();
    }
  }
  std::vector<double> half_max(length), left_max(length), right_max(length), credits(length);
  for (int id = 0; id < static_cast<int>(pairs.size()); ++id) if (allowed[id] && unary[id] > 0) {
    const auto& pair = pairs[id];
    half_max[pair.left] = std::max(half_max[pair.left], unary[id]);
    half_max[pair.right] = std::max(half_max[pair.right], unary[id]);
    left_max[pair.left] = std::max(left_max[pair.left], unary[id]);
    right_max[pair.right] = std::max(right_max[pair.right], unary[id]);
  }
  double best = std::numeric_limits<double>::infinity();
  for (int start = 0; start < 3; ++start) {
    if (start == 0) {
      for (int i = 0; i < length; ++i)
        credits[i] = half_max[i] > 0
            ? std::nextafter(half_max[i] / 2, std::numeric_limits<double>::infinity()) : 0;
    } else credits = start == 1 ? left_max : right_max;
    for (bool reverse : {false, true}) for (int position = 0; position < length; ++position) {
      const int i = reverse ? length - position - 1 : position;
      double needed = 0;
      for (int id : incident[i]) {
        const auto& pair = pairs[id];
        const int other = pair.left == i ? pair.right : pair.left;
        if (unary[id] > credits[other]) needed = std::max(needed, subtract_up(unary[id], credits[other]));
      }
      credits[i] = std::min(credits[i], needed);
    }
    double result = 0;
    for (double credit : credits) result = add_up(result, credit);
    best = std::min(best, result);
  }
  return best;
}
