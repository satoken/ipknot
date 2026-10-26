#include "dd_joint_bound.h"
#include "dd_recovery.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <unordered_map>

namespace {
double up_add(double a, double b) {
  if (std::isinf(a) && a > 0) return a;
  if (std::isinf(b) && b > 0) return b;
  if (a == 0) return b;
  if (b == 0) return a;
  return std::nextafter(a + b, std::numeric_limits<double>::infinity());
}
double up_subtract(double a, double b) {
  if (b == 0) return a;
  if (a == b) return 0;
  return std::nextafter(a - b, std::numeric_limits<double>::infinity());
}
double up_half(double value) {
  return value > 0 ? std::nextafter(value / 2, std::numeric_limits<double>::infinity()) : 0;
}
std::uint64_t pair_key(int left, int right) {
  return (std::uint64_t(static_cast<unsigned>(left)) << 32) | static_cast<unsigned>(right);
}
bool crossing(const DDPair& a, const DDPair& b) {
  return (a.left < b.left && b.left < a.right && a.right < b.right) ||
         (b.left < a.left && a.left < b.right && b.right < a.right);
}
struct LocalRow {
  std::vector<int> witnesses;
  bool external = false;
};
}

struct DDJointBound::Workspace {
  int length, levels;
  bool no_lonely, matching_only, overflow = false, ready = false;
  std::size_t state_budget;
  const std::vector<DDPair>& pairs;
  std::vector<unsigned char> allowed, selected, borrow_stack;
  std::vector<int> block_of, local_position, pair_block, external, chosen;
  std::vector<std::vector<int>> positions, internal, by_left, external_incident, row_ids;
  std::vector<std::array<int, 2>> stack_neighbors;
  std::vector<std::vector<std::pair<int, double>>> internal_contacts;
  std::vector<LocalRow> local_rows;
  std::vector<double> unary, internal_bonus;
  std::array<std::vector<double>, 3> credits, adjusted;
  DDJointBoundResult cached;

  Workspace(int n, const std::vector<DDPair>& p, int l, bool lonely,
            const std::vector<unsigned char>& mask, const std::vector<DDRecoveryRow>& rows,
            int width, int offset, bool matching, std::size_t budget, bool contact_clusters)
      : length(n), levels(l), no_lonely(lonely), matching_only(matching), state_budget(budget),
        pairs(p), allowed(mask), selected(p.size()), borrow_stack(p.size()), block_of(n), local_position(n),
        pair_block(p.size(), -1), by_left(n), external_incident(n), row_ids(p.size()),
        stack_neighbors(p.size(), {-1, -1}), internal_contacts(p.size()),
        unary(p.size()), internal_bonus(p.size()) {
    std::vector<std::unordered_map<std::uint64_t, int>> lookup(levels);
    for (int id = 0; id < static_cast<int>(p.size()); ++id) {
      const auto& pair = p[id];
      if (pair.left < 0 || pair.right >= n || pair.left >= pair.right || pair.level < 0 || pair.level >= levels ||
          !std::isfinite(pair.weight)) throw std::invalid_argument("Invalid joint-bound pair");
      if (!lookup[pair.level].emplace(pair_key(pair.left, pair.right), id).second)
        throw std::invalid_argument("Duplicate physical pair within a joint-bound level");
      unary[id] = pair.weight;
    }
    if (contact_clusters) {
      // Direct assignments plus explicit small member lists avoid an
      // unbounded union-find/search cost. Every rewrite touches <=width
      // bases; contact/root deduplication touches at most four roots.
      std::vector<int> root(n);
      std::vector<std::vector<int>> members(n);
      for (int i = 0; i < n; ++i) { root[i] = i; members[i].push_back(i); }
      const auto merge = [&](const std::array<int, 4>& endpoints, int count) {
        std::array<int, 4> roots{};
        int unique = 0, size = 0, target = n;
        for (int position = 0; position < count; ++position) {
          const int r = root[endpoints[position]];
          bool seen = false;
          for (int previous = 0; previous < unique; ++previous) if (roots[previous] == r) seen = true;
          if (seen) continue;
          roots[unique++] = r;
          size += static_cast<int>(members[r].size());
          target = std::min(target, r);
        }
        if (size > width || unique < 2) return;
        for (int position = 0; position < unique; ++position) if (roots[position] != target) {
          auto& source = members[roots[position]];
          for (int base : source) { root[base] = target; members[target].push_back(base); }
          source.clear();
        }
      };
      for (bool scored : {true, false}) for (const auto& row : rows) {
        if (row.upper < 0 || row.upper >= static_cast<int>(p.size()) || p[row.upper].level < 1)
          throw std::invalid_argument("Invalid joint-bound upper support row");
        for (const auto& [lower, score] : row.contacts) {
          if (lower < 0 || lower >= static_cast<int>(p.size()) || p[lower].level >= p[row.upper].level ||
              !std::isfinite(score) || !crossing(p[row.upper], p[lower]))
            throw std::invalid_argument("Invalid joint-bound support contact");
          if (allowed[row.upper] && allowed[lower] && (score != 0) == scored)
            merge({p[row.upper].left, p[row.upper].right, p[lower].left, p[lower].right}, 4);
        }
      }
      for (int id = 0; id < static_cast<int>(p.size()); ++id) if (allowed[id])
        merge({p[id].left, p[id].right, 0, 0}, 2);
      for (int i = 0; i < n; ++i) if (!members[i].empty()) {
        // Each sort is bounded by width, never a global n/E sort.
        std::sort(members[i].begin(), members[i].end());
        positions.push_back(std::move(members[i]));
      }
    } else {
      for (int begin = 0; begin < n;) {
        const int size = begin == 0 && offset > 0 ? std::min(offset, n) : std::min(width, n - begin);
        positions.emplace_back();
        for (int i = 0; i < size; ++i) positions.back().push_back(begin + i);
        begin += size;
      }
    }
    internal.resize(positions.size());
    for (int block = 0; block < static_cast<int>(positions.size()); ++block)
      for (int position = 0; position < static_cast<int>(positions[block].size()); ++position) {
        block_of[positions[block][position]] = block;
        local_position[positions[block][position]] = position;
      }
    for (int id = 0; id < static_cast<int>(p.size()); ++id) {
      const auto& pair = p[id];
      if (!allowed[id]) continue;
      if (block_of[pair.left] == block_of[pair.right]) {
        pair_block[id] = block_of[pair.left];
        internal[pair_block[id]].push_back(id);
        by_left[pair.left].push_back(id);
      } else {
        external.push_back(id);
        external_incident[pair.left].push_back(id);
        external_incident[pair.right].push_back(id);
      }
    }
    if (no_lonely && !matching_only) for (int id = 0; id < static_cast<int>(p.size()); ++id) {
      if (pair_block[id] < 0) continue;
      const auto& pair = p[id];
      const std::array<std::pair<int, int>, 2> adjacent{{{pair.left + 1, pair.right - 1},
                                                     {pair.left - 1, pair.right + 1}}};
      for (int side = 0; side < 2; ++side) {
        const auto [left, right] = adjacent[side];
        if (left < 0 || right >= n || left >= right) continue;
        const auto found = lookup[pair.level].find(pair_key(left, right));
        if (found == lookup[pair.level].end() || !allowed[found->second]) continue;
        const int neighbor = found->second;
        if (pair_block[neighbor] == pair_block[id]) stack_neighbors[id][side] = neighbor;
        else borrow_stack[id] = 1;
      }
    }
    std::vector<std::vector<unsigned char>> lower_seen(p.size());
    std::vector<unsigned char> empty_row(p.size());
    for (int id = 0; id < static_cast<int>(p.size()); ++id) lower_seen[id].resize(p[id].level);
    for (const auto& row : rows) {
      if (row.upper < 0 || row.upper >= static_cast<int>(p.size()) || p[row.upper].level < 1)
        throw std::invalid_argument("Invalid joint-bound upper support row");
      int lower_level = -1;
      LocalRow local;
      for (const auto& [lower, score] : row.contacts) {
        if (lower < 0 || lower >= static_cast<int>(p.size()) || p[lower].level >= p[row.upper].level ||
            !std::isfinite(score) || !crossing(p[row.upper], p[lower]))
          throw std::invalid_argument("Invalid joint-bound support contact");
        if (lower_level < 0) lower_level = p[lower].level;
        else if (lower_level != p[lower].level) throw std::invalid_argument("Mixed levels in joint-bound support row");
        if (!allowed[row.upper] || !allowed[lower]) continue;
        const bool inside = pair_block[row.upper] >= 0 && pair_block[row.upper] == pair_block[lower];
        if (pair_block[row.upper] >= 0) {
          if (inside) local.witnesses.push_back(lower);
          else local.external = true;
        }
        if (inside && !matching_only) {
          internal_contacts[row.upper].push_back({lower, score});
          internal_contacts[lower].push_back({row.upper, score});
          if (score > 0) {
            internal_bonus[row.upper] = up_add(internal_bonus[row.upper], up_half(score));
            internal_bonus[lower] = up_add(internal_bonus[lower], up_half(score));
          }
        } else if (score > 0) {
          unary[row.upper] = up_add(unary[row.upper], up_half(score));
          unary[lower] = up_add(unary[lower], up_half(score));
        }
      }
      if (lower_level < 0) empty_row[row.upper] = 1;
      else if (lower_seen[row.upper][lower_level]++)
        throw std::invalid_argument("Duplicate joint-bound lower-level support row");
      if (allowed[row.upper] && pair_block[row.upper] >= 0 && !matching_only) {
        row_ids[row.upper].push_back(static_cast<int>(local_rows.size()));
        local_rows.push_back(std::move(local));
      }
    }
    if (!matching_only) for (int id = 0; id < static_cast<int>(p.size()); ++id)
      if (allowed[id] && !empty_row[id]) for (unsigned char present : lower_seen[id])
        if (!present) throw std::invalid_argument("Missing joint-bound lower-level support row");
    if (matching_only) {
      // Once every contact has been majorized, level labels impose no rule.
      // At the same physical coordinates only the largest unary can improve
      // any of the three credit scenarios (their endpoint debit is shared).
      std::unordered_map<std::uint64_t, int> physical_best;
      for (const auto& block : internal) for (int id : block) {
        const auto token = pair_key(p[id].left, p[id].right);
        const auto found = physical_best.find(token);
        if (found == physical_best.end()) physical_best.emplace(token, id);
        else if (unary[id] > unary[found->second]) found->second = id;
      }
      for (auto& row : by_left) row.clear();
      // Preserve deterministic candidate input order rather than hash order.
      for (const auto& block : internal) for (int id : block)
        if (physical_best.at(pair_key(p[id].left, p[id].right)) == id) by_left[p[id].left].push_back(id);
    }
    for (int id = 0; id < static_cast<int>(p.size()); ++id)
      if (allowed[id] && (!std::isfinite(unary[id]) || !std::isfinite(internal_bonus[id]))) overflow = true;
    for (int start = 0; start < 3; ++start) {
      credits[start].resize(n);
      adjusted[start].resize(p.size());
    }
  }

  void make_credits() {
    for (int id : external) if (unary[id] > 0) {
      const auto& pair = pairs[id];
      for (int endpoint : {pair.left, pair.right}) credits[0][endpoint] = std::max(credits[0][endpoint], unary[id]);
      credits[1][pair.left] = std::max(credits[1][pair.left], unary[id]);
      credits[2][pair.right] = std::max(credits[2][pair.right], unary[id]);
    }
    for (double& value : credits[0]) value = up_half(value);
    for (int start = 0; start < 3; ++start) {
      auto& credit = credits[start];
      for (bool reverse : {false, true}) for (int position = 0; position < length; ++position) {
        const int i = reverse ? length - position - 1 : position;
        double needed = 0;
        for (int id : external_incident[i]) {
          const auto& pair = pairs[id];
          const int other = pair.left == i ? pair.right : pair.left;
          if (unary[id] > credit[other]) needed = std::max(needed, up_subtract(unary[id], credit[other]));
        }
        credit[i] = std::min(credit[i], needed);
      }
      for (const auto& block : internal) for (int id : block) {
        const auto& pair = pairs[id];
        adjusted[start][id] = up_subtract(up_subtract(unary[id], credit[pair.left]), credit[pair.right]);
      }
    }
  }

  bool valid_leaf() const {
    if (matching_only) return true;
    for (int id : chosen) {
      for (int row : row_ids[id]) {
        if (local_rows[row].external) continue;
        bool supported = false;
        for (int lower : local_rows[row].witnesses) if (selected[lower]) { supported = true; break; }
        if (!supported) return false;
      }
      if (no_lonely && !borrow_stack[id]) {
        bool stacked = false;
        for (int neighbor : stack_neighbors[id]) if (neighbor >= 0 && selected[neighbor]) stacked = true;
        if (!stacked) return false;
      }
    }
    return true;
  }

  std::array<double, 3> fallback(int block) const {
    std::array<double, 3> value{};
    const int size = positions[block].size();
    for (int start = 0; start < 3; ++start) {
      std::array<double, 12> maxima{};
      for (int id : internal[block]) {
        const double unary_cover = up_add(adjusted[start][id], internal_bonus[id]);
        maxima[local_position[pairs[id].left]] = std::max(maxima[local_position[pairs[id].left]], unary_cover);
        maxima[local_position[pairs[id].right]] = std::max(maxima[local_position[pairs[id].right]], unary_cover);
      }
      for (int i = 0; i < size; ++i) value[start] = up_add(value[start], up_half(maxima[i]));
    }
    return value;
  }

  std::array<double, 3> enumerate_window(int block, bool& exhausted, std::size_t& states) {
    const int size = positions[block].size();
    std::array<double, 3> best{};
    if (internal[block].empty()) { states = 1; return best; }
    std::function<void(int, unsigned, const std::array<double, 3>&)> visit;
    visit = [&](int position, unsigned used, const std::array<double, 3>& score) {
      if (state_budget && states >= state_budget) { exhausted = true; return; }
      ++states;
      while (position < size && (used & (1u << position))) ++position;
      if (position == size) {
        if (valid_leaf()) for (int start = 0; start < 3; ++start) best[start] = std::max(best[start], score[start]);
        return;
      }
      visit(position + 1, used, score);
      if (exhausted) return;
      for (int id : by_left[positions[block][position]]) {
        const auto& pair = pairs[id];
        const unsigned endpoint = 1u << local_position[pair.right];
        if (used & endpoint) continue;
        bool invalid = false;
        if (!matching_only) for (int other : chosen)
          if (pair.level == pairs[other].level && crossing(pair, pairs[other])) { invalid = true; break; }
        if (invalid) continue;
        std::array<double, 3> next = score;
        for (int start = 0; start < 3; ++start) next[start] = up_add(next[start], adjusted[start][id]);
        if (!matching_only) for (const auto& [other, contact_score] : internal_contacts[id]) if (selected[other])
          for (int start = 0; start < 3; ++start) next[start] = up_add(next[start], contact_score);
        selected[id] = 1; chosen.push_back(id);
        visit(position + 1, used | endpoint | (1u << position), next);
        chosen.pop_back(); selected[id] = 0;
        if (exhausted) return;
      }
    };
    visit(0, 0, {});
    return exhausted ? fallback(block) : best;
  }

  DDJointBoundResult evaluate() {
    if (ready) return cached;
    ready = true;
    if (overflow) return cached;
    make_credits();
    std::array<double, 3> total{};
    for (int start = 0; start < 3; ++start)
      for (double value : credits[start]) total[start] = up_add(total[start], value);
    for (int block = 0; block < static_cast<int>(internal.size()); ++block) {
      bool exhausted = false;
      std::size_t states = 0;
      const auto local = enumerate_window(block, exhausted, states);
      cached.states += states;
      if (exhausted) ++cached.fallback_windows;
      else ++cached.exact_windows;
      for (int start = 0; start < 3; ++start) total[start] = up_add(total[start], local[start]);
    }
    cached.upper_bound = std::min({total[0], total[1], total[2]});
    return cached;
  }
};

DDJointBound::DDJointBound(int length, const std::vector<DDPair>& pairs,
                         int levels, bool no_lonely_pairs,
                         const std::vector<unsigned char>& allowed,
                         const std::vector<DDRecoveryRow>& rows,
                         int width, int offset, bool matching_only, std::size_t state_budget, bool contact_clusters) {
  if (length < 0 || levels < 1 || width < 1 || width > 12 || offset < 0 || offset >= width ||
      (contact_clusters && offset != 0) || allowed.size() != pairs.size() ||
      pairs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("Invalid joint-bound dimensions");
  workspace_ = std::make_unique<Workspace>(length, pairs, levels, no_lonely_pairs, allowed, rows,
                                          width, offset, matching_only, state_budget, contact_clusters);
}
DDJointBound::~DDJointBound() = default;
DDJointBoundResult DDJointBound::evaluate() { return workspace_->evaluate(); }
