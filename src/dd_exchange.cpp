#include "dd_exchange.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace {
bool crosses(const DDPair& a, const DDPair& b) {
  return (a.left < b.left && b.left < a.right && a.right < b.right) ||
         (b.left < a.left && a.left < b.right && b.right < a.right);
}
std::uint64_t pair_key(int left, int right) {
  return (std::uint64_t(static_cast<unsigned>(left)) << 32) |
         static_cast<unsigned>(right);
}
}

struct DDLocalExchange::Workspace {
  struct Factor { int upper, lower; double score; };
  int length, levels;
  bool no_lonely;
  const std::vector<DDPair>& pairs;
  std::vector<unsigned char> allowed;
  std::vector<DDRecoveryRow> rows;
  std::vector<std::vector<int>> by_left, own_rows, incoming_rows, incident;
  std::vector<std::array<int, 2>> neighbors;
  std::vector<Factor> factors;
  // Reusable indices are reset only for entries touched by a local window.
  std::vector<int> local_slot, row_slot;
  std::vector<unsigned char> selected;
  std::vector<int> occupied, support;

  Workspace(int n, const std::vector<DDPair>& ps, int ls, bool lonely,
            const std::vector<unsigned char>& mask,
            const std::vector<DDRecoveryRow>& rs)
      : length(n), levels(ls), no_lonely(lonely), pairs(ps), allowed(mask),
        rows(rs) {
    if (n < 0 || ls < 1 || mask.size() != ps.size())
      throw std::invalid_argument("invalid exchange graph dimensions");
    const int m = static_cast<int>(ps.size());
    by_left.resize(n); own_rows.resize(m); incoming_rows.resize(m);
    incident.resize(m); neighbors.assign(m, {{-1, -1}});
    local_slot.assign(m, -1); row_slot.assign(rows.size(), -1);
    std::vector<std::unordered_map<std::uint64_t, int>> lookup(ls);
    for (int id = 0; id < m; ++id) {
      const auto& p = ps[id];
      if (p.left < 0 || p.left >= p.right || p.right >= n ||
          p.level < 0 || p.level >= ls || !std::isfinite(p.weight) || mask[id] > 1)
        throw std::invalid_argument("invalid exchange candidate");
      if (!lookup[p.level].emplace(pair_key(p.left, p.right), id).second)
        throw std::invalid_argument("duplicate exchange candidate");
      by_left[p.left].push_back(id);
    }
    for (int id = 0; id < m; ++id) {
      const auto& p = ps[id];
      const auto& table = lookup[p.level];
      auto inner = table.find(pair_key(p.left + 1, p.right - 1));
      auto outer = table.find(pair_key(p.left - 1, p.right + 1));
      if (inner != table.end()) neighbors[id][0] = inner->second;
      if (outer != table.end()) neighbors[id][1] = outer->second;
    }
    std::vector<std::vector<unsigned char>> seen_level(m);
    std::vector<int> seen_contact(m, -1);
    for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
      const auto& row = rows[r];
      if (row.upper < 0 || row.upper >= m || ps[row.upper].level == 0)
        throw std::invalid_argument("invalid exchange witness row");
      own_rows[row.upper].push_back(r);
      int lower_level = -1;
      for (const auto& contact : row.contacts) {
        const int lower = contact.first;
        if (lower < 0 || lower >= m || !std::isfinite(contact.second) ||
            ps[lower].level >= ps[row.upper].level ||
            !crosses(ps[row.upper], ps[lower]) || seen_contact[lower] == r)
          throw std::invalid_argument("invalid exchange witness contact");
        seen_contact[lower] = r;
        if (lower_level < 0) lower_level = ps[lower].level;
        if (lower_level != ps[lower].level)
          throw std::invalid_argument("mixed exchange witness levels");
        incoming_rows[lower].push_back(r);
        int f = static_cast<int>(factors.size());
        factors.push_back({row.upper, lower, contact.second});
        incident[row.upper].push_back(f); incident[lower].push_back(f);
      }
      if (lower_level >= 0) {
        auto& seen = seen_level[row.upper];
        if (seen.empty()) seen.assign(ps[row.upper].level, 0);
        if (seen[lower_level])
          throw std::invalid_argument("duplicate exchange witness level");
        seen[lower_level] = 1;
      } else if (allowed[row.upper]) {
        throw std::invalid_argument("empty exchange witness row");
      }
    }
    for (int id = 0; id < m; ++id)
      if (allowed[id] && ps[id].level > 0) {
        const auto& seen = seen_level[id];
        if (seen.size() != static_cast<std::size_t>(ps[id].level) ||
            std::find(seen.begin(), seen.end(), 0) != seen.end())
          throw std::invalid_argument("missing exchange witness level");
      }
  }

  long double objective() const {
    long double value = 0;
    for (int id = 0; id < static_cast<int>(pairs.size()); ++id)
      if (selected[id]) value += pairs[id].weight;
    for (const auto& f : factors)
      if (selected[f.upper] && selected[f.lower]) value += f.score;
    return value;
  }
  void initialize(const std::vector<unsigned char>& incumbent) {
    if (incumbent.size() != pairs.size())
      throw std::invalid_argument("invalid exchange incumbent size");
    selected = incumbent; occupied.assign(length, -1);
    support.assign(rows.size(), 0);
    for (int id = 0; id < static_cast<int>(pairs.size()); ++id) {
      if (selected[id] > 1 || (selected[id] && !allowed[id]))
        throw std::invalid_argument("invalid exchange incumbent mask");
      if (!selected[id]) continue;
      const auto& p = pairs[id];
      if (occupied[p.left] >= 0 || occupied[p.right] >= 0)
        throw std::invalid_argument("exchange incumbent reuses a base");
      occupied[p.left] = occupied[p.right] = id;
      for (int r : incoming_rows[id]) ++support[r];
    }
    std::vector<std::vector<int>> ends(levels);
    for (int base = 0; base < length; ++base) {
      int id = occupied[base];
      if (id < 0 || pairs[id].left != base) continue;
      const auto& p = pairs[id];
      auto& stack = ends[p.level];
      while (!stack.empty() && stack.back() < base) stack.pop_back();
      if (!stack.empty() && p.right > stack.back())
        throw std::invalid_argument("crossing exchange incumbent level");
      stack.push_back(p.right);
    }
    for (int r = 0; r < static_cast<int>(rows.size()); ++r)
      if (selected[rows[r].upper] && !support[r])
        throw std::invalid_argument("unsupported exchange incumbent");
    if (no_lonely)
      for (int id = 0; id < static_cast<int>(pairs.size()); ++id)
        if (selected[id] && !stacked(id))
          throw std::invalid_argument("lonely exchange incumbent pair");
  }
  bool stacked(int id) const {
    for (int neighbor : neighbors[id])
      if (neighbor >= 0 && selected[neighbor]) return true;
    return false;
  }
  void change_global(int id, bool take) {
    selected[id] = take;
    const auto& p = pairs[id];
    occupied[p.left] = occupied[p.right] = take ? id : -1;
    for (int r : incoming_rows[id]) support[r] += take ? 1 : -1;
  }

  void window(int left, int right, std::size_t budget, DDExchangeResult& result) {
    std::vector<int> old, local;
    for (int base = left; base < right; ++base)
      for (int id : by_left[base])
        if (pairs[id].right < right && selected[id]) old.push_back(id);
    for (int id : old) change_global(id, false);
    // Only frozen arcs with an endpoint in this window can cross an internal
    // arc. A remote frozen arc is either nested around it or disjoint.
    std::vector<int> boundary;
    for (int base = left; base < right; ++base)
      if (occupied[base] >= 0 &&
          std::find(boundary.begin(), boundary.end(), occupied[base]) == boundary.end())
        boundary.push_back(occupied[base]);
    for (int base = left; base < right; ++base)
      for (int id : by_left[base]) {
        const auto& p = pairs[id];
        if (!allowed[id] || p.right >= right ||
            occupied[p.left] >= 0 || occupied[p.right] >= 0) continue;
        bool compatible = true;
        for (int frozen : boundary)
          if (p.level == pairs[frozen].level && crosses(p, pairs[frozen])) {
            compatible = false; break;
          }
        if (compatible) {
          local_slot[id] = static_cast<int>(local.size()); local.push_back(id);
        }
      }
    const int count = static_cast<int>(local.size());
    std::vector<long double> linear(count);
    std::vector<std::vector<std::pair<int, double>>> interactions(count);
    std::vector<std::vector<int>> choices(right - left);
    for (int k = 0; k < count; ++k) {
      int id = local[k]; linear[k] = pairs[id].weight;
      choices[pairs[id].left - left].push_back(k);
      for (int f : incident[id]) {
        const auto& factor = factors[f];
        int other = factor.upper == id ? factor.lower : factor.upper;
        int slot = local_slot[other];
        if (slot >= 0) interactions[k].push_back({slot, factor.score});
        else if (selected[other]) linear[k] += factor.score;
      }
    }
    std::vector<unsigned char> was_old(count, 0);
    for (int id : old) {
      if (local_slot[id] < 0)
        throw std::logic_error("feasible exchange incumbent excluded from window");
      was_old[local_slot[id]] = 1;
    }
    long double old_value = 0;
    for (int k = 0; k < count; ++k) if (was_old[k]) {
      old_value += linear[k];
      for (const auto& c : interactions[k])
        if (c.first > k && was_old[c.first]) old_value += c.second;
    }
    long double best_value = old_value;
    std::vector<int> best = old, chosen;
    for (auto& alternatives : choices)
      std::sort(alternatives.begin(), alternatives.end(), [&](int a, int b) {
        if (was_old[a] != was_old[b]) return was_old[a] > was_old[b];
        if (linear[a] != linear[b]) return linear[a] > linear[b];
        return local[a] < local[b];
      });

    std::vector<int> relevant_rows;
    auto add_row = [&](int r) {
      if (row_slot[r] < 0) {
        row_slot[r] = static_cast<int>(relevant_rows.size());
        relevant_rows.push_back(r);
      }
    };
    for (int id : local) {
      for (int r : own_rows[id]) add_row(r);
      for (int r : incoming_rows[id])
        if (selected[rows[r].upper] || local_slot[rows[r].upper] >= 0) add_row(r);
    }
    std::vector<int> local_support(relevant_rows.size());
    std::vector<std::vector<int>> lower_rows(count), witnesses(relevant_rows.size());
    for (int s = 0; s < static_cast<int>(relevant_rows.size()); ++s)
      local_support[s] = support[relevant_rows[s]];
    for (int k = 0; k < count; ++k)
      for (int r : incoming_rows[local[k]]) if (row_slot[r] >= 0) {
        lower_rows[k].push_back(row_slot[r]);
        witnesses[row_slot[r]].push_back(k);
      }
    std::vector<int> stack_affected = local;
    if (no_lonely)
      for (int id : local)
        for (int neighbor : neighbors[id])
          if (neighbor >= 0 && selected[neighbor] &&
              std::find(stack_affected.begin(), stack_affected.end(), neighbor) == stack_affected.end())
            stack_affected.push_back(neighbor);
    auto compatible_chosen = [&](int id) {
      for (int other : chosen)
        if (pairs[id].level == pairs[other].level && crosses(pairs[id], pairs[other]))
          return false;
      return true;
    };
    auto available = [&](int k, int next) {
      const auto& p = pairs[local[k]];
      return p.left >= next && occupied[p.left] < 0 && occupied[p.right] < 0 &&
             compatible_chosen(local[k]);
    };
    // These prunes use only necessary discrete constraints, not a rounded
    // objective bound. A complete search therefore remains an integer proof.
    auto can_complete = [&](int next) {
      for (int s = 0; s < static_cast<int>(relevant_rows.size()); ++s) {
        int upper = rows[relevant_rows[s]].upper;
        if (!selected[upper] || local_support[s]) continue;
        bool possible = false;
        for (int k : witnesses[s])
          if (available(k, next)) { possible = true; break; }
        if (!possible) return false;
      }
      if (no_lonely)
        for (int id : stack_affected) {
          if (!selected[id] || stacked(id)) continue;
          bool possible = false;
          for (int neighbor : neighbors[id])
            if (neighbor >= 0 && local_slot[neighbor] >= 0 &&
                available(local_slot[neighbor], next)) { possible = true; break; }
          if (!possible) return false;
        }
      return true;
    };
    std::size_t visited = 0;
    bool interrupted = false;
    std::function<void(int, long double)> dfs = [&](int next, long double value) {
      if (budget && visited >= budget) { interrupted = true; return; }
      ++visited;
      while (next < right && occupied[next] >= 0) ++next;
      if (!can_complete(next)) return;
      if (next == right) {
        if (value > best_value) { best_value = value; best = chosen; }
        return;
      }
      for (int k : choices[next - left]) {
        int id = local[k]; const auto& p = pairs[id];
        if (occupied[p.right] >= 0 || !compatible_chosen(id)) continue;
        long double delta = linear[k];
        for (const auto& c : interactions[k])
          if (selected[local[c.first]]) delta += c.second;
        selected[id] = 1; occupied[p.left] = occupied[p.right] = id;
        for (int s : lower_rows[k]) ++local_support[s];
        chosen.push_back(id);
        dfs(next + 1, value + delta);
        chosen.pop_back();
        for (int s : lower_rows[k]) --local_support[s];
        selected[id] = 0; occupied[p.left] = occupied[p.right] = -1;
        if (interrupted) break;
      }
      if (!interrupted) dfs(next + 1, value);
    };
    dfs(left, 0);
    for (int id : best) change_global(id, true);
    ++result.windows; result.states_visited += visited;
    if (best_value > old_value) ++result.improved_windows;
    if (interrupted) ++result.budget_windows;
    if (left == 0 && right == length && !interrupted) result.globally_exact = true;
    for (int id : local) local_slot[id] = -1;
    for (int r : relevant_rows) row_slot[r] = -1;
  }
};

DDLocalExchange::DDLocalExchange(int length, const std::vector<DDPair>& pairs,
    int levels, bool no_lonely_pairs, const std::vector<unsigned char>& allowed,
    const std::vector<DDRecoveryRow>& rows)
    : workspace_(new Workspace(length, pairs, levels, no_lonely_pairs, allowed, rows)) {}
DDLocalExchange::~DDLocalExchange() = default;

DDExchangeResult DDLocalExchange::improve(const std::vector<unsigned char>& incumbent,
    int width, int passes, std::size_t state_budget) {
  if (width < 1 || width > 12 || passes < 1 || passes > 4)
    throw std::invalid_argument("exchange width must be 1..12 and passes 1..4");
  auto& w = *workspace_;
  w.initialize(incumbent);
  DDExchangeResult result;
  if (!w.length) result.globally_exact = true;
  for (int pass = 0; pass < passes; ++pass) {
    // Search a complete small graph as one window on every pass. For longer
    // sequences, the shifted tiling also visits its clipped edge windows.
    int offset = w.length <= width ? 0 : (pass % 2) * (width / 2);
    for (int start = -offset; start < w.length; start += width)
      w.window(std::max(0, start), std::min(w.length, start + width),
               state_budget, result);
  }
  result.selected = w.selected;
  result.objective = static_cast<double>(w.objective());
  return result;
}
