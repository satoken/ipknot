// SPDX-License-Identifier: GPL-3.0-or-later
// Structural constraints extracted from IPknot, Copyright (C) 2010 Kengo Sato.
#include "ilp_model.h"

namespace ipknot::detail {
void add_level_constraints(IP &ip, const PairVariables &v_l) {
  // constraint 2: disallow pseudoknots in x[lv]
  for (auto lv = 0; lv != v_l.size(); ++lv)
    for (auto i = 0; i < v_l[lv].size(); ++i)
      for (auto [j, v_ij] : v_l[lv][i])
        for (auto k = i + 1; k < j; ++k)
          for (auto [l, v_kl] : v_l[lv][k])
            if (j < l) {
              auto row = ip.make_constraint(IP::UP, 0, 1);
              ip.add_constraint(row, v_ij, 1);
              ip.add_constraint(row, v_kl, 1);
            }

  // constraint 3: any x[t]_kl must be pseudoknotted with x[u]_ij for t>u
  for (auto lv = 1; lv != v_l.size(); ++lv)
    for (auto k = 0; k < v_l[lv].size(); ++k)
      for (auto [l, v_kl] : v_l[lv][k])
        for (auto plv = 0; plv != lv; ++plv) {
          int row = ip.make_constraint(IP::LO, 0, 0);
          ip.add_constraint(row, v_kl, -1);
          for (auto i = 0; i < k; ++i)
            for (auto [j, v_ij] : v_l[plv][i])
              if (k < j && j < l)
                ip.add_constraint(row, v_ij, 1);

          for (auto i = k + 1; i < l; ++i)
            for (auto [j, v_ij] : v_l[plv][i])
              if (l < j)
                ip.add_constraint(row, v_ij, 1);
        }
}

void add_stacking_constraints(IP &ip, const PairVariables &v_l, const PairVariables &v_r,
                              int max_neighbor_distance) {
  const auto L = v_l.front().size();
  for (auto lv = 0; lv != v_l.size(); ++lv) {
    // upstream
    for (auto i = 0; i < L; ++i) {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (auto [j, v_ji] : v_r[lv][i])
        ip.add_constraint(row, v_ji, -1);
      for (int d = 1; d <= max_neighbor_distance; ++d) {
        if (i >= static_cast<unsigned int>(d))
          for (auto [j, v_ji] : v_r[lv][i - d])
            ip.add_constraint(row, v_ji, 1);
        if (i + d < L)
          for (auto [j, v_ji] : v_r[lv][i + d])
            ip.add_constraint(row, v_ji, 1);
      }
    }

    // downstream
    for (auto i = 0; i < L; ++i) {
      auto row = ip.make_constraint(IP::LO, 0, 0);
      for (auto [j, v_ij] : v_l[lv][i])
        ip.add_constraint(row, v_ij, -1);
      for (int d = 1; d <= max_neighbor_distance; ++d) {
        if (i >= static_cast<unsigned int>(d))
          for (auto [j, v_ij] : v_l[lv][i - d])
            ip.add_constraint(row, v_ij, 1);
        if (i + d < L)
          for (auto [j, v_ij] : v_l[lv][i + d])
            ip.add_constraint(row, v_ij, 1);
      }
    }
  }
}

DecodeResult decode_ilp(int length, const std::vector<DDPair> &pairs, int levels,
                        bool no_lonely_pairs, int threads) {
  DecodeResult result;
  result.bpseq.assign(length, -1);
  result.levels.assign(length, -1);
  if (pairs.empty())
    return result;
  IP ip(IP::MAX, threads);
  PairVariables left(levels, std::vector<std::vector<std::pair<unsigned int, int>>>(length));
  PairVariables right = left;
  std::vector<int> columns;
  for (const auto &pair : pairs) {
    const int col = ip.make_variable(pair.weight);
    columns.push_back(col);
    left[pair.level][pair.left].emplace_back(pair.right, col);
    right[pair.level][pair.right].emplace_back(pair.left, col);
  }
  ip.update();
  for (int i = 0; i < length; ++i) {
    const int row = ip.make_constraint(IP::UP, 0, 1);
    for (int level = 0; level < levels; ++level) {
      for (const auto &[j, col] : left[level][i])
        ip.add_constraint(row, col, 1);
      for (const auto &[j, col] : right[level][i])
        ip.add_constraint(row, col, 1);
    }
  }
  add_level_constraints(ip, left);
  if (no_lonely_pairs)
    add_stacking_constraints(ip, left, right);
  ip.update();
  ip.solve();
  for (std::size_t id = 0; id < pairs.size(); ++id) {
    if (ip.get_value(columns[id]) <= 0.5)
      continue;
    const auto &pair = pairs[id];
    result.bpseq[pair.left] = pair.right;
    result.bpseq[pair.right] = pair.left;
    result.levels[pair.left] = result.levels[pair.right] = pair.level;
    result.objective += pair.weight;
  }
  return result;
}
} // namespace ipknot::detail
