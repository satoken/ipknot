#include "dual_decomposition.h"
#include "dd_parallel.h"
#include "dd_constrained.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

void check(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

// An independent dense interval recurrence, without candidate sparsification,
// suffix queries, beams, or trace arrows. Covers arbitrary signed weights.
double reference(int n, const std::vector<DDPair>& pairs,
    const std::vector<double>& weights, const std::vector<unsigned char>& mask, bool stacked) {
  std::vector<std::vector<double>> f(n, std::vector<double>(n));
  std::vector<std::vector<double>> p(n, std::vector<double>(n, -INFINITY));
  std::vector<std::vector<int>> id(n, std::vector<int>(n, -1));
  for (std::size_t k = 0; k < pairs.size(); ++k) if (mask[k]) id[pairs[k].left][pairs[k].right] = k;
  for (int span = 1; span < n; ++span) for (int i = 0; i + span < n; ++i) {
    const int j = i + span;
    f[i][j] = f[i][j - 1];
    for (int k = i; k < j; ++k) if (id[k][j] >= 0) {
      const double inside = k + 1 < j ? f[k + 1][j - 1] : 0;
      const double child = k + 1 < j - 1 ? p[k + 1][j - 1] : -INFINITY;
      p[k][j] = weights[id[k][j]] + (stacked ? std::max(inside, child) : inside);
      const double complete = stacked ? weights[id[k][j]] + child : p[k][j];
      f[i][j] = std::max(f[i][j], (k > i ? f[i][k - 1] : 0) + complete);
    }
  }
  return f[0][n - 1];
}

int main() {
  std::mt19937 rng(38433200);
  for (int trial = 0; trial < 80; ++trial) {
    const int n = 10 + rng() % 26;
    std::vector<DDPair> pairs;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j)
      if (rng() % 3) pairs.push_back({i, j, 0, 0});
    std::shuffle(pairs.begin(), pairs.end(), rng);
    for (bool stacked : {false, true}) {
      DDNussinov decoder(n, pairs, 0, 0, stacked);
      for (int reuse = 0; reuse < 4; ++reuse) {
        std::vector<double> weights;
        std::vector<unsigned char> mask;
        for (std::size_t k = 0; k < pairs.size(); ++k) {
          weights.push_back((int(rng() % 25) - 12) / 7.0);
          mask.push_back(rng() % 4 != 0);
        }
        std::vector<int> selected;
        const double actual = decoder.decode(weights, mask, selected);
        const double expected = reference(n, pairs, weights, mask, stacked);
        check(std::abs(actual - expected) < 1e-9, "Sparsified exact DP differs from dense reference");
        check(decoder.pruned_states() == 0, "Exact sparsification counted as beam pruning");
        double score = 0;
        for (int id : selected) { check(mask[id], "Masked pair in traceback"); score += weights[id]; }
        check(std::abs(actual - score) < 1e-9, "Sparsified traceback has wrong score");
      }
    }
  }
  std::vector<DDPair> dense;
  for (int i = 0; i < 30; ++i) for (int j = i + 1; j < 30; ++j) dense.push_back({i, j, 0, 1});
  DDNussinov decoder(30, dense, 0, 0, false);
  std::vector<int> selected;
  std::vector<double> weights(dense.size(), 1);
  std::vector<unsigned char> mask(dense.size(), 1);
  check(decoder.decode(weights, mask, selected) == 15 && decoder.sparsified_candidates() > 100,
        "Decomposable candidates were not sparsified");
  std::fill(mask.begin(), mask.end(), 0);
  check(decoder.decode(weights, mask, selected) == 0 && decoder.sparsified_candidates() == 0,
        "Sparsification state leaked between decodes");
  const std::vector<DDPair> long_pairs{{0, 20004, 0, 2}};
  DDNussinov long_sparse(20005, long_pairs, 0, 0, false);
  check(long_sparse.decode({2}, {1}, selected) == 2 && selected == std::vector<int>{0},
        "Long sparse exact chart lost a distant pair");
  check(long_sparse.decode({-2}, {1}, selected) == 0 && selected.empty(),
        "Long reused exact chart retained a negative component");

  std::vector<DDPair> pairs;
  for (int i = 0; i < 24; ++i) for (int j = i + 1; j < 24; ++j)
    if (rng() % 4 == 0) for (int level = 0; level < 3; ++level)
      pairs.push_back({i, j, level, (int(rng() % 23) - 5) / 9.0});
  std::vector<int> thread_counts{0, 1};
#ifdef _OPENMP
  thread_counts = {0, 1, 2, 8};
#else
  {
    DDOptions invalid; invalid.threads = 2;
    bool rejected = false;
    try { invalid.validate(); } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Non-OpenMP build silently accepted parallel execution");
  }
#endif
  for (int mode : {0, 1, 2}) for (bool stacked : {false, true}) {
    DDOptions options;
    options.beam = mode == 0 ? 0 : 5;
    options.improved_beam = mode == 2;
    options.max_iterations = 25;
    options.recovery_every = 3;
    options.bound_block = 8;
    const auto serial = solve_dual_decomposition(24, pairs, 3, stacked, options);
    for (int threads : thread_counts) {
      options.threads = threads;
      const auto parallel = solve_dual_decomposition(24, pairs, 3, stacked, options);
      check(serial.bpseq == parallel.bpseq && serial.levels == parallel.levels &&
            serial.objective == parallel.objective && serial.upper_bound == parallel.upper_bound &&
            serial.iterations == parallel.iterations && serial.stop_reason == parallel.stop_reason,
            "Parallel DD differs from serial execution");
    }
  }

  IPModel model; IP ip(model);
  std::vector<int> columns;
  for (const auto& pair : pairs) columns.push_back(ip.make_variable(pair.weight));
  for (int base = 0; base < 24; ++base) {
    const auto row = ip.make_constraint(IP::UP, 0, 1);
    for (std::size_t id = 0; id < pairs.size(); ++id)
      if (pairs[id].left == base || pairs[id].right == base) ip.add_constraint(row, columns[id], 1);
  }
  for (int mode : {0, 1, 2}) {
    DDOptions options; options.max_iterations = 15; options.constraint_states = 100;
    options.beam = mode == 0 ? 0 : 5; options.improved_beam = mode == 2;
    auto serial_model = model;
    const auto serial = solve_constrained_dd(24, pairs, columns, 3, serial_model, options);
    for (int threads : thread_counts) {
      options.threads = threads;
      auto parallel_model = model;
      const auto parallel = solve_constrained_dd(24, pairs, columns, 3, parallel_model, options);
      check(serial_model.solution == parallel_model.solution && serial.objective == parallel.objective &&
            serial.upper_bound == parallel.upper_bound && serial.iterations == parallel.iterations &&
            serial.dp_pruned_states == parallel.dp_pruned_states,
            "Parallel constrained DD differs from serial execution");
    }
  }
  // A worker failure is rethrown safely after all other level tasks finish.
  std::vector<std::unique_ptr<DDNussinov>> decoders;
  for (int level = 0; level < 3; ++level)
    decoders.push_back(std::make_unique<DDNussinov>(24, pairs, level, 5, false));
  std::vector<std::vector<int>> outputs;
  std::vector<double> scores;
  bool rejected = false;
  try { dd_decode_levels(decoders, {}, {}, outputs, scores, thread_counts.back()); }
  catch (const std::invalid_argument&) { rejected = true; }
  check(rejected, "Worker exception did not reach caller");
  std::cout << "Sparse exact reference and deterministic parallel DD checks passed\n";
}
