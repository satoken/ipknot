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

void check_trace(int n, const std::vector<DDPair>& pairs, const std::vector<double>& weights,
    const std::vector<unsigned char>& mask, const std::vector<int>& selected, double score, bool stacked) {
  std::vector<int> partner(n, -1);
  double sum = 0;
  for (int id : selected) {
    const auto& p = pairs.at(id);
    check(mask[id] && partner[p.left] < 0 && partner[p.right] < 0, "Invalid traceback endpoint or mask");
    partner[p.left] = p.right; partner[p.right] = p.left;
    sum += weights[id];
  }
  for (int id : selected) {
    const auto& p = pairs[id];
    for (int k = p.left + 1; k < p.right; ++k)
      check(partner[k] < 0 || (partner[k] > p.left && partner[k] < p.right), "Crossing traceback");
    if (stacked) check((p.left && p.right + 1 < n && partner[p.left - 1] == p.right + 1) ||
        (p.left + 1 < p.right && partner[p.left + 1] == p.right - 1), "Lonely pair in stacked traceback");
  }
  check(std::abs(sum - score) < 1e-9, "Traceback score mismatch");
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
  const std::vector<DDPair> empty;
  for (int n : {0, 1}) {
    DDNussinov tiny(n, empty, 0, 0, true, false, 1, true);
    std::vector<int> selected{42};
    check(tiny.decode({}, {}, selected) == 0 && selected.empty(), "Empty wavefront boundary failed");
  }
  std::mt19937 rng(38433200);
  for (int trial = 0; trial < 80; ++trial) {
    const int n = 10 + rng() % 26;
    std::vector<DDPair> pairs;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j)
      if (rng() % 3) pairs.push_back({i, j, 0, 0});
    std::shuffle(pairs.begin(), pairs.end(), rng);
    for (bool stacked : {false, true}) {
      DDNussinov decoder(n, pairs, 0, 0, stacked);
      DDNussinov wavefront(n, pairs, 0, 0, stacked, false, 1, true);
      int workers = 1;
#ifdef _OPENMP
      workers = 2;
#endif
      DDNussinov parallel_wavefront(n, pairs, 0, 0, stacked, false, workers, true);
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
        std::vector<int> wave_selected, parallel_selected;
        const double wave_score = wavefront.decode(weights, mask, wave_selected);
        const double parallel_score = parallel_wavefront.decode(weights, mask, parallel_selected);
        check(std::abs(wave_score - expected) < 1e-9, "Wavefront differs from independent reference");
        check(wave_score == parallel_score && wave_selected == parallel_selected,
              "Wavefront changed with thread count");
        check(wavefront.pruned_states() == 0 && parallel_wavefront.used_workers() <= workers,
              "Wavefront pruning or worker count mismatch");
        check_trace(n, pairs, weights, mask, wave_selected, wave_score, stacked);
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

  // Wide, dense columns exercise the actual private-scatter/merge path, with
  // shuffled IDs, tied/signed weights, masks, stacking and workspace reuse.
  std::vector<DDPair> wide;
  constexpr int wide_n = 384;
  for (int i = 0; i < wide_n; ++i) for (int j = i + 1; j < wide_n; ++j)
    wide.push_back({i, j, 0, 0});
  std::shuffle(wide.begin(), wide.end(), rng);
  for (bool improved : {false, true}) for (bool stacked : {false, true}) {
    DDNussinov serial(wide_n, wide, 0, 256, stacked, improved);
    int workers = 1;
#ifdef _OPENMP
    workers = 4;
#endif
    DDNussinov parallel(wide_n, wide, 0, 256, stacked, improved, workers);
    for (int reuse = 0; reuse < 3; ++reuse) {
      std::vector<double> wide_weights(wide.size());
      std::vector<unsigned char> wide_mask(wide.size());
      for (std::size_t id = 0; id < wide.size(); ++id) {
        wide_weights[id] = (int(rng() % 13) - 6) / 4.;
        wide_mask[id] = reuse != 2 && rng() % 4;
      }
      std::vector<int> expected, actual;
      const double score = serial.decode(wide_weights, wide_mask, expected);
      check(score == parallel.decode(wide_weights, wide_mask, actual) && expected == actual &&
            serial.pruned_states() == parallel.pruned_states(), "Column parallel beam changed serial results");
      check(parallel.used_workers() <= workers && (reuse != 2 || parallel.used_workers() == 1),
            "Column parallel worker/fallback mismatch");
      check_trace(wide_n, wide, wide_weights, wide_mask, actual, score, stacked);
    }
  }

  std::vector<DDPair> pairs;
  for (int i = 0; i < 24; ++i) for (int j = i + 1; j < 24; ++j)
    if (rng() % 4 == 0) for (int level = 0; level < 3; ++level)
      pairs.push_back({i, j, level, (int(rng() % 23) - 5) / 9.0});
  std::vector<int> thread_counts{0, 1};
#ifdef _OPENMP
  thread_counts = {0, 1, 2, 8};
#else
  for (bool inner : {false, true}) {
    DDOptions invalid; invalid.threads = inner ? 1 : 2; invalid.dp_threads = inner ? 2 : 1;
    bool rejected = false;
    try { invalid.validate(); } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Non-OpenMP build silently accepted parallel execution");
  }
#endif
  for (int mode : {0, 1, 2, 3}) for (bool stacked : {false, true}) {
    DDOptions options;
    options.beam = mode == 0 ? 0 : 5;
    options.improved_beam = mode == 2;
    options.wavefront_dp = mode == 3;
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
      options.dp_threads = threads;
      const auto inner_parallel = solve_dual_decomposition(24, pairs, 3, stacked, options);
      check(serial.bpseq == inner_parallel.bpseq && serial.levels == inner_parallel.levels &&
            serial.objective == inner_parallel.objective && serial.upper_bound == inner_parallel.upper_bound &&
            serial.iterations == inner_parallel.iterations && serial.stop_reason == inner_parallel.stop_reason,
            "Inner-thread DD differs from serial execution");
      // Also propagate the settings to recovery when workspaces are unshared.
      options.recovery_share = false;
      const auto unshared = solve_dual_decomposition(24, pairs, 3, stacked, options);
      check(inner_parallel.bpseq == unshared.bpseq && inner_parallel.upper_bound == unshared.upper_bound,
            "Unshared parallel recovery changed results");
      options.recovery_share = true;
      options.dp_threads = 1;
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
  for (int mode : {0, 1, 2, 3}) {
    DDOptions options; options.max_iterations = 15; options.constraint_states = 100;
    options.beam = mode == 0 ? 0 : 5; options.improved_beam = mode == 2;
    options.wavefront_dp = mode == 3;
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
      options.dp_threads = threads;
      parallel_model = model;
      const auto inner = solve_constrained_dd(24, pairs, columns, 3, parallel_model, options);
      check(serial_model.solution == parallel_model.solution && serial.objective == inner.objective &&
            serial.upper_bound == inner.upper_bound && serial.iterations == inner.iterations &&
            serial.dp_pruned_states == inner.dp_pruned_states, "Inner-thread constrained DD changed results");
      options.dp_threads = 1;
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
