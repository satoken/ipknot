#include <ipknot/decoder.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <random>
#include <stdexcept>

using namespace ipknot;
void check(bool condition) {
  if (!condition)
    throw std::runtime_error("Decoder API test failed");
}
template <class F> void invalid(F operation) {
  try {
    operation();
  } catch (const std::invalid_argument &) {
    return;
  }
  throw std::runtime_error("Invalid input was accepted");
}

bool crosses(const PairProbability &a, const PairProbability &b) {
  return (a.left < b.left && b.left < a.right && a.right < b.right) ||
         (b.left < a.left && a.left < b.right && b.right < a.right);
}

// Independent enumeration of physical matchings and their level assignments.
double exhaustive(int n, const std::vector<PairProbability> &pairs, const DecoderOptions &options) {
  std::vector<int> used(n), assignment(pairs.size(), -1);
  double best = 0;
  std::function<void(std::size_t, double)> visit = [&](std::size_t id, double score) {
    if (id == pairs.size()) {
      for (std::size_t i = 0; i < pairs.size(); ++i)
        if (assignment[i] >= 0) {
          for (int lower = 0; lower < assignment[i]; ++lower) {
            bool support = false;
            for (std::size_t j = 0; j < pairs.size(); ++j)
              if (assignment[j] == lower && crosses(pairs[i], pairs[j]))
                support = true;
            if (!support)
              return;
          }
        }
      best = std::max(best, score);
      return;
    }
    visit(id + 1, score);
    const auto &pair = pairs[id];
    if (used[pair.left] || used[pair.right])
      return;
    for (int level = 0; level < static_cast<int>(options.thresholds.size()); ++level) {
      if (pair.probability <= options.thresholds[level])
        continue;
      bool planar = true;
      for (std::size_t j = 0; j < id; ++j)
        if (assignment[j] == level && crosses(pair, pairs[j]))
          planar = false;
      if (!planar)
        continue;
      used[pair.left] = used[pair.right] = 1;
      assignment[id] = level;
      visit(id + 1,
            score + (pair.probability - options.thresholds[level]) * options.weights[level]);
      assignment[id] = -1;
      used[pair.left] = used[pair.right] = 0;
    }
  };
  visit(0, 0);
  return best;
}

int main() {
  DecoderOptions options;
  options.thresholds = {0.2f, 0.1f};
  const std::vector<PairProbability> pairs{{0, 7, .9f}, {1, 6, .9f}, {3, 10, .9f}, {4, 9, .9f}};
  SparsePosterior upper(12), symmetric(12);
  for (const auto &pair : pairs) {
    upper[pair.left + 1].emplace_back(pair.right + 1, pair.probability);
    symmetric[pair.left + 1].emplace_back(pair.right + 1, pair.probability);
    symmetric[pair.right + 1].emplace_back(pair.left + 1, pair.probability);
  }
  for (int mode = 0; mode < 4; ++mode) {
    if (mode == 3 && !Decoder::ilp_available())
      continue;
    options.backend = mode == 3 ? Backend::ILP : Backend::DD;
    options.dd.nussinov_dp = mode == 1 || mode == 2;
    options.dd.crossing_beam = mode == 2 ? 0 : 100;
    options.dd.witnesses = mode == 2 ? 0 : 16;
    Decoder decoder(options);
    const auto result = decoder.decode(11, pairs);
    check(result.bpseq == std::vector<int>({7, 6, -1, 10, 9, -1, 1, 0, -1, 4, 3}));
    check(std::abs(result.objective - 1.5) < 1e-6);
    check(result.dd.has_value() == (mode != 3));
    for (int i = 0; i < 11; ++i)
      if (result.bpseq[i] >= 0)
        check(result.levels[i] == result.levels[result.bpseq[i]]);
    check(decoder.decode(upper).bpseq == result.bpseq);
    check(decoder.decode(symmetric).bpseq == result.bpseq);
    check(decoder.decode(11, pairs).bpseq == result.bpseq); // repeated calls
    check(decoder.decode(0, {}).bpseq.empty());
    check(decoder.decode(SparsePosterior(1)).bpseq.empty());
    check(decoder.decode(11, {}).bpseq == std::vector<int>(11, -1));
  }

  Decoder decoder;
  // A single-level decoder uses the same API, and isolates the stacking toggle
  // and the strict p > threshold candidate rule.
  DecoderOptions single;
  single.thresholds = {0.5f};
  single.no_lonely_pairs = false;
  check(Decoder(single).decode(4, {{0, 3, 0.8f}}).bpseq == std::vector<int>({3, -1, -1, 0}));
  check(Decoder(single).decode(4, {{0, 3, 0.5f}}).bpseq == std::vector<int>(4, -1));
  single.no_lonely_pairs = true;
  check(Decoder(single).decode(4, {{0, 3, 0.8f}}).bpseq == std::vector<int>(4, -1));
  invalid([&] { decoder.decode(-1, {}); });
  invalid([&] { decoder.decode(5, {{0, 5, .5f}}); });
  invalid([&] { decoder.decode(5, {{2, 2, .5f}}); });
  invalid([&] { decoder.decode(5, {{3, 2, .5f}}); });
  invalid([&] { decoder.decode(5, {{-1, 2, .5f}}); });
  invalid([&] { decoder.decode(5, {{0, 2, 1.1f}}); });
  invalid([&] { decoder.decode(5, {{0, 2, -.1f}}); });
  invalid([&] { decoder.decode(5, {{0, 2, std::numeric_limits<float>::quiet_NaN()}}); });
  invalid([&] { decoder.decode(5, {{0, 2, .1f}, {0, 2, .1f}}); }); // below threshold too
  invalid([&] { decoder.decode(SparsePosterior{}); });
  invalid([&] {
    auto bad = upper;
    bad[0].emplace_back(1, .5f);
    decoder.decode(bad);
  });
  invalid([&] {
    auto bad = upper;
    bad[1].emplace_back(8, .5f);
    decoder.decode(bad);
  });
  invalid([&] {
    auto bad = upper;
    bad[1].emplace_back(0, .5f);
    decoder.decode(bad);
  });
  invalid([&] {
    auto bad = upper;
    bad[1].emplace_back(12, .5f);
    decoder.decode(bad);
  });
  invalid([&] {
    auto bad = upper;
    bad[1].emplace_back(1, .5f);
    decoder.decode(bad);
  });
  invalid([&] {
    auto bad = upper;
    bad[8].emplace_back(1, 2.f);
    decoder.decode(bad);
  });
  invalid([&] {
    auto bad = options;
    bad.thresholds.clear();
    Decoder{bad};
  });
  invalid([&] {
    auto bad = options;
    bad.weights = {1.f};
    Decoder{bad};
  });
  invalid([&] {
    auto bad = options;
    bad.thresholds[0] = -1;
    Decoder{bad};
  });
  invalid([&] {
    auto bad = options;
    bad.weights = {1.f, 0.f};
    Decoder{bad};
  });
  invalid([&] {
    auto bad = options;
    bad.backend = Backend::DD;
    bad.dd.max_iterations = 0;
    Decoder{bad};
  });
  if (!Decoder::ilp_available()) {
    options.backend = Backend::ILP;
    bool unavailable = false;
    try {
      Decoder(options).decode(11, pairs);
    } catch (const std::runtime_error &) {
      unavailable = true;
    }
    check(unavailable);
    unavailable = false;
    try {
      Decoder(options).decode(0, {});
    } catch (const std::runtime_error &) {
      unavailable = true;
    }
    check(unavailable);
  }

  // Exact level DP bounds the full model; finite-iteration DD may have a gap.
  options.no_lonely_pairs = false;
  options.weights = {.7f, .3f};
  options.dd.nussinov_dp = true;
  options.dd.crossing_beam = options.dd.witnesses = 0;
  std::mt19937 random(42);
  for (int trial = 0; trial < 25; ++trial) {
    std::vector<PairProbability> sample;
    for (int i = 0; i < 8; ++i)
      for (int j = i + 1; j < 8; ++j)
        if (random() % 4 == 0)
          sample.push_back({i, j, (random() % 10 + 1) / 10.f});
    const double best = exhaustive(8, sample, options);
    options.backend = Backend::DD;
    const auto dd = Decoder(options).decode(8, sample);
    check(dd.objective <= best + 1e-6);
    check(dd.dd->upper_bound >= best - 1e-6);
    if (Decoder::ilp_available()) {
      options.backend = Backend::ILP;
      check(std::abs(Decoder(options).decode(8, sample).objective - best) < 1e-6);
    }
  }
}
