// Build against either revision's decoder library. Define
// IPKNOT_PERFORMANCE_THREADS for revisions providing DDOptions::threads.
#include "dual_decomposition.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <string>

int main(int argc, char** argv) {
  try {
    if (argc != 9) throw std::runtime_error("mode(dp|dd) shape(sparse|dense) length degree beam levels repeats threads");
    const std::string mode = argv[1], shape = argv[2];
    const int n = std::stoi(argv[3]), degree = std::stoi(argv[4]), beam = std::stoi(argv[5]);
    const int levels = std::stoi(argv[6]), repeats = std::stoi(argv[7]), threads = std::stoi(argv[8]);
    if ((mode != "dp" && mode != "dd") || (shape != "sparse" && shape != "dense") ||
        n < 2 || degree < 1 || beam < 0 || levels < 1 || repeats < 1 || threads < 0)
      throw std::runtime_error("Invalid benchmark argument");
    std::mt19937 rng(38433200);
    std::vector<DDPair> pairs;
    for (int j = 1; j < n; ++j) {
      std::set<int> starts;
      if (shape == "dense") for (int i = 0; i < j; ++i) starts.insert(i);
      else while (int(starts.size()) < std::min(degree, j)) starts.insert(rng() % j);
      for (int i : starts) for (int level = 0; level < (mode == "dp" ? 1 : levels); ++level)
        pairs.push_back({i, j, level, (int(rng() % 101) - 20) / 100.0});
    }
    std::vector<double> weights;
    for (const auto& pair : pairs) weights.push_back(pair.weight);
    std::vector<unsigned char> allowed(pairs.size(), 1);
    std::vector<int> selected;
    DDNussinov dp(n, pairs, 0, beam, false);
    DDOptions options;
    options.beam = beam;
    options.max_iterations = 20;
#ifdef IPKNOT_PERFORMANCE_THREADS
    options.threads = threads;
#else
    if (threads != 1) throw std::runtime_error("Baseline supports threads=1 only");
#endif
    std::vector<double> times;
    double objective = 0, upper = 0;
    int iterations = 0;
    for (int repeat = -1; repeat < repeats; ++repeat) {
      const auto start = std::chrono::steady_clock::now();
      if (mode == "dp") objective = upper = dp.decode(weights, allowed, selected);
      else {
        const auto result = solve_dual_decomposition(n, pairs, levels, false, options);
        objective = result.objective; upper = result.upper_bound; iterations = result.iterations;
      }
      const auto end = std::chrono::steady_clock::now();
      if (repeat >= 0) times.push_back(std::chrono::duration<double>(end - start).count());
    }
    std::sort(times.begin(), times.end());
    const double median = (times[(times.size()-1)/2] + times[times.size()/2]) / 2;
    std::cout << std::setprecision(17) << "{\"mode\":\"" << mode << "\",\"shape\":\"" << shape
              << "\",\"length\":" << n << ",\"pairs\":" << pairs.size()
              << ",\"beam\":" << beam << ",\"levels\":" << levels << ",\"threads\":" << threads
              << ",\"median_seconds\":" << median << ",\"objective\":" << objective
              << ",\"upper_bound\":" << upper << ",\"iterations\":" << iterations << "}\n";
  } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
