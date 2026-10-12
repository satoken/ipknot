// Compile with IPKNOT_INTRA_THREADS against the new library; without it,
// compile against 1b12107 for a baseline with the identical fixture generator.
#include "dual_decomposition.h"
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <string>

int main(int argc, char** argv) {
  try {
    if (argc != 9) throw std::runtime_error("engine(sparse|wavefront|beam|improved) shape n degree beam repeats threads stacked");
    const std::string engine = argv[1], shape = argv[2];
    const int n = std::stoi(argv[3]), degree = std::stoi(argv[4]), beam = std::stoi(argv[5]);
    const int repeats = std::stoi(argv[6]), threads = std::stoi(argv[7]);
    const bool stacked = std::stoi(argv[8]);
    if (n < 2 || degree < 1 || repeats < 1 || threads < 0 || beam < 0 ||
        (shape != "sparse" && shape != "dense") ||
        (engine != "sparse" && engine != "wavefront" && engine != "beam" && engine != "improved"))
      throw std::runtime_error("Invalid argument");
    std::mt19937 rng(38433200);
    std::vector<DDPair> pairs;
    for (int j = 1; j < n; ++j) {
      std::set<int> starts;
      if (shape == "dense") for (int i = 0; i < j; ++i) starts.insert(i);
      else while (int(starts.size()) < std::min(degree, j)) starts.insert(rng() % j);
      for (int i : starts) pairs.push_back({i, j, 0, (int(rng() % 101) - 20) / 100.});
    }
    const int width = engine == "sparse" || engine == "wavefront" ? 0 : beam;
#ifdef IPKNOT_INTRA_THREADS
    DDNussinov dp(n, pairs, 0, width, stacked, engine == "improved", threads, engine == "wavefront");
#else
    if (threads != 1 || engine == "wavefront") throw std::runtime_error("Unsupported baseline mode");
    DDNussinov dp(n, pairs, 0, width, stacked, engine == "improved");
#endif
    std::vector<double> weights;
    for (const auto& pair : pairs) weights.push_back(pair.weight);
    std::vector<unsigned char> allowed(pairs.size(), 1);
    std::vector<int> selected;
    std::vector<double> times;
    double score = 0;
    for (int repeat = -1; repeat < repeats; ++repeat) {
      const auto begin = std::chrono::steady_clock::now();
      score = dp.decode(weights, allowed, selected);
      const auto end = std::chrono::steady_clock::now();
      if (repeat >= 0) times.push_back(std::chrono::duration<double>(end - begin).count());
    }
    std::sort(times.begin(), times.end());
    std::sort(selected.begin(), selected.end());
    std::uint64_t hash = 14695981039346656037ULL;
    for (int id : selected) { hash ^= id; hash *= 1099511628211ULL; }
    std::cout << std::setprecision(17) << "{\"engine\":\"" << engine << "\",\"shape\":\"" << shape
        << "\",\"n\":" << n << ",\"pairs\":" << pairs.size() << ",\"beam\":" << width
        << ",\"stacked\":" << stacked << ",\"threads\":" << threads
        << ",\"median_seconds\":" << (times[(times.size()-1)/2] + times[times.size()/2]) / 2
        << ",\"objective\":" << score << ",\"trace_hash\":\"" << hash
        << "\",\"pruned\":" << dp.pruned_states();
#ifdef IPKNOT_INTRA_THREADS
    std::cout << ",\"workers\":" << dp.used_workers() << ",\"parallel_steps\":" << dp.parallel_steps();
#endif
    std::cout << "}\n";
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
