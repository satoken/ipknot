#ifndef IPKNOT_PK_ENSEMBLE_H
#define IPKNOT_PK_ENSEMBLE_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <vector>

// Sequence-independent DP loop component approximation. A connected component
// of crossing maximal stems is treated as one PK loop. This is exact for an
// isolated simple H core; complex components are NOT a full DP decomposition.
struct PKComponentEnergy { double kcal = 0; int components = 0; };
inline PKComponentEnergy pk_component_energy(const std::vector<int>& pairs) {
  const int n = pairs.size();
  struct Stem { int left, right, size; };
  std::vector<Stem> stems;
  for (int i = 0; i < n; ++i) {
    const int j = pairs[i];
    if (j < -1 || j >= n || j == i || (j >= 0 && pairs[j] != i))
      throw std::invalid_argument("Invalid PK ensemble structure");
    if (j <= i || (i > 0 && j + 1 < n && pairs[i-1] == j+1)) continue;
    int length = 1;
    while (i+length < j-length && pairs[i+length] == j-length) ++length;
    stems.push_back({i,j,length});
  }
  std::vector<int> parent(stems.size());
  for (std::size_t k = 0; k < parent.size(); ++k) parent[k] = k;
  auto root = [&](int k) { while (parent[k] != k) k = parent[k]; return k; };
  std::vector<bool> crossed(stems.size(), false);
  for (std::size_t a = 0; a < stems.size(); ++a)
    for (std::size_t b = a+1; b < stems.size(); ++b) {
      const auto& s = stems[a]; const auto& t = stems[b];
      if (s.left < t.left && t.left < s.right && s.right < t.right) {
        parent[root(b)] = root(a); crossed[a] = crossed[b] = true;
      }
    }
  std::map<int, std::vector<int>> groups;
  for (std::size_t k = 0; k < stems.size(); ++k) if (crossed[k]) groups[root(k)].push_back(k);
  PKComponentEnergy result;
  for (const auto& [_, group] : groups) {
    int left = n, right = -1;
    for (int k : group) { left = std::min(left, stems[k].left); right = std::max(right, stems[k].right); }
    int unpaired = 0, boundaries = group.size();
    bool enclosed = false;
    for (int i = left; i <= right; ++i) if (pairs[i] < 0) ++unpaired;
    for (std::size_t k = 0; k < stems.size(); ++k) {
      const auto& s = stems[k];
      if (s.left < left && s.right > right) enclosed = true;
      // Noncrossing child helices border this approximated PK loop.
      if (!crossed[k] && left < s.left && s.right < right) {
        ++boundaries;
        // Their enclosed unpaired bases belong to their own secondary loops.
        for (int i = s.left; i <= s.right; ++i) if (pairs[i] < 0) --unpaired;
        // Only the outermost such child is a boundary of this PK loop.
        for (std::size_t h = 0; h < k; ++h)
          if (!crossed[h] && left < stems[h].left && stems[h].left < s.left &&
              s.right < stems[h].right && stems[h].right < right) {
            --boundaries;
            for (int i = s.left; i <= s.right; ++i) if (pairs[i] < 0) ++unpaired;
            break;
          }
      }
    }
    result.kcal += (enclosed ? 15.0 : 9.6) + .1 * boundaries + .1 * unpaired;
    ++result.components;
  }
  return result;
}

struct PKEnsembleCandidate {
  std::vector<int> pairs, levels;
  double utility = 0, penalty = 0;
  PKComponentEnergy energy;
};
struct PKEnsembleResult {
  std::vector<PKEnsembleCandidate> candidates;
  std::vector<double> probabilities, expected_gain;
  std::map<std::pair<int,int>, double> marginals;
  std::size_t selected = 0, visits = 0;
};
inline PKEnsembleResult pk_finite_ensemble(const std::vector<PKEnsembleCandidate>& pool,
    double temperature, double energy_scale, double intercept, double rt, double threshold) {
  if (pool.empty() || !std::isfinite(temperature) || temperature <= 0 ||
      !std::isfinite(energy_scale) || energy_scale < 0 || !std::isfinite(intercept) ||
      !std::isfinite(rt) || rt <= 0 || !std::isfinite(threshold) || threshold < 0 || threshold > 1)
    throw std::invalid_argument("Invalid finite PK ensemble parameters");
  PKEnsembleResult result; result.visits = pool.size();
  std::map<std::vector<int>, PKEnsembleCandidate> unique;
  for (const auto& c : pool) {
    if (!std::isfinite(c.utility) || !std::isfinite(c.penalty) || c.penalty < 0)
      throw std::invalid_argument("Invalid PK ensemble candidate utility");
    auto [it, inserted] = unique.emplace(c.pairs, c);
    if (!inserted && c.utility - c.penalty > it->second.utility - it->second.penalty) it->second = c;
  }
  double max_log = -INFINITY;
  for (const auto& [_, c] : unique) {
    auto copy = c; copy.energy = pk_component_energy(c.pairs);
    if (copy.levels.size() != copy.pairs.size()) throw std::invalid_argument("Invalid PK ensemble levels");
    result.candidates.push_back(copy);
    const double log_weight = (copy.utility - copy.penalty) / temperature +
        intercept * copy.energy.components - energy_scale * copy.energy.kcal / rt;
    if (!std::isfinite(log_weight)) throw std::invalid_argument("PK ensemble weight overflow");
    result.probabilities.push_back(log_weight); max_log = std::max(max_log, log_weight);
  }
  double norm = 0;
  for (auto& w : result.probabilities) { w = std::exp(w - max_log); norm += w; }
  for (std::size_t k = 0; k < result.candidates.size(); ++k) {
    const auto& c = result.candidates[k]; result.probabilities[k] /= norm;
    for (std::size_t i = 0; i < c.pairs.size(); ++i) if (c.pairs[i] > static_cast<int>(i))
      result.marginals[{i, c.pairs[i]}] += result.probabilities[k];
  }
  double best = -INFINITY;
  for (std::size_t k = 0; k < result.candidates.size(); ++k) {
    const auto& c = result.candidates[k]; double gain = -c.penalty;
    for (std::size_t i = 0; i < c.pairs.size(); ++i) if (c.pairs[i] > static_cast<int>(i))
      gain += result.marginals.at({i,c.pairs[i]}) - threshold;
    result.expected_gain.push_back(gain);
    if (gain > best) { best = gain; result.selected = k; }
  }
  return result;
}

// Optional compact diagnostics outside timed runs. Contains no sequence data.
inline void write_pk_ensemble(const std::string& filename, const PKEnsembleResult& result) {
  std::ofstream out(filename, std::ios::app);
  if (!out) throw std::runtime_error("Cannot write PK ensemble diagnostics");
  out << std::setprecision(17) << "{\"visits\":" << result.visits << ",\"selected\":" << result.selected << ",\"candidates\":[";
  for (std::size_t k = 0; k < result.candidates.size(); ++k) {
    if (k) out << ',';
    const auto& c = result.candidates[k];
    out << "{\"utility\":" << c.utility << ",\"penalty\":" << c.penalty << ",\"energy\":" << c.energy.kcal
        << ",\"components\":" << c.energy.components << ",\"probability\":" << result.probabilities[k]
        << ",\"gain\":" << result.expected_gain[k] << ",\"pairs\":[";
    bool first = true;
    for (std::size_t i = 0; i < c.pairs.size(); ++i) if (c.pairs[i] > static_cast<int>(i)) {
      if (!first) out << ','; first = false; out << '[' << i << ',' << c.pairs[i] << ']';
    }
    out << "]}";
  }
  out << "]}\n";
  if (!out) throw std::runtime_error("Failed writing PK ensemble diagnostics");
}
#endif
