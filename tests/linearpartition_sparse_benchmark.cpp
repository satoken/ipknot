// Measure the sparse fixed-beam posterior path independently of exact MIP.
#include "linearpartition/LinearPartition.h"
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  const int length = argc > 1 ? std::atoi(argv[1]) : 1000;
  const bool vienna = argc > 2 && std::string(argv[2]) == "lpv";
  if (length < 5) return 2;
  std::mt19937 random(19);
  std::string sequence;
  sequence.reserve(length);
  for (int i = 0; i < length; ++i) sequence.push_back("ACGU"[random() % 4]);
  LinearPartition::BeamCKYParser parser(100, true, false, false, 0.001, false);
  const auto start = std::chrono::steady_clock::now();
  if (vienna) parser.parse<true, int>(sequence);
  else parser.parse<false, float>(sequence);
  std::vector<std::vector<std::pair<unsigned int, float>>> posterior;
  parser.get_posterior(posterior);
  size_t entries = 0;
  for (const auto& row : posterior) entries += row.size();
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  std::cout << "engine=" << (vienna ? "lpv" : "lpc") << " length=" << length
            << " beam=100 entries=" << entries << " seconds=" << seconds << '\n';
}
