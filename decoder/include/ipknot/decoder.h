// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef IPKNOT_DECODER_H
#define IPKNOT_DECODER_H

#include "dual_decomposition.h"
#include <optional>

namespace ipknot {

enum class Backend { DD, ILP };

// Each physical pair occurs once, with zero-based left < right.
struct PairProbability {
  int left, right;
  float probability;
};

// Native IPknot sparse format: n+1 rows, unused row 0, one-based partners.
// Upper-triangular and symmetric matrices are accepted; lower entries are ignored.
using SparsePosterior = std::vector<std::vector<std::pair<unsigned int, float>>>;

struct DecoderOptions {
  Backend backend = Backend::DD;
  std::vector<float> thresholds{0.25f, 0.125f}; // one threshold per level
  std::vector<float> weights;                   // empty: 1 / number of levels
  bool no_lonely_pairs = true;
  int threads = 1; // ILP only; 0: backend default
  DDOptions dd;    // default: bounded beam DD; nussinov_dp=true: exact level DP
};

struct DecodeResult {
  // Zero-based partners and levels; -1 for unpaired bases, on both endpoints.
  std::vector<int> bpseq, levels;
  double objective = 0;
  // Present only for DD. Its upper_bound refers to the retained witness graph.
  std::optional<DDResult> dd;
};

// Decodes supplied probabilities without a sequence or a probability engine.
// Invalid inputs throw std::invalid_argument; unavailable ILP throws std::runtime_error.
class Decoder {
public:
  explicit Decoder(DecoderOptions options = {});
  DecodeResult decode(int length, const std::vector<PairProbability> &pairs) const;
  DecodeResult decode(const SparsePosterior &posterior) const;
  static bool ilp_available();

private:
  DecoderOptions options_;
};

} // namespace ipknot
#endif
