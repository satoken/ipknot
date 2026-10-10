// SPDX-License-Identifier: GPL-3.0-or-later
#include "ipknot/decoder.h"
#include "decoder_internal.h"
#include "ilp_model.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace ipknot {
namespace {
void validate_pair(int length, int left, int right, float probability) {
  if (left < 0 || left >= right || right >= length || !std::isfinite(probability) ||
      probability < 0 || probability > 1)
    throw std::invalid_argument("Invalid base-pair probability or coordinates");
}

void append_pair(std::vector<DDPair> &pairs, int left, int right, float p,
                 const DecoderOptions &options) {
  for (int level = 0; level < static_cast<int>(options.thresholds.size()); ++level)
    if (p > options.thresholds[level])
      pairs.push_back(
          {left, right, level, (p - options.thresholds[level]) * options.weights[level]});
}

DecoderOptions normalize_options(DecoderOptions options) {
  if (options.backend != Backend::DD && options.backend != Backend::ILP)
    throw std::invalid_argument("Unknown decoder backend");
  if (options.thresholds.empty() ||
      options.thresholds.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("At least one decoding level is required");
  if (options.weights.empty())
    options.weights.assign(options.thresholds.size(), 1.0 / options.thresholds.size());
  if (options.weights.size() != options.thresholds.size() || options.threads < 0)
    throw std::invalid_argument("Invalid level weights or thread count");
  for (std::size_t level = 0; level < options.thresholds.size(); ++level) {
    const float t = options.thresholds[level], w = options.weights[level];
    if (!std::isfinite(t) || t < 0 || t > 1 || !std::isfinite(w) || w <= 0)
      throw std::invalid_argument("Thresholds must be in [0,1] and weights positive and finite");
  }
  if (options.backend == Backend::DD)
    options.dd.validate();
  return options;
}
} // namespace

Decoder::Decoder(DecoderOptions options) : options_(normalize_options(std::move(options))) {}

bool Decoder::ilp_available() { return IP::available(); }

namespace {
std::vector<DDPair> posterior_candidates(const SparsePosterior &posterior,
                                         const DecoderOptions &options, bool probabilities) {
  if (posterior.empty() || !posterior[0].empty() ||
      posterior.size() - 1 > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("Sparse posterior requires n+1 rows and an empty row 0");
  const int length = static_cast<int>(posterior.size() - 1);
  std::vector<DDPair> pairs;
  for (int i = 0; i < length; ++i) {
    std::unordered_set<unsigned int> seen;
    for (const auto &[j, p] : posterior[i + 1]) {
      if (j == 0 || j > static_cast<unsigned int>(length) ||
          j == static_cast<unsigned int>(i + 1) || !std::isfinite(p) || p < 0 ||
          (probabilities && p > 1))
        throw std::invalid_argument("Invalid sparse posterior entry");
      if (j <= static_cast<unsigned int>(i + 1))
        continue;
      if (!seen.insert(j).second)
        throw std::invalid_argument("Duplicate upper-triangular posterior entry");
      append_pair(pairs, i, j - 1, p, options);
    }
  }
  return pairs;
}

DecodeResult decode_candidates(int length, const std::vector<DDPair> &pairs,
                               const DecoderOptions &options) {
  if (options.backend == Backend::ILP) {
    if (!Decoder::ilp_available())
      throw std::runtime_error("No ILP solver is linked to ipknot_decoder");
    return detail::decode_ilp(length, pairs, options.thresholds.size(), options.no_lonely_pairs,
                              options.threads);
  }
  auto dd = solve_dual_decomposition(length, pairs, options.thresholds.size(),
                                     options.no_lonely_pairs, options.dd);
  DecodeResult result{dd.bpseq, dd.levels, dd.objective, std::move(dd)};
  return result;
}
} // namespace

DecodeResult Decoder::decode(int length, const std::vector<PairProbability> &probabilities) const {
  if (length < 0)
    throw std::invalid_argument("Negative sequence length");
  std::vector<DDPair> pairs;
  std::unordered_set<std::uint64_t> seen;
  for (const auto &pair : probabilities) {
    validate_pair(length, pair.left, pair.right, pair.probability);
    const auto key = (std::uint64_t(pair.left) << 32) | std::uint32_t(pair.right);
    if (!seen.insert(key).second)
      throw std::invalid_argument("Duplicate base-pair probability");
    append_pair(pairs, pair.left, pair.right, pair.probability, options_);
  }
  return decode_candidates(length, pairs, options_);
}

DecodeResult Decoder::decode(const SparsePosterior &posterior) const {
  const auto pairs = posterior_candidates(posterior, options_, true);
  return decode_candidates(static_cast<int>(posterior.size() - 1), pairs, options_);
}

DecodeResult detail::decode_posterior_scores(const SparsePosterior &posterior,
                                             DecoderOptions options) {
  options = normalize_options(std::move(options));
  const auto pairs = posterior_candidates(posterior, options, false);
  return decode_candidates(static_cast<int>(posterior.size() - 1), pairs, options);
}
} // namespace ipknot
