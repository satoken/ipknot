#ifndef IPKNOT_PK_LEARNED_H
#define IPKNOT_PK_LEARNED_H

#include <array>
#include <cstddef>
#include <map>
#include <unordered_map>
#include <string>
#include <utility>
#include <vector>

// Input coordinates are one-based, as in IPknot's sparse VSVF. The upper
// triangle is authoritative, matching the decoder; validated reverse entries
// are ignored rather than summed. Context/stem coordinates are zero-based.
using PKPosteriorPairs = std::vector<std::vector<std::pair<unsigned int, float>>>;
constexpr std::size_t PK_LEARNED_FEATURE_COUNT = 12;
using PKLearnedFeatureArray = std::array<double, PK_LEARNED_FEATURE_COUNT>;
struct PKLearnedStem;
struct PKLearnedFeatures;

class PKPosteriorContext {
public:
  explicit PKPosteriorContext(const PKPosteriorPairs& posterior);
  int length() const { return length_; }
  bool contains(int left, int right) const;
  double evidence(int left, int right) const;
  double support(int left, int right) const;
  double specificity(int left, int right) const;
  double bounded_margin(int left, int right) const;
  // Proxy for a whole-stem event, not its thermodynamic joint probability.
  double mean_probability(const PKLearnedStem& stem) const;

private:
  using Pair = std::pair<int, int>;
  struct Alternative {
    double evidence = -1;
    int partner = -1;
  };
  int length_;
  struct PairHash {
    std::size_t operator()(const Pair& p) const {
      return std::hash<unsigned long long>()((static_cast<unsigned long long>(static_cast<unsigned>(p.first)) << 32) | static_cast<unsigned>(p.second));
    }
  };
  struct StemHash {
    std::size_t operator()(const std::array<int, 3>& s) const {
      return PairHash()({s[0], s[1]}) ^ (std::hash<int>()(s[2]) * 0x9e3779b9U);
    }
  };
  std::unordered_map<Pair, double, PairHash> evidence_;
  std::vector<std::array<Alternative, 2>> top_two_;
  struct StemEvidence {
    double support = 0;
    double min_support = 1;
    double specificity = 0;
  };
  // A context belongs to one BPP matrix and its sequential threshold search.
  // This logical-const local cache computes each physical stem summary once.
  mutable std::unordered_map<std::array<int, 3>, StemEvidence, StemHash> stem_cache_;
  mutable std::unordered_map<std::array<int, 3>, std::array<double, 2>, StemHash> bounded_stem_cache_;
  mutable std::unordered_map<std::array<int, 3>, double, StemHash> probability_stem_cache_;
  StemEvidence stem_evidence(const PKLearnedStem& stem) const;
  std::array<double, 2> bounded_stem_evidence(const PKLearnedStem& stem) const;
  double strongest_alternative(int base, int partner) const;
  friend PKLearnedFeatures pk_learned_features(const PKPosteriorContext&,
      const PKLearnedStem&, const PKLearnedStem&, const std::array<int, 5>&);
  friend PKLearnedFeatures pk_bounded_features(const PKPosteriorContext&,
      const PKLearnedStem&, const PKLearnedStem&, const std::array<int, 5>&);
};

struct PKLearnedStem {
  int left;
  int right;
  int length;
};

struct PKLearnedFeatures {
  PKLearnedFeatureArray values{};
  bool anchor_first = true;
};

const std::array<const char*, PK_LEARNED_FEATURE_COUNT>& pk_learned_feature_names();

// The physical anchor is selected from bounded mean evidence before
// optimization, independently of the level assigned to either stem.
PKLearnedFeatures pk_learned_features(const PKPosteriorContext& context,
                                     const PKLearnedStem& first,
                                     const PKLearnedStem& second,
                                     const std::array<int, 5>& geometry);
// Four signed features: -1, min mean bounded support, min mean bounded
// competitor margin, -sum(log1p(loop length)). No sequence input.
PKLearnedFeatures pk_bounded_features(const PKPosteriorContext& context,
                                     const PKLearnedStem& first,
                                     const PKLearnedStem& second,
                                     const std::array<int, 5>& geometry);

class PKLearnedModel {
public:
  PKLearnedFeatureArray weights{};
  bool loaded() const { return loaded_; }
  bool bounded() const { return kind_ == Kind::Bounded; }
  bool dp_conversion() const { return kind_ == Kind::DPLinear || kind_ == Kind::DPLocal; }
  bool exclusion() const { return kind_ == Kind::Exclusion; }
  bool block_normalized() const { return kind_ != Kind::Linear; }
  bool active() const;
  PKLearnedFeatures features(const PKPosteriorContext& context,
      const PKLearnedStem& first, const PKLearnedStem& second,
      const std::array<int, 5>& geometry) const;
  void load(const std::string& filename);
  double score(const PKLearnedFeatureArray& features) const;
  double score(const PKLearnedFeatures& features) const {
    return score(features.values);
  }
  static const std::array<const char*, PK_LEARNED_FEATURE_COUNT>& feature_names() {
    return pk_learned_feature_names();
  }

private:
  enum class Kind { Linear, LinearBlock, Bounded, DPLinear, DPLocal, Exclusion };
  Kind kind_ = Kind::Linear;
  double cap_ = 0;
  bool bounded_cc_cost_ = false;
  bool loaded_ = false;
};

#endif
