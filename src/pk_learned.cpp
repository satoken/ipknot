#include "pk_learned.h"
#include "pk_energy.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace {
constexpr double specificity_epsilon = 1e-6;

double bounded_length(int value, double denominator) {
  return std::min(1.0, std::log1p(static_cast<double>(value)) / denominator);
}

} // namespace

double PKPosteriorContext::mean_probability(const PKLearnedStem& stem) const {
  const std::array<int, 3> key{stem.left, stem.right, stem.length};
  auto found = probability_stem_cache_.find(key);
  if (found != probability_stem_cache_.end()) return found->second;
  if (stem.length < 1 || stem.left < 0 || stem.right >= length() ||
      stem.length > (stem.right - stem.left + 1) / 2)
    throw std::invalid_argument("Invalid DP stem coordinates");
  double result = 0;
  for (int d = 0; d < stem.length; ++d)
    result += std::min(1.0, evidence(stem.left + d, stem.right - d));
  result /= stem.length;
  probability_stem_cache_.emplace(key, result);
  return result;
}

PKPosteriorContext::StemEvidence PKPosteriorContext::stem_evidence(
    const PKLearnedStem& stem) const {
  const std::array<int, 3> key{stem.left, stem.right, stem.length};
  const auto cached = stem_cache_.find(key);
  if (cached != stem_cache_.end()) return cached->second;
  if (stem.left < 0 || stem.right >= length() || stem.length < 1 ||
      stem.left >= stem.right ||
      stem.length > (stem.right - stem.left + 1) / 2)
    throw std::invalid_argument("Invalid learned PK stem coordinates or length");
  StemEvidence result;
  for (int d = 0; d < stem.length; ++d) {
    const double value = support(stem.left + d, stem.right - d);
    result.support += value;
    result.min_support = std::min(result.min_support, value);
    result.specificity += specificity(stem.left + d, stem.right - d);
  }
  result.support /= stem.length;
  result.specificity /= stem.length;
  stem_cache_.emplace(key, result);
  return result;
}

std::array<double, 2> PKPosteriorContext::bounded_stem_evidence(
    const PKLearnedStem& stem) const {
  const std::array<int, 3> key{stem.left, stem.right, stem.length};
  const auto cached = bounded_stem_cache_.find(key);
  if (cached != bounded_stem_cache_.end()) return cached->second;
  if (stem.left < 0 || stem.right >= length() || stem.length < 1 ||
      stem.left >= stem.right || stem.length > (stem.right - stem.left + 1) / 2)
    throw std::invalid_argument("Invalid bounded PK stem coordinates or length");
  std::array<double, 2> result{};
  for (int d = 0; d < stem.length; ++d) {
    const int left = stem.left + d, right = stem.right - d;
    const double own = support(left, right);
    const double alternative = std::max(strongest_alternative(left, right),
                                        strongest_alternative(right, left));
    result[0] += own;
    result[1] += own - alternative / (1 + alternative);
  }
  for (auto& value : result) value /= stem.length;
  bounded_stem_cache_.emplace(key, result);
  return result;
}

PKPosteriorContext::PKPosteriorContext(const PKPosteriorPairs& posterior)
    : length_(0) {
  if (posterior.empty() || posterior.size() - 1 >
      static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("Learned PK posterior requires one-based sparse rows");
  if (!posterior.front().empty())
    throw std::invalid_argument("Learned PK posterior row zero must be empty");
  length_ = static_cast<int>(posterior.size() - 1);
  top_two_.resize(length_);
  auto insert_alternative = [&](int base, int partner, double value) {
    const Alternative candidate{value, partner};
    auto better = [](const Alternative& a, const Alternative& b) {
      return a.evidence > b.evidence ||
          (a.evidence == b.evidence && a.partner < b.partner);
    };
    auto& top = top_two_[base];
    if (better(candidate, top[0])) {
      top[1] = top[0];
      top[0] = candidate;
    } else if (better(candidate, top[1])) {
      top[1] = candidate;
    }
  };
  for (int i = 1; i <= length_; ++i) {
    for (const auto& [j, value] : posterior[i]) {
      if (j == 0 || j == static_cast<unsigned int>(i) ||
          j > static_cast<unsigned int>(length_) ||
          !std::isfinite(value) || value < 0)
        throw std::invalid_argument("Invalid learned PK sparse posterior pair or evidence");
      // FASTA engines store symmetric entries. The ILP creates physical
      // variables from i<j only, so use that same authoritative evidence;
      // reverse values may differ for other engines and are never added.
      if (j < static_cast<unsigned int>(i)) continue;
      const Pair pair{i - 1, static_cast<int>(j) - 1};
      if (!evidence_.emplace(pair, value).second)
        throw std::invalid_argument("Duplicate learned PK sparse posterior pair");
      insert_alternative(pair.first, pair.second, value);
      insert_alternative(pair.second, pair.first, value);
    }
  }
}

bool PKPosteriorContext::contains(int left, int right) const {
  return evidence_.find({left, right}) != evidence_.end();
}

double PKPosteriorContext::evidence(int left, int right) const {
  const auto found = evidence_.find({left, right});
  if (found == evidence_.end())
    throw std::invalid_argument("Missing posterior evidence for learned PK pair");
  return found->second;
}

double PKPosteriorContext::support(int left, int right) const {
  const double value = evidence(left, right);
  // Sparse refinement adds conditional evidence to the existing matrix and
  // can exceed one. This transform is bounded without treating it as a
  // calibrated probability or taking logit(value).
  return value / (1 + value);
}

double PKPosteriorContext::strongest_alternative(int base, int partner) const {
  const auto& top = top_two_.at(base);
  if (top[0].partner != partner && top[0].partner >= 0) return top[0].evidence;
  return top[1].partner >= 0 ? top[1].evidence : 0;
}

double PKPosteriorContext::bounded_margin(int left, int right) const {
  const double p=evidence(left,right);
  const double q=std::max(strongest_alternative(left,right),strongest_alternative(right,left));
  return p/(1+p)-q/(1+q);
}

double PKPosteriorContext::specificity(int left, int right) const {
  const double value = evidence(left, right);
  const double alternative = std::max(strongest_alternative(left, right),
                                      strongest_alternative(right, left));
  // Zero here denotes absent sparse alternative evidence. It is a decoder
  // feature, not an assertion of exact zero physical probability.
  const double ratio = std::log((value + specificity_epsilon) /
                               (alternative + specificity_epsilon)) / 4;
  return std::clamp(ratio, -1.0, 1.0);
}

const std::array<const char*, PK_LEARNED_FEATURE_COUNT>& pk_learned_feature_names() {
  static const std::array<const char*, PK_LEARNED_FEATURE_COUNT> names{
      "bias", "anchor_support", "weak_support", "weak_min_support",
      "support_product", "support_gap", "weak_specificity",
      "anchor_specificity_product", "min_length", "weak_length",
      "outer_loop", "middle_loop"};
  return names;
}

PKLearnedFeatures pk_learned_features(const PKPosteriorContext& context,
                                     const PKLearnedStem& first,
                                     const PKLearnedStem& second,
                                     const std::array<int, 5>& geometry) {
  if (geometry[0] != first.length || geometry[1] != second.length ||
      geometry[2] < 0 || geometry[3] < 0 || geometry[4] < 0)
    throw std::invalid_argument("Invalid learned PK feature geometry");
  const auto a = context.stem_evidence(first);
  const auto b = context.stem_evidence(second);
  // Prefer longer support at an exact evidence tie, then smaller physical
  // coordinates. Thus swapping the input stems leaves the same anchor.
  const bool anchor_first = a.support > b.support ||
      (a.support == b.support && (first.length > second.length ||
       (first.length == second.length &&
        std::tie(first.left, first.right) <= std::tie(second.left, second.right))));
  const auto& anchor = anchor_first ? a : b;
  const auto& weak = anchor_first ? b : a;
  const auto& weak_stem = anchor_first ? second : first;
  const double length_scale = std::log(13.0);
  const double loop_scale = std::log(31.0);
  PKLearnedFeatures result;
  result.anchor_first = anchor_first;
  result.values = {
      1, anchor.support, weak.support, weak.min_support,
      anchor.support * weak.support, anchor.support - weak.support,
      weak.specificity, anchor.support * weak.specificity,
      bounded_length(std::min(first.length, second.length), length_scale),
      bounded_length(weak_stem.length, length_scale),
      std::min(1.0, (std::log1p(static_cast<double>(geometry[2])) +
                     std::log1p(static_cast<double>(geometry[4]))) / (2 * loop_scale)),
      bounded_length(geometry[3], loop_scale)};
  return result;
}

PKLearnedFeatures pk_bounded_features(const PKPosteriorContext& context,
                                     const PKLearnedStem& first,
                                     const PKLearnedStem& second,
                                     const std::array<int, 5>& geometry) {
  if (geometry[0] != first.length || geometry[1] != second.length ||
      geometry[2] < 0 || geometry[3] < 0 || geometry[4] < 0)
    throw std::invalid_argument("Invalid bounded PK feature geometry");
  const auto a = context.bounded_stem_evidence(first);
  const auto b = context.bounded_stem_evidence(second);
  PKLearnedFeatures result;
  result.values[0] = -1;
  result.values[1] = std::min(a[0], b[0]);
  result.values[2] = std::min(a[1], b[1]);
  result.values[3] = -(std::log1p(geometry[2]) + std::log1p(geometry[3]) +
                       std::log1p(geometry[4]));
  return result;
}

bool PKLearnedModel::active() const {
  if (exclusion()) return loaded_; // J=0 still restores the excluded 11 state
  if (dp_conversion()) return loaded_ && (weights[0] != 0 || weights[1] != 0);
  return loaded_ && (!bounded() || cap_ != 0) &&
      std::any_of(weights.begin(), weights.end(), [](double value) { return value != 0; });
}

PKLearnedFeatures PKLearnedModel::features(const PKPosteriorContext& context,
    const PKLearnedStem& first, const PKLearnedStem& second,
    const std::array<int, 5>& geometry) const {
  if (dp_conversion() || exclusion()) {
    if (geometry[0] != first.length || geometry[1] != second.length)
      throw std::invalid_argument("Invalid DP conversion geometry");
    PKLoopEnergy energy;
    energy.model = PKLoopEnergyModel::DP;
    energy.temperature_celsius = weights[2];
    PKLearnedFeatures result;
    result.values[0] = context.mean_probability(first);
    result.values[1] = context.mean_probability(second);
    result.values[2] = first.length;
    result.values[3] = second.length;
    result.values[4] = energy.evaluate(geometry).q;
    return result;
  }
  if (!bounded()) return pk_learned_features(context, first, second, geometry);
  auto result = pk_bounded_features(context, first, second, geometry);
  if (bounded_cc_cost_) {
    PKLoopEnergy energy;
    energy.model = PKLoopEnergyModel::CC;
    result.values[3] = -std::max(0.0, energy.evaluate(geometry).q);
  }
  return result;
}

void PKLearnedModel::load(const std::string& filename) {
  std::ifstream source(filename);
  if (!source) throw std::runtime_error("Cannot read learned PK model: " + filename);
  PKLearnedFeatureArray candidate{};
  std::array<bool, PK_LEARNED_FEATURE_COUNT> seen{};
  bool header = false;
  Kind kind = Kind::Linear;
  double cap = 0;
  bool seen_cap = false;
  bool seen_loop_model = false, cc_cost = false;
  std::string line;
  std::size_t line_number = 0;
  const auto& names = feature_names();
  while (std::getline(source, line)) {
    ++line_number;
    const auto comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    std::istringstream row(line);
    std::string name, extra;
    if (!(row >> name)) continue;
    if (!header) {
      if (name == "IPKNOT_PK_BOUNDED_V1") kind = Kind::Bounded;
      else if (name == "IPKNOT_PK_DP_LINEAR_V1") kind = Kind::DPLinear;
      else if (name == "IPKNOT_PK_DP_LOCAL_V1") kind = Kind::DPLocal;
      else if (name == "IPKNOT_PK_EXCLUSION_V1") kind = Kind::Exclusion;
      else if (name == "IPKNOT_PK_LINEAR_AB_V1") kind = Kind::LinearBlock;
      else if (name != "IPKNOT_PK_LINEAR_V1")
        throw std::invalid_argument("Invalid learned PK model header at line " +
                                    std::to_string(line_number));
      if (row >> extra)
        throw std::invalid_argument("Invalid learned PK model header at line " +
                                    std::to_string(line_number));
      header = true;
      continue;
    }
    if (kind == Kind::Bounded && name == "loop_model") {
      std::string value;
      if (seen_loop_model || !(row >> value) || (row >> extra) ||
          (value != "log" && value != "cc"))
        throw std::invalid_argument("Invalid or duplicate bounded PK loop model");
      cc_cost = value == "cc"; seen_loop_model = true; continue;
    }
    double value;
    if (!(row >> value) || !std::isfinite(value) || (row >> extra))
      throw std::invalid_argument("Invalid learned PK model row at line " +
                                  std::to_string(line_number));
    if (kind == Kind::DPLinear || kind == Kind::DPLocal || kind == Kind::Exclusion) {
      const std::array<const char*, 3> dp_names{"intercept", "energy_scale", "temperature"};
      const auto found = std::find_if(dp_names.begin(), dp_names.end(),
          [&](const char* item) { return name == item; });
      if (found == dp_names.end()) throw std::invalid_argument("Unknown DP conversion parameter");
      const auto index = std::distance(dp_names.begin(), found);
      if (seen[index] || (index == 1 && value < 0) || (index == 2 && value <= -273.15))
        throw std::invalid_argument("Invalid or duplicate DP conversion parameter");
      candidate[index] = value; seen[index] = true; continue;
    }
    if (kind == Kind::Bounded) {
      if (value < 0) throw std::invalid_argument("Bounded PK coefficients and cap must be nonnegative");
      if (name == "cap") {
        if (seen_cap) throw std::invalid_argument("Duplicate bounded PK cap");
        cap = value; seen_cap = true; continue;
      }
      static const std::array<const char*, 4> bounded_names{
          "bias", "support", "competition", "loop_cost"};
      const auto found = std::find_if(bounded_names.begin(), bounded_names.end(),
          [&](const char* item) { return name == item; });
      if (found == bounded_names.end()) throw std::invalid_argument("Unknown bounded PK coefficient");
      const std::size_t index = std::distance(bounded_names.begin(), found);
      if (seen[index]) throw std::invalid_argument("Duplicate bounded PK coefficient");
      candidate[index] = value; seen[index] = true; continue;
    }
    const auto found = std::find_if(names.begin(), names.end(),
                                  [&](const char* item) { return name == item; });
    if (found == names.end())
      throw std::invalid_argument("Unknown learned PK feature at line " +
                                  std::to_string(line_number));
    const std::size_t index = std::distance(names.begin(), found);
    if (seen[index])
      throw std::invalid_argument("Duplicate learned PK feature at line " +
                                  std::to_string(line_number));
    candidate[index] = value;
    seen[index] = true;
  }
  if (!source.eof()) throw std::runtime_error("Error reading learned PK model: " + filename);
  const auto required_end = kind == Kind::Bounded ? seen.begin()+4 :
      (kind == Kind::DPLinear || kind == Kind::DPLocal || kind == Kind::Exclusion) ? seen.begin()+3 : seen.end();
  if (!header || (kind == Kind::Bounded && !seen_cap) ||
      !std::all_of(seen.begin(), required_end, [](bool item) { return item; }))
    throw std::invalid_argument("Learned PK model is missing its header or required features");
  weights = candidate;
  kind_ = kind;
  cap_ = cap;
  bounded_cc_cost_ = cc_cost;
  loaded_ = true;
}

double PKLearnedModel::score(const PKLearnedFeatureArray& features) const {
  if (dp_conversion() || exclusion()) {
    for (int i = 0; i < 5; ++i)
      if (!std::isfinite(features[i])) throw std::invalid_argument("Nonfinite DP feature");
    const double a = features[0], b = features[1], A = features[2], B = features[3];
    if (a < 0 || a > 1 || b < 0 || b > 1 || A < 1 || B < 1 || features[4] < 0)
      throw std::invalid_argument("Invalid DP conversion feature");
    const double J = weights[0] - weights[1] * features[4];
    if (!std::isfinite(J)) throw std::invalid_argument("DP conversion overflow");
    if (exclusion()) {
      // The two observed marginals belong to a distribution excluding 11.
      // Means are whole-core proxies. Abstain if approximate BPP is outside
      // the three-state simplex or nearly saturates it; never divide by a
      // tiny residual mass. With valid mass, recover a latent independent
      // prior and optionally tilt 11. This is not a global PK posterior.
      const double empty = 1 - a - b;
      if (a == 0 || b == 0 || empty < .01) return 0;
      const double logit = std::log(a) + std::log(b) - std::log(empty) + J;
      const double joint = logit >= 0 ? 1 / (1 + std::exp(-logit)) :
          std::exp(logit) / (1 + std::exp(logit));
      return joint * (A * (1 - a) + B * (1 - b));
    }
    if (kind_ == Kind::DPLinear) return (A + B) * J;
    const double prior = a * b;
    if (J == 0 || prior == 0 || prior == 1) return 0;
    const double logit = std::log(prior) - std::log1p(-prior) + J;
    const double joint = logit >= 0 ? 1 / (1 + std::exp(-logit)) :
        std::exp(logit) / (1 + std::exp(logit));
    // Algebraically identical to summing the normalized four-state table.
    return (joint - prior) * (A * (1 - a) + B * (1 - b)) / (1 - prior);
  }
  double result = 0;
  for (std::size_t i = 0; i < weights.size(); ++i) {
    if (!std::isfinite(weights[i]) || !std::isfinite(features[i]))
      throw std::invalid_argument("Nonfinite learned PK weight or feature");
    result += weights[i] * features[i];
  }
  if (!std::isfinite(result)) throw std::invalid_argument("Learned PK score overflow");
  return bounded() ? std::clamp(result, -cap_, cap_) : result;
}
