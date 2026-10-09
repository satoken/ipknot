#include "pk_learned.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
void close(double actual, double expected, const char* message) {
  require(std::isfinite(actual) && std::abs(actual - expected) < 1e-10, message);
}
template<class F> void invalid(F function, const char* message) {
  try { function(); } catch (const std::invalid_argument&) { return; }
  throw std::runtime_error(message);
}
void same(const PKLearnedFeatureArray& a, const PKLearnedFeatureArray& b,
          const char* message) {
  for (std::size_t i = 0; i < a.size(); ++i) close(a[i], b[i], message);
}
void add(PKPosteriorPairs& p, int left, int right, float value) {
  p[left + 1].emplace_back(right + 1, value);
}
struct TemporaryFile {
  std::filesystem::path path;
  ~TemporaryFile() { std::error_code error; std::filesystem::remove(path, error); }
  void write(const std::string& text) const {
    std::ofstream output(path);
    output << text;
    require(static_cast<bool>(output), "Cannot write learned model fixture");
  }
};
std::string model_text(std::size_t count = PK_LEARNED_FEATURE_COUNT,
                       bool reverse = false) {
  std::ostringstream text;
  text << "# portable test coefficients\n\nIPKNOT_PK_LINEAR_V1 # version\n";
  const auto& names = pk_learned_feature_names();
  for (std::size_t item = 0; item < count; ++item) {
    const auto i = reverse ? count - item - 1 : item;
    text << names[i] << ' ' << (static_cast<int>(i) - 5) * .125 << " # coefficient\n";
  }
  return text.str();
}
} // namespace

int main() {
  try {
    PKPosteriorPairs posterior(21);
    add(posterior, 0, 12, 2); // accumulated refined evidence greater than one
    add(posterior, 1, 11, 1);
    add(posterior, 2, 10, .5);
    add(posterior, 5, 18, .25);
    add(posterior, 6, 17, .125);
    add(posterior, 0, 14, .5);
    add(posterior, 8, 12, 4); // stronger competitor at the second endpoint
    add(posterior, 5, 15, .5);
    PKPosteriorContext context(posterior);
    require(context.length() == 20, "One-based posterior length conversion failed");
    require(context.contains(0, 12) && !context.contains(1, 12),
            "Physical posterior coordinates were shifted incorrectly");
    close(context.evidence(0, 12), 2, "Evidence above one was lost");
    close(context.support(0, 12), 2.0 / 3, "Accumulated evidence was not bounded");
    close(context.specificity(0, 12), std::log((2 + 1e-6) / (4 + 1e-6)) / 4,
          "Specificity missed the stronger alternative at one endpoint");
    close(context.specificity(1, 11), 1,
          "Target pair was not excluded from its own competitors");
    close(context.specificity(5, 18), std::log((.25 + 1e-6) / (.5 + 1e-6)) / 4,
          "Specificity used the wrong endpoint competition");
    invalid([&] { context.evidence(3, 9); }, "Missing posterior evidence was accepted");
    invalid([&] { context.specificity(-1, 9); }, "Invalid posterior coordinates were accepted");

    const PKLearnedStem a{0, 12, 3}, b{5, 18, 2};
    const std::array<int, 5> geometry{3, 2, 2, 3, 4};
    const auto features = pk_learned_features(context, a, b, geometry);
    auto symmetric = posterior;
    for (std::size_t i = 1; i < posterior.size(); ++i)
      for (const auto& [j, value] : posterior[i])
        symmetric[j].emplace_back(static_cast<unsigned int>(i), value);
    same(features.values, pk_learned_features(PKPosteriorContext(symmetric), a, b, geometry).values,
         "Symmetric posterior evidence was duplicated or shifted");
    symmetric[13].front().second = .01f;
    same(features.values, pk_learned_features(PKPosteriorContext(symmetric), a, b, geometry).values,
         "Reverse evidence overrode the decoder's upper-triangle value");
    require(features.anchor_first, "Higher supported stem was not selected as anchor");
    const double weak_support = (.2 + 1.0 / 9) / 2;
    const double weak_specificity = (context.specificity(5, 18) + 1) / 2;
    close(features.values[0], 1, "Bias feature is wrong");
    close(features.values[1], .5, "Anchor mean support is wrong");
    close(features.values[2], weak_support, "Weak mean support is wrong");
    close(features.values[3], 1.0 / 9, "Weak minimum support is wrong");
    close(features.values[4], .5 * weak_support, "Support interaction is wrong");
    close(features.values[5], .5 - weak_support, "Support gap is wrong");
    close(features.values[6], weak_specificity, "Weak mean specificity is wrong");
    close(features.values[7], .5 * weak_specificity, "Anchor/specificity interaction is wrong");
    close(features.values[8], std::log(3.0) / std::log(13.0), "Minimum length feature is wrong");
    close(features.values[9], std::log(3.0) / std::log(13.0), "Weak length feature is wrong");
    close(features.values[10], (std::log(3.0) + std::log(5.0)) / (2 * std::log(31.0)),
          "Outer-loop feature is wrong");
    close(features.values[11], std::log(4.0) / std::log(31.0), "Middle-loop feature is wrong");
    for (double value : features.values)
      require(std::isfinite(value) && value >= -1 && value <= 1, "Feature is not finite/bounded");
    const auto swapped = pk_learned_features(context, b, a, {2, 3, 4, 3, 2});
    require(!swapped.anchor_first, "Swapping stems changed the physical anchor");
    same(features.values, swapped.values, "Swapping stems changed feature values");
    same(features.values, pk_learned_features(context, a, b, geometry).values,
         "Cached stem evidence changed feature values");
    const auto bounded = pk_bounded_features(context, a, b, geometry);
    close(bounded.values[0], -1, "Bounded intercept feature is wrong");
    close(bounded.values[1], weak_support, "Bounded joint support is wrong");
    // A margins: 2/3-4/5, 1/2, 1/3; B margins: 1/5-1/3, 1/9.
    close(bounded.values[2], ((.2-1.0/3)+1.0/9)/2,
          "Bounded margin used logarithmic specificity or included its own pair");
    close(bounded.values[3], -(std::log(3.0)+std::log(4.0)+std::log(5.0)),
          "Bounded loop cost is wrong");
    same(bounded.values,pk_bounded_features(context,b,a,{2,3,4,3,2}).values,
         "Bounded score is not symmetric between physical stems");
    same(bounded.values,pk_bounded_features(context,a,b,geometry).values,
         "Bounded stem cache changed evidence");
    invalid([&] { pk_learned_features(context, {0, 12, 4}, b, {4, 2, 2, 3, 4}); },
            "Missing required stem pair did not throw");
    invalid([&] { pk_learned_features(context, {-1, 12, 3}, b, geometry); },
            "Negative stem position did not throw");
    invalid([&] { pk_learned_features(context, a, b, {2, 2, 2, 3, 4}); },
            "Geometry/stem length mismatch did not throw");
    invalid([&] { pk_learned_features(context, a, b, {3, 2, -1, 3, 4}); },
            "Negative loop length did not throw");

    PKPosteriorPairs tied(21);
    for (int d = 0; d < 3; ++d) add(tied, a.left + d, a.right - d, .5);
    for (int d = 0; d < 2; ++d) add(tied, b.left + d, b.right - d, .5);
    PKPosteriorContext tied_context(tied);
    require(pk_learned_features(tied_context, a, b, geometry).anchor_first,
            "Support tie did not prefer the longer stem");
    require(!pk_learned_features(tied_context, b, a, {2, 3, 4, 3, 2}).anchor_first,
            "Support/length tie handling is not deterministic");
    add(tied, 7, 16, .5);
    PKPosteriorContext coordinate_tie(tied);
    const PKLearnedStem same_length_b{5, 18, 3};
    require(pk_learned_features(coordinate_tie, a, same_length_b, {3, 3, 2, 2, 3}).anchor_first,
            "Coordinate tie did not prefer the smaller physical stem");
    require(!pk_learned_features(coordinate_tie, same_length_b, a, {3, 3, 3, 2, 2}).anchor_first,
            "Coordinate tie changed after swapping stems");

    PKPosteriorPairs long_posterior(91);
    for (int d = 0; d < 13; ++d) {
      add(long_posterior, d, 49 - d, 1e12f);
      add(long_posterior, 20 + d, 89 - d, 1e12f);
    }
    const auto large = pk_learned_features(PKPosteriorContext(long_posterior),
        {0, 49, 13}, {20, 89, 13}, {13, 13, 2000000000, 2000000000, 2000000000});
    for (double value : large.values)
      require(std::isfinite(value) && value >= -1 && value <= 1, "Large input feature overflowed");
    for (int i : {8, 9, 10, 11}) close(large.values[i], 1, "Long geometry was not clipped");

    invalid([] { PKPosteriorContext empty(PKPosteriorPairs{}); }, "Empty sparse layout was accepted");
    for (auto bad_pair : {std::pair<unsigned int, float>{0, .5f}, {1, .5f}, {4, .5f},
                          {3, -.5f}, {3, std::numeric_limits<float>::infinity()},
                          {3, std::numeric_limits<float>::quiet_NaN()}}) {
      PKPosteriorPairs bad(4);
      bad[1].push_back(bad_pair);
      invalid([&] { PKPosteriorContext rejected(bad); }, "Malformed sparse evidence was accepted");
    }
    PKPosteriorPairs duplicate(4);
    duplicate[1] = {{3, .5f}, {3, .25f}};
    invalid([&] { PKPosteriorContext rejected(duplicate); }, "Duplicate evidence was accepted");
    duplicate[0] = {{3, .5f}};
    invalid([&] { PKPosteriorContext rejected(duplicate); }, "Nonempty zero row was accepted");

    PKLearnedModel model;
    require(!model.loaded(), "Default model claims to be loaded");
    close(model.score(features), 0, "Default model is not zero");
    TemporaryFile file{std::filesystem::temp_directory_path() /
        ("ipknot-learned-test-" + std::to_string(reinterpret_cast<std::size_t>(&model)) + ".txt")};
    file.write(model_text(PK_LEARNED_FEATURE_COUNT, true));
    model.load(file.path.string());
    require(model.loaded(), "Valid model did not load");
    double expected = 0;
    for (std::size_t i = 0; i < model.weights.size(); ++i) {
      close(model.weights[i], (static_cast<int>(i) - 5) * .125, "Named model coefficient is wrong");
      expected += model.weights[i] * features.values[i];
    }
    close(model.score(features), expected, "Learned linear dot product is wrong");
    auto doubled = model;
    for (auto& value : doubled.weights) value *= 2;
    close(doubled.score(features), 2 * expected, "Model score is not linear");
    const auto old_weights = model.weights;
    for (const std::string& text : {
        std::string("# empty\n"), std::string("WRONG\n"),
        std::string("IPKNOT_PK_LINEAR_V1 trailing\n"), model_text(11),
        model_text() + "bias 0\n", model_text() + "unknown 0\n",
        std::string("IPKNOT_PK_LINEAR_V1\nbias nan\n"),
        std::string("IPKNOT_PK_LINEAR_V1\nbias inf\n"),
        std::string("IPKNOT_PK_LINEAR_V1\nbias 1 extra\n"),
        std::string("IPKNOT_PK_LINEAR_V1\nbias\n")}) {
      file.write(text);
      invalid([&] { model.load(file.path.string()); }, "Invalid model text was accepted");
      require(model.loaded(), "Failed reload cleared the loaded model");
    same(model.weights, old_weights, "Failed model reload corrupted old coefficients");
    }
    file.write("IPKNOT_PK_BOUNDED_V1\nbias .01\nsupport 1\ncompetition .5\nloop_cost .001\ncap .025\n");
    PKLearnedModel capped;
    capped.load(file.path.string());
    require(capped.bounded() && capped.block_normalized() && capped.active(),
            "Bounded model did not select whole-block allocation");
    close(capped.score(bounded),.025,"Bounded positive score exceeded its cap");
    auto negative_bounded=bounded.values;
    negative_bounded[1]=0;negative_bounded[2]=-1;
    close(capped.score(negative_bounded),-.025,"Bounded negative score was clipped to zero");
    for (const std::string& text : {
        std::string("IPKNOT_PK_BOUNDED_V1\nbias -1\n"),
        std::string("IPKNOT_PK_BOUNDED_V1\nbias 0\nsupport 1\ncompetition 0\nloop_cost 0\n"),
        std::string("IPKNOT_PK_BOUNDED_V1\nbias 0\nsupport 1\ncompetition 0\nloop_cost 0\ncap -1\n"),
        std::string("IPKNOT_PK_BOUNDED_V1\nbias 0\nsupport 1\ncompetition 0\nloop_cost 0\ncap 1\ncap 1\n")}) {
      file.write(text);
      invalid([&] { capped.load(file.path.string()); },"Invalid bounded model was accepted");
      close(capped.score(bounded),.025,"Failed bounded reload changed score or cap");
    }
    file.write("IPKNOT_PK_BOUNDED_V1\nbias 0\nsupport 1\ncompetition 0\nloop_cost 0\ncap 0\n");
    capped.load(file.path.string());
    require(!capped.active(),"Zero cap did not disable bounded scoring");
    close(capped.score(bounded),0,"Zero cap changed objective");
    auto invalid_features = features.values;
    invalid_features[0] = std::numeric_limits<double>::quiet_NaN();
    invalid([&] { model.score(invalid_features); }, "Nonfinite feature was accepted");
    doubled.weights[0] = std::numeric_limits<double>::infinity();
    invalid([&] { doubled.score(features); }, "Nonfinite coefficient was accepted");
    doubled.weights.fill(std::numeric_limits<double>::max());
    invalid([&] { doubled.score(features); }, "Overflowing dot product was accepted");
    std::cout << "Verified posterior coordinates, accumulated evidence, competition, anchor symmetry, bounded features and atomic learned model import\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
