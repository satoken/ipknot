#include "pk_energy.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
void close(double actual, double expected, const char* message) {
  require(std::abs(actual - expected) < 1e-10, message);
}
template<class F> void invalid(F function, const char* message) {
  try { function(); } catch (const std::invalid_argument&) { return; }
  throw std::runtime_error(message);
}
struct TemporaryFile {
  std::filesystem::path path;
  ~TemporaryFile() { std::filesystem::remove(path); }
  void write(const std::string& text) const {
    std::ofstream output(path);
    output << text;
    require(static_cast<bool>(output), "Cannot write energy table fixture");
  }
};
} // namespace

int main() {
  try {
    PKLoopEnergy energy;
    close(energy.evaluate({3, 5, 1, 0, 2}).q, 0, "Disabled energy is nonzero");
    require(energy.evaluate({3, 5, 1, 0, 2}).source == PKLoopEnergySource::None,
            "Disabled energy has a source");
    energy.model = PKLoopEnergyModel::DP;
    close(energy.rt_kcal_per_mol(), 0.6163314008174533, "Incorrect kcal/mol RT conversion");
    close(energy.evaluate({3, 5, 2, 0, 4}).q, 16.874038846968144,
          "DP loop energy does not match 10.4 kcal/mol / RT");
    close(energy.evaluate({9, 11, 2, 0, 4}).q, energy.evaluate({3, 5, 2, 0, 4}).q,
          "DP incorrectly counts all stem pairs as pseudoloop boundaries");
    auto normal_dp = energy.evaluate({3, 5, 2, 0, 4}).q;
    energy.temperature_celsius = 0;
    close(energy.evaluate({3, 5, 2, 0, 4}).q / normal_dp, 310.15 / 273.15,
          "DP temperature scaling is incorrect");
    energy.temperature_celsius = -273.15;
    invalid([&] { energy.validate(); }, "Absolute zero was accepted");
    energy.temperature_celsius = std::numeric_limits<double>::quiet_NaN();
    invalid([&] { energy.validate(); }, "Nonfinite temperature was accepted");
    energy.temperature_celsius = 37;
    invalid([&] { energy.evaluate({1, 5, 2, 0, 4}); }, "One-pair stem was accepted");
    invalid([&] { energy.evaluate({3, 5, -1, 0, 4}); }, "Negative loop was accepted");

    energy.model = PKLoopEnergyModel::CC;
    auto value = energy.evaluate({3, 5, 1, 0, 2});
    require(value.source == PKLoopEnergySource::CC06, "Known CC06 shape fell back");
    close(value.q, 2.3 + 6.5 + std::log(9.0), "CC06 stem/loop numbering or assembly term is wrong");
    close(energy.evaluate({3, 5, 1, 1, 2}).q, value.q,
          "CC06 incorrectly charged the short middle loop");
    require(energy.evaluate({5, 3, 1, 0, 2}).source == PKLoopEnergySource::DPFallback,
            "CC06 major/minor groove orientation was lost");
    close(energy.evaluate({4, 4, 1, 0, 3}).q, 4.4 + 9.2 + std::log(9.0),
          "Finite starred CC06 table values were rejected");
    require(energy.evaluate({6, 5, 1, 0, 7}).q > energy.evaluate({6, 5, 1, 0, 12}).q,
            "CC06 minor-groove nonmonotonic loop behavior was lost");
    require(energy.evaluate({3, 5, 12, 0, 2}).source == PKLoopEnergySource::DPFallback,
            "An unspecified CC06 table cell was treated as zero or extrapolated");
    require(energy.evaluate({13, 5, 1, 0, 2}).source == PKLoopEnergySource::DPFallback,
            "CC06 extrapolated beyond its stem-length domain");
    close(energy.evaluate({3, 5, 13, 1, 14}).q, 18.810628356513146,
          "Published CC06 long-loop fits are incorrect");
    require(energy.evaluate({3, 5, 13, 1, 14}).source == PKLoopEnergySource::CC06,
            "CC06 long-loop fit was not used");
    energy.temperature_celsius = 0;
    close(energy.evaluate({3, 5, 1, 0, 2}).q, value.q,
          "Dimensionless CC entropy should not change with temperature");
    energy.temperature_celsius = 37;

    TemporaryFile table{std::filesystem::temp_directory_path() /
        ("ipknot-cc09-test-" + std::to_string(reinterpret_cast<size_t>(&energy)) + ".tsv")};
    table.write("# CC09 direct + long-loop fixtures, q=-DeltaS/R\n"
                "Q 3 3 1 2 1 6.8137\n"
                "L1 3 3 2 1 1.25 4.5\n"
                "L3 3 3 1 2 2.0 5.0\n"
                "L3 3 3 9 2 3.0 6.0\n");
    energy.load_cc09_table(table.path.string());
    require(energy.cc09_table_size() == 4, "Incorrect CC09 record count");
    close(energy.evaluate({3, 3, 1, 2, 1}).q, 6.8137,
          "CC09 direct joint entropy or loop numbering is wrong");
    require(energy.evaluate({3, 3, 1, 2, 1}).source == PKLoopEnergySource::CC09,
            "Known CC09 direct shape fell back");
    close(energy.evaluate({3, 3, 9, 2, 1}).q, 1.25 * std::log(9.0) + 4.5,
          "CC09 long-L1 fit key or long-loop length is wrong");
    close(energy.evaluate({3, 3, 1, 2, 11}).q, 2.0 * std::log(11.0) + 5.0,
          "CC09 long-outer-loop fit key or loop numbering is wrong");
    close(energy.evaluate({3, 3, 9, 2, 11}).q, 3.0 * std::log(11.0) + 6.0,
          "CC09 two-long-loop case did not use its actual-L1 keyed outer-loop fit");
    auto copy = energy;
    require(copy.cc09_table_size() == 4, "CC09 data did not survive option copying");
    for (auto geometry : {std::array<int, 5>{3, 3, 2, 2, 1}, {2, 3, 1, 2, 1},
                         {3, 3, 1, 7, 1}, {3, 3, 100, 2, 11}, {3, 3, 0, 2, 1}}) {
      auto fallback = energy.evaluate(geometry);
      require(fallback.source == PKLoopEnergySource::DPFallback,
              "Unknown CC09 geometry did not fall back to DP");
      auto dp = energy;
      dp.model = PKLoopEnergyModel::DP;
      close(fallback.q, dp.evaluate(geometry).q, "DP fallback uses a different energy");
    }
    for (const char* invalid_row : {
             "Q 3 3 1 2 1 6.8\nQ 3 3 1 2 1 6.8\n",
             "Q 3 3 1 2 1 -1\n", "Q 3 3 1 2 1 nan\n",
             "Q 3 3 8 2 1 6.8\n", "Q 3 3 1 7 1 6.8\n",
             "L3 3 3 100 2 1 4\n", "L1 3 3 2 8 1 4\n",
             "L1 3 3 2 1 -1 4\n", "L1 3 3 2 1 0 -1\n",
             "Q 3 3 1 2 1 6.8 trailing\n", "BAD 3 3 1 2 1 6.8\n", "# empty\n"}) {
      table.write(invalid_row);
      invalid([&] { energy.load_cc09_table(table.path.string()); }, "Invalid CC09 record was accepted");
      require(energy.cc09_table_size() == 4, "Failed table reload corrupted the old data");
    }
    std::cout << "Verified DP conversion, CC06 tables/fits, CC09 joint/fitted lookups, fallback and import validation\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
