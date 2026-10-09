#include "pk_energy.h"

#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>

// Numerical entropy facts from Cao & Chen (2006), Tables 1 and 2:
// https://pmc.ncbi.nlm.nih.gov/articles/PMC1463895/ (doi:10.1093/nar/gkl346).
// Verified against HotKnots 2.0 bin/params/pkmodelCC2006.dat, SHA256
// a85819050e56ab29f955d74fd8493484a0bf8a6cd874a79134441239f23745fe.
// Unspecified cells are represented by -1; finite starred approximations in
// the original paper are retained. We do not extrapolate across stem lengths.
namespace {
constexpr double gas_constant = 0.00198720425864083; // kcal/(mol K)
constexpr double unknown = -1.0;
constexpr double major[11][12] = {
    {unknown, unknown, unknown, 6.2, 6.4, 6.4, 6.6, 6.8, 6.9, 7.1, 7.2, unknown},
    {unknown, 6.4, 6.4, 6.4, 6.6, 6.6, 6.8, 6.9, 7.1, 7.3, 7.5, unknown},
    {4.4, 4.4, 4.5, 5.4, 5.6, 6.0, 6.3, 6.6, 6.9, 7.1, 7.3, unknown},
    {2.3, 4.4, 4.6, 5.7, 6.0, 6.5, 6.9, 7.2, 7.5, 7.8, 8.0, unknown},
    {2.3, 4.4, 4.8, 5.8, 6.0, 6.5, 6.8, 7.1, 7.4, 7.6, 7.8, unknown},
    {2.3, 4.4, 5.0, 5.9, 6.2, 6.8, 7.0, 7.3, 7.6, 7.8, 8.0, unknown},
    {unknown, 4.4, 5.2, 5.7, 6.4, 6.7, 7.1, 7.3, 7.5, 7.7, 7.9, unknown},
    {unknown, 5.5, 5.5, 6.4, 6.7, 7.2, 7.5, 7.9, 8.1, 8.3, 8.5, unknown},
    {unknown, 6.9, 6.9, 6.9, 7.5, 7.7, 8.1, 8.3, 8.6, 8.8, 8.9, unknown},
    {unknown, unknown, unknown, unknown, 8.7, 8.8, 8.9, 9.1, 9.2, 9.3, 9.3, unknown},
    {unknown, unknown, unknown, unknown, 9.8, 9.2, 9.5, 9.6, 9.7, 9.8, 9.8, unknown},
};

constexpr double minor[11][12] = {
    {unknown, unknown, unknown, 7.6, 7.0, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.7},
    {unknown, 6.5, 6.5, 6.5, 6.6, 6.7, 6.9, 7.1, 7.2, 7.4, 7.6, 7.7},
    {unknown, unknown, 9.2, 9.2, 9.2, 8.9, 8.9, 8.9, 9.0, 9.0, 9.1, 9.2},
    {unknown, unknown, unknown, 9.8, 9.8, 9.8, 9.1, 8.9, 8.8, 8.8, 8.8, 8.8},
    {unknown, unknown, unknown, 11.9, 11.9, 11.9, 11.9, 11.0, 10.4, 10.1, 9.9, 9.8},
    {unknown, unknown, unknown, unknown, 12.4, 12.4, 12.4, 12.4, 11.4, 11.0, 10.7, 10.5},
    {unknown, unknown, unknown, unknown, 12.1, 12.1, 12.1, 12.1, 11.6, 11.4, 11.2, 11.1},
    {unknown, unknown, unknown, unknown, unknown, 13.7, 13.7, 13.7, 13.7, 12.6, 12.0, 11.5},
    {unknown, unknown, unknown, unknown, unknown, 13.7, 13.7, 13.7, 12.7, 12.2, 11.8, 11.5},
    {unknown, unknown, unknown, unknown, unknown, unknown, unknown, unknown, 15.9, 14.1, 13.0, 12.4},
    {unknown, unknown, unknown, unknown, unknown, unknown, unknown, unknown, 18.7, 15.8, 14.2, 13.2},
};

struct EntropyFit { int minimum; double a; double b; double c; };
constexpr EntropyFit major_fit[11] = {
    {4, 0.12, 1.96, 0.52},
    {2, 0.39, 1.92, -3.89},
    {1, -2.14, 2.15, -2.09},
    {1, -2.22, 2.11, -2.25},
    {1, -2.4, 2.18, -2.33},
    {1, -2.61, 2.21, -2.32},
    {2, -1.17, 2.03, -1.96},
    {2, -1.66, 2.09, -1.98},
    {2, -1.43, 2.09, -2.93},
    {5, -0.14, 2.06, 0.15},
    {5, 0.77, 1.84, -0.65},
};

constexpr EntropyFit minor_fit[11] = {
    {4, 0.95, 1.84, -0.67},
    {2, 0.32, 1.92, -3.9},
    {3, 1.77, 1.82, -5.76},
    {4, 3.99, 1.55, -5.86},
    {4, 7.73, 1.29, -12.67},
    {5, 8.38, 1.16, -11.45},
    {5, 4.52, 1.61, -7.58},
    {6, 9.05, 1.15, -11.45},
    {6, 4.77, 1.68, -6.78},
    {9, 2.74, 2.05, 1.38},
    {9, 4.69, 1.8, -1.11},
};

void validate_geometry(const std::array<int, 5>& g) {
  if (g[0] < 2 || g[1] < 2 || g[2] < 0 || g[3] < 0 || g[4] < 0) {
    throw std::invalid_argument("PK loop energy requires stems >= 2 and loops >= 0");
  }
}

double entropy(int stem, int loop, const double table[11][12],
               const EntropyFit fit[11]) {
  if (stem < 2 || stem > 12 || loop <= 0) return unknown;
  if (loop <= 12) return table[stem - 2][loop - 1];
  const auto& f = fit[stem - 2];
  const double extension = static_cast<double>(loop) - f.minimum + 1.0;
  const double result = 2.14 * loop + 0.10 - f.a * std::log(extension)
      - f.b * extension - f.c;
  // The fitted model is intended for short and moderate loops. Never turn a
  // nonphysical extrapolated entropy into a stabilizing energy bonus.
  return std::isfinite(result) && result > 0 ? result : unknown;
}

bool valid_cc09_stems(int a, int b) {
  return a >= 3 && a <= 10 && b >= 3 && b <= 10;
}

bool valid_cc09_middle(int middle) { return middle >= 2 && middle <= 6; }

} // namespace

struct PKLoopEnergy::CC09Data {
  std::map<std::array<int, 5>, double> direct;
  // Long L1 keys: S1,S2,middle,fixed_outer. Long outer L3 keys:
  // S1,S2,fixed_L1,middle. These match the author distribution's fits.
  std::map<std::array<int, 4>, std::pair<double, double>> long_l1;
  std::map<std::array<int, 4>, std::pair<double, double>> long_l3;
};

const char* pk_loop_energy_source_name(PKLoopEnergySource source) {
  switch (source) {
  case PKLoopEnergySource::None: return "none";
  case PKLoopEnergySource::DP: return "dp";
  case PKLoopEnergySource::CC06: return "cc06";
  case PKLoopEnergySource::CC09: return "cc09";
  case PKLoopEnergySource::DPFallback: return "dp-fallback";
  }
  throw std::invalid_argument("Unknown PK loop energy source");
}

double PKLoopEnergy::rt_kcal_per_mol() const {
  const double kelvin = temperature_celsius + 273.15;
  if (!std::isfinite(temperature_celsius) || !std::isfinite(kelvin) || kelvin <= 0) {
    throw std::invalid_argument("PK energy temperature must be finite and above absolute zero");
  }
  return gas_constant * kelvin;
}

void PKLoopEnergy::validate() const {
  if (model != PKLoopEnergyModel::None && model != PKLoopEnergyModel::DP &&
      model != PKLoopEnergyModel::CC) {
    throw std::invalid_argument("Unknown PK loop energy model");
  }
  rt_kcal_per_mol();
}

PKLoopEnergyValue PKLoopEnergy::evaluate(const std::array<int, 5>& g) const {
  validate_geometry(g);
  if (model == PKLoopEnergyModel::None) return {};
  // DP-style loop component, with IPknot's bundled NUPACK parameter values
  // beta1=9.6, beta2=0.1, beta3=0.1 kcal/mol. A simple two-band H motif has
  // two interior boundary pairs, hence beta1+2*beta2+beta3*unpaired.
  // Nested pseudoknots, multiloop contexts and ordinary helix terms are not
  // included. The H-domain geometry uses loop spans when scored in the ILP.
  const double unpaired = static_cast<double>(g[2]) + g[3] + g[4];
  const double dp = (9.8 + 0.1 * unpaired) / rt_kcal_per_mol();
  if (model == PKLoopEnergyModel::DP) return {dp, PKLoopEnergySource::DP};
  if (model != PKLoopEnergyModel::CC) {
    throw std::invalid_argument("Unknown PK loop energy model");
  }

  if (g[3] <= 1) {
    // CC06 L1 crosses the major groove of S2; the outer IPknot L3 crosses
    // the minor groove of S1. Its short middle segment contributes no term.
    const double first = entropy(g[1], g[2], major, major_fit);
    const double outer = entropy(g[0], g[4], minor, minor_fit);
    if (first > 0 && outer > 0) {
      return {first + outer + std::log(9.0), PKLoopEnergySource::CC06};
    }
  } else if (cc09_ && valid_cc09_stems(g[0], g[1]) &&
             valid_cc09_middle(g[3]) && g[2] > 0 && g[4] > 0) {
    // CC09 joint entropy already includes all three loop correlations; no
    // independent CC06 assembly term is added (DotKnot's CC09 convention).
    if (g[2] <= 7 && g[4] <= 7) {
      auto it = cc09_->direct.find(g);
      if (it != cc09_->direct.end()) return {it->second, PKLoopEnergySource::CC09};
    } else {
      const bool long_outer = g[4] > 7;
      const auto& fits = long_outer ? cc09_->long_l3 : cc09_->long_l1;
      std::array<int, 4> key = long_outer
          ? std::array<int, 4>{g[0], g[1], g[2], g[3]}
          : std::array<int, 4>{g[0], g[1], g[3], g[4]};
      auto it = fits.find(key);
      if (it != fits.end()) {
        const double q = it->second.first * std::log(long_outer ? g[4] : g[2])
            + it->second.second;
        if (std::isfinite(q) && q > 0) return {q, PKLoopEnergySource::CC09};
      }
    }
  }
  return {dp, PKLoopEnergySource::DPFallback};
}

size_t PKLoopEnergy::cc09_table_size() const {
  return cc09_ ? cc09_->direct.size() + cc09_->long_l1.size() + cc09_->long_l3.size() : 0;
}

void PKLoopEnergy::load_cc09_table(const std::string& filename) {
  std::ifstream input(filename);
  if (!input) throw std::invalid_argument("Cannot read CC09 entropy table: " + filename);
  auto data = std::make_shared<CC09Data>();
  std::string line;
  int number = 0;
  while (std::getline(input, line)) {
    ++number;
    const auto comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
    std::istringstream row(line);
    std::string type, extra;
    bool valid = static_cast<bool>(row >> type);
    if (valid && type == "Q") {
      std::array<int, 5> g{};
      double q = 0;
      valid = static_cast<bool>(row >> g[0] >> g[1] >> g[2] >> g[3] >> g[4] >> q)
          && !(row >> extra) && valid_cc09_stems(g[0], g[1])
          && valid_cc09_middle(g[3]) && g[2] >= 1 && g[2] <= 7
          && g[4] >= 1 && g[4] <= 7 && std::isfinite(q) && q > 0;
      if (valid) valid = data->direct.emplace(g, q).second;
    } else if (valid && (type == "L1" || type == "L3")) {
      std::array<int, 4> key{};
      double a = 0, b = 0;
      valid = static_cast<bool>(row >> key[0] >> key[1] >> key[2] >> key[3] >> a >> b)
          && !(row >> extra) && valid_cc09_stems(key[0], key[1])
          && std::isfinite(a) && std::isfinite(b) && a >= 0;
      if (type == "L1") {
        valid = valid && valid_cc09_middle(key[2]) && key[3] >= 1 && key[3] <= 7;
      } else {
        valid = valid && key[2] >= 1 && key[2] <= 99 && valid_cc09_middle(key[3]);
      }
      valid = valid && a * std::log(8.0) + b > 0;
      if (valid) {
        auto& fits = type == "L1" ? data->long_l1 : data->long_l3;
        valid = fits.emplace(key, std::make_pair(a, b)).second;
      }
    } else {
      valid = false;
    }
    if (!valid) {
      throw std::invalid_argument("Invalid or duplicate CC09 entropy table row "
                                  + std::to_string(number));
    }
  }
  if (input.bad()) throw std::invalid_argument("Error reading CC09 entropy table: " + filename);
  if (data->direct.empty() && data->long_l1.empty() && data->long_l3.empty()) {
    throw std::invalid_argument("CC09 entropy table is empty");
  }
  // Publish only after every row passed validation: failed reloads preserve
  // the previous immutable data, which is shared across option copies.
  cc09_ = std::move(data);
}
