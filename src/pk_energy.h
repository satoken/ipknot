#ifndef IPKNOT_PK_ENERGY_H
#define IPKNOT_PK_ENERGY_H

#include <array>
#include <cstddef>
#include <memory>
#include <string>

enum class PKLoopEnergyModel { None, DP, CC };
enum class PKLoopEnergySource { None, DP, CC06, CC09, DPFallback };

struct PKLoopEnergyValue {
  // Dimensionless loop cost G_loop/(RT), not an IPknot objective coefficient.
  double q = 0.0;
  PKLoopEnergySource source = PKLoopEnergySource::None;
};

// Pure H motif geometry: (S1, S2, L1, middle L2, outer L3), in sequence order
// S1_left, L1, S2_left, L2, S1_right, L3, S2_right. The CC papers instead call
// the outer loop L2 and the middle loop L3. No ordinary helix energy is added.
// CC06 tables/fits are built in. CC09 joint entropies/fits can be imported once
// and are shared by copies of this object; absent CC geometries use DP.
class PKLoopEnergy {
public:
  PKLoopEnergyModel model = PKLoopEnergyModel::None;
  // This sets the conversion scale for the decoder correction only; it does
  // not change the temperature used by the base-pair probability calculation.
  double temperature_celsius = 37.0;

  PKLoopEnergyValue evaluate(const std::array<int, 5>& geometry) const;
  double rt_kcal_per_mol() const;
  void validate() const;

  // Compact CC09 facts, with # comments and one record per line:
  // Q  S1 S2 L1 middle outer q
  // L1 S1 S2 middle fixed_outer q_a q_b
  // L3 S1 S2 fixed_L1 middle q_a q_b
  // A fit gives q=q_a*ln(long_loop_length)+q_b. Q contains direct entries.
  // Unknown shapes remain eligible and use the DP loop cost.
  void load_cc09_table(const std::string& filename);
  size_t cc09_table_size() const;

private:
  struct CC09Data;
  std::shared_ptr<const CC09Data> cc09_;
};

const char* pk_loop_energy_source_name(PKLoopEnergySource source);

#endif
