// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef IPKNOT_DECODER_ILP_MODEL_H
#define IPKNOT_DECODER_ILP_MODEL_H

#include "ip.h"
#include "ipknot/decoder.h"

namespace ipknot::detail {
using PairVariables = std::vector<std::vector<std::vector<std::pair<unsigned int, int>>>>;

// Shared by the reusable decoder and IPknot's constrained/NMR formulation.
void add_level_constraints(IP &ip, const PairVariables &left);
void add_stacking_constraints(IP &ip, const PairVariables &left, const PairVariables &right,
                              int neighbor_distance = 1);
DecodeResult decode_ilp(int length, const std::vector<DDPair> &pairs, int levels,
                        bool no_lonely_pairs, int threads);
} // namespace ipknot::detail
#endif
