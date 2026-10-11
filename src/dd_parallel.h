#ifndef IPKNOT_DD_PARALLEL_H
#define IPKNOT_DD_PARALLEL_H

#include "dual_decomposition.h"

// Each level owns its workspace and output. Aggregate scores and update shared
// primal/dual state in level order after joining, preserving serial arithmetic.
void dd_decode_levels(const std::vector<std::unique_ptr<DDNussinov>>& decoders,
    const std::vector<double>& weights, const std::vector<unsigned char>& allowed,
    std::vector<std::vector<int>>& selected, std::vector<double>& scores, int threads);

#endif
