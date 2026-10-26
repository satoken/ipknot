#ifndef IPKNOT_DD_CONSTRAINT_BUILDER_H
#define IPKNOT_DD_CONSTRAINT_BUILDER_H
#include "ipknot.h"

void decode_linear_constraints(
    const std::string &sequence, const VSVF &posterior, const VF &thresholds,
    const VF &alpha, int levels, bool stacking,
    const DDOptions &options, VI &bpseq, VI &plevel, bool fixed);
#endif
