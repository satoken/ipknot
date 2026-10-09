#ifndef IPKNOT_DD_CONSTRAINT_BUILDER_H
#define IPKNOT_DD_CONSTRAINT_BUILDER_H
#include "ipknot.h"

double decode_linear_constraints(
    const std::string &sequence, const VSVF &posterior, const VF &thresholds,
    const VF &alpha, int levels, bool stacking, bool canonical_neighbor,
    bool coaxial, const NMRConstraintOptions &nmr,
    const DDOptions &options, VI &bpseq, VI &plevel, bool fixed,
    const BPConstraints &counts, const StackConstraints &stacks);
#endif
