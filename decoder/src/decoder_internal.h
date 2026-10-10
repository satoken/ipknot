// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef IPKNOT_DECODER_INTERNAL_H
#define IPKNOT_DECODER_INTERNAL_H
#include "ipknot/decoder.h"

namespace ipknot::detail {
// IPknot refinement adds posterior contributions. Preserve those nonnegative
// scores without applying the public probability API's upper limit of one.
DecodeResult decode_posterior_scores(const SparsePosterior &posterior, DecoderOptions options);
} // namespace ipknot::detail
#endif
