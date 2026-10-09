#ifndef IPKNOT_NMR_SEQUENCE_INDEX_H
#define IPKNOT_NMR_SEQUENCE_INDEX_H

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Type matching is indexed by the sequence alphabet, not by all L*L pairs.
// Preprocessing/storage are O(L + alphabet_size^2); enumeration emits exactly
// the requested types, in the former i-then-j order. With a fixed alphabet its
// cost is O(L + emitted_pairs), which can itself be quadratic for dense types.
class NMRSequenceIndex {
public:
  static uint16_t type_code(const std::string& type) {
    return static_cast<uint16_t>((static_cast<unsigned char>(type[0]) << 8) |
                                static_cast<unsigned char>(type[1]));
  }

  template<class NormalizeBase, class NormalizePair>
  NMRSequenceIndex(const std::string& sequence, NormalizeBase normalize_base,
                   NormalizePair normalize_pair) {
    std::array<int, 256> ids;
    ids.fill(-1);
    std::vector<char> alphabet;
    base_ids_.reserve(sequence.size());
    for (char base : sequence) {
      base = normalize_base(base);
      int& id = ids[static_cast<unsigned char>(base)];
      if (id < 0) {
        id = static_cast<int>(alphabet.size());
        alphabet.push_back(base);
      }
      base_ids_.push_back(static_cast<size_t>(id));
    }
    alphabet_size_ = alphabet.size();
    pair_types_.resize(alphabet_size_ * alphabet_size_);
    for (size_t a = 0; a < alphabet_size_; ++a)
      for (size_t b = 0; b < alphabet_size_; ++b)
        pair_types_[a * alphabet_size_ + b] = type_code(normalize_pair(alphabet[a], alphabet[b]));
    positions_.resize(alphabet_size_);
    for (size_t i = 0; i < base_ids_.size(); ++i) positions_[base_ids_[i]].push_back(i);
  }

  uint16_t pair_type(size_t left, size_t right) const {
    return pair_types_[base_ids_[left] * alphabet_size_ + base_ids_[right]];
  }

  template<class Emit>
  void for_each_pair(std::vector<uint16_t> types, size_t minimum_gap, Emit emit) const {
    std::sort(types.begin(), types.end());
    std::vector<std::vector<size_t>> partners(alphabet_size_);
    for (size_t a = 0; a < alphabet_size_; ++a)
      for (size_t b = 0; b < alphabet_size_; ++b)
        if (std::binary_search(types.begin(), types.end(), pair_types_[a * alphabet_size_ + b]))
          partners[a].push_back(b);
    std::vector<size_t> offsets(alphabet_size_, 0), cursors(alphabet_size_);
    const size_t length = base_ids_.size();
    if (minimum_gap >= length) return;
    for (size_t i = 0; i < length - minimum_gap; ++i) {
      const auto& ids = partners[base_ids_[i]];
      for (size_t id : ids) {
        auto& offset = offsets[id];
        while (offset < positions_[id].size() && positions_[id][offset] < i + minimum_gap) ++offset;
        cursors[id] = offset;
      }
      if (ids.size() == 1) {
        const size_t id = ids.front();
        for (size_t k = cursors[id]; k < positions_[id].size(); ++k) emit(i, positions_[id][k]);
      } else if (!ids.empty()) {
        // A fixed-alphabet merge preserves j order without sorting each row.
        while (true) {
          size_t next = length, next_id = 0;
          for (size_t id : ids)
            if (cursors[id] < positions_[id].size() && positions_[id][cursors[id]] < next) {
              next = positions_[id][cursors[id]];
              next_id = id;
            }
          if (next == length) break;
          emit(i, next);
          ++cursors[next_id];
        }
      }
    }
  }

private:
  size_t alphabet_size_;
  std::vector<size_t> base_ids_;
  std::vector<uint16_t> pair_types_;
  std::vector<std::vector<size_t>> positions_;
};

#endif
