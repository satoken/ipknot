#include "nmr_sequence_index.h"
#include <cassert>
#include <cctype>
#include <random>
#include <utility>

static char normalize_base(char c) {
  c = std::toupper(static_cast<unsigned char>(c));
  return c == 'T' ? 'U' : c;
}
static std::string normalize_pair(char a, char b) {
  a = normalize_base(a); b = normalize_base(b);
  if ((a == 'G' && b == 'C') || (a == 'C' && b == 'G')) return "GC";
  if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G')) return "GU";
  if (a > b) std::swap(a, b);
  return std::string{a, b};
}

int main() {
  std::mt19937 random(12);
  const std::string alphabet = "ACGUTacgutNX";
  for (size_t length = 0; length < 120; ++length) {
    std::string sequence;
    for (size_t i = 0; i < length; ++i) sequence.push_back(alphabet[random() % alphabet.size()]);
    NMRSequenceIndex index(sequence, normalize_base, normalize_pair);
    for (size_t gap = 1; gap < 10; ++gap) {
      std::vector<uint16_t> types;
      for (size_t i = 0; i < alphabet.size(); ++i)
        if (random() % 3 == 0)
          types.push_back(NMRSequenceIndex::type_code(normalize_pair(alphabet[i], alphabet[random() % alphabet.size()])));
      std::vector<std::pair<size_t, size_t>> actual, expected;
      index.for_each_pair(types, gap, [&](size_t i, size_t j) { actual.emplace_back(i, j); });
      for (size_t i = 0; i < length; ++i)
        for (size_t j = i + gap; j < length; ++j) {
          const uint16_t code = NMRSequenceIndex::type_code(normalize_pair(sequence[i], sequence[j]));
          assert(index.pair_type(i, j) == code);
          if (std::find(types.begin(), types.end(), code) != types.end()) expected.emplace_back(i, j);
        }
      assert(actual == expected);
    }
  }
  // An absent type on a long input should require only linear preprocessing /
  // row enumeration, not 100,000^2 unsuccessful pair comparisons.
  const std::string sequence(100000, 'A');
  NMRSequenceIndex index(sequence, normalize_base, normalize_pair);
  size_t count = 0;
  index.for_each_pair({NMRSequenceIndex::type_code("GC")}, 4,
      [&](size_t, size_t) { ++count; });
  assert(count == 0);
}
