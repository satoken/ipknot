#include "dd_constraint_builder.h"
#include "bpseq.h"
#include "dd_constrained.h"
#include "ip.h"
#include "spdlog/spdlog.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

namespace {
using Pair = std::pair<int, int>;
std::uint64_t key(Pair p) {
  return (std::uint64_t(unsigned(p.first)) << 32) | unsigned(p.second);
}
bool crosses(Pair a, Pair b) {
  return (a.first < b.first && b.first < a.second && a.second < b.second) ||
         (b.first < a.first && a.first < b.second && b.second < a.second);
}
bool encloses(Pair a, Pair b) {
  return a.first < b.first && b.second < a.second;
}
bool shares(Pair a, Pair b) {
  return a.first == b.first || a.first == b.second || a.second == b.first ||
         a.second == b.second;
}
struct Physical {
  Pair pair;
  std::string type;
  double probability;
  std::vector<int> columns;
};
struct Witness {
  std::vector<Pair> required;
  double score = 0;
  bool bulged = false, coaxial = false;
  CoaxialKind kind = CoaxialKind::CLOSING_FIRST_CHILD;
  HelixFace face1 = HelixFace::INNER, face2 = HelixFace::OUTER;
};
struct Builder {
  const std::string &seq;
  const VSVF &posterior;
  const VF &th;
  const VF &alpha;
  int levels, length;
  const DDOptions &options;
  const VI &fixed;
  bool constrained;
  IPModel model;
  IP ip;
  std::vector<Physical> physical;
  std::unordered_map<std::uint64_t, int> lookup;
  std::unordered_map<std::uint64_t, double> probabilities;
  std::array<std::vector<int>, 256> positions;
  std::vector<unsigned char> alphabet;
  std::size_t pair_drops = 0, witness_drops = 0, seed_visits = 0;

  Builder(const std::string &s, const VSVF &bp, const VF &t, const VF &a, int k,
          const DDOptions &o, const VI &f, bool c)
      : seq(s), posterior(bp), th(t), alpha(a), levels(k), length(s.size()),
        options(o), fixed(f), constrained(c), ip(model) {
    for (int i = 0; i < length; ++i) {
      const auto base = static_cast<unsigned char>(normalize_base(seq[i]));
      if (positions[base].empty())
        alphabet.push_back(base);
      positions[base].push_back(i);
    }
    for (int i = 1; i <= length; ++i)
      for (const auto &[j, p] : posterior[i]) {
        if (j == 0 || j >= posterior.size() || j == static_cast<unsigned>(i) ||
            !std::isfinite(p) || p < 0)
          throw std::invalid_argument("Invalid sparse posterior pair for DD");
        if (i < static_cast<int>(j))
          probabilities.emplace(key({i - 1, int(j) - 1}), p);
      }
  }
  bool compatible(Pair p) const {
    if (p.first < 0 || p.second >= length || p.first >= p.second)
      return false;
    if (!constrained)
      return true;
    for (auto end : {p.first, p.second}) {
      const int other = end == p.first ? p.second : p.first, value = fixed[end];
      if (value == BPSEQ::U || (value >= 0 && value != other) ||
          (value == BPSEQ::L && other < end) ||
          (value == BPSEQ::R && other > end))
        return false;
    }
    return true;
  }
  double probability(Pair p) const {
    auto it = probabilities.find(key(p));
    return it == probabilities.end() ? 0 : it->second;
  }
  int add(Pair p, double probability, bool all) {
    if (!compatible(p))
      return -1;
    auto inserted = lookup.emplace(key(p), physical.size());
    if (inserted.second)
      physical.push_back({p,
                          normalize_base_pair_type(seq[p.first], seq[p.second]),
                          probability, std::vector<int>(levels, -1)});
    auto &pair = physical[inserted.first->second];
    for (int lv = 0; lv < levels; ++lv) {
      if (pair.columns[lv] < 0 && (all || probability > th[lv]))
        pair.columns[lv] = ip.make_variable((probability - th[lv]) * alpha[lv]);
    }
    return inserted.first->second;
  }
  void add_pair_sum(int row, int id, double coefficient) {
    for (int col : physical[id].columns)
      if (col >= 0)
        ip.add_constraint(row, col, coefficient);
  }
  template <class F> void sampled_pairs(const std::string &type, F visit) {
    std::array<std::size_t, 256> cursors{};
    for (int i = 0; i < length; ++i)
      for (unsigned char base : alphabet) {
        if (normalize_base_pair_type(seq[i], char(base)) != type)
          continue;
        const auto &list = positions[base];
        auto &start = cursors[base];
        while (start < list.size() && list[start] < i + 4)
          ++start;
        const std::size_t available = list.size() - start,
                          retained = std::min<std::size_t>(
                              available, options.nmr_pair_beam);
        pair_drops += available - retained;
        for (std::size_t t = 0; t < retained; ++t) {
          const std::size_t offset =
              retained > 1 ? t * (available - 1) / (retained - 1) : 0;
          visit(Pair{i, list[start + offset]});
          ++seed_visits;
        }
      }
  }
  void keep(std::vector<Witness> &beam, Witness candidate) {
    for (const auto &old : beam)
      if (old.required == candidate.required &&
          old.coaxial == candidate.coaxial &&
          (!old.coaxial || old.kind == candidate.kind))
        return;
    beam.push_back(std::move(candidate));
    if (beam.size() > static_cast<std::size_t>(options.nmr_witnesses)) {
      auto better = [](const Witness &a, const Witness &b) {
        return a.score != b.score ? a.score > b.score : a.required < b.required;
      };
      beam.erase(std::max_element(beam.begin(), beam.end(), better));
      ++witness_drops;
    }
  }
  // Direct matching needs one complementary right-strand string for each
  // left window. This determines fallback globally without a quadratic scan.
  bool direct_exists(const std::vector<std::string> &types) const {
    const int width = types.size();
    if (width < 2 || length < 2 * width + 3)
      return false;
    struct Ends {
      std::vector<int> positions;
      std::size_t cursor = 0;
    };
    std::unordered_map<std::string, Ends> index;
    for (int end = width - 1; end < length; ++end) {
      std::string word;
      word.reserve(width);
      for (int d = width - 1; d >= 0; --d)
        word.push_back(normalize_base(seq[end - d]));
      index[word].positions.push_back(end);
    }
    for (int i = 0; i + 2 * width + 3 <= length; ++i) {
      std::string word(width, ' ');
      bool valid = true;
      for (int d = 0; d < width; ++d) {
        const char base = normalize_base(seq[i + d]);
        const auto &type = types[d];
        if (base != type[0] && base != type[1]) {
          valid = false;
          break;
        }
        word[width - 1 - d] = base == type[0] ? type[1] : type[0];
      }
      if (!valid)
        continue;
      auto found = index.find(word);
      if (found == index.end())
        continue;
      auto &ends = found->second;
      while (ends.cursor < ends.positions.size() &&
             ends.positions[ends.cursor] < i + 4 + 2 * (width - 1))
        ++ends.cursor;
      if (ends.cursor < ends.positions.size())
        return true;
    }
    return false;
  }
  std::vector<Witness> stack_witnesses(const StackConstraint &pattern,
                                       int observation) {
    std::vector<Witness> direct, bulged;
    const int width = pattern.size();
    if (width < 2 || width > (length - 3) / 2)
      return {};
    for (int direction = 0; direction < 2; ++direction) {
      auto types = pattern.bp_types;
      if (direction) {
        std::reverse(types.begin(), types.end());
        if (types == pattern.bp_types)
          continue;
      }
      auto extend = [&](Pair seed) {
        if (seed.second - seed.first < 4 + 2 * (width - 1) || !compatible(seed))
          return;
        std::vector<Witness> frontier{
            {{seed}, probability(seed), false, false}};
        for (int depth = 1; depth < width && !frontier.empty(); ++depth) {
          std::vector<Witness> next;
          for (const auto &prefix : frontier) {
            const auto last = prefix.required.back();
            for (Pair step : {Pair{1, 1}, Pair{2, 1}, Pair{1, 2}}) {
              if (options.nmr_bulge_mode == 0 && step != Pair{1, 1})
                continue;
              const Pair p{last.first + step.first, last.second - step.second};
              if (p.second - p.first < 4 + 2 * (width - depth - 1) ||
                  !compatible(p) ||
                  normalize_base_pair_type(seq[p.first], seq[p.second]) !=
                      types[depth])
                continue;
              auto candidate = prefix;
              candidate.required.push_back(p);
              candidate.score += probability(p);
              candidate.bulged |= step != Pair{1, 1};
              next.push_back(std::move(candidate));
              if (next.size() >
                  static_cast<std::size_t>(options.nmr_pattern_beam)) {
                auto better = [](const Witness &a, const Witness &b) {
                  if (a.bulged != b.bulged)
                    return !a.bulged;
                  return a.score != b.score ? a.score > b.score
                                            : a.required < b.required;
                };
                next.erase(std::max_element(next.begin(), next.end(), better));
                ++witness_drops;
              }
            }
          }
          frontier = std::move(next);
        }
        for (auto &witness : frontier)
          keep(witness.bulged ? bulged : direct, std::move(witness));
      };
      for (const auto &p : physical)
        if (p.type == types.front())
          extend(p.pair);
      sampled_pairs(types.front(), [&](Pair p) {
        if (!lookup.count(key(p)))
          extend(p);
      });
    }
    if (options.nmr_bulge_mode == 1) {
      auto reverse = pattern.bp_types;
      std::reverse(reverse.begin(), reverse.end());
      if (direct_exists(pattern.bp_types) || direct_exists(reverse)) {
        spdlog::info("NMR bulge fallback: constraint {} uses {} direct "
                     "instance(s); discarded {} bulged instance(s)",
                     observation + 1, direct.size(), bulged.size());
        return direct;
      }
      spdlog::info("NMR bulge fallback: constraint {} has no direct instance; "
                   "retained {} bulged instance(s)",
                   observation + 1, bulged.size());
      return bulged;
    }
    for (auto &w : bulged)
      keep(direct, std::move(w));
    return direct;
  }
};
std::array<Pair, 2> stem_continuations(const Witness &w) {
  auto continuation = [](Pair pair, HelixFace face) {
    return face == HelixFace::INNER
        ? Pair{pair.first - 1, pair.second + 1}
        : Pair{pair.first + 1, pair.second - 1};
  };
  return {continuation(w.required[0], w.face1),
          continuation(w.required[1], w.face2)};
}
bool blocks(const Witness &w, Pair p) {
  for (auto required : w.required)
    if (shares(p, required))
      return false;
  const auto a = w.required[0], b = w.required[1], s = w.required[2];
  if (crosses(p, a) || crosses(p, b) || crosses(p, s))
    return true;
  if (w.kind == CoaxialKind::CLOSING_FIRST_CHILD)
    return encloses(a, p) && (encloses(p, b) || encloses(p, s));
  if (w.kind == CoaxialKind::ADJACENT_CHILDREN)
    return encloses(s, p) && (encloses(p, a) || encloses(p, b));
  return encloses(b, p) && (encloses(p, a) || encloses(p, s));
}
} // namespace

namespace {
void coaxial_witnesses(Builder &b, const StackConstraints &stacks,
                       const std::vector<int> &ordinary,
                       std::vector<std::vector<Witness>> &witnesses) {
  std::unordered_set<std::string> types;
  for (const auto &pattern : stacks.constraints)
    if (pattern.size() == 2)
      types.insert(pattern.bp_types.begin(), pattern.bp_types.end());
  struct Terminal {
    Pair pair;
    std::string type;
    double score;
  };
  std::vector<Terminal> terminals;
  std::unordered_set<std::uint64_t> seen;
  std::vector<std::vector<int>> by_left(b.length), by_right(b.length);
  auto add = [&](Pair p) {
    if (!b.compatible(p) || p.second - p.first < 4 ||
        !seen.insert(key(p)).second)
      return;
    const int id = terminals.size();
    terminals.push_back(
        {p, normalize_base_pair_type(b.seq[p.first], b.seq[p.second]),
         b.probability(p)});
    auto offer = [&](std::vector<int> &row) {
      row.push_back(id);
      if (row.size() > static_cast<std::size_t>(b.options.nmr_pair_beam)) {
        auto better = [&](int a, int c) {
          return terminals[a].score != terminals[c].score
                     ? terminals[a].score > terminals[c].score
                     : terminals[a].pair < terminals[c].pair;
        };
        row.erase(std::max_element(row.begin(), row.end(), better));
        ++b.pair_drops;
      }
    };
    offer(by_left[p.first]);
    offer(by_right[p.second]);
  };
  for (const auto &p : b.physical)
    if (types.count(p.type))
      add(p.pair);
  for (const auto &type : types)
    b.sampled_pairs(type, add);
  for (int observation = 0;
       observation < static_cast<int>(stacks.constraints.size());
       ++observation) {
    const auto &pattern = stacks.constraints[observation];
    if (pattern.size() != 2)
      continue;
    std::vector<Witness> contexts;
    for (int direction = 0; direction < 2; ++direction) {
      const auto &type1 = pattern.bp_types[direction ? 1 : 0];
      const auto &type2 = pattern.bp_types[direction ? 0 : 1];
      if (direction && type1 == type2)
        continue;
      for (const auto &first : terminals)
        if (first.type == type1) {
          const auto a = first.pair;
          auto offer = [&](int other, CoaxialKind kind, HelixFace face1,
                           HelixFace face2) {
            const auto &second = terminals[other];
            if (second.type != type2)
              return;
            const auto c = second.pair;
            if (kind == CoaxialKind::CLOSING_FIRST_CHILD && !encloses(a, c))
              return;
            if (kind == CoaxialKind::LAST_CHILD_CLOSING && !encloses(c, a))
              return;
            Witness context;
            context.required = {a, c};
            context.score = first.score + second.score;
            context.coaxial = true;
            context.kind = kind;
            context.face1 = face1;
            context.face2 = face2;
            b.keep(contexts, std::move(context));
          };
          if (a.first + 1 < b.length)
            for (int id : by_left[a.first + 1])
              offer(id, CoaxialKind::CLOSING_FIRST_CHILD, HelixFace::INNER,
                    HelixFace::OUTER);
          if (a.second + 1 < b.length) {
            for (int id : by_left[a.second + 1])
              offer(id, CoaxialKind::ADJACENT_CHILDREN, HelixFace::OUTER,
                    HelixFace::OUTER);
            for (int id : by_right[a.second + 1])
              offer(id, CoaxialKind::LAST_CHILD_CLOSING, HelixFace::OUTER,
                    HelixFace::INNER);
          }
        }
    }
    for (const auto &context : contexts)
      for (int support_id : ordinary) {
        const auto s = b.physical[support_id].pair, a = context.required[0],
                   c = context.required[1];
        bool valid = context.kind == CoaxialKind::CLOSING_FIRST_CHILD
                         ? c.second < s.first && s.second < a.second
                     : context.kind == CoaxialKind::ADJACENT_CHILDREN
                         ? encloses(s, a) && encloses(s, c)
                         : c.first < s.first && s.second < a.first;
        if (!valid)
          continue;
        auto witness = context;
        witness.required.push_back(s);
        witness.score += b.probability(s);
        b.keep(witnesses[observation], std::move(witness));
      }
  }
}
} // namespace

double decode_linear_constraints(
    const std::string &sequence, const VSVF &posterior, const VF &thresholds,
    const VF &alpha, int levels, bool stacking, bool canonical_neighbor,
    bool coaxial, const NMRConstraintOptions &nmr,
    const DDOptions &options, VI &bpseq, VI &plevel, bool fixed,
    const BPConstraints &counts, const StackConstraints &stacks) {
  options.validate();
  const int length = sequence.size();
  if (posterior.size() != sequence.size() + 1 ||
      thresholds.size() != static_cast<std::size_t>(levels) ||
      alpha.size() != thresholds.size() ||
      (fixed && bpseq.size() != sequence.size()))
    throw std::invalid_argument("Invalid constrained DD input dimensions");
  if (options.crossing_beam == 0 ||
      options.witnesses == 0 || options.constraint_states == 0)
    throw std::invalid_argument(
        "Linear constrained DD needs positive crossing beams and repair budget; use "
        "--dd-constraints full for unlimited diagnostics");
  if (fixed)
    for (int i = 0; i < length; ++i)
      if (bpseq[i] >= length || bpseq[i] == i || bpseq[i] < BPSEQ::LR ||
          (bpseq[i] >= 0 && bpseq[bpseq[i]] != i))
        throw std::invalid_argument("Invalid fixed structure for DD");
  const VI specification = fixed ? bpseq : VI(length, BPSEQ::DOT);
  Builder b(sequence, posterior, thresholds, alpha, levels, options,
            specification, fixed);
  for (int i = 1; i <= length; ++i)
    for (const auto &[j, p] : posterior[i])
      if (i < static_cast<int>(j)) {
        bool eligible = false;
        for (int lv = 0; lv < levels; ++lv)
          eligible |= p > thresholds[lv];
        const Pair pair{i - 1, int(j) - 1};
        if (eligible || (fixed && specification[i - 1] == int(j) - 1))
          b.add(pair, p, fixed && specification[i - 1] == int(j) - 1);
      }
  if (fixed)
    for (int i = 0; i < length; ++i)
      if (specification[i] > i)
        b.add({i, specification[i]}, b.probability({i, specification[i]}), true);
  std::vector<int> ordinary;
  for (int id = 0; id < static_cast<int>(b.physical.size()); ++id)
    ordinary.push_back(id);
  std::unordered_set<std::string> noncanonical;
  for (const auto &[type, count] : counts.constraints)
    if (count >= 0 && !is_canonical_base_pair(type))
      noncanonical.insert(type);
  for (const auto &pattern : stacks.constraints)
    for (const auto &type : pattern.bp_types)
      if (!is_canonical_base_pair(type))
        noncanonical.insert(type);
  if (canonical_neighbor) {
    for (int id : ordinary) {
      const auto p = b.physical[id].pair;
      // Match the existing level-0 candidate-neighbor admission rule.
      if (b.physical[id].columns[0] < 0)
        continue;
      for (Pair neighbor :
           {Pair{p.first - 1, p.second + 1}, Pair{p.first + 1, p.second - 1}})
        if (neighbor.first >= 0 && neighbor.second < length &&
            neighbor.second - neighbor.first >= 4 &&
            noncanonical.count(normalize_base_pair_type(
                sequence[neighbor.first], sequence[neighbor.second])))
          b.add(neighbor, 0, true);
    }
  } else
    for (const auto &type : noncanonical)
      b.sampled_pairs(type, [&](Pair p) { b.add(p, 0, true); });
  std::vector<std::vector<Witness>> witnesses(stacks.constraints.size());
  for (int observation = 0;
       observation < static_cast<int>(stacks.constraints.size()); ++observation)
    witnesses[observation] =
        b.stack_witnesses(stacks.constraints[observation], observation);
  if (coaxial)
    coaxial_witnesses(b, stacks, ordinary, witnesses);
  bool bulged = false;
  std::size_t retained = 0, coaxial_count = 0;
  for (auto &group : witnesses)
    for (auto &witness : group) {
      for (const auto &p : witness.required)
        b.add(p, 0, true);
      bulged |= witness.bulged;
      coaxial_count += witness.coaxial;
      ++retained;
    }
  if (nmr.coaxial_no_adjacent_bulge) {
    std::size_t removed = 0;
    for (auto &group : witnesses) {
      auto end = std::remove_if(group.begin(), group.end(), [&](const Witness &w) {
        if (!w.coaxial) return false;
        for (const auto &stem : stem_continuations(w)) {
          if (!b.lookup.count(key(stem)) || blocks(w, stem)) return true;
          for (const auto &required : w.required)
            if (shares(stem, required) || crosses(stem, required)) return true;
        }
        return false;
      });
      removed += std::distance(end, group.end());
      group.erase(end, group.end());
    }
    retained -= removed;
    coaxial_count -= removed;
    spdlog::info("NMR coaxial no-adjacent-bulge prior: {} retained witnesses removed", removed);
  }
  spdlog::info("Found {} stack instances in sequence",
               retained - coaxial_count);
  if (coaxial)
    spdlog::info("Found {} flush coaxial-stacking instances", coaxial_count);
  std::vector<DDPair> pairs;
  std::vector<int> columns;
  std::vector<std::vector<int>> left(levels * length), right(levels * length),
      endpoints(length);
  std::unordered_map<std::string, std::vector<int>> by_type;
  for (const auto &physical : b.physical)
    for (int lv = 0; lv < levels; ++lv)
      if (physical.columns[lv] >= 0) {
        const int col = physical.columns[lv], i = physical.pair.first,
                  j = physical.pair.second;
        pairs.push_back({i, j, lv, b.model.variables[col].coefficient});
        columns.push_back(col);
        left[lv * length + i].push_back(col);
        right[lv * length + j].push_back(col);
        endpoints[i].push_back(col);
        endpoints[j].push_back(col);
        by_type[physical.type].push_back(col);
      }
  for (int i = 0; i < length; ++i) {
    const int row =
        b.ip.make_constraint(IP::UP, 0, specification[i] == BPSEQ::U ? 0 : 1);
    for (int col : endpoints[i])
      b.ip.add_constraint(row, col, 1);
    const int value = specification[i];
    if (value == BPSEQ::L || value == BPSEQ::R || value == BPSEQ::LR) {
      if (endpoints[i].empty())
        throw DDInfeasible("Fixed paired-base constraint has no candidate "
                           "partner in retained DD graph");
      const int required = b.ip.make_constraint(IP::FX, 1, 1);
      for (int col : endpoints[i])
        b.ip.add_constraint(required, col, 1);
    } else if (value > i) {
      const int required = b.ip.make_constraint(IP::FX, 1, 1);
      b.add_pair_sum(required, b.lookup.at(key({i, value})), 1);
    }
  }
  if (stacking) {
    const int distance = bulged ? 2 : 1;
    spdlog::debug("Stacking neighbor distance: {}", distance);
    for (int lv = 0; lv < levels; ++lv)
      for (int i = 0; i < length; ++i)
        for (auto *side : {&left, &right}) {
          const int row = b.ip.make_constraint(IP::LO, 0, 0);
          for (int col : (*side)[lv * length + i])
            b.ip.add_constraint(row, col, -1);
          for (int d = 1; d <= distance; ++d)
            for (int neighbor : {i - d, i + d})
              if (neighbor >= 0 && neighbor < length)
                for (int col : (*side)[lv * length + neighbor])
                  b.ip.add_constraint(row, col, 1);
        }
  }
  const auto graph = dd_bounded_graph(length, pairs, levels, options);
  for (const auto &support : graph.rows) {
    const int upper = columns[support.upper],
              row = b.ip.make_constraint(IP::LO, 0, 0);
    b.ip.add_constraint(row, upper, -1);
    for (int id : support.contacts)
      b.ip.add_constraint(row, columns[id], 1);
  }
  struct CountSlack {
    std::string type;
    int requested, missing, excess;
  };
  std::vector<CountSlack> slacks;
  for (const auto &[type, expected] : counts.constraints)
    if (expected >= 0) {
      const bool lower = nmr.count_mode == NMRCountMode::LOWER_BOUND;
      const int row =
          b.ip.make_constraint(lower ? IP::LO : IP::FX, expected, expected);
      for (int col : by_type[type])
        b.ip.add_constraint(row, col, 1);
      if (nmr.soft) {
        const int missing = b.ip.make_variable(-nmr.count_penalty, 0,
                                               std::max(length, expected));
        b.ip.add_constraint(row, missing, 1);
        int excess = -1;
        if (!lower) {
          excess = b.ip.make_variable(-nmr.count_penalty, 0,
                                      std::max(length, expected));
          b.ip.add_constraint(row, excess, -1);
        }
        slacks.push_back({type, expected, missing, excess});
      }
    }
  std::vector<std::vector<int>> pair_witnesses(b.physical.size());
  std::unordered_map<std::uint64_t, std::vector<int>> faces;
  std::vector<std::vector<int>> observation_columns(witnesses.size());
  std::vector<int> violations(witnesses.size(), -1);
  for (std::size_t observation = 0; observation < witnesses.size();
       ++observation) {
    std::vector<std::vector<int>> blockers(b.physical.size());
    for (const auto &witness : witnesses[observation]) {
      const int w = b.ip.make_variable(0);
      b.ip.mark_noe_variable(w);
      observation_columns[observation].push_back(w);
      for (const auto &p : witness.required) {
        const int id = b.lookup.at(key(p));
        if (!nmr.allow_shared_stack_pairs)
          pair_witnesses[id].push_back(w);
        else {
          const int link = b.ip.make_constraint(IP::UP, 0, 0);
          b.ip.add_constraint(link, w, 1);
          b.add_pair_sum(link, id, -1);
        }
      }
      if (witness.coaxial) {
        if (nmr.coaxial_no_adjacent_bulge) {
          // Continuation pairs support this context but do not consume an
          // observation's pair-sharing capacity, matching the full formulation.
          for (const auto &stem : stem_continuations(witness)) {
            const int row = b.ip.make_constraint(IP::UP, 0, 0);
            b.ip.add_constraint(row, w, 1);
            b.add_pair_sum(row, b.lookup.at(key(stem)), -1);
          }
        }
        for (int id = 0; id < static_cast<int>(b.physical.size()); ++id)
          if (blocks(witness, b.physical[id].pair))
            blockers[id].push_back(w);
        if (nmr.allow_shared_stack_pairs) {
          faces[std::uint64_t(b.lookup.at(key(witness.required[0]))) * 2 +
                unsigned(witness.face1)]
              .push_back(w);
          faces[std::uint64_t(b.lookup.at(key(witness.required[1]))) * 2 +
                unsigned(witness.face2)]
              .push_back(w);
        }
      }
    }
    if (observation_columns[observation].empty() && !nmr.soft)
      throw DDInfeasible(
          "Stack constraint cannot be satisfied: no valid instances found in "
          "retained DD graph; raise --dd-nmr-witnesses/--dd-nmr-pair-beam or "
          "use --dd-constraints full");
    const int row = b.ip.make_constraint(IP::FX, 1, 1);
    for (int col : observation_columns[observation])
      b.ip.add_constraint(row, col, 1);
    if (nmr.soft) {
      violations[observation] = b.ip.make_variable(-nmr.stack_penalty);
      b.ip.mark_noe_variable(violations[observation]);
      b.ip.add_constraint(row, violations[observation], 1);
    }
    for (int id = 0; id < static_cast<int>(blockers.size()); ++id)
      if (!blockers[id].empty()) {
        const int link = b.ip.make_constraint(IP::UP, 0, 1);
        for (int w : blockers[id])
          b.ip.add_constraint(link, w, 1);
        b.add_pair_sum(link, id, 1);
      }
  }
  if (!nmr.allow_shared_stack_pairs)
    for (int id = 0; id < static_cast<int>(pair_witnesses.size()); ++id)
      if (!pair_witnesses[id].empty()) {
        const int row = b.ip.make_constraint(IP::UP, 0, 0);
        for (int w : pair_witnesses[id])
          b.ip.add_constraint(row, w, 1);
        b.add_pair_sum(row, id, -1);
      }
  for (const auto &entry : faces)
    if (entry.second.size() > 1) {
      const int row = b.ip.make_constraint(IP::UP, 0, 1);
      for (int w : entry.second)
        b.ip.add_constraint(row, w, 1);
    }
  // The empty structure is always a useful seed when all NMR observations
  // are soft and there are no mandatory fixed pairs. The decoder independently
  // validates this seed, including every topology and capacity row.
  if (nmr.soft) {
    b.model.solution.assign(b.model.variables.size(), 0);
    for (const auto &slack : slacks)
      b.model.solution[slack.missing] = slack.requested;
    for (int violation : violations)
      if (violation >= 0)
        b.model.solution[violation] = 1;
  }
  const auto result =
      solve_constrained_dd(length, pairs, columns, levels, b.model, options);
  if (options.noe_ilp)
    spdlog::info("DD NOE ILP: variables={}, rows={}, calls={}, cache_hits={}, time={:.6f}s, "
                 "primal_calls={}, primal_feasible={}, primal_time={:.6f}s",
                 result.noe_ilp_variables, result.noe_ilp_rows, result.noe_ilp_calls,
                 result.noe_ilp_cache_hits, result.noe_ilp_seconds,
                 result.noe_primal_calls, result.noe_primal_feasible, result.noe_primal_seconds);
  bpseq.assign(length, -1);
  plevel.assign(length, -1);
  double penalty = 0;
  for (std::size_t id = 0; id < pairs.size(); ++id)
    if (b.model.solution[columns[id]] > .5) {
      const auto &p = pairs[id];
      bpseq[p.left] = p.right;
      bpseq[p.right] = p.left;
      plevel[p.left] = plevel[p.right] = p.level;
    }
  for (const auto &slack : slacks) {
    const double missing = b.model.solution[slack.missing],
                 excess =
                     slack.excess >= 0 ? b.model.solution[slack.excess] : 0;
    double actual = 0;
    for (int col : by_type[slack.type])
      actual += b.model.solution[col];
    const double cost = nmr.count_penalty * (missing + excess);
    penalty += cost;
    if (missing > .5 || excess > .5)
      spdlog::info("Soft NMR base-pair constraint {} violated: expected {}, "
                   "actual {:.0f}, missing {:.0f}, excess {:.0f}, penalty {}",
                   slack.type, slack.requested, actual, missing, excess, cost);
    else
      spdlog::info("Soft NMR base-pair constraint {}={} satisfied", slack.type,
                   slack.requested);
  }
  bool satisfied = true;
  for (std::size_t observation = 0; observation < witnesses.size();
       ++observation) {
    if (violations[observation] >= 0 &&
        b.model.solution[violations[observation]] > .5) {
      penalty += nmr.stack_penalty;
      satisfied = false;
      spdlog::info("Soft NMR stack constraint {} violated, penalty {}",
                   observation + 1, nmr.stack_penalty);
    } else {
      if (nmr.soft)
        spdlog::info("Soft NMR stack constraint {} satisfied", observation + 1);
      for (std::size_t id = 0; id < witnesses[observation].size(); ++id)
        if (b.model.solution[observation_columns[observation][id]] > .5) {
          if (witnesses[observation][id].coaxial)
            spdlog::info(
                "Stack constraint {} satisfied by flush coaxial stacking",
                observation + 1);
          else
            spdlog::info("Stack constraint {} satisfied by instance",
                         observation + 1);
        }
    }
  }
  if (stacks.has_constraints())
    spdlog::info(satisfied ? "All stack constraints are satisfied."
                           : "Some stack constraints are NOT satisfied "
                             "(allowed in soft NMR mode).");
  spdlog::info("DD: iterations={}, pairs={}, constraint_rows={}, nonzeros={}, "
               "objective={:.12g}, graph_upper_bound={:.12g}, "
               "repair_states={}, propagation_work={}, structural_work={}, "
               "nmr_seed_visits={}, nmr_pair_drops={}, nmr_witness_drops={}, "
               "crossing_beam_drops={}, witness_drops={}, linear=1, stop={}",
               result.iterations, pairs.size(), b.model.rows.size(),
               result.nonzeros, result.objective, result.upper_bound,
               result.repair_states, result.propagation_work,
               result.structural_work, b.seed_visits, b.pair_drops,
               b.witness_drops, graph.crossing_drops, graph.witness_drops,
               result.stop_reason);
  return penalty;
}
