/*
 * $Id$
 * 
 * Copyright (C) 2010 Kengo Sato
 *
 * This file is part of IPknot.
 *
 * IPknot is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * IPknot is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with IPknot.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef IPKNOT_H
#define IPKNOT_H

#include <vector>
#include <utility>
#include <string>
#include <list>

class IP;

using uint = unsigned int;
using VF = std::vector<float>;
using VVF = std::vector<VF>;
using VI = std::vector<int>;
using VVI = std::vector<VI>;
using VVVI = std::vector<VVI>;
using SVI = std::vector<std::pair<uint,int>>;
using VSVI = std::vector<SVI>; 
using VVSVI = std::vector<VSVI>;
using SVF = std::vector<std::pair<uint,float>>;
using VSVF = std::vector<SVF>; 

// Structure for base pair constraints
struct BPConstraints {
  int GC = -1;  // -1 means no constraint
  int AU = -1;
  int GU = -1;
  int UU = -1;

  bool has_constraints() const {
    return GC >= 0 || AU >= 0 || GU >= 0 || UU >= 0;
  }
};

// Structure for stacking constraints
// Represents a pattern of consecutive stacked base pairs specified by base types
struct StackConstraint {
  std::vector<std::string> bp_types;  // e.g., ["GC", "AU", "GU"] for a 3-bp stack

  StackConstraint() = default;

  // Add a base pair type to the stack
  void add_bp_type(const std::string& bp_type) {
    bp_types.push_back(bp_type);
  }

  // Check if valid (at least 2 base pairs)
  bool is_valid() const {
    return bp_types.size() >= 2;
  }

  size_t size() const {
    return bp_types.size();
  }
};

// Represents a concrete instance of a stack constraint in the sequence
struct StackInstance {
  std::vector<std::pair<int, int>> pairs;  // Actual (i, j) positions (0-indexed)
  int constraint_id;  // Which constraint this instance belongs to

  StackInstance(int id) : constraint_id(id) {}

  void add_pair(int i, int j) {
    if (i > j) std::swap(i, j);
    pairs.emplace_back(i, j);
  }

  size_t size() const {
    return pairs.size();
  }
};

// Container for multiple stack constraints
struct StackConstraints {
  std::vector<StackConstraint> constraints;
  std::vector<StackInstance> instances;  // All found instances in the sequence

  bool has_constraints() const {
    return !constraints.empty();
  }

  void add_constraint(const StackConstraint& constraint) {
    if (constraint.is_valid()) {
      constraints.push_back(constraint);
    }
  }

  void add_instance(const StackInstance& instance) {
    instances.push_back(instance);
  }

  void clear_instances() {
    instances.clear();
  }
};

class IPknot
{
public:
  template < class T > class EnumParam;

public:
  IPknot(uint pk_level, const float* alpha,
         bool levelwise, bool stacking_constraints, int n_th);

public:
  void solve(const std::string& seq, const VF& bp, const VI& offset,
             const VF& th, VI& bpseq, VI& plevel, bool constraint,
             const BPConstraints& bp_constraints = BPConstraints(),
             const StackConstraints& stack_constraints = StackConstraints()) const;

  void solve(const std::string& seq, const VSVF& bp,
             const VF& th, VI& bpseq, VI& plevel, bool constraint,
             const BPConstraints& bp_constraints = BPConstraints(),
             const StackConstraints& stack_constraints = StackConstraints()) const;

  auto solve(const std::string& seq, const VSVF& bp,
             EnumParam<float>& ep, VI& bpseq, VI& plevel, bool constraint,
             const BPConstraints& bp_constraints = BPConstraints(),
             const StackConstraints& stack_constraints = StackConstraints()) const -> std::pair<float,float>;

  static int decompose_plevel(const std::vector<int>& bpseq, std::vector<int>& plevel);

  auto check_pseudoknots(const VI& bpseq) -> VI;

  static uint length(const std::string& seq);
  static uint length(const std::list<std::string>& aln);

private:
  void solve(const std::string& seq, IP& ip, const VVSVI& v_l, const VVSVI& v_r, const VI& c_l, const VI& c_r,
             const VF& th, VI& bpseq, VI& plevel, bool constraint,
             const BPConstraints& bp_constraints,
             const StackConstraints& stack_constraints) const;

  static auto compute_expected_accuracy(float etp, float etn, float efp, float efn) -> std::tuple<float,float,float,float>;
  static auto compute_expected_accuracy(const VI& bpseq, const VF& bp, const VI& offset) -> std::tuple<float,float,float,float>;
  static auto compute_expected_accuracy(const VI& bpseq, const VSVF& bp) -> std::tuple<float,float,float,float>;
  static auto compute_expected_accuracy_pk(const VI& bpseq, const VSVF& bp) -> std::tuple<float,float,float,float>;
  static auto compute_sump_pk(const VSVF& bp) -> float;
  static auto compute_expected_accuracy_pk(const VI& bpseq, const VSVF& bp, float sump) -> std::tuple<float,float,float,float>;

private:
  // options
  uint pk_level_;
  std::vector<float> alpha_;
  bool levelwise_;
  bool stacking_constraints_;
  int n_th_;
};

template < class T >
class IPknot::EnumParam
{
public:
  EnumParam(const std::vector<std::vector<T> >& p);

  uint size() const;
  void get(std::vector<T>& q) const;
  bool succ();

private:
  static bool succ(int n, const int* m, int* v);

private:
  const std::vector<std::vector<T> >& p_;
  std::vector<int> m_;
  std::vector<int> v_;
};

#endif // IPKNOT_H