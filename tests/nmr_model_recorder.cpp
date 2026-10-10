// Test-only solver substitute: fingerprint every model without optimizing it.
// Linked with the normal CLI/engines to compare complete formulations.
#ifdef NMR_RECORDER_LEGACY_IP_API
// Retain the pre-DD ABI when linking immutable baseline object files.
#include "nmr_legacy_ip.h"
#else
#include "ip.h"
#endif
#include <cstdint>
#include <cstring>
#include <iostream>
#include <vector>

class IPimpl {
public:
  uint64_t hash = 14695981039346656037ULL;
  size_t columns = 0, rows = 0, coefficients = 0;
  void value(double value) {
#ifndef NMR_RECORD_COUNTS_ONLY
    uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    for (int i = 0; i < 8; ++i) { hash ^= (bits >> (8 * i)) & 255; hash *= 1099511628211ULL; }
#endif
  }
};
#ifdef NMR_RECORDER_LEGACY_IP_API
IP::IP(DirType direction, int threads) : impl_(new IPimpl) { impl_->value(direction); }
#else
IP::IP(DirType direction, int threads, bool exact) : impl_(new IPimpl) { impl_->value(direction); }
IP::IP(IPModel& model) : impl_(nullptr), model_(&model) {}
bool IP::available() { return false; } // A recorder is not an optimization backend.
int IP::make_continuous_variable(double coefficient, double lower, double upper) {
  if (model_) {
    model_->variables.push_back({coefficient, lower, upper, false});
    return static_cast<int>(model_->variables.size()) - 1;
  }
  impl_->value(4); impl_->value(coefficient); impl_->value(lower); impl_->value(upper);
  return impl_->columns++;
}
void IP::add_objective_coefficient(int column, double coefficient) {
  if (model_) { model_->variables.at(column).coefficient += coefficient; return; }
  impl_->value(5); impl_->value(column); impl_->value(coefficient);
}
void IP::mark_noe_variable(int column) {
  if (model_) {
    if (column < 0 || column >= static_cast<int>(model_->variables.size()))
      throw std::invalid_argument("Invalid NOE variable tag");
    model_->noe_columns.push_back(column);
  }
}
#endif
IP::~IP() { delete impl_; }
int IP::make_variable(double coefficient) { return make_variable(coefficient, 0, 1); }
int IP::make_variable(double coefficient, int lower, int upper) {
#ifndef NMR_RECORDER_LEGACY_IP_API
  if (model_) {
    model_->variables.push_back({coefficient, double(lower), double(upper), true});
    return static_cast<int>(model_->variables.size()) - 1;
  }
#endif
  impl_->value(1); impl_->value(coefficient); impl_->value(lower); impl_->value(upper);
  return impl_->columns++;
}
int IP::make_constraint(BoundType bound, double lower, double upper) {
#ifndef NMR_RECORDER_LEGACY_IP_API
  if (model_) {
    model_->rows.push_back({bound, lower, upper, {}});
    return static_cast<int>(model_->rows.size()) - 1;
  }
#endif
  impl_->value(2); impl_->value(bound); impl_->value(lower); impl_->value(upper);
  return impl_->rows++;
}
void IP::add_constraint(int row, int column, double value) {
#ifndef NMR_RECORDER_LEGACY_IP_API
  if (model_) { model_->rows.at(row).terms.emplace_back(column, value); return; }
#endif
  impl_->value(3); impl_->value(row); impl_->value(column); impl_->value(value);
  ++impl_->coefficients;
}
void IP::update() {}
double IP::solve() {
#ifndef NMR_RECORDER_LEGACY_IP_API
  if (model_) {
    if (!model_->optimize) throw std::logic_error("Missing recorded-model decoder");
    return model_->optimize();
  }
#endif
  std::cout << "MODEL " << impl_->columns << ' ' << impl_->rows << ' '
            << impl_->coefficients << ' ' << impl_->hash << '\n';
  return 0;
}
double IP::get_value(int column) const {
#ifndef NMR_RECORDER_LEGACY_IP_API
  if (model_) return model_->solution.at(column);
#endif
  return 0;
}

#ifndef NMR_RECORDER_LEGACY_IP_API
struct IPModelSolver::Impl {};
IPModelSolver::IPModelSolver(const IPModel&, IP::DirType, int) {
  throw std::logic_error("The model recorder cannot optimize IP models");
}
IPModelSolver::~IPModelSolver() = default;
IPModelSolver::Result IPModelSolver::solve(const std::vector<double>&) {
  throw std::logic_error("The model recorder cannot optimize IP models");
}
double IPModelSolver::get_value(int) const {
  throw std::logic_error("The model recorder cannot optimize IP models");
}
#endif
