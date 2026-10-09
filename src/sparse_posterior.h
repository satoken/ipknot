#ifndef IPKNOT_SPARSE_POSTERIOR_H
#define IPKNOT_SPARSE_POSTERIOR_H

#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

// Append and accumulate in the original visitation order. The hash tables
// store vector indices, so reallocation cannot invalidate the lookup and
// changing the lookup never changes floating-point addition order.
class SparsePosteriorAccumulator {
public:
  using Matrix = std::vector<std::vector<std::pair<unsigned int, float>>>;
  explicit SparsePosteriorAccumulator(Matrix& matrix)
      : matrix_(matrix), index_(matrix.size()) {
    for (std::size_t i = 0; i < matrix.size(); ++i) {
      index_[i].reserve(matrix[i].size());
      for (std::size_t k = 0; k < matrix[i].size(); ++k)
        index_[i].emplace(matrix[i][k].first, k);
    }
  }
  void add(std::size_t row, unsigned int partner, float contribution) {
    auto& values = matrix_[row];
    auto [it, inserted] = index_[row].emplace(partner, values.size());
    if (inserted) values.emplace_back(partner, contribution);
    else values[it->second].second += contribution;
  }
private:
  Matrix& matrix_;
  std::vector<std::unordered_map<unsigned int, std::size_t>> index_;
};

#endif
