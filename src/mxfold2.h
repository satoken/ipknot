#pragma once

#include <vector>
#include <string>
#include <pybind11/embed.h>
#include "fold.h"
namespace py = pybind11;
#pragma GCC visibility push(hidden)

class MXfold2Model : public BPEngineSeq
{
public:
  MXfold2Model(uint n_th=1, const std::string& config="", int gpu=-1);

  void calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::string& seq, float th=DEFAULT_THRESHOLD) const -> std::vector<std::vector<std::pair<uint, float>>>;
  void calculate_posterior(const std::string& seq, const std::string& paren, std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::string& seq, const std::string& paren, float th=DEFAULT_THRESHOLD) const -> std::vector<std::vector<std::pair<uint, float>>>;

private:
  py::scoped_interpreter guard_;
  py::object model_;
};
