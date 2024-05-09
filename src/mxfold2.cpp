#include <iostream>
#include "mxfold2.h"
#include <pybind11/embed.h>
namespace py = pybind11;
using namespace py::literals;


static
auto parse_structure(const std::string& stru)
{
  std::vector<int> stack;
  std::vector<int> ret(stru.size()+1, 0);
  for (uint i=0; i!=stru.size(); ++i)
  {
    if (stru[i] == '(')
      stack.push_back(i);
    else if (stru[i] == ')')
    {
      ret[i+1] = stack.back()+1;
      ret[stack.back()+1] = i+1;
      stack.pop_back();
    }
  }
  return ret;
}

static
auto make_sbp(const py::list& bpp, float th)
{
  std::vector<std::vector<std::pair<uint, float>>> sbp(bpp.size());
  for (auto i=1; i!=bpp.size(); ++i)
  {
    for (auto bp: py::cast<py::list>(bpp[i]))
    {
      auto p = py::cast<py::tuple>(bp);
      auto j = py::cast<uint>(p[0]);
      auto pr = py::cast<float>(p[1]);
      if (pr >= th)
      {
        sbp[i].emplace_back(j, pr);
        sbp[j].emplace_back(i, pr);
      }
    }
  }
  return sbp;
}

static
auto
make_bp_offset(const std::vector<std::vector<std::pair<uint, float>>>& sbp,
              std::vector<float>& bp, std::vector<int>& offset)
{
  uint L=sbp.size()-1;
  bp.resize((L+1)*(L+2)/2);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;

  for (auto i=1; i!=sbp.size(); ++i)
    for (auto [j, pr]: sbp[i])
      if (i<j)
        bp[offset[i]+j] = pr;
}

MXfold2Model::MXfold2Model(uint n_th, const std::string& config, int gpu) : BPEngineSeq(), guard_(), model_()
{
  auto argparse = py::module_::import("argparse");
  auto __main__ = py::module_::import("mxfold2.__main__");
  auto predict = py::module_::import("mxfold2.predict");
  auto interface = py::module_::import("mxfold2.interface");
  auto torch = py::module_::import("torch");
  auto pathlib = py::module_::import("pathlib");

  auto parser = argparse.attr("ArgumentParser")("fromfile_prefix_chars"_a="@");
  auto subparser = parser.attr("add_subparsers")("title"_a="Sub-commands");
  predict.attr("Predict").attr("add_args")(subparser);
  auto l = py::list();
  l.append("predict");
  if (!config.empty())
    l.append(py::str("@") + py::str(config));
  else
    l.append(py::str("@")+__main__.attr("default_conf"));
  l.append("input");
  auto args = parser.attr("parse_args")("args"_a=l);
  auto model_conf = predict.attr("Predict")().attr("build_model")(args);
  model_ = py::cast<py::tuple>(model_conf)[0];

  auto param = pathlib.attr("Path")(args.attr("param"));
  if (!py::cast<bool>(param.attr("exists")()))
      param = pathlib.attr("Path")(__main__.attr("default_conf")).attr("parent").attr("joinpath")(param);
  auto p = torch.attr("load")(param, "map_location"_a="cpu");
  if (py::isinstance<py::dict>(p) && py::cast<py::dict>(p).contains("model_state_dict"))
    p = py::cast<py::dict>(p)["model_state_dict"];
  model_.attr("load_state_dict")(p);
  if (gpu >= 0)
    model_.attr("cuda")(gpu);

  torch.attr("set_grad_enabled")(false);
  model_.attr("eval")();

  torch.attr("set_num_threads")(n_th);
  interface.attr("set_num_threads")(n_th);
}

auto MXfold2Model::calculate_posterior(const std::string& seq, float th) const -> std::vector<std::vector<std::pair<uint, float>>>
{
  auto seq_ = py::list();
  seq_.append(seq);
  auto ret = py::cast<py::tuple>(model_(seq_, "return_partfunc"_a=true));
  auto bpps = py::cast<py::list>(ret[4]);
  auto bpp = py::cast<py::list>(bpps[0]);

  return make_sbp(bpp, th);
}

void MXfold2Model::calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const
{
  auto sbp = calculate_posterior(seq);
  make_bp_offset(sbp, bp, offset);
}

auto MXfold2Model::calculate_posterior(const std::string& seq, const std::string& paren, float th) const -> std::vector<std::vector<std::pair<uint, float>>>
{
  char symbol[] = " ";
  auto seq_ = py::list();
  seq_.append(seq);

  auto p = parse_structure(paren);
  auto paren_ = py::list();
  paren_.append(0);
  for (auto i=0; i!=paren.size(); ++i)
  {
    switch (paren[i])
    {
      case '(':
      case ')': paren_.append(p[i+1]); break;
      case 'x': paren_.append(0); break;
      case '.': paren_.append(-1); break;
      case '<': paren_.append(-2); break;
      case '>': paren_.append(-3); break;
      case '|': paren_.append(-4); break;
      default:
        throw std::runtime_error("Invalid structure");
    }
  }

  auto ret = py::cast<py::tuple>(model_(seq_, "constraint"_a=paren_, "return_partfunc"_a=true));
  auto bpps = py::cast<py::list>(ret[4]);
  auto bpp = py::cast<py::list>(bpps[0]);

  return make_sbp(bpp, th);
}

void MXfold2Model::calculate_posterior(const std::string& seq, const std::string& paren, std::vector<float>& bp, std::vector<int>& offset) const
{
  auto sbp = calculate_posterior(seq, paren);
  make_bp_offset(sbp, bp, offset);
}