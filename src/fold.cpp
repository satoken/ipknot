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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "fold.h"

#include <cstring>
#include <iostream>
#include <sstream>

#include "contrafold/SStruct.hpp"
#include "contrafold/InferenceEngine.hpp"
#include "contrafold/Defaults.hpp"

#if defined(HAVE_VIENNA18) || defined(HAVE_VIENNA20)
namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/aln_util.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/PS_dot.h>
  extern void read_parameter_file(const char fname[]);
};
};

#ifndef FLT_OR_DBL
typedef Vienna::FLT_OR_DBL FLT_OR_DBL;
#endif
#endif

#include "nupack/nupack.h"

extern "C" {
#include "boltzmann_param.h"
};

#include "linearpartition/LinearPartition.h"

typedef unsigned int uint;

// CONTRAfold model

void
CONTRAfoldModel::
calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const
{
  SStruct ss("unknown", seq);
  ParameterManager<float> pm;
  InferenceEngine<float> en(false);
  std::vector<float> w = GetDefaultComplementaryValues<float>();
  bp.resize((seq.size()+1)*(seq.size()+2)/2, 0.0);
  en.RegisterParameters(pm);
  en.LoadValues(w);
  en.LoadSequence(ss);
  en.ComputeInside();
  en.ComputeOutside();
  en.ComputePosterior();
  en.GetPosterior(0, bp, offset);
}

auto
CONTRAfoldModel::
calculate_posterior(const std::string& seq, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  std::vector<float> bp;
  std::vector<int> offset;
  this->calculate_posterior(seq, bp, offset);

  std::vector<std::vector<std::pair<uint, float>>> sbp(seq.size()+1);
  for (auto i=0; i!=seq.size(); i++)
    for (auto j=i+1; j!=seq.size(); j++)
    {
      const auto p = bp[offset[i+1]+(j+1)];
      if (p>=th)
      {
        sbp[i+1].emplace_back(j+1, p);
        sbp[j+1].emplace_back(i+1, p);
      }
    }

  return sbp;
}

void
CONTRAfoldModel::
calculate_posterior(const std::string& seq, const std::string& paren,
                    std::vector<float>& bp, std::vector<int>& offset) const
{
  SStruct ss("unknown", seq, paren);
  ParameterManager<float> pm;
  InferenceEngine<float> en(false);
  std::vector<float> w = GetDefaultComplementaryValues<float>();
  bp.resize((seq.size()+1)*(seq.size()+2)/2, 0.0);
  en.RegisterParameters(pm);
  en.LoadValues(w);
  en.LoadSequence(ss);
  en.UseConstraints(ss.GetMapping());
  en.ComputeInside();
  en.ComputeOutside();
  en.ComputePosterior();
  en.GetPosterior(0, bp, offset);
}

auto 
CONTRAfoldModel::
calculate_posterior(const std::string& seq, const std::string& paren, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  std::vector<float> bp;
  std::vector<int> offset;
  this->calculate_posterior(seq, paren, bp, offset);

  std::vector<std::vector<std::pair<uint, float>>> sbp(seq.size()+1);
  for (auto i=0; i!=seq.size(); i++)
    for (auto j=i+1; j!=seq.size(); j++)
    {
      const auto p = bp[offset[i+1]+(j+1)];
      if (p>=th)
      {
        sbp[i+1].emplace_back(j+1, p);
        sbp[j+1].emplace_back(i+1, p);
      }
    }

  return sbp;
}

// RNAfold model
#if defined(HAVE_VIENNA18) || defined(HAVE_VIENNA20)
RNAfoldModel::
RNAfoldModel(const char* param)
{
  if (param)
  {
    if (strcmp(param, "default")!=0)
      Vienna::read_parameter_file(param);
  }
  else
    copy_boltzmann_parameters();
}

void
RNAfoldModel::
calculate_posterior(const std::string& seq, const std::string& paren,
                    std::vector<float>& bp, std::vector<int>& offset) const
{
  std::string p(paren);
  std::replace(p.begin(), p.end(), '.', 'x');
  std::replace(p.begin(), p.end(), '?', '.');

  int bk = Vienna::fold_constrained;
  Vienna::fold_constrained = 1;

  uint L=seq.size();
  bp.resize((L+1)*(L+2)/2);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
#if 0
  std::string str(seq.size()+1, '.');
  float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &str[0]);
  float sfact = 1.07;
  float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
#else
  Vienna::pf_scale = -1;
#endif
#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(L);
#endif
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), &p[0]);
#ifdef HAVE_VIENNA20
  FLT_OR_DBL* pr = Vienna::export_bppm();
  int* iindx = Vienna::get_iindx(seq.size());
#else
  FLT_OR_DBL* pr = Vienna::pr;
  int* iindx = Vienna::iindx;
#endif
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
      bp[offset[i+1]+(j+1)] = pr[iindx[i+1]-(j+1)];
  Vienna::free_pf_arrays();
  Vienna::fold_constrained = bk;
}

auto
RNAfoldModel::
calculate_posterior(const std::string& seq, const std::string& paren, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  std::string p(paren);
  std::replace(p.begin(), p.end(), '.', 'x');
  std::replace(p.begin(), p.end(), '?', '.');

  int bk = Vienna::fold_constrained;
  Vienna::fold_constrained = 1;

  uint L=seq.size();
  std::vector<std::vector<std::pair<uint, float>>> sbp(L+1);
#if 0
  std::string str(seq.size()+1, '.');
  float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &str[0]);
  float sfact = 1.07;
  float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
#else
  Vienna::pf_scale = -1;
#endif
#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(L);
#endif
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), &p[0]);
#ifdef HAVE_VIENNA20
  FLT_OR_DBL* pr = Vienna::export_bppm();
  int* iindx = Vienna::get_iindx(seq.size());
#else
  FLT_OR_DBL* pr = Vienna::pr;
  int* iindx = Vienna::iindx;
#endif
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
    {
      const auto p = pr[iindx[i+1]-(j+1)];
      if (p>=th)
      {
        sbp[i+1].emplace_back(j+1, p);
        sbp[j+1].emplace_back(i+1, p);
      }
    }
  Vienna::free_pf_arrays();
  Vienna::fold_constrained = bk;

  return sbp;
}
    
void
RNAfoldModel::
calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const
{
  uint L=seq.size();
  bp.resize((L+1)*(L+2)/2);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
#if 0
  std::string str(seq.size()+1, '.');
  float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &str[0]);
  float sfact = 1.07;
  float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
#else
  Vienna::pf_scale = -1;
#endif
#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(L);
#endif
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
#ifdef HAVE_VIENNA20
  FLT_OR_DBL* pr = Vienna::export_bppm();
  int* iindx = Vienna::get_iindx(seq.size());
#else
  FLT_OR_DBL* pr = Vienna::pr;
  int* iindx = Vienna::iindx;
#endif
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
      bp[offset[i+1]+(j+1)] = pr[iindx[i+1]-(j+1)];
  Vienna::free_pf_arrays();
}

auto
RNAfoldModel::
calculate_posterior(const std::string& seq, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  uint L=seq.size();
  std::vector<std::vector<std::pair<uint, float>>> sbp(L+1);
#if 0
  std::string str(seq.size()+1, '.');
  float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &str[0]);
  float sfact = 1.07;
  float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
#else
  Vienna::pf_scale = -1;
#endif
#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(L);
#endif
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
#ifdef HAVE_VIENNA20
  FLT_OR_DBL* pr = Vienna::export_bppm();
  int* iindx = Vienna::get_iindx(seq.size());
#else
  FLT_OR_DBL* pr = Vienna::pr;
  int* iindx = Vienna::iindx;
#endif
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
    {
      const auto p = pr[iindx[i+1]-(j+1)];
      if (p>=th)
      {
        sbp[i+1].emplace_back(j+1, p);
        sbp[j+1].emplace_back(i+1, p);
      }
    }
  Vienna::free_pf_arrays();

  return sbp;
}

#else

RNAfoldModel::
RNAfoldModel(const char* param)
{
  throw "ViennaRNA package is not linked.";
}

void
RNAfoldModel::
calculate_posterior(const std::string& seq, const std::string& paren,
                    std::vector<float>& bp, std::vector<int>& offset) const
{
  throw "ViennaRNA package is not linked.";
}

auto
RNAfoldModel::
calculate_posterior(const std::string& seq, const std::string& paren, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  throw "ViennaRNA package is not linked.";
  uint L=seq.size();
  std::vector<std::vector<std::pair<uint, float>>> sbp(L+1);
  return sbp;
}
    
void
RNAfoldModel::
calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const
{
  throw "ViennaRNA package is not linked.";
}

auto
RNAfoldModel::
calculate_posterior(const std::string& seq, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  throw "ViennaRNA package is not linked.";
  uint L=seq.size();
  std::vector<std::vector<std::pair<uint, float>>> sbp(L+1);
  return sbp;
}
#endif

// Nupack model

void
NupackModel::
calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const
{
  Nupack<long double> nu;
  if (param_)
    nu.load_parameters(param_);
  else
    nu.load_default_parameters(/*model_*/);
  //nu.dump_parameters(std::cout);
  nu.load_sequence(seq);
  nu.calculate_partition_function();
  nu.calculate_posterior();
  nu.get_posterior(bp, offset);
}

auto
NupackModel::
calculate_posterior(const std::string& seq, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  std::vector<float> bp;
  std::vector<int> offset;
  this->calculate_posterior(seq, bp, offset);

  std::vector<std::vector<std::pair<uint, float>>> sbp(seq.size()+1);
  for (auto i=0; i!=seq.size(); i++)
    for (auto j=i+1; j!=seq.size(); j++)
    {
      const auto p = bp[offset[i+1]+(j+1)];
      if (p>=th)
      {
        sbp[i+1].emplace_back(j+1, p);
        sbp[j+1].emplace_back(i+1, p);
      }
    }

  return sbp;
}

// LinearPartition model

void
LinearPartitionModel::
calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const
{
  auto seq2 = seq;
  // convert to uppercase
  transform(seq2.begin(), seq2.end(), seq2.begin(), ::toupper);
  // convert T to U
  replace(seq2.begin(), seq2.end(), 'T', 'U');

  LinearPartition::BeamCKYParser parser(beam_size_, true, false, false,  0.0, true);
  if (use_vienna_)
    parser.parse<true,int>(seq2);
  else
    parser.parse<false,float>(seq2);
  parser.get_posterior(bp, offset);
}

auto
LinearPartitionModel::
calculate_posterior(const std::string& seq, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  auto seq2 = seq;
  // convert to uppercase
  transform(seq2.begin(), seq2.end(), seq2.begin(), ::toupper);
  // convert T to U
  replace(seq2.begin(), seq2.end(), 'T', 'U');

  LinearPartition::BeamCKYParser parser(beam_size_, true, false, false, th, false);
  if (use_vienna_)
    parser.parse<true,int>(seq2);
  else
    parser.parse<false,float>(seq2);
  std::vector<std::vector<std::pair<uint, float>>> bp;
  parser.get_posterior(bp);
  return bp;
}

void
LinearPartitionModel::
calculate_posterior(const std::string& seq, const std::string& paren, std::vector<float>& bp, std::vector<int>& offset) const
{
  auto seq2 = seq;
  // convert to uppercase
  std::transform(seq2.begin(), seq2.end(), seq2.begin(), ::toupper);
  // convert T to U
  std::replace(seq2.begin(), seq2.end(), 'T', 'U');

  LinearPartition::BeamCKYParser parser(beam_size_, true, false, false, 0.0, true);
  if (use_vienna_)
    parser.parse<true,int>(seq2, paren);
  else
    parser.parse<false,float>(seq2, paren);
  parser.get_posterior(bp, offset);
}

auto
LinearPartitionModel::
calculate_posterior(const std::string& seq, const std::string& paren, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  auto seq2 = seq;
  // convert to uppercase
  std::transform(seq2.begin(), seq2.end(), seq2.begin(), ::toupper);
  // convert T to U
  std::replace(seq2.begin(), seq2.end(), 'T', 'U');

  LinearPartition::BeamCKYParser parser(beam_size_, true, false, false, th, true);
  if (use_vienna_)
    parser.parse<true,int>(seq2, paren);
  else
    parser.parse<false,float>(seq2, paren);
  std::vector<std::vector<std::pair<uint, float>>> bp;
  parser.get_posterior(bp);
  return bp;
}



// Alifold model
#if defined(HAVE_VIENNA18) || defined(HAVE_VIENNA20)
AlifoldModel::
AlifoldModel(const char* param)
{
  if (param)
  {
    if (strcmp(param, "default")!=0)
      Vienna::read_parameter_file(param);
  }
  else
    copy_boltzmann_parameters();
}

static char** alloc_aln(const std::list<std::string>& aln)
{
  uint L=aln.front().size();
  char** seqs = new char*[aln.size()+1];
  seqs[aln.size()]=NULL;
  std::list<std::string>::const_iterator x;
  uint i=0;
  for (x=aln.begin(); x!=aln.end(); ++x)
  {
    seqs[i] = new char[L+1];
    strcpy(seqs[i++], x->c_str());
  }
  return seqs;
}

static void free_aln(char** seqs)
{
  for (uint i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];
  delete[] seqs;
}

void
AlifoldModel::
calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                    std::vector<float>& bp, std::vector<int>& offset) const
{
  std::string p(paren);
  std::replace(p.begin(), p.end(), '.', 'x');
  std::replace(p.begin(), p.end(), '?', '.');

  int bk = Vienna::fold_constrained;
  Vienna::fold_constrained = 1;

  //uint N=aln.size();
  uint L=aln.front().size();
  bp.resize((L+1)*(L+2)/2, 0.0);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;

  char** seqs=alloc_aln(aln);
  std::string res(p);
  // scaling parameters to avoid overflow
#ifdef HAVE_VIENNA20
  double min_en = Vienna::alifold((const char**)seqs, &res[0]);  
#else
  double min_en = Vienna::alifold(seqs, &res[0]);
#endif
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/L);
  Vienna::free_alifold_arrays();

#ifdef HAVE_VIENNA18
  Vienna::plist* pi;
#else
  Vienna::pair_info* pi;
#endif
#ifdef HAVE_VIENNA20
  Vienna::alipf_fold((const char**)seqs, &p[0], &pi);
#else
  Vienna::alipf_fold(seqs, &p[0], &pi);
#endif
  for (uint k=0; pi[k].i!=0; ++k)
    bp[offset[pi[k].i]+pi[k].j]=pi[k].p;
  free(pi);

  Vienna::free_alipf_arrays();
  free_aln(seqs);
  Vienna::fold_constrained = bk;
}

auto
AlifoldModel::
calculate_posterior(const std::list<std::string>& aln, const std::string& paren, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  std::string p(paren);
  std::replace(p.begin(), p.end(), '.', 'x');
  std::replace(p.begin(), p.end(), '?', '.');

  int bk = Vienna::fold_constrained;
  Vienna::fold_constrained = 1;

  //uint N=aln.size();
  uint L=aln.front().size();
  std::vector<std::vector<std::pair<uint, float>>> bp(L+1);

  char** seqs=alloc_aln(aln);
  std::string res(p);
  // scaling parameters to avoid overflow
#ifdef HAVE_VIENNA20
  double min_en = Vienna::alifold((const char**)seqs, &res[0]);  
#else
  double min_en = Vienna::alifold(seqs, &res[0]);
#endif
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/L);
  Vienna::free_alifold_arrays();

#ifdef HAVE_VIENNA18
  Vienna::plist* pi;
#else
  Vienna::pair_info* pi;
#endif
#ifdef HAVE_VIENNA20
  Vienna::alipf_fold((const char**)seqs, &p[0], &pi);
#else
  Vienna::alipf_fold(seqs, &p[0], &pi);
#endif
  for (uint k=0; pi[k].i!=0; ++k)
    if (pi[k].p>=th && pi[k].i<pi[k].j) 
    {
      bp[pi[k].i].emplace_back(pi[k].j, pi[k].p);
      bp[pi[k].j].emplace_back(pi[k].i, pi[k].p);
    }
  free(pi);

  Vienna::free_alipf_arrays();
  free_aln(seqs);
  Vienna::fold_constrained = bk;

  return bp;
}

void
AlifoldModel::
calculate_posterior(const std::list<std::string>& aln,
                    std::vector<float>& bp, std::vector<int>& offset) const
{

  //uint N=aln.size();
  uint L=aln.front().size();
  bp.resize((L+1)*(L+2)/2, 0.0);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;

  char** seqs=alloc_aln(aln);
  std::string res(L+1, ' ');
  // scaling parameters to avoid overflow
#ifdef HAVE_VIENNA20
  double min_en = Vienna::alifold((const char**)seqs, &res[0]);  
#else
  double min_en = Vienna::alifold(seqs, &res[0]);
#endif
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/L);
  Vienna::free_alifold_arrays();

#ifdef HAVE_VIENNA18
  Vienna::plist* pi;
#else
  Vienna::pair_info* pi;
#endif
#ifdef HAVE_VIENNA20
  Vienna::alipf_fold((const char**)seqs, NULL, &pi);
#else
  Vienna::alipf_fold(seqs, NULL, &pi);
#endif
  for (uint k=0; pi[k].i!=0; ++k)
    bp[offset[pi[k].i]+pi[k].j]=pi[k].p;
  free(pi);

  Vienna::free_alipf_arrays();
  free_aln(seqs);
}

auto
AlifoldModel::
calculate_posterior(const std::list<std::string>& aln, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{

  //uint N=aln.size();
  uint L=aln.front().size();
  std::vector<std::vector<std::pair<uint, float>>> bp(L+1);

  char** seqs=alloc_aln(aln);
  std::string res(L+1, ' ');
  // scaling parameters to avoid overflow
#ifdef HAVE_VIENNA20
  double min_en = Vienna::alifold((const char**)seqs, &res[0]);  
#else
  double min_en = Vienna::alifold(seqs, &res[0]);
#endif
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/L);
  Vienna::free_alifold_arrays();

#ifdef HAVE_VIENNA18
  Vienna::plist* pi;
#else
  Vienna::pair_info* pi;
#endif
#ifdef HAVE_VIENNA20
  Vienna::alipf_fold((const char**)seqs, NULL, &pi);
#else
  Vienna::alipf_fold(seqs, NULL, &pi);
#endif
  for (uint k=0; pi[k].i!=0; ++k)
    if (pi[k].p>=th && pi[k].i<pi[k].j) 
    {
      bp[pi[k].i].emplace_back(pi[k].j, pi[k].p);
      bp[pi[k].j].emplace_back(pi[k].i, pi[k].p);
    }
  free(pi);

  Vienna::free_alipf_arrays();
  free_aln(seqs);

  return bp;
}

#else

AlifoldModel::
AlifoldModel(const char* param)
{
  throw "ViennaRNA package is not linked.";
}

void
AlifoldModel::
calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                    std::vector<float>& bp, std::vector<int>& offset) const
{
  throw "ViennaRNA package is not linked.";
}

auto
AlifoldModel::
calculate_posterior(const std::list<std::string>& aln, const std::string& paren, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  throw "ViennaRNA package is not linked.";
  uint L=aln.front().size();
  std::vector<std::vector<std::pair<uint, float>>> bp(L+1);
  return bp;
}

void
AlifoldModel::
calculate_posterior(const std::list<std::string>& aln,
                    std::vector<float>& bp, std::vector<int>& offset) const
{
  throw "ViennaRNA package is not linked.";
}

auto
AlifoldModel::
calculate_posterior(const std::list<std::string>& aln, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  throw "ViennaRNA package is not linked.";
  //uint N=aln.size();
  uint L=aln.front().size();
  std::vector<std::vector<std::pair<uint, float>>> bp(L+1);
  return bp;
}
#endif

// Averaged model
void
AveragedModel::
calculate_posterior(const std::list<std::string>& aln,
                    std::vector<float>& bp, std::vector<int>& offset) const
{
  throw "ViennaRNA package is not linked.";
  uint N=aln.size();
  uint L=aln.front().size();
  bp.resize((L+1)*(L+2)/2, 0.0);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
  for (std::list<std::string>::const_iterator s=aln.begin(); s!=aln.end(); ++s)
  {
    std::vector<float> lbp;
    std::vector<int> loffset;
    std::string seq;
    std::vector<int> idx;
    for (uint i=0; i!=s->size(); ++i)
    {
      if ((*s)[i]!='-')
      {
        seq.push_back((*s)[i]);
        idx.push_back(i);
      }
    }
    en_->calculate_posterior(seq, lbp, loffset);
    for (uint i=0; i!=seq.size()-1; ++i)
      for (uint j=i+1; j!=seq.size(); ++j)
        bp[offset[idx[i]+1]+(idx[j]+1)] += lbp[loffset[i+1]+(j+1)]/N;
  }
}

auto
AveragedModel::
calculate_posterior(const std::list<std::string>& aln, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  uint N=aln.size();
  uint L=aln.front().size();
  std::vector<std::vector<std::pair<uint, float>>> bp(L+1);

  for (const auto& s: aln)
  {
    std::string seq;
    std::vector<int> idx;
    for (uint i=0; i!=s.size(); ++i)
    {
      if (s[i]!='-')
      {
        seq.push_back(s[i]);
        idx.push_back(i);
      }
    }
    auto lbp = en_->calculate_posterior(seq);
    for (auto i=1; i!=lbp.size(); ++i)
      for (const auto [j, p]: lbp[i])
      {
        auto res = std::find_if(std::begin(bp[idx[i-1]+1]), std::end(bp[idx[i-1]+1]),
                        [&, &j=j](const auto& x) { return x.first==idx[j-1]+1; });
        if (res != std::end(bp[idx[i-1]+1]))
          res->second += p/N;
        else
          bp[idx[i-1]+1].emplace_back(idx[j-1]+1, p/N);
      }
  }

  return bp;
}

static
std::vector<int>
bpseq(const std::string& paren)
{
  std::vector<int> ret(paren.size(), -1);
  std::stack<int> st;
  for (uint i=0; i!=paren.size(); ++i)
  {
    if (paren[i]=='(')
    {
      st.push(i);
    }
    else if (paren[i]==')')
    {
      ret[st.top()]=i;
      ret[i]=st.top();
      st.pop();
    }
  }
  return ret;
}

void
AveragedModel::
calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                    std::vector<float>& bp, std::vector<int>& offset) const
{
  uint N=aln.size();
  uint L=aln.front().size();
  bp.resize((L+1)*(L+2)/2, 0.0);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
  std::vector<int> p = bpseq(paren);
  for (std::list<std::string>::const_iterator s=aln.begin(); s!=aln.end(); ++s)
  {
    std::vector<float> lbp;
    std::vector<int> loffset;
    std::string seq;
    std::vector<int> idx;
    std::vector<int> rev(s->size(), -1);
    for (uint i=0; i!=s->size(); ++i)
    {
      if ((*s)[i]!='-')
      {
        seq.push_back((*s)[i]);
        idx.push_back(i);
        rev[i]=seq.size()-1;
      }
    }
    std::string lparen(seq.size(), '.');
    for (uint i=0; i!=p.size(); ++i)
    {
      if (rev[i]>=0)
      {
        if (p[i]<0 || rev[p[i]]>=0)
          lparen[rev[i]] = paren[i];
        else
          lparen[rev[i]] = '.';
      }
    }
      
    en_->calculate_posterior(seq, lparen, lbp, loffset);
    for (uint i=0; i!=seq.size()-1; ++i)
      for (uint j=i+1; j!=seq.size(); ++j)
        bp[offset[idx[i]+1]+(idx[j]+1)] += lbp[loffset[i+1]+(j+1)]/N;
  }
}

auto
AveragedModel::
calculate_posterior(const std::list<std::string>& aln, const std::string& paren, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  uint N=aln.size();
  uint L=aln.front().size();
  std::vector<std::vector<std::pair<uint, float>>> bp(L+1);

  std::vector<int> p = bpseq(paren);
  for (const auto& s: aln)
  {
    std::string seq;
    std::vector<int> idx;
    std::vector<int> rev(s.size(), -1);
    for (uint i=0; i!=s.size(); ++i)
    {
      if (s[i]!='-')
      {
        seq.push_back(s[i]);
        idx.push_back(i);
        rev[i]=seq.size()-1;
      }
    }
    std::string lparen(seq.size(), '.');
    for (uint i=0; i!=p.size(); ++i)
    {
      if (rev[i]>=0)
      {
        if (p[i]<0 || rev[p[i]]>=0)
          lparen[rev[i]] = paren[i];
        else
          lparen[rev[i]] = '.';
      }
    }
      
    auto lbp = en_->calculate_posterior(seq, lparen);
    for (auto i=1; i!=lbp.size(); ++i)
      for (const auto [j, p]: lbp[i])
      {
        auto res = std::find_if(std::begin(bp[idx[i-1]+1]), std::end(bp[idx[i-1]+1]),
                        [&, &j=j](const auto& x) { return x.first==idx[j-1]+1; });
        if (res != std::end(bp[idx[i-1]+1]))
          res->second += p/N;
        else
          bp[idx[i-1]+1].emplace_back(idx[j-1]+1, p/N);
      }
  }

  return bp;
}

// MixtureModel
void
MixtureModel::
calculate_posterior(const std::list<std::string>& aln,
                    std::vector<float>& bp, std::vector<int>& offset) const
{
  uint L=aln.front().size();
  bp.resize((L+1)*(L+2)/2);
  std::fill(bp.begin(), bp.end(), 0.0);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
  assert(en_.size()==w_.size());
  for (uint i=0; i!=en_.size(); ++i)
  {
    std::vector<float> lbp;
    std::vector<int> loffset;
    en_[i]->calculate_posterior(aln, lbp, loffset);
    for (uint j=0; j!=bp.size(); ++j)
      bp[j] += lbp[j]*w_[i];
  }
}

auto
MixtureModel::
calculate_posterior(const std::list<std::string>& aln, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  uint L=aln.front().size();
  std::vector<std::vector<std::pair<uint, float>>> bp(L+1);
  assert(en_.size()==w_.size());
  for (uint k=0; k!=en_.size(); ++k)
  {
    auto lbp = en_[k]->calculate_posterior(aln);
    for (auto i=1; i!=lbp.size(); ++i)
    {
      for (const auto [j, p]: lbp[i])
      {
        auto res = std::find_if(std::begin(bp[i]), std::end(bp[i]),
                          [&, &j=j](const auto& x) { return x.first==j; });
        if (res != std::end(bp[i]))
          res->second += p * w_[k];
        else
          bp[i].emplace_back(j, p * w_[k]);
      }
    }
  }

  return bp;
}

void
MixtureModel::
calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                    std::vector<float>& bp, std::vector<int>& offset) const
{
  uint L=aln.front().size();
  bp.resize((L+1)*(L+2)/2);
  std::fill(bp.begin(), bp.end(), 0.0);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
  assert(en_.size()==w_.size());
  for (uint i=0; i!=en_.size(); ++i)
  {
    std::vector<float> lbp;
    std::vector<int> loffset;
    en_[i]->calculate_posterior(aln, paren, lbp, loffset);
    for (uint j=0; j!=bp.size(); ++j)
      bp[j] += lbp[j]*w_[i];
  }
}

auto
MixtureModel::
calculate_posterior(const std::list<std::string>& aln, const std::string& paren, float th) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  uint L=aln.front().size();
  std::vector<std::vector<std::pair<uint, float>>> bp(L+1);
  assert(en_.size()==w_.size());
  for (uint k=0; k!=en_.size(); ++k)
  {
    auto lbp = en_[k]->calculate_posterior(aln, paren);
    for (auto i=1; i!=lbp.size(); ++i)
    {
      for (const auto [j, p]: lbp[i])
      {
        auto res = std::find_if(std::begin(bp[i]), std::end(bp[i]),
                          [&, &j=j](const auto& x) { return x.first==j; });
        if (res != std::end(bp[i]))
          res->second += p * w_[k];
        else
          bp[i].emplace_back(j, p * w_[k]);
      }
    }
  }

  return bp;
}

// Aux model
bool
AuxModel::
calculate_posterior(const char* filename, std::string& seq,
                    std::vector<float>& bp, std::vector<int>& offset) const
{
  std::string l;
  uint L=0;
  std::ifstream in(filename);
  if (!in) return false;
  while (std::getline(in, l)) ++L;
  bp.resize((L+1)*(L+2)/2, 0.0);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
  seq.resize(L);
  in.clear();
  in.seekg(0, std::ios::beg);
  while (std::getline(in, l))
  {
    std::vector<std::string> v;
    std::istringstream ss(l);
    std::string s;
    while (ss >> s) v.push_back(s);
    uint up = atoi(v[0].c_str());
    seq[up-1] = v[1][0];
    for (uint i=2; i!=v.size(); ++i)
    {
      uint down;
      float p;
      if (sscanf(v[i].c_str(), "%u:%f", &down, &p)==2)
        bp[offset[up]+down]=p;
    }
  }
  return true;
}

auto
AuxModel::
calculate_posterior(const char* filename, std::string& seq) const
  -> std::vector<std::vector<std::pair<uint, float>>>
{
  std::string l;
  uint L=0;
  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);    
  while (std::getline(in, l)) ++L;
  std::vector<std::vector<std::pair<uint, float>>> bp(L+1);
  seq.resize(L);
  in.clear();
  in.seekg(0, std::ios::beg);
  while (std::getline(in, l))
  {
    std::vector<std::string> v;
    std::istringstream ss(l);
    std::string s;
    while (ss >> s) v.push_back(s);
    uint up = atoi(v[0].c_str());
    seq[up-1] = v[1][0];
    for (uint i=2; i!=v.size(); ++i)
    {
      uint down;
      float p;
      if (sscanf(v[i].c_str(), "%u:%f", &down, &p)==2)
        bp[up].emplace_back(down, p);
    }
  }
  return bp;
}