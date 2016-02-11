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

#include "config.h"
#include "fold.h"

#include <cstring>
#include <iostream>
#include <sstream>

#include "contrafold/SStruct.hpp"
#include "contrafold/InferenceEngine.hpp"
#include "contrafold/Defaults.ipp"

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

#include "nupack/nupack.h"

extern "C" {
#include "boltzmann_param.h"
};

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

// RNAfold model

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


// Alifold model

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

// Averaged model
void
AveragedModel::
calculate_posterior(const std::list<std::string>& aln,
                    std::vector<float>& bp, std::vector<int>& offset) const
{
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
