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

#ifndef __INC_FOLD_H__
#define __INC_FOLD_H__

#include <string>
#include <vector>
#include <list>
#include <stdexcept>

// The base class for calculating base-pairing probabilities of an indivisual sequence
class BPEngineSeq
{
public:
  BPEngineSeq() { }
  virtual ~BPEngineSeq() {}

  virtual void calculate_posterior(const std::string& seq,
                                   std::vector<float>& bp, std::vector<int>& offset) const = 0;

  virtual void calculate_posterior(const std::string& seq, const std::string& paren,
                                   std::vector<float>& bp, std::vector<int>& offset) const = 0;
};

// The base class for calculating base-pairing probabilities of aligned sequences
class BPEngineAln
{
public:
  BPEngineAln() { }
  virtual ~BPEngineAln() {}

  virtual void calculate_posterior(const std::list<std::string>& aln,
                                   std::vector<float>& bp, std::vector<int>& offset) const = 0;
  virtual void calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                                   std::vector<float>& bp, std::vector<int>& offset) const = 0;
};

class CONTRAfoldModel : public BPEngineSeq
{
public:
  CONTRAfoldModel() : BPEngineSeq() { }

  void calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const;

  void calculate_posterior(const std::string& seq, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const;
};

class RNAfoldModel : public BPEngineSeq
{
public:
  RNAfoldModel(const char* param);
  
  void calculate_posterior(const std::string& seq, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const;
    
  void calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const;
};

class NupackModel : public BPEngineSeq
{
public:
  //NupackModel(int model) : BPEngineSeq(), model_(model), param_(NULL) { }
  NupackModel(const char* param) : BPEngineSeq(), param_(param) { }
  
  void calculate_posterior(const std::string& seq, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const
  {
    throw std::runtime_error("not supported yet");
  }
    
  void calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const;

private:
  //int model_;
  const char* param_;
};

class AlifoldModel : public BPEngineAln
{
public:
  AlifoldModel(const char* param);

  void calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const;

  void calculate_posterior(const std::list<std::string>& aln,
                           std::vector<float>& bp, std::vector<int>& offset) const;
};

class AveragedModel : public BPEngineAln
{
public:
  AveragedModel(BPEngineSeq* en) : en_(en) { }

  void calculate_posterior(const std::list<std::string>& aln,
                           std::vector<float>& bp, std::vector<int>& offset) const;

  void calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const;

private:
  BPEngineSeq* en_;
};

class MixtureModel : public BPEngineAln
{
public:
  MixtureModel(std::vector<BPEngineAln*> en)
    : en_(en),
      w_(en_.size(), 1.0/en_.size())
  { }

  void calculate_posterior(const std::list<std::string>& aln,
                           std::vector<float>& bp, std::vector<int>& offset) const;

  void calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const;

private:
  std::vector<BPEngineAln*> en_;
  std::vector<float> w_;
};

class AuxModel
{
public:
  AuxModel() { }

  bool calculate_posterior(const char* filename, std::string& seq,
                           std::vector<float>& bp, std::vector<int>& offset) const;
};

#endif  // __INC_FOLD_H__

// Local Variables:
// mode: C++
// End:
