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
#include <memory>
#include <stdexcept>

// The base class for calculating base-pairing probabilities of an indivisual sequence
class BPEngineSeq
{
public:
  BPEngineSeq() { }
  virtual ~BPEngineSeq() {}

  virtual void calculate_posterior(const std::string& seq,
                                   std::vector<float>& bp, std::vector<int>& offset) const = 0;
  virtual auto calculate_posterior(const std::string& seq, float th=0.0) const 
    -> std::vector<std::vector<std::pair<uint, float>>>
  {
    throw std::runtime_error("not supported yet");
  }

  virtual void calculate_posterior(const std::string& seq, const std::string& paren,
                                   std::vector<float>& bp, std::vector<int>& offset) const = 0;
  virtual auto calculate_posterior(const std::string& seq, const std::string& paren, float th=0.0) const
    -> std::vector<std::vector<std::pair<uint, float>>>
  {
    throw std::runtime_error("not supported yet");
  }
};

// The base class for calculating base-pairing probabilities of aligned sequences
class BPEngineAln
{
public:
  BPEngineAln() { }
  virtual ~BPEngineAln() {}

  virtual void calculate_posterior(const std::list<std::string>& aln,
                                   std::vector<float>& bp, std::vector<int>& offset) const = 0;
  virtual auto calculate_posterior(const std::list<std::string>& aln, float th=0.0) const 
    -> std::vector<std::vector<std::pair<uint, float>>>
  {
    throw std::runtime_error("not supported yet");
  }

  virtual void calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                                   std::vector<float>& bp, std::vector<int>& offset) const = 0;
  virtual auto calculate_posterior(const std::list<std::string>& aln, const std::string& paren, float th=0.0) const
    -> std::vector<std::vector<std::pair<uint, float>>>
  {
    throw std::runtime_error("not supported yet");
  }
};

class CONTRAfoldModel : public BPEngineSeq
{
public:
  CONTRAfoldModel() : BPEngineSeq() { }

  void calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::string& seq, float th=0.0) const 
    -> std::vector<std::vector<std::pair<uint, float>>>;

  void calculate_posterior(const std::string& seq, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::string& seq, const std::string& paren, float th=0.0) const
    -> std::vector<std::vector<std::pair<uint, float>>>;
};

class RNAfoldModel : public BPEngineSeq
{
public:
  RNAfoldModel(const char* param);
  
  void calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::string& seq, float th=0.0) const 
    -> std::vector<std::vector<std::pair<uint, float>>>;
    
  void calculate_posterior(const std::string& seq, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::string& seq, const std::string& paren, float th=0.0) const
    -> std::vector<std::vector<std::pair<uint, float>>>;
};

class NupackModel : public BPEngineSeq
{
public:
  //NupackModel(int model) : BPEngineSeq(), model_(model), param_(NULL) { }
  NupackModel(const char* param) : BPEngineSeq(), param_(param) { }
  
  void calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::string& seq, float th=0.0) const 
    -> std::vector<std::vector<std::pair<uint, float>>>;

  void calculate_posterior(const std::string& seq, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const
  {
    throw std::runtime_error("not supported yet");
  }
//  auto calculate_posterior(const std::string& seq, const std::string& paren, float th=0.0) const
//    -> std::vector<std::vector<std::pair<uint, float>>>;
    
private:
  //int model_;
  const char* param_;
};

class LinearPartitionModel : public BPEngineSeq
{
public:
  LinearPartitionModel(bool use_vienna=false) : BPEngineSeq(), use_vienna_(use_vienna) { }
  
  void calculate_posterior(const std::string& seq, std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::string& seq, float th=0.0) const -> std::vector<std::vector<std::pair<uint, float>>>;
  void calculate_posterior(const std::string& seq, const std::string& paren, std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::string& seq, const std::string& paren, float th=0.0) const -> std::vector<std::vector<std::pair<uint, float>>>;

private:
  bool use_vienna_;
};

class AlifoldModel : public BPEngineAln
{
public:
  AlifoldModel(const char* param);

  void calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::list<std::string>& aln, const std::string& paren, float th=0.0) const
    -> std::vector<std::vector<std::pair<uint, float>>>;

  void calculate_posterior(const std::list<std::string>& aln,
                           std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::list<std::string>& aln, float th=0.0) const 
    -> std::vector<std::vector<std::pair<uint, float>>>;
};

class AveragedModel : public BPEngineAln
{
public:
  AveragedModel(std::unique_ptr<BPEngineSeq>&& en) : en_(std::move(en)) { }

  void calculate_posterior(const std::list<std::string>& aln,
                           std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::list<std::string>& aln, float th=0.0) const 
    -> std::vector<std::vector<std::pair<uint, float>>>;

  void calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::list<std::string>& aln, const std::string& paren, float th=0.0) const
    -> std::vector<std::vector<std::pair<uint, float>>>;

private:
  std::unique_ptr<BPEngineSeq> en_;
};

class MixtureModel : public BPEngineAln
{
public:
  MixtureModel(std::vector<std::unique_ptr<BPEngineAln>>&& en)
    : en_(std::move(en)),
      w_(en_.size(), 1.0/en_.size())
  { }

  void calculate_posterior(const std::list<std::string>& aln,
                           std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::list<std::string>& aln, float th=0.0) const 
    -> std::vector<std::vector<std::pair<uint, float>>>;

  void calculate_posterior(const std::list<std::string>& aln, const std::string& paren,
                           std::vector<float>& bp, std::vector<int>& offset) const;
  auto calculate_posterior(const std::list<std::string>& aln, const std::string& paren, float th=0.0) const
    -> std::vector<std::vector<std::pair<uint, float>>>;

private:
  std::vector<std::unique_ptr<BPEngineAln>> en_;
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
