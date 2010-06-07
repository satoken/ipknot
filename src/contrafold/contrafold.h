/*
 * wrapper routines for CONTRAfold
 *
 * Copyright (C) 2008-2010 Kengo Sato
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifndef __INC_CONTRAFOLD_H__
#define __INC_CONTRAFOLD_H__

#include <string>

template < class T >
class CONTRAfold
{
private:
  struct Impl;

public:
  CONTRAfold(bool canonical_only = true, int max_bp_dist = 0);
  ~CONTRAfold();

  void SetParameters(const std::string& params);
  void SetConstraint(const std::string& paren);
  
  const T* ComputePosterior(const std::string& seq);
  const T* ComputePosterior(const std::string& seq, std::vector<T>& p);
  T ComputeInside(const std::string& seq);
  T ComputeLogPartitionCoefficient() const;
  T ComputeViterbi(const std::string& seq);

  void PrepareStochasticTraceback(const std::string& seq);
  std::vector<int> StochasticTraceback() const;

  int max_bp_dist() const { return max_bp_dist_; }

  void init_rand(unsigned long seed);

private:
  Impl* impl_;
  int max_bp_dist_;
};

template < class T >
class CONTRAfoldM
{
private:
  struct Impl;

public:
  CONTRAfoldM(bool canonical_only = true, int max_bp_dist = 0);
  ~CONTRAfoldM();

  void SetParameters(const std::string& params);
  //void SetConstraint(const std::string& paren);
  
  //const T* ComputePosterior(const std::string& seq);
  //const T* ComputePosterior(const std::string& seq, std::vector<T>& p);
  //T ComputeInside(const std::string& seq);
  //T ComputeLogPartitionCoefficient() const;
  //T ComputeViterbi(const std::string& seq);

  void PrepareStochasticTraceback(const std::vector<std::string>& aln);
  std::vector<int> StochasticTraceback() const;

  int max_bp_dist() const { return max_bp_dist_; }

  void init_rand(unsigned long seed);
  
private:
  Impl* impl_;
  int max_bp_dist_;
};

#endif	// __INC_CONTRAFOLD_H__

// Local Variables:
// mode: C++
// End:
