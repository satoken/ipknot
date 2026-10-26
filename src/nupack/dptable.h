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

#ifndef __INC_DP_TABLE_H__
#define __INC_DP_TABLE_H__

#include <cassert>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <limits>
#include <vector>

template < class T >
class DPTable2
{
public:
  DPTable2() : V_(), N_(0)
  {
  }

  int size() const { return N_; }

  void resize(int n)
  {
    N_=n;
    V_.resize(std::size_t(N_)*(N_+1)/2+(N_+1));
  }

  void fill(const T& v)
  {
    std::fill(V_.begin(), V_.end(), v);
  }

  T& operator()(int i, int j)
  {
    return V_[index(i, j)];
  }

  const T& operator()(int i, int j) const
  {
    return V_[index(i, j)];
  }

private:
  std::size_t index(int i, int j) const
  {
    //assert(i<=j);
    assert(j<=N_);

    return j==i-1 ? std::size_t(N_)*(N_+1)/2 + i :
      std::size_t(i)*N_+j-std::size_t(i)*(1+i)/2;
  }

private:
  std::vector<T> V_;
  int N_;
};

template < class T >
class DPTable4
{
public:
  DPTable4() : V_(), N_(0)
  {
  }

  void resize(int n)
  {
    if (n < 0) throw std::invalid_argument("negative DP table size");
    if (n == N_ && !V_.empty()) return;
    N_=n;
    const auto count = choose(n, 4);
    // The degenerate one-pair state has its own slot, even for n < 4.
    V_.resize(static_cast<std::size_t>(count)+1);
    left_.resize(n);
    inner_left_.resize(n);
    inner_right_.resize(n);
    for (int i=0; i<n; ++i)
    {
      left_[i] = count-choose(n-i,4)+choose(n-i-1,3);
      inner_left_[i] = -choose(n-i,3)+choose(n-i-1,2);
      inner_right_[i] = -choose(n-i,2)-i-1;
    }
  }

  void fill(const T& v)
  {
    std::fill(V_.begin(), V_.end(), v);
  }

  T& operator()(int i, int d, int e, int j)
  {
    return V_[index(i, d, e, j)];
  }

  const T& operator()(int i, int d, int e, int j) const
  {
    return V_[index(i, d, e, j)];
  }

private:
  static std::int64_t choose(int n, int k)
  {
    if (n < k) return 0;
    std::int64_t result=1;
    for (int i=1; i<=k; ++i)
    {
      if (result>std::numeric_limits<std::int64_t>::max()/(n-i+1))
        throw std::length_error("NUPACK DP table is too large");
      result=result*(n-i+1)/i;
    }
    return result;
  }

  std::size_t index(int h, int r, int m, int s) const
  {
    assert(h>=0);
    assert(h<=r);
    assert(r<=m);
    assert(m<=s);
    assert(s<N_);
    if (h==r && m==s) return V_.size()-1;
    assert(h<r && r<m && m<s);
    return static_cast<std::size_t>(left_[h]+inner_left_[r]+inner_right_[m]+s);
  }

private:
  std::vector<T> V_;
  int N_;
  std::vector<std::int64_t> left_, inner_left_, inner_right_;
};

template < class T >
class DPTableX
{
public:
  DPTableX() : V_(), N_(0), D_(0)
  {
  }

  void resize(int d, int n)
  {
    N_=n;
    D_=d;
    std::size_t max_sz=0;
    for (int i=d; i<d+3; ++i)
      if (i>5 && i<N_)
        max_sz = std::max(max_sz, std::size_t(N_-i)*(i-5)*(i-1)*(i-2)/2);
    V_.resize(max_sz);
  }

  void fill(const T& v)
  {
    std::fill(V_.begin(), V_.end(), v);
  }

  T& operator()(int i, int d, int e, int s)
  {
    return V_[index(i, d, e, s)];
  }

  const T& operator()(int i, int d, int e, int s) const
  {
    return V_[index(i, d, e, s)];
  }

  void swap(DPTableX& x)
  {
    std::swap(V_, x.V_);
    std::swap(N_, x.N_);
    std::swap(D_, x.D_);
  }

private:
  std::size_t index(int i, int h1, int m1, int s) const
  {
    int d=D_;
    std::size_t d1d2 = std::size_t(d-1)*(d-2);
    int d5 = d-5;
    int h1_i_1 = h1-i-1;
    assert(i+d<N_);
    assert(d-6>=s);
    assert(i<h1);
    return std::size_t(i)*d5*d1d2/2 + s*d1d2/2 +
      std::size_t(h1_i_1)*(d-1) - std::size_t(h1_i_1)*(h1-i)/2 + m1 - h1 - 1;
  }

private:
  std::vector<T> V_;
  int N_;
  int D_;
};

#endif  //  __INC_DP_TABLE_H__

// Local Variables:
// mode: C++
// End:
