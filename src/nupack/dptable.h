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

template < class T >
class DPTable2
{
public:
  DPTable2() : V_(), N_(0)
  {
  }

  void resize(int n)
  {
    N_=n;
    V_.resize(N_*(N_+1)/2+(N_+1));
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
  int index(int i, int j) const
  {
    //assert(i<=j);
    assert(j<=N_);

    return j==i-1 ? N_*(N_+1)/2 + i : i*N_+j-i*(1+i)/2;
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
    N_=n;
    V_.resize(N_*(N_-1)*(N_-2)*(N_-3)/2/3/4);
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
  int index(int h, int r, int m, int s) const
  {
    int n = N_;
    int h2 = h*h;
    int h3 = h2*h;
    int h4 = h3*h;
    int m2 = m*m;
    int n2 = n*n;
    int n3 = n2*n;
    int r2 = r*r;
    int r3 = r2*r;

    assert(h<=r);
    assert(r<=m);
    assert(m<=s);
    assert(s<=N_);

    return (h==r && m==s) ? V_.size()-1 :
      (-24 - 50*h - 35*h2 - 10*h3 - h4 - 36*m -12*m2 +
       12*n + 70*h*n + 30*h2*n + 4*h3*n + 24*m*n - 12*n2 -30*h*n2 -
       6*h2*n2 + 4*h*n3 + 44*r - 48*n*r + 12*n2*r + 
       24*r2 - 12*n*r2 +  4*r3 + 24*s)/24 ;
  }

private:
  std::vector<T> V_;
  int N_;
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
    int max_sz=0;
    for (int i=d; i<d+3; ++i)
      max_sz = std::max(max_sz, (N_-i)*(i-5)*(i-1)*(i-2)/2);
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
  int index(int i, int h1, int m1, int s) const
  {
    int d=D_;
    int d1d2 = (d-1)*(d-2);
    int d5 = d-5;
    int h1_i_1 = h1-i-1;
    assert(i+d<N_);
    assert(d-6>=s);
    assert(i<h1);
    return i*d5*d1d2/2 + s*d1d2/2 + 
      h1_i_1*(d-1) - h1_i_1*(h1-i)/2 + m1 - h1 - 1;
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
