/*
 * $Id$
 *
 * PHMMTS -- an implementation of Pair Hidden Markov Models on Tree Structures,
 *           which is based on "Pair Hidden Markov Models on Tree Structures",
 *           Yasubumi Sakakibara, ISMB 2003.
 *
 * Copyright (C) 2003-2006 Kengo Sato
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

#ifndef __INC_LOG_VALUE_H__
#define __INC_LOG_VALUE_H__

#define FAST_LOG1EXP0      

#ifdef HAVE_NUMERIC_LIMITS
#include <limits>
#else
#include <climits>
#endif

/*
  Example:

  LogValue<double> a = 2.0;
  LogValue<double> b = 10.0;
  LogValue<double> c = a * b;  

  In this example,
  a's internal expression is log(2.0),
  b's internal expression is log(10.0),
  and c is calculated as log(2.0)+log(10.0)
*/

#include <cmath>
#include <iosfwd>


template <class T>
class LogValue
{
public:
  typedef T value_type;

  class ExpOf {
    value_type val_;
  public:
    explicit ExpOf(value_type val) : val_(val) { }
    value_type log() const { return val_; }
  };

  static value_type log(const ExpOf& v) { return v.log(); }

private:
  T val_;
  
public:
  // default constructor
  LogValue() : val_() {}

  // constructor
  LogValue(const ExpOf& v) : val_(v.log()) { }
    
  // constructor
  LogValue(double v) : val_(my_log(v)) { }

  // constructor
#if 0
  template <class U>
  LogValue(U v) : val_(static_cast<T>(my_log(v))) { }
#endif
  
  // copy constructor
  LogValue(const LogValue& x) : val_(x.val_) { }

  // cast to double
  operator double() const { return exp(val_); }

  // cast to float
  //operator float() const { return exp(val_); }

  // assign U
#if 0
  template <class U>
  LogValue& operator=(const U& v)
  {
    if (v==0.0) val_=-2e20;
    else val_=static_cast<T>(my_log(v));
    return *this;
  }
#endif

  // assign LogValue
  LogValue& operator=(const LogValue& x)
  {
    val_ = x.val_;
    return *this;
  }

  // assign logged value
  void set_log_value(const T& v)
  {
    val_ = v;
  }

  // assign logged value
  const T& log() const
  {
    return val_;
  }

  // multiplied by T
  template <class U>
  LogValue operator*(const U& v) const
  {
    LogValue r(*this);
    r.val_ += static_cast<T>(my_log(v));
    return r;
  }

  // multiplied by LogValue
  LogValue operator*(const LogValue& x) const
  {
    LogValue r(*this);
    r.val_ += x.val_;
    return r;
  }

  // multiplied and assigned by T
  template <class U>
  LogValue& operator*=(const U& v)
  {
    val_ += static_cast<T>(my_log(v));
    return *this;
  }

  // multiplied and assigned by LogValue
  LogValue& operator*=(const LogValue& x)
  {
    val_ += x.val_;
    return *this;
  }

  // divided by T
  template <class U>
  LogValue operator/(const U& v) const
  {
    LogValue r(*this);
    r.val_ -= static_cast<T>(my_log(v));
    return r;
  }

  // divided by LogValue
  LogValue operator/(const LogValue& x) const
  {
    LogValue r(*this);
    r.val_ -= x.val_;
    return r;
  }

  // divided and assigned by T
  template <class U>
  LogValue& operator/=(const U& v)
  {
    val_ -= static_cast<U>(my_log(v));
    return *this;
  }

  // divided and assigned by LogValue
  LogValue& operator/=(const LogValue& x)
  {
    val_ -= x.val_;
    return *this;
  }

  // added by LogValue
  LogValue operator+(const LogValue& x) const
  {
    if (zerop(x)) return *this;
    if (zerop(*this)) return x;
    if (x.val_>val_) {
#if 0
      if (x.val_>val_+30) {
	return x;
      } else {
#endif
	LogValue r;
	r.val_ = val_ + log1exp0(x.val_-val_);
	return r;
#if 0
      }
#endif
    } else {
#if 0
      if (val_>x.val_+30) {
	return *this;
      } else {
#endif
	LogValue r;
	r.val_ = x.val_ + log1exp0(val_-x.val_);
	return r;
#if 0
      }
#endif
    }
  }

  // added by T
  template <class U>
  LogValue operator+(const U& v) const
  {
    return *this + LogValue(v);
  }

  // added and assinged by LogValue
  LogValue& operator+=(const LogValue& x)
  {
    if (zerop(x)) return *this;
    if (zerop(*this)) {
      val_ = x.val_;
      return *this;
    }
    if (x.val_>val_) {
#if 0
      if (x.val_>val_+30) {
	val_ = x.val_;
      } else {
#endif
	val_ = val_ + log1exp0(x.val_-val_);
#if 0
      }
#endif
    } else {
#if 0
      if (val_>x.val_+30) {
	
      } else {
#endif
	val_ = x.val_ + log1exp0(val_-x.val_);
#if 0
      }
#endif
    }
    return *this;
  }

  // added and assinged by T
  template <class U>
  LogValue& operator+=(const U& v)
  {
    return operator+=(LogValue(v));
  }

  // is equal to LogValue
  bool operator==(const LogValue& x) const
  {
    return val_ == x.val_;
  }

  // is equal to T
  template <class U>
  bool operator==(const U& v) const
  {
    return val_ == my_log(v);
  }

  // is not equal to LogValue
  bool operator!=(const LogValue& x) const
  {
    return val_ != x.val_;
  }

  // is not equal to T
  template <class U>
  bool operator!=(const U& v) const
  {
    return val_ != my_log(v);
  }

  // is less than LogValue
  bool operator<(const LogValue& x) const
  {
    return val_ < x.val_;
  }

  // is less than T
  template <class U>
  bool operator<(const U& v) const
  {
    return val_ < my_log(v);
  }

  // is equal or less than LogValue
  bool operator<=(const LogValue& x) const
  {
    return val_ <= x.val_;
  }

  // is equal or less than T
  template <class U>
  bool operator<=(const U& v) const
  {
    return val_ <= my_log(v);
  }

  // is greater than LogValue
  bool operator>(const LogValue& x) const
  {
    return val_ > x.val_;
  }

  // is greater than T
  template <class U>
  bool operator>(const U& v) const
  {
    return val_ > my_log(v);
  }

  // is equal or greater than LogValue
  bool operator>=(const LogValue& x) const
  {
    return val_ >= x.val_;
  }

  // is equal or greater than T
  template <class U>
  bool operator>=(const U& v) const
  {
    return val_ >= my_log(v);
  }

private:
  T log1exp0(double x) const
  {
#ifdef FAST_LOG1EXP0
    // These values are taken from probcons_v1_10
    if (x > 10.0f) return x;
    if (x >= 0.00f) {
      if (x <= 1.00f)
	return static_cast<T>(((-0.009350833524763f * x
				+ 0.130659527668286f) * x
			       + 0.498799810682272f) * x
			      + 0.693203116424741f);
      if (x <= 2.50f)
	return static_cast<T>(((-0.014532321752540f * x
				+ 0.139942324101744f) * x
			       + 0.495635523139337f) * x
			      + 0.692140569840976f);
      if (x <= 4.50f)
	return static_cast<T>(((-0.004605031767994f * x
				+ 0.063427417320019f) * x
			       + 0.695956496475118f) * x
			      + 0.514272634594009f);
      if (x <= 7.50f)
	return static_cast<T>(((-0.000458661602210f * x
				+ 0.009695946122598f) * x
			       + 0.930734667215156f) * x
			      + 0.168037164329057f);
      //if (x <= 10.00f)
      return static_cast<T>((((0.00000051726300753785 * x
			       - 0.00002720671238876090) * x
			      + 0.00053403733818413500) * x
			     + 0.99536021775747900000) * x
			    + 0.01507065715532010000);
    }
#endif
    return static_cast<T>(::log(exp(x)+1.0));
  }

  static T my_log(T v)
  {
    return v==0.0 ? -2e20 : static_cast<T>(::log(v));
  }
};

template <class T>
T zero()
{
  return static_cast<T>(0.0);
}

template < class T >
inline LogValue<T> zero()
{
  LogValue<T> r;
  r.set_log_value(-2e20);
  return r;
}

#if 0
template < >
inline LogValue<short int> zero()
{
  LogValue<short int> r;
#ifdef HAVE_NUMERIC_LIMITS 
  r.set_log_value(std::numeric_limits<short int>::min());
#else
  r.set_log_value(SHRT_MIN);
#endif
  return r;
}
#endif

template < class T >
bool zerop(const T& x)
{
  return x==0.0;
}

template < class T >
bool zerop(const LogValue<T>& x)
{
  //return std::isinf(x.log())<0;
  return x.log()<=-2e20;
}

template < >
inline bool zerop(const LogValue<short int>& x)
{
#ifdef HAVE_NUMERIC_LIMITS 
  return x.log()==std::numeric_limits<short int>::min();
#else
  return x.log()==SHRT_MIN;
#endif
}


// return logged value
template <class T>
T std::log(const LogValue<T>& x)
{
  return x.log();
}

template <class T>
LogValue<T> std::pow(const LogValue<T>& x, double p)
{
  LogValue<T> r;
  r.set_log_value(static_cast<T>(log(x) * p));
  return r;
}

// input from a stream
template <class T>
std::istream& operator>>(std::istream& in, LogValue<T>& x)
{
  T t;
  in >> t;
  x=t;
  return in;
}

#endif // __INC_LOG_VALUE_H__

// Local Variables:
// mode: C++
// End:
