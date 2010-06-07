// $Id:$

#ifndef __INC_RAND_H__
#define __INC_RAND_H__

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

#ifdef STD_RAND
#include <cstdlib>
class Die
{
public:
  Die(unsigned int seed)
  {
    srand(seed);
  }

  double operator()()
  {
    return rand()/(RAND_MAX+1.0);
  }
};

#else
#include <boost/random.hpp>

class Die
{
  typedef boost::mt19937 RNG;
  typedef boost::uniform_real<> DST;
  typedef boost::variate_generator<RNG&, DST> DIE;

public:
  Die(unsigned int seed)
    : rng_(seed), dst_(0, 1), die_(rng_, dst_)
  { }

  double operator()() { return die_(); }

private:
  RNG rng_;
  DST dst_;
  DIE die_;
};
#endif

#endif  //  __INC_RAND_H__
