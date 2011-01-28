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

#ifndef __INC_NUPACK_H__
#define __INC_NUPACK_H__

#include <vector>
#include <string>
#include <boost/multi_array.hpp>
#include "dptable.h"

#define kB 0.00198717 // Boltzmann constant in kcal/mol/K
#define ZERO_C_IN_KELVIN 273.15 // Zero degrees C in Kelvin
#define AVOGADRO 6.022e23 // Avogadro's number

template <class T> class DPTable2;
template <class T> class DPTable4;
template <class T> class DPTableX;

template < class PF_TYPE >
class Nupack
{
public:
  //enum { DP03, DP09 };
  typedef PF_TYPE pf_type;
  typedef double DBL_TYPE;
  typedef float energy_t;

public:
  Nupack();
  void load_sequence(const std::string& s);
  void load_parameters_fm363(const std::vector<float>& v);
  void load_default_parameters(/*int which*/);
  bool load_parameters(const char* filename);
  void dump_parameters(std::ostream& os) const;
  pf_type calculate_partition_function();
  void calculate_posterior();
  void get_posterior(std::vector<float>& bp, std::vector<int>& offset) const;
  void get_posterior(std::vector<float>& bp1, std::vector<float>& bp2, std::vector<int>& offset) const;
  
private:
  void fastiloops(int i, int j, DPTable4<PF_TYPE>& Qg, DPTableX<PF_TYPE>& Qx, DPTableX<PF_TYPE>& Qx2);
  void fastiloops_pr(int i, int j,
#ifdef ENABLE_RECALCULATE
                     DPTableX<float>& Precx, DPTableX<float>& Precx2,
#endif
                     DPTable4<PF_TYPE>& Qg, DPTableX<PF_TYPE>& Qx, DPTableX<PF_TYPE>& Qx2,
                     DPTable4<DBL_TYPE>& Pg, DPTableX<DBL_TYPE>& Px, DPTableX<DBL_TYPE>& Px2);

  energy_t score_hairpin(int i, int j) const;
  energy_t score_loop(int l) const;
  energy_t score_interior(int i, int d, int e, int j, bool pk) const;
  energy_t score_interior_mismatch(int i, int j) const;
  energy_t score_interior_mismatch(int i, int j, int k, int l) const;
  energy_t score_interior_asymmetry(int l1, int l2) const;
  energy_t score_multiloop(bool pk) const;
  energy_t score_multiloop_paired(int n, bool pk) const;
  energy_t score_multiloop_unpaired(int n, bool pk) const;
  energy_t score_at_penalty(int i, int j) const;
  energy_t score_dangle(int i, int j) const;
  energy_t score_pk() const;  
  energy_t score_pk_multiloop() const;
  energy_t score_pk_pk() const;  
  energy_t score_pk_paired(int n) const;
  energy_t score_pk_unpaired(int n) const;
  energy_t score_pk_band(int n) const;

  int base(char x) const;
  bool allow_paired(int i, int j) const;
  bool wc_pair(int i, int j) const;
  int pair_type(int i, int j) const;
  int pair_type(int i) const;

private:
  std::vector<int> base_map;
  boost::multi_array<int,2> pair_map;
  std::vector<int> seq;
  int N;
  float RT;
  DPTable2<PF_TYPE> Q;
  DPTable2<PF_TYPE> Qb;
  DPTable2<PF_TYPE> Qm;
  DPTable2<PF_TYPE> Qp;
  DPTable2<PF_TYPE> Qz;
  DPTable4<PF_TYPE> Qg;
  DPTable4<PF_TYPE> Qgl;
  DPTable4<PF_TYPE> Qgr;
  DPTable4<PF_TYPE> Qgls;
  DPTable4<PF_TYPE> Qgrs;
  DPTable2<DBL_TYPE> P;
  DPTable2<DBL_TYPE> Pb;
  DPTable2<DBL_TYPE> Pm;
  DPTable2<DBL_TYPE> Pp;
  DPTable2<DBL_TYPE> Pz;
  DPTable2<DBL_TYPE> Pbg;
  DPTable4<DBL_TYPE> Pg;
  DPTable4<DBL_TYPE> Pgl;
  DPTable4<DBL_TYPE> Pgr;
  DPTable4<DBL_TYPE> Pgls;
  DPTable4<DBL_TYPE> Pgrs;
  
  // energy parameters
#if 0
  std::vector<energy_t> hairpin37;
  std::vector<energy_t> bulge37;
  std::vector<energy_t> interior37;
  boost::multi_array<energy_t,2> stack37;
  boost::multi_array<energy_t,4> int11_37;
  boost::multi_array<energy_t,5> int21_37;
  boost::multi_array<energy_t,6> int22_37;
  boost::multi_array<energy_t,2> dangle3_37;
  boost::multi_array<energy_t,2> dangle5_37;
  boost::multi_array<energy_t,5> triloop37;
  boost::multi_array<energy_t,6> tloop37;
  boost::multi_array<energy_t,3> mismatch_hairpin37;
  boost::multi_array<energy_t,3> mismatch_interior37;
  std::vector<energy_t> asymmetry_penalty;
#else
  energy_t hairpin37[30];
  energy_t bulge37[30];
  energy_t interior37[30];
  energy_t stack37[6][6];
  energy_t int11_37[6][6][4][4];
  energy_t int21_37[6][4][4][6][4];
  energy_t int22_37[6][6][4][4][4][4];
  energy_t dangle3_37[6][4];
  energy_t dangle5_37[6][4];
  energy_t triloop37[4][4][4][4][4];
  energy_t tloop37[4][4][4][4][4][4];
  energy_t mismatch_hairpin37[4][4][6];
  energy_t mismatch_interior37[4][4][6];
  energy_t asymmetry_penalty[4];
#endif
  energy_t polyC_penalty, polyC_slope, polyC_int;
  energy_t at_penalty;
  energy_t multiloop_penalty; // alpha1
  energy_t multiloop_paired_penalty; // alpha2
  energy_t multiloop_unpaired_penalty; // alpha3
  energy_t pk_penalty; // beta1
  energy_t pk_multiloop_penalty; // beta1m
  energy_t pk_pk_penalty; // beta1p
  energy_t pk_paired_penalty; // beta2
  energy_t pk_unpaired_penalty; // beta3
  energy_t pk_band_penalty;
  energy_t pk_stack_span;
  energy_t pk_interior_span;
  energy_t multiloop_penalty_pk;
  energy_t multiloop_paired_penalty_pk;
  energy_t multiloop_unpaired_penalty_pk;

  energy_t max_asymmetry;
  energy_t SALT_CORRECTION;
  energy_t loop_greater30;
  energy_t hairpin_GGG;
  float intermolecular_initiation;
};

#endif // __INC_NUPACK_H__

// Local Variables:
// mode: C++
// End:
