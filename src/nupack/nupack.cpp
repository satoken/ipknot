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

#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/array.hpp>
#include "dptable.h"

typedef float energy_t;

#include "nupack.h"

enum { PAIR_AU=0, PAIR_CG, PAIR_GC, PAIR_UA, PAIR_GU, PAIR_UG };
enum { BASE_N=0, BASE_A, BASE_C, BASE_G, BASE_U };
enum { A=BASE_A-1, C=BASE_C-1, G=BASE_G-1, U=BASE_U-1 };
enum { AU=PAIR_AU, CG=PAIR_CG, GC=PAIR_GC, UA=PAIR_UA, GU=PAIR_GU, UG=PAIR_UG };

#define EXP expl
#define LOG logl

template < class PF_TYPE >
Nupack<PF_TYPE>::
Nupack()
  : base_map('z'-'a'+1),
    pair_map(boost::extents[5][5]),
    RT(kB*(ZERO_C_IN_KELVIN+37)),
#if 0
    hairpin37(30),
    bulge37(30),
    interior37(30),
    stack37(boost::extents[6][6]),
    int11_37(boost::extents[6][6][4][4]),
    int21_37(boost::extents[6][4][4][6][4]),
    int22_37(boost::extents[6][6][4][4][4][4]),
    dangle3_37(boost::extents[6][4]),
    dangle5_37(boost::extents[6][4]),
    triloop37(boost::extents[4][4][4][4][4]),
    tloop37(boost::extents[4][4][4][4][4][4]),
    mismatch_hairpin37(boost::extents[4][4][6]),
    mismatch_interior37(boost::extents[4][4][6]),
    asymmetry_penalty(4),
#endif
    SALT_CORRECTION(0),
    loop_greater30(1.079 /*=1.75*RT*/),
    hairpin_GGG(0.0)
{
  std::fill(base_map.begin(), base_map.end(), (int)BASE_N);
  base_map['a'-'a'] = BASE_A;
  base_map['c'-'a'] = BASE_C;
  base_map['g'-'a'] = BASE_G;
  base_map['u'-'a'] = BASE_U;
  base_map['t'-'a'] = BASE_U;

  std::fill(pair_map.data(), pair_map.data()+pair_map.num_elements(), -1);
  pair_map[BASE_A][BASE_U] = PAIR_AU;
  pair_map[BASE_U][BASE_A] = PAIR_UA;
  pair_map[BASE_C][BASE_G] = PAIR_CG;
  pair_map[BASE_G][BASE_C] = PAIR_GC;
  pair_map[BASE_G][BASE_U] = PAIR_GU;
  pair_map[BASE_U][BASE_G] = PAIR_UG;
#if 0
  std::fill(hairpin37.begin(), hairpin37.end(), 0);
  std::fill(bulge37.begin(), bulge37.end(), 0);
  std::fill(interior37.begin(), interior37.end(), 0);
  std::fill(stack37.data(), stack37.data()+stack37.num_elements(), 0);
  std::fill(int11_37.data(), int11_37.data()+int11_37.num_elements(), 0);
  std::fill(int21_37.data(), int21_37.data()+int21_37.num_elements(), 0);
  std::fill(int22_37.data(), int22_37.data()+int22_37.num_elements(), 0);
  std::fill(dangle3_37.data(), dangle3_37.data()+dangle3_37.num_elements(), 0);
  std::fill(dangle5_37.data(), dangle5_37.data()+dangle5_37.num_elements(), 0);
  std::fill(triloop37.data(), triloop37.data()+triloop37.num_elements(), 0);
  std::fill(tloop37.data(), tloop37.data()+tloop37.num_elements(), 0);
  std::fill(mismatch_hairpin37.data(), mismatch_hairpin37.data()+mismatch_hairpin37.num_elements(), 0);
  std::fill(mismatch_interior37.data(), mismatch_interior37.data()+mismatch_interior37.num_elements(), 0);
  std::fill(asymmetry_penalty.begin(), asymmetry_penalty.end(), 0);
#endif
}

template < class PF_TYPE >
int
Nupack<PF_TYPE>::
base(char x) const
{
  if (x>='a' && x<='z')
    return base_map[x-'a'];
  else if (x>='A' && x<='Z')
    return base_map[x-'A'];
  else
    return BASE_N;
}

template < class PF_TYPE >
int
Nupack<PF_TYPE>::
pair_type(int i, int j) const
{
  return pair_map[seq[i]][seq[j]];
}

template < class PF_TYPE >
int
Nupack<PF_TYPE>::
pair_type(int i) const
{
  // assume Watson-Crick pairs
  switch (seq[i])
  {
    case BASE_A: return PAIR_AU; break;
    case BASE_C: return PAIR_CG; break;
    case BASE_G: return PAIR_GC; break;
    case BASE_U: return PAIR_UA; break;
  }
  return -1;
}

template < class PF_TYPE >
bool
Nupack<PF_TYPE>::
wc_pair(int i, int j) const
{
  return pair_type(i, j)!=PAIR_GU && pair_type(i, j)!=PAIR_UG;
}

template < class PF_TYPE >
bool
Nupack<PF_TYPE>::
allow_paired(int i, int j) const
{
  return j-i-1>=3 && pair_type(i,j)>=0;
}

template < class PF_TYPE >
void
Nupack<PF_TYPE>::
load_sequence(const std::string& s)
{
  N=s.size();
  seq.resize(N);
  for (int i=0; i!=N; ++i) seq[i] = base(s[i]);
}

int
check_stability_and_size (int k, int l, int o, int p)
// helper function, to detect which delta we need for the int22 parameters
{
  // having at least one AC mismatch is the simplest, test first
  if ((k==A && l==C) || (k==C && l==A) || (o==A && p==C) || (o==C && p==A))
    return 4;
    
  // combination of all mismatches of equal size (purine-purine, purine-pyrimidine, and pyrimidine-pyrimidine are different sizes)
  // purine =  A, G
  // pyrimidine = C, U
  // if all purine-purines
  if ((k==A || k==G) && (l==A || l==G) && (o==A || o==G) && (p==A || p==G))
    return 1;
  // if all pyrimidine-pyrimidine
  if ((k==C || k==U) && (l==C || l==U) && (o==C || o==U) && (p==C || p==U))
    return 1;
  // if both  purine-pyrimidine
  // assume the A-C pairs have been found above
  if ( (((k==A || k==G) && (l==C || l==U)) || ((k==C || k==U) && (l==A || l==G))) &&
       (((o==A || o==G) && (p==C || p==U)) || ((o==C || o==U) && (p==A || p==G))) )
    return 1;
  // or any combination of 2 unstable mismatches except AC: AA, CC, CU, GG
  if ( ((k==A && l==A) || (k==C && l==C) || (k==C && l==U) || (k==U && l==C) || (k==G && l==G)) &&
       ((o==A && p==A) || (o==C && p==C) || (o==C && p==U) || (o==U && p==C) || (o==G && p==G)) )
    return 1;
                 
  // two stabilizing mismatches (GU, GA, UU) of different sizes  (purine-purine, purine-pyrimidine, and pyrimidine-pyrimidine are different sizes)
  if ( (((k==G && l==U) || (k==U && l==G))   &&   ((o==G && p==A) || (o==A && p==G) || (o==U && p==U))) ||
       (((k==G && l==A) || (k==A && l==G))   &&   ((o==G && p==U) || (o==U && p==G) || (o==U && p==U))) ||
       ((k==U && l==U)                       &&   ((o==G && p==A) || (o==A && p==G) || (o==G && p==U) || (o==U && p==G))) )
    return 2;
        
  // one stable (GU, GA, UU) and one unstable mismatch (excluding AC) (AA, CC, CU, GG) of different sizes
  // GU
  if ( ((k==G && l==U) || (k==U && l==G)) && 
       ((o==A && p==A) || (o==C && p==C) || (o==C && p==U) || (o==U && p==C) || (o==G && p==G)) )
    return 3;    
  if ( ((o==G && p==U) || (o==U && p==G)) &&
       ((k==A && l==A) || (k==C && l==C) || (k==C && l==U) || (k==U && l==C) || (k==G && l==G)) )
    return 3;
  // GA        
  if ( ((k==G && l==A) || (k==A && l==G)) &&
       ((o==C && p==C) || (o==C && p==U) || (o==U && p==C)) )
    return 3;    
  if ( ((o==G && p==A) || (o==A && p==G)) &&
       ((k==C && l==C) || (k==C && l==U) || (k==U && l==C)) )
    return 3;
  // UU        
  if ( (k==U && l==U) &&
       ((o==A && p==A) || (o==G && p==G)) )
    return 3;    
  if ( (o==U && p==U) &&
       ((k==A && l==A) || (k==G && l==G)) )
    return 3;    
        
  return -1;            
}

#include "def_param.h"

#if 0
static float bl[] = {
  -0.70, -1.30, -1.39, -0.14, -0.85, -0.81, -1.32, -2.08, -1.33, -0.38,
  -1.47, -1.23, -2.05, -0.92, -1.51, -0.58, -0.68, -0.23, -0.69, -0.03,
  -0.38, 0.42, 0.77, 0.65, 1.18, -0.03, 0.09, 0.45, 0.43, -0.32,
  0.23, 0.16, 0.60, 0.30, -0.07, 0.32, -0.10, -0.35, -0.14, 0.04,
  0.13, -0.13, -0.17, -0.98, -0.23, -0.51, -0.24, -0.02, 0.03, -0.16,
  -0.58, -0.97, -0.39, 0.02, -0.23, -0.65, 1.10, 0.12, -0.27, -0.18,
  0.00, -0.67, 0.02, -0.46, 0.37, 0.08, -0.40, -0.56, -0.92, 0.60,
  0.76, 0.43, 1.81, 0.51, 0.45, 0.71, -0.02, -0.93, 0.16, 0.40,
  1.14, 0.25, 0.07, 0.67, 0.37, 0.18, 0.52, 0.76, 0.69, 0.54,
  0.37, 0.38, 0.11, -0.39, 0.47, 0.17, 0.69, 0.50, -0.16, 0.46,
  0.29, 0.43, 0.53, 0.88, 1.10, 0.52, 0.36, -0.10, 0.40, -0.26,
  -0.28, 0.06, 0.93, 0.43, -0.38, -0.15, -0.31, 0.63, -0.51, -0.46,
  0.63, -0.18, 0.54, 0.46, 0.77, 1.58, 0.19, 0.75, 0.06, 0.98,
  0.47, 0.85, 0.89, 1.24, 0.79, 0.73, -0.27, 0.59, 0.52, 0.36,
  -0.20, -0.13, 1.14, 0.12, 0.39, 0.43, -0.11, -0.52, 0.69, -1.61,
  1.25, 0.35, -0.21, 1.72, 1.94, 1.67, 1.46, 1.57, 1.32, 0.05,
  1.67, 1.21, 1.23, 1.53, 1.71, 1.50, 1.55, 1.48, 1.67, 1.84,
  1.77, 1.43, 1.93, 0.78, 2.30, 1.59, 1.08, 1.49, 0.77, 1.37,
  1.57, 1.19, 1.45, 1.86, 1.55, 2.50, 2.05, 0.95, 0.46, 1.71,
  0.86, 1.84, 1.65, 1.27, 1.95, 1.60, 1.58, 0.55, 0.72, 0.92,
  -0.24, 1.70, 1.51, 2.33, 1.62, 1.71, 0.67, 1.61, 1.06, 1.87,
  2.42, 1.28, 0.90, 0.51, 1.26, 2.05, 1.65, 2.21, 0.05, 1.00,
  1.23, -0.01, 1.59, 1.52, 1.51, -0.16, 1.23, 1.44, 0.97, 0.22,
  -0.36, 1.42, 1.16, -0.09, 1.09, 1.21, 0.72, -0.51, 0.72, -0.33,
  1.23, 0.84, -0.51, 2.64, 1.15, 1.82, 1.64, 1.85, 2.40, 0.30,
  2.04, 2.66, 1.97, 2.32, -0.08, 0.21, 1.44, 0.69, 0.48, 2.43,
  -0.11, -0.30, -0.44, -0.08, -0.42, 0.00, -0.46, -0.35, -0.11, -0.09,
  -0.52, -0.15, -0.11, -0.40, -0.61, -0.00, -0.11, -0.13, -0.25, -0.08,
  -0.11, 0.00, -0.51, -0.03, -0.04, 0.00, 0.00, -0.00, -0.11, 0.00,
  0.00, 0.00, -0.09, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
  -0.11, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.84, 1.18,
  0.91, 2.82, 1.57, 2.01, 2.88, 2.98, 2.73, 3.65, 2.82, 2.97,
  2.87, 2.60, 2.60, 2.73, 0.57, 0.06, 0.27, -1.34, 0.05, 3.16,
  0.15, -0.02, -1.49, -0.34, -2.28, -1.61, -0.85, -2.09, -1.46, -1.47,
  -1.91, -2.32, -2.20, -1.57, -1.54, -1.84, -2.20, -0.55, -1.40, -1.25,
  -1.26, -1.14, -0.83, -0.79, -0.74, -1.18, -0.75, -0.72, -1.08, -1.20
  -0.30, -0.41, 0.01, 9.6, 15.0, 15.0, 0.2, 0.1, 0.1, 0.83,
  0.83, 3.4, 0.4, 0.0, 3.4, 0.4, 0.0
};

// from Hotknot 2.0
static float dp03[] = {
  -0.90, -2.20, -2.10, -0.60, -1.10, -1.40, -2.10, -3.30, -2.40, -1.40,
  -2.10, -2.40, -3.40, -1.50, -2.50, -1.30, -0.50, 1.30, -1.30, -1.00,
  0.30, -0.30, -0.50, -0.30, -0.30, -0.10, -0.20, -1.50, -0.20, -1.10,
  -1.20, -0.20, 0.20, -0.30, -0.30, -0.60, -1.10, -1.50, -1.50, -1.40,
  -1.80, -1.00, -0.90, -2.90, -0.80, -2.20, -2.00, -1.60, -1.10, -1.70,
  -1.40, -1.80, -2.00, -1.10, -1.50, -1.30, -2.10, -1.10, -0.70, -2.40,
  -0.50, -2.40, -2.90, -1.40, -1.20, -1.90, -1.00, -2.20, -1.50, 0.20,
  -0.50, -0.30, -0.30, -0.10, -0.20, -1.50, -0.20, -0.90, -1.10, -0.30,
  0.00, -0.30, -0.30, -0.40, -1.10, -0.50, -0.30, -0.60, -0.50, -0.20,
  -0.10, -1.20, 0.00, -1.40, -1.20, -0.70, -0.20, -0.30, -0.10, -0.50,
  -0.80, -0.50, -0.30, -0.60, -0.50, -0.20, -0.10, -1.70, 0.00, -0.80,
  -1.20, -0.30, -0.70, -0.60, -0.10, -0.60, -0.80, 0.65, -1.10, -0.70,
  1.50, 1.00, 1.10, 1.20, 0.40, 1.10, -0.40, 0.40, 0.40, 0.40,
  0.30, 0.50, 0.40, 0.50, 0.40, -0.10, -1.70, -1.40, 0.00, 1.10,
  -0.30, 0.40, 0.80, 0.40, 0.40, 0.40, 0.40, -2.10, 1.10, -0.70,
  1.80, 0.40, -1.70, 2.30, 2.10, 0.80, 2.20, 1.70, 0.60, 1.10,
  1.60, 0.40, 2.30, 2.20, 2.50, 2.20, 2.50, 1.90, 2.20, 2.20,
  2.20, 1.70, 0.80, 0.80, 2.20, 2.20, 1.70, 1.50, 1.20, 2.50,
  2.10, 1.20, 2.20, 1.70, 0.60, 2.10, 1.60, 0.40, 2.30, 2.20,
  2.50, 2.20, 2.50, 1.90, 2.20, 2.20, 2.20, 1.70, 0.80, 0.80,
  2.20, 2.20, 1.70, 1.20, 1.20, 4.00, 0.75, 2.80, 2.50, 0.30,
  2.30, 2.20, 2.20, 0.30, 1.40, -0.10, 2.20, -2.10, 0.60, 1.30,
  2.00, -0.70, 1.10, 1.70, 1.40, -0.70, 0.80, -1.10, 1.40, -4.20,
  -0.40, 1.50, 0.90, -1.30, 1.00, 1.00, 1.10, -2.60, 0.80, -4.10,
  -1.00, -4.90, -0.50, 2.80, 2.80, 0.70, 1.90, 2.80, 2.20, 0.70,
  1.50, -0.30, 2.80, -2.90, 1.10, 0.00, 1.80, 1.00, 0.00, 2.00,
  -0.80, -0.50, -0.80, -0.60, -1.70, -0.80, -1.70, -1.20, -1.10, -0.40,
  -1.30, -0.60, -0.80, -0.50, -0.80, -0.60, -0.70, -0.10, -0.70, -0.10,
  -0.70, -0.10, -0.70, -0.10, -0.30, -0.10, -0.20, -0.10, -0.20, -0.10,
  0.00, 0.00, -0.50, -0.10, -0.20, -0.10, -0.30, -0.10, -0.20, -0.10,
  -0.30, -0.10, -0.40, -0.10, -0.30, -0.10, -0.40, -0.10, 1.70, 1.80,
  2.00, 3.80, 2.80, 3.20, 3.60, 4.00, 4.40, 5.70, 5.60, 5.60,
  5.40, 5.90, 5.60, 6.40, 0.50, -2.20, 0.30, 1.60, 1.40, 3.40,
  0.40, 0.00, 4.10, -3.00, -3.00, -3.00, -3.00, -3.00, -3.00, -3.00,
  -3.00, -3.00, -2.50, -2.50, -2.50, -2.50, -2.50, -2.00, -2.00, -2.00,
  -2.00, -2.00, -1.50, -1.50, -1.50, -1.50, -1.50, -1.50, -1.50, -1.50,
  -1.50, -1.50, -1.50,
  9.6, 15.0, 15.0, 0.2, 0.1, 0.1, 0.83,
  0.83, 3.4, 0.4, 0.0, 3.4, 0.4, 0.0
};

// from Hotknot 2.0
static float dp09[] = {
  -0.71, -1.32, -1.54, -0.12, -0.74, -0.91, -1.49, -2.16, -1.52, -0.58,
  -1.68, -1.39, -2.16, -0.97, -1.57, -0.69, -0.73, -0.29, -0.77, -0.12,
  -0.58,  0.27,  0.87,  0.45,  1.18, -0.09,  0.26,  0.42,  0.42, -0.30,
  0.31,  0.16,  1.18,  0.28, -0.03,  0.86,  0.02, -0.44, -0.34, -0.00,
  -0.09, -0.27, -0.23, -1.14, -0.37, -0.68, -0.06, -0.04, 0.03, -0.32,
  -0.59, -0.95, -0.29, 0.05, -0.34, -0.96, 1.01, -0.28, 0.32, -0.32,
  0.30, -0.68, -0.14, -0.28, 0.30, 0.07, -0.40, -0.37, -0.91, 0.60,
  0.76, 0.74, 1.88, 0.35, 0.84, 0.58, 0.10, -1.17, 0.13, 0.36,
  1.00, 0.07, 0.41, 0.63, 0.43, 0.29, 0.53, 0.92, 0.64, 0.72,
  0.42, 0.12, 0.69, -0.44, 0.41, 0.17, 1.18, 0.57, -0.25, 0.30,
  0.30, 0.31, 0.52, 1.09, 1.04, 0.57, 0.40, -0.03, 0.56, -0.27,
  0.04, -0.22, 0.44, 0.42, -0.34, 0.09, -0.33, 0.77, -0.57, -0.46,
  0.77, -0.38, 0.48, 0.37, 1.00, 1.79, 0.02, 0.42, -0.04, 0.64,
  0.22, 0.65, 0.97, 1.26, 0.65, 0.37, -0.38, 0.79, 0.54, 0.21,
  -0.19, -0.53, 0.89, -0.13, 0.12, 0.68, -0.38, -0.71, 0.40, -1.93,
  1.03, 0.05, -0.52, 1.56, 1.72, 1.16, 1.06, 1.12, 0.83, 0.50,
  1.13, 0.60, 0.61, 1.66, 1.38, 1.74, 1.32, 1.02, 1.63, 1.28,
  1.62, 1.13, 1.47, -0.15, 2.22, 1.10, 0.57, 1.25, 0.53, 1.03,
  1.36, 1.08, 1.12, 1.45, 1.25, 1.44, 1.59, 0.62, -0.23, 1.66,
  0.39, 1.34, 1.43, 1.30, 1.70, 1.24, 1.53, 0.12, 0.12, 0.65,
  -0.82, 1.34, 1.35, 2.52, 1.48, 1.54, 0.76, 1.48, 0.70, 1.55,
  2.09, 1.37, 1.14, 0.11, 1.36, 1.70, 2.06, 1.71, 0.31, 0.93,
  0.58, -0.52, 0.56, 1.25, 1.25, -0.43, 0.71, 0.62, 0.64, -0.58,
  -0.76, 1.25, 0.56, 0.52, -0.01, 1.12, 1.16, -0.82, 0.35, -1.05,
  1.02, -0.04, -0.94, 2.26, 2.04, 1.43, 1.71, 2.07, 2.28, 0.27,
  1.97, 2.17, 1.70, 1.67, 0.22, 0.22, 1.68, 0.79, 0.97, 2.08
  -0.32, -0.70, -0.90, -0.35, -0.87, -0.24, -0.98, -0.86, -0.46, -0.37,
  -0.95, -0.54, -0.32, -0.87, -1.08, -0.11, -0.32, -0.42, -0.68, -0.45,
  -0.32, -0.19, -0.86, -0.11, -0.14, 0.00, 0.00, -0.07, -0.32, -0.19,
  -0.13, -0.05, -0.09, -0.17, 0.00, -0.01, 0.00, -0.01, 0.00, -0.11,
  -0.32, -0.07, -0.06, -0.06, 0.00, 0.00, 0.00, 0.00, 0.44, 0.78,
  0.47, 2.81, 1.52, 2.03, 3.12, 2.67, 2.78, 3.69, 2.83, 3.23,
  2.87, 2.59, 2.65, 2.77, 0.50, 0.05, 0.38, -1.44, 0.31, 3.39,
  0.03, 0.02, -0.46, -0.65, -2.03, -1.30, -0.62, -1.66, -1.68, -1.58,
  -1.84, -1.96, -2.33, -1.86, -1.39, -1.42, -2.16, 0.20, -1.35, -1.02,
  -1.33, -0.97, -1.17, -0.80, -0.67, -0.84, -0.52, -0.69, -0.68, -0.18,
  -0.63, -0.33, 0.06,
  // -1.38, 10.07, 15.00, 2.46, 0.06, 0.96, 0.89,
  // 0.74, 3.39, 0.03, 0.02, 3.41, 0.56, 0.12
  9.6, 15.0, 15.0, 0.2, 0.1, 0.1, 0.83,
  0.83, 3.4, 0.4, 0.0, 3.4, 0.4, 0.0
};

template < class PF_TYPE >
void
Nupack<PF_TYPE>::
load_parameters_fm363(const std::vector<float>& v)
{
  std::vector<float>::const_iterator x=v.begin();
  //int b[] = { BASE_A-1, BASE_C-1, BASE_G-1, BASE_U-1 };
  int bp[4][4] = 
    {
      { -1, -1, -1, AU },
      { -1, -1, CG, -1 },
      { -1, GC, -1, GU },
      { UA, -1, UG, -1 }
    };
  int wc[4][4] = 
    {
      { 0, 0, 0, 1 },
      { 0, 0, 1, 0 },
      { 0, 1, 0, 0 },
      { 1, 0, 0, 0 }
    };
  int gc[4][4] = 
  {
      { 0, 0, 0, 0 },
      { 0, 0, 1, 0 },
      { 0, 1, 0, 0 },
      { 0, 0, 0, 0 }
  };

  // stack
  stack37[AU][AU] = stack37[UA][UA] = *x++;
  stack37[AU][CG] = stack37[GC][UA] = *x++;
  stack37[AU][GC] = stack37[CG][UA] = *x++;
  stack37[AU][GU] = stack37[UG][UA] = *x++;
  stack37[AU][UA] = stack37[AU][UA] = *x++;
  stack37[AU][UG] = stack37[GU][UA] = *x++;
  stack37[CG][AU] = stack37[UA][GC] = *x++;
  stack37[CG][CG] = stack37[GC][GC] = *x++;
  stack37[CG][GC] = stack37[CG][GC] = *x++;
  stack37[CG][GU] = stack37[UG][GC] = *x++;
  stack37[CG][UG] = stack37[GU][GC] = *x++;
  stack37[GC][AU] = stack37[UA][CG] = *x++;
  stack37[GC][CG] = stack37[GC][CG] = *x++;
  stack37[GC][GU] = stack37[UG][CG] = *x++;
  stack37[GC][UG] = stack37[GU][CG] = *x++;
  stack37[GU][AU] = stack37[UA][UG] = *x++;
  stack37[GU][GU] = stack37[UG][UG] = *x++;
  stack37[GU][UG] = stack37[GU][UG] = *x++;
  stack37[UA][AU] = stack37[UA][AU] = *x++;
  stack37[UA][GU] = stack37[UG][AU] = *x++;
  stack37[UG][GU] = stack37[UG][GU] = *x++;

  // terminal mismatch for hairpin loops
  mismatch_hairpin37[A][A][AU] = *x++;
  mismatch_hairpin37[A][C][AU] = *x++;
  mismatch_hairpin37[A][G][AU] = *x++;
  mismatch_hairpin37[A][U][AU] = *x++;
  mismatch_hairpin37[C][A][AU] = *x++;
  mismatch_hairpin37[C][C][AU] = *x++;
  mismatch_hairpin37[C][G][AU] = *x++;
  mismatch_hairpin37[C][U][AU] = *x++;
  mismatch_hairpin37[G][A][AU] = *x++;
  mismatch_hairpin37[G][C][AU] = *x++;
  mismatch_hairpin37[G][G][AU] = *x++;
  mismatch_hairpin37[G][U][AU] = *x++;
  mismatch_hairpin37[U][A][AU] = *x++;
  mismatch_hairpin37[U][C][AU] = *x++;
  mismatch_hairpin37[U][G][AU] = *x++;
  mismatch_hairpin37[U][U][AU] = *x++;
  mismatch_hairpin37[A][A][CG] = *x++;
  mismatch_hairpin37[A][C][CG] = *x++;
  mismatch_hairpin37[A][G][CG] = *x++;
  mismatch_hairpin37[A][U][CG] = *x++;
  mismatch_hairpin37[C][A][CG] = *x++;
  mismatch_hairpin37[C][C][CG] = *x++;
  mismatch_hairpin37[C][G][CG] = *x++;
  mismatch_hairpin37[C][U][CG] = *x++;
  mismatch_hairpin37[G][A][CG] = *x++;
  mismatch_hairpin37[G][C][CG] = *x++;
  mismatch_hairpin37[G][G][CG] = *x++;
  mismatch_hairpin37[G][U][CG] = *x++;
  mismatch_hairpin37[U][A][CG] = *x++;
  mismatch_hairpin37[U][C][CG] = *x++;
  mismatch_hairpin37[U][G][CG] = *x++;
  mismatch_hairpin37[U][U][CG] = *x++;
  mismatch_hairpin37[A][A][GC] = *x++;
  mismatch_hairpin37[A][C][GC] = *x++;
  mismatch_hairpin37[A][G][GC] = *x++;
  mismatch_hairpin37[A][U][GC] = *x++;
  mismatch_hairpin37[C][A][GC] = *x++;
  mismatch_hairpin37[C][C][GC] = *x++;
  mismatch_hairpin37[C][G][GC] = *x++;
  mismatch_hairpin37[C][U][GC] = *x++;
  mismatch_hairpin37[G][A][GC] = *x++;
  mismatch_hairpin37[G][C][GC] = *x++;
  mismatch_hairpin37[G][G][GC] = *x++;
  mismatch_hairpin37[G][U][GC] = *x++;
  mismatch_hairpin37[U][A][GC] = *x++;
  mismatch_hairpin37[U][C][GC] = *x++;
  mismatch_hairpin37[U][G][GC] = *x++;
  mismatch_hairpin37[U][U][GC] = *x++;
  mismatch_hairpin37[A][A][GU] = *x++;
  mismatch_hairpin37[A][C][GU] = *x++;
  mismatch_hairpin37[A][G][GU] = *x++;
  mismatch_hairpin37[A][U][GU] = *x++;
  mismatch_hairpin37[C][A][GU] = *x++;
  mismatch_hairpin37[C][C][GU] = *x++;
  mismatch_hairpin37[C][G][GU] = *x++;
  mismatch_hairpin37[C][U][GU] = *x++;
  mismatch_hairpin37[G][A][GU] = *x++;
  mismatch_hairpin37[G][C][GU] = *x++;
  mismatch_hairpin37[G][G][GU] = *x++;
  mismatch_hairpin37[G][U][GU] = *x++;
  mismatch_hairpin37[U][A][GU] = *x++;
  mismatch_hairpin37[U][C][GU] = *x++;
  mismatch_hairpin37[U][G][GU] = *x++;
  mismatch_hairpin37[U][U][GU] = *x++;
  mismatch_hairpin37[A][A][UA] = *x++;
  mismatch_hairpin37[A][C][UA] = *x++;
  mismatch_hairpin37[A][G][UA] = *x++;
  mismatch_hairpin37[A][U][UA] = *x++;
  mismatch_hairpin37[C][A][UA] = *x++;
  mismatch_hairpin37[C][C][UA] = *x++;
  mismatch_hairpin37[C][G][UA] = *x++;
  mismatch_hairpin37[C][U][UA] = *x++;
  mismatch_hairpin37[G][A][UA] = *x++;
  mismatch_hairpin37[G][C][UA] = *x++;
  mismatch_hairpin37[G][G][UA] = *x++;
  mismatch_hairpin37[G][U][UA] = *x++;
  mismatch_hairpin37[U][A][UA] = *x++;
  mismatch_hairpin37[U][C][UA] = *x++;
  mismatch_hairpin37[U][G][UA] = *x++;
  mismatch_hairpin37[U][U][UA] = *x++;
  mismatch_hairpin37[A][A][UG] = *x++;
  mismatch_hairpin37[A][C][UG] = *x++;
  mismatch_hairpin37[A][G][UG] = *x++;
  mismatch_hairpin37[A][U][UG] = *x++;
  mismatch_hairpin37[C][A][UG] = *x++;
  mismatch_hairpin37[C][C][UG] = *x++;
  mismatch_hairpin37[C][G][UG] = *x++;
  mismatch_hairpin37[C][U][UG] = *x++;
  mismatch_hairpin37[G][A][UG] = *x++;
  mismatch_hairpin37[G][C][UG] = *x++;
  mismatch_hairpin37[G][G][UG] = *x++;
  mismatch_hairpin37[G][U][UG] = *x++;
  mismatch_hairpin37[U][A][UG] = *x++;
  mismatch_hairpin37[U][C][UG] = *x++;
  mismatch_hairpin37[U][G][UG] = *x++;
  mismatch_hairpin37[U][U][UG] = *x++;

  // terminal mismatch for interior loops
  float internal_AU_closure = *x++;
  float internal_AG_mismatch = *x++;
  float internal_UU_mismatch = *x++;

  // interior loops 1x1
  int11_37[AU][AU][U][U] = int11_37[UA][UA][U][U] = *x++;
  int11_37[AU][CG][U][U] = int11_37[GC][UA][U][U] = *x++;
  int11_37[AU][GC][U][U] = int11_37[CG][UA][U][U] = *x++;
  int11_37[AU][UA][U][U] = int11_37[AU][UA][U][U] = *x++;
  int11_37[CG][CG][A][A] = int11_37[GC][GC][A][A] = *x++;
  int11_37[CG][GC][A][A] = int11_37[CG][GC][A][A] = *x++;
  int11_37[CG][CG][A][C] = int11_37[GC][GC][C][A] = *x++;
  int11_37[CG][GC][A][C] = int11_37[CG][GC][C][A] = *x++;
  int11_37[CG][CG][A][G] = int11_37[GC][GC][G][A] = *x++;
  int11_37[CG][GC][A][G] = int11_37[CG][GC][G][A] = *x++;
  int11_37[CG][CG][C][A] = int11_37[GC][GC][A][C] = *x++;
  int11_37[CG][CG][C][C] = int11_37[GC][GC][C][C] = *x++;
  int11_37[CG][GC][C][C] = int11_37[CG][GC][C][C] = *x++;
  int11_37[CG][CG][C][U] = int11_37[GC][GC][U][C] = *x++;
  int11_37[CG][GC][C][U] = int11_37[CG][GC][U][C] = *x++;
  int11_37[CG][CG][G][A] = int11_37[GC][GC][A][G] = *x++;
  int11_37[CG][CG][G][G] = int11_37[GC][GC][G][G] = *x++;
  int11_37[CG][GC][G][G] = int11_37[CG][GC][G][G] = *x++;
  int11_37[CG][CG][U][C] = int11_37[GC][GC][C][U] = *x++;
  int11_37[CG][AU][U][U] = int11_37[UA][GC][U][U] = *x++;
  int11_37[CG][CG][U][U] = int11_37[GC][GC][U][U] = *x++;
  int11_37[CG][GC][U][U] = int11_37[CG][GC][U][U] = *x++;
  int11_37[GC][CG][A][A] = int11_37[GC][CG][A][A] = *x++;
  int11_37[GC][CG][A][C] = int11_37[GC][CG][C][A] = *x++;
  int11_37[GC][CG][A][G] = int11_37[GC][CG][G][A] = *x++;
  int11_37[GC][CG][C][C] = int11_37[GC][CG][C][C] = *x++;
  int11_37[GC][CG][C][U] = int11_37[GC][CG][U][C] = *x++;
  int11_37[GC][CG][G][G] = int11_37[GC][CG][G][G] = *x++;
  int11_37[GC][AU][U][U] = int11_37[UA][CG][U][U] = *x++;
  int11_37[GC][CG][U][U] = int11_37[GC][CG][U][U] = *x++;
  int11_37[UA][AU][U][U] = int11_37[UA][AU][U][U] = *x++;
  float internal11_basic_mismatch = *x++;
  float internal11_GG_mismatch = *x++;
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      if (bp[i][j]>=0)
        for (int k=0; k!=4; k++)
          for (int l=0; l!=4; l++)
            for (int m=0; m!=4; m++)
              for (int n=0; n!=4; n++)
                if (bp[m][n]>=0)
                  if ( ((i==C && j==G) || (i==G && j==C)) &&
                       ((m==C && n==G) || (m==G && n==C)))
                  {
                    if (bp[k][l]>=0)
                      int11_37[bp[i][j]][bp[m][n]][k][l] = internal11_basic_mismatch;
                  }
                  else
                  {
                    if (!(wc[i][j] && wc[m][n] && k==U && l==U))
                    {
                      int11_37[bp[i][j]][bp[m][n]][k][l] = (k==G && l==G) ?
                        internal11_GG_mismatch : internal11_basic_mismatch;
                      if (!gc[i][j])
                        int11_37[bp[i][j]][bp[m][n]][k][l] += internal_AU_closure;
                      if (!gc[m][n])
                        int11_37[bp[i][j]][bp[m][n]][k][l] += internal_AU_closure;
                    }
                  }
                      
  // interoir loops 2x1
  int21_37[CG][A][A][CG][A] = *x++;
  int21_37[CG][A][A][CG][C] = *x++;
  int21_37[CG][A][A][CG][G] = *x++;
  int21_37[CG][C][A][CG][A] = *x++;
  int21_37[CG][C][A][CG][C] = *x++;
  int21_37[CG][C][A][CG][G] = *x++;
  int21_37[CG][G][A][CG][A] = *x++;
  int21_37[CG][G][A][CG][C] = *x++;
  int21_37[CG][G][A][CG][G] = *x++;
  int21_37[CG][A][C][CG][A] = *x++;
  int21_37[CG][A][C][CG][C] = *x++;
  int21_37[CG][A][C][CG][U] = *x++;
  int21_37[CG][C][C][CG][A] = *x++;
  int21_37[CG][C][C][CG][C] = *x++;
  int21_37[CG][C][C][CG][U] = *x++;
  int21_37[CG][U][C][CG][A] = *x++;
  int21_37[CG][U][C][CG][C] = *x++;
  int21_37[CG][U][C][CG][U] = *x++;
  int21_37[CG][A][G][CG][A] = *x++;
  int21_37[CG][A][G][CG][G] = *x++;
  int21_37[CG][G][G][CG][A] = *x++;
  int21_37[CG][G][G][CG][G] = *x++;
  int21_37[CG][C][U][CG][C] = *x++;
  int21_37[CG][C][U][CG][U] = *x++;
  int21_37[CG][U][U][CG][C] = *x++;
  int21_37[CG][U][U][CG][U] = *x++;
  int21_37[GC][A][A][GC][A] = *x++;
  int21_37[GC][A][A][GC][C] = *x++;
  int21_37[GC][A][A][GC][G] = *x++;
  int21_37[GC][C][A][GC][A] = *x++;
  int21_37[GC][C][A][GC][C] = *x++;
  int21_37[GC][C][A][GC][G] = *x++;
  int21_37[GC][G][A][GC][A] = *x++;
  int21_37[GC][G][A][GC][C] = *x++;
  int21_37[GC][G][A][GC][G] = *x++;
  int21_37[GC][A][C][GC][A] = *x++;
  int21_37[GC][A][C][GC][C] = *x++;
  int21_37[GC][A][C][GC][U] = *x++;
  int21_37[GC][C][C][GC][A] = *x++;
  int21_37[GC][C][C][GC][C] = *x++;
  int21_37[GC][C][C][GC][U] = *x++;
  int21_37[GC][U][C][GC][A] = *x++;
  int21_37[GC][U][C][GC][C] = *x++;
  int21_37[GC][U][C][GC][U] = *x++;
  int21_37[GC][A][G][GC][A] = *x++;
  int21_37[GC][A][G][GC][G] = *x++;
  int21_37[GC][G][G][GC][A] = *x++;
  int21_37[GC][G][G][GC][G] = *x++;
  int21_37[GC][C][U][GC][C] = *x++;
  int21_37[GC][C][U][GC][U] = *x++;
  int21_37[GC][U][U][GC][C] = *x++;
  int21_37[GC][U][U][GC][U] = *x++;
  float internal21_match = *x++;
  float internal21_AU_closure = *x++;

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      if (bp[i][j]>=0)
        for (int k=0; k<4; k++)
          for (int l=0; l<4; l++)
            for (int m=0; m<4; m++)
              for (int n=0; n<4; n++)
                if (bp[m][n]>=0)
                  for(int o=0; o<4; o++)
                    if ((i==C && j==G && m==C && n==G) ||
                        (i==G && j==C && m==G && n==C))
                    {
                      if (bp[k][l]>=0 || bp[k][o]>=0)
                        int21_37[bp[i][j]][l][k][bp[m][n]][o] = internal21_match;
                    }
                    else
                    {
                      if (bp[k][l]>=0 || bp[k][o]>=0)
                        int21_37[bp[i][j]][l][k][bp[m][n]][o] = internal21_match;
                      else
                        int21_37[bp[i][j]][l][k][bp[m][n]][o] =
                          int21_37[CG][l][k][CG][o]/2.0 + int21_37[GC][l][k][GC][o]/2.0;
                      if (!gc[i][j])
                        int21_37[bp[i][j]][l][k][bp[m][n]][o] += internal21_AU_closure;
                      if (!gc[m][n])    
                        int21_37[bp[i][j]][l][k][bp[m][n]][o] += internal21_AU_closure;
                    }

  // interior loops 2x2
  int22_37[AU][UA][A][A][A][A] = int22_37[AU][UA][A][A][A][A] = *x++;
  int22_37[AU][UA][A][C][C][A] = int22_37[AU][UA][A][C][C][A] = *x++;
  int22_37[AU][UA][A][G][G][A] = int22_37[AU][UA][A][G][G][A] = *x++;
  int22_37[AU][UA][C][A][A][C] = int22_37[AU][UA][C][A][A][C] = *x++;
  int22_37[AU][UA][C][C][C][C] = int22_37[AU][UA][C][C][C][C] = *x++;
  int22_37[AU][UA][C][U][U][C] = int22_37[AU][UA][C][U][U][C] = *x++;
  int22_37[AU][UA][G][A][A][G] = int22_37[AU][UA][G][A][A][G] = *x++;
  int22_37[AU][UA][G][G][G][G] = int22_37[AU][UA][G][G][G][G] = *x++;
  int22_37[AU][UA][G][U][U][G] = int22_37[AU][UA][G][U][U][G] = *x++;
  int22_37[AU][UA][U][C][C][U] = int22_37[AU][UA][U][C][C][U] = *x++;
  int22_37[AU][UA][U][G][G][U] = int22_37[AU][UA][U][G][G][U] = *x++;
  int22_37[AU][UA][U][U][U][U] = int22_37[AU][UA][U][U][U][U] = *x++;
  int22_37[CG][GC][A][A][A][A] = int22_37[CG][GC][A][A][A][A] = *x++;
  int22_37[CG][GC][A][C][C][A] = int22_37[CG][GC][A][C][C][A] = *x++;
  int22_37[CG][GC][A][G][G][A] = int22_37[CG][GC][A][G][G][A] = *x++;
  int22_37[CG][GC][C][A][A][C] = int22_37[CG][GC][C][A][A][C] = *x++;
  int22_37[CG][GC][C][C][C][C] = int22_37[CG][GC][C][C][C][C] = *x++;
  int22_37[CG][GC][C][U][U][C] = int22_37[CG][GC][C][U][U][C] = *x++;
  int22_37[CG][GC][G][A][A][G] = int22_37[CG][GC][G][A][A][G] = *x++;
  int22_37[CG][GC][G][G][G][G] = int22_37[CG][GC][G][G][G][G] = *x++;
  int22_37[CG][GC][G][U][U][G] = int22_37[CG][GC][G][U][U][G] = *x++;
  int22_37[CG][GC][U][C][C][U] = int22_37[CG][GC][U][C][C][U] = *x++;
  int22_37[CG][GC][U][G][G][U] = int22_37[CG][GC][U][G][G][U] = *x++;
  int22_37[CG][GC][U][U][U][U] = int22_37[CG][GC][U][U][U][U] = *x++;
  int22_37[GC][CG][A][A][A][A] = int22_37[GC][CG][A][A][A][A] = *x++;
  int22_37[GC][CG][A][C][C][A] = int22_37[GC][CG][A][C][C][A] = *x++;
  int22_37[GC][CG][A][G][G][A] = int22_37[GC][CG][A][G][G][A] = *x++;
  int22_37[GC][CG][C][A][A][C] = int22_37[GC][CG][C][A][A][C] = *x++;
  int22_37[GC][CG][C][C][C][C] = int22_37[GC][CG][C][C][C][C] = *x++;
  int22_37[GC][CG][C][U][U][C] = int22_37[GC][CG][C][U][U][C] = *x++;
  int22_37[GC][CG][G][A][A][G] = int22_37[GC][CG][G][A][A][G] = *x++;
  int22_37[GC][CG][G][G][G][G] = int22_37[GC][CG][G][G][G][G] = *x++;
  int22_37[GC][CG][G][U][U][G] = int22_37[GC][CG][G][U][U][G] = *x++;
  int22_37[GC][CG][U][C][C][U] = int22_37[GC][CG][U][C][C][U] = *x++;
  int22_37[GC][CG][U][G][G][U] = int22_37[GC][CG][U][G][G][U] = *x++;
  int22_37[GC][CG][U][U][U][U] = int22_37[GC][CG][U][U][U][U] = *x++;
  int22_37[UA][AU][A][A][A][A] = int22_37[UA][AU][A][A][A][A] = *x++;
  int22_37[UA][AU][A][C][C][A] = int22_37[UA][AU][A][C][C][A] = *x++;
  int22_37[UA][AU][A][G][G][A] = int22_37[UA][AU][A][G][G][A] = *x++;
  int22_37[UA][AU][C][A][A][C] = int22_37[UA][AU][C][A][A][C] = *x++;
  int22_37[UA][AU][C][C][C][C] = int22_37[UA][AU][C][C][C][C] = *x++;
  int22_37[UA][AU][C][U][U][C] = int22_37[UA][AU][C][U][U][C] = *x++;
  int22_37[UA][AU][G][A][A][G] = int22_37[UA][AU][G][A][A][G] = *x++;
  int22_37[UA][AU][G][G][G][G] = int22_37[UA][AU][G][G][G][G] = *x++;
  int22_37[UA][AU][G][U][U][G] = int22_37[UA][AU][G][U][U][G] = *x++;
  int22_37[UA][AU][U][C][C][U] = int22_37[UA][AU][U][C][C][U] = *x++;
  int22_37[UA][AU][U][G][G][U] = int22_37[UA][AU][U][G][G][U] = *x++;
  int22_37[UA][AU][U][U][U][U] = int22_37[UA][AU][U][U][U][U] = *x++;
  float internal22_delta_same_size = *x++;
  float internal22_delta_different_size = *x++;
  float internal22_delta_1stable_1unstable = *x++;
  float internal22_delta_AC = *x++;
  float internal22_match = *x++;

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      if (bp[i][j]>=0)
        for (int k=0; k<4; k++)
          for (int l=0; l<4; l++)
            for (int m=0; m<4; m++)
              for (int n=0; n<4; n++)
                if (bp[m][n]>=0)
                  for(int o=0; o<4; o++)
                    for (int p=0; p<4; p++)
                    {
                      int ii, jj, mm, nn;
                      if (i==G && j==U)   ii = A;     else ii = i;
                      if (i==U && j==G)   jj = A;     else jj = j;
                      if (m==G && n==U)   mm = A;     else mm = m;
                      if (m==U && n==G)   nn = A;     else nn = n;
                      if (wc[k][l] || wc[o][p])
                      {
                        int22_37[bp[i][j]][bp[m][n]][k][l][o][p] = internal22_match;
                      }
                      else if (nn==ii && mm==jj && p==k && o==l)  // the UG closing pairs are the same as UA
                      {
                        int22_37[bp[i][j]][bp[m][n]][k][l][o][p] = int22_37[bp[ii][jj]][bp[mm][nn]][k][l][o][p];
                      }
                      else //if (!(n==i && m==j && p==k && o==l))   // was already filled above
                      {
                        int result = check_stability_and_size (k, l, o, p);
                        // I approximate it to int, so that we obtain the same as by counte_each_structure_type and summing up
                        float temp = int22_37[bp[ii][jj]][bp[jj][ii]][k][l][l][k]/2.0 +
                          int22_37[bp[nn][mm]][bp[mm][nn]][p][o][o][p]/2.0;
                        // rounf it to match Turner parameters
                        //if (temp%10 == 5) temp -= 5; if (temp%10 == -5) temp += 5;
                        int22_37[bp[i][j]][bp[m][n]][k][l][o][p] =  temp;
                        switch (result)
                        {
                          case 1: int22_37[bp[i][j]][bp[m][n]][k][l][o][p] += internal22_delta_same_size; break;
                          case 2: int22_37[bp[i][j]][bp[m][n]][k][l][o][p] += internal22_delta_different_size; break;
                          case 3: int22_37[bp[i][j]][bp[m][n]][k][l][o][p] += internal22_delta_1stable_1unstable; break;
                          case 4: int22_37[bp[i][j]][bp[m][n]][k][l][o][p] += internal22_delta_AC; break;
                          default: assert(!"unreachable");
                        }                                                
                      }                                        
                    }

  // dangle5
  dangle5_37[AU][A] = *x++;
  dangle5_37[AU][C] = *x++;
  dangle5_37[AU][G] = *x++;
  dangle5_37[AU][U] = *x++;
  dangle5_37[CG][A] = *x++;
  dangle5_37[CG][C] = *x++;
  dangle5_37[CG][G] = *x++;
  dangle5_37[CG][U] = *x++;
  dangle5_37[GC][A] = *x++;
  dangle5_37[GC][C] = *x++;
  dangle5_37[GC][G] = *x++;
  dangle5_37[GC][U] = *x++;
  dangle5_37[GU][A] = *x++;
  dangle5_37[GU][C] = *x++;
  dangle5_37[GU][G] = *x++;
  dangle5_37[GU][U] = *x++;
  dangle5_37[UA][A] = *x++;
  dangle5_37[UA][C] = *x++;
  dangle5_37[UA][G] = *x++;
  dangle5_37[UA][U] = *x++;
  dangle5_37[UG][A] = *x++;
  dangle5_37[UG][C] = *x++;
  dangle5_37[UG][G] = *x++;
  dangle5_37[UG][U] = *x++;

  // dangle3
  dangle3_37[AU][A] = *x++;
  dangle3_37[AU][C] = *x++;
  dangle3_37[AU][G] = *x++;
  dangle3_37[AU][U] = *x++;
  dangle3_37[CG][A] = *x++;
  dangle3_37[CG][C] = *x++;
  dangle3_37[CG][G] = *x++;
  dangle3_37[CG][U] = *x++;
  dangle3_37[GC][A] = *x++;
  dangle3_37[GC][C] = *x++;
  dangle3_37[GC][G] = *x++;
  dangle3_37[GC][U] = *x++;
  dangle3_37[GU][A] = *x++;
  dangle3_37[GU][C] = *x++;
  dangle3_37[GU][G] = *x++;
  dangle3_37[GU][U] = *x++;
  dangle3_37[UA][A] = *x++;
  dangle3_37[UA][C] = *x++;
  dangle3_37[UA][G] = *x++;
  dangle3_37[UA][U] = *x++;
  dangle3_37[UG][A] = *x++;
  dangle3_37[UG][C] = *x++;
  dangle3_37[UG][G] = *x++;
  dangle3_37[UG][U] = *x++;

  // loop length for interior loops
  for (int i=4; i!=7; ++i) interior37[i-1] = *x++;
  for (int i=7; i<=30; ++i)
    interior37[i-1] = interior37[6-1] + loop_greater30 * LOG(1.0*i/7);

  // loop length for bulge loops
  for (int i=1; i!=7; ++i) bulge37[i-1] = *x++;
  for (int i=7; i<=30; ++i)
    bulge37[i-1] = bulge37[6-1] + loop_greater30 * LOG(1.0*i/7);

  // loop length for hairpin loops
  for (int i=3; i!=10; ++i) hairpin37[i-1] = *x++;
  for (int i=10; i<=30; ++i)
    hairpin37[i-1] = hairpin37[9-1] + loop_greater30 * LOG(1.0*i/10);

  // misc loops
  at_penalty = *x++;
  hairpin_GGG = *x++;
  polyC_slope = *x++;
  polyC_int = *x++;
  polyC_penalty = *x++;

  // terminal mismatch for interior loops
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      if (bp[i][j]>=0)
        for (int k=0; k<4; k++)
          for (int l=0; l<4; l++)
          {
            mismatch_interior37[k][l][bp[i][j]] = 0;
            if (((i == A || i == G) && j == U) || ((j == A || j == G) && i == U))
              mismatch_interior37[k][l][bp[i][j]] += internal_AU_closure;
            if ((k == A && l == G) || (l == A && k == G))
              mismatch_interior37[k][l][bp[i][j]] += internal_AG_mismatch;
            if (k == U && l == U)
              mismatch_interior37[k][l][bp[i][j]] += internal_UU_mismatch;
          }

  // multiloop penalty
  multiloop_penalty = *x++;
  multiloop_paired_penalty = *x++;
  multiloop_unpaired_penalty = *x++;
  intermolecular_initiation = *x++;

  // triloop
  //std::fill(triloop37.data(), triloop37.data()+triloop37.num_elements(), 0.0);
  std::fill(&triloop37[0][0][0][0][0], &triloop37[0][0][0][0][0]+4*4*4*4*4, 0.0);

  // tloops
  //std::fill(tloop37.data(), tloop37.data()+tloop37.num_elements(), 0.0);
  std::fill(&tloop37[0][0][0][0][0][0], &tloop37[0][0][0][0][0][0]+4*4*4*4*4*4, 0.0);
  tloop37[G][G][G][G][A][C] = *x++;
  tloop37[G][G][U][G][A][C] = *x++;
  tloop37[C][G][A][A][A][G] = *x++;
  tloop37[G][G][A][G][A][C] = *x++;
  tloop37[C][G][C][A][A][G] = *x++;
  tloop37[G][G][A][A][A][C] = *x++;
  tloop37[C][G][G][A][A][G] = *x++;
  tloop37[C][U][U][C][G][G] = *x++;
  tloop37[C][G][U][G][A][G] = *x++;
  tloop37[C][G][A][A][G][G] = *x++;
  tloop37[C][U][A][C][G][G] = *x++;
  tloop37[G][G][C][A][A][C] = *x++;
  tloop37[C][G][C][G][A][G] = *x++;
  tloop37[U][G][A][G][A][G] = *x++;
  tloop37[C][G][A][G][A][G] = *x++;
  tloop37[A][G][A][A][A][U] = *x++;
  tloop37[C][G][U][A][A][G] = *x++;
  tloop37[C][U][A][A][C][G] = *x++;
  tloop37[U][G][A][A][A][G] = *x++;
  tloop37[G][G][A][A][G][C] = *x++;
  tloop37[G][G][G][A][A][C] = *x++;
  tloop37[U][G][A][A][A][A] = *x++;
  tloop37[A][G][C][A][A][U] = *x++;
  tloop37[A][G][U][A][A][U] = *x++;
  tloop37[C][G][G][G][A][G] = *x++;
  tloop37[A][G][U][G][A][U] = *x++;
  tloop37[G][G][C][G][A][C] = *x++;
  tloop37[G][G][G][A][G][C] = *x++;
  tloop37[G][U][G][A][A][C] = *x++;
  tloop37[U][G][G][A][A][A] = *x++;

  // asymmetry penalty
  max_asymmetry = 3.00;
  asymmetry_penalty[0] = 0.50;
  asymmetry_penalty[1] = 0.50;
  asymmetry_penalty[2] = 0.50;
  asymmetry_penalty[3] = 0.50;

  // peudoknot parameters
  pk_penalty = *x++;
  pk_multiloop_penalty = *x++;
  pk_pk_penalty = *x++;
  pk_band_penalty = *x++;
  pk_unpaired_penalty = *x++;
  pk_paired_penalty = *x++;
  pk_stack_span = *x++;
  pk_interior_span = *x++;
  multiloop_penalty_pk = *x++;
  multiloop_paired_penalty_pk = *x++;
  multiloop_unpaired_penalty_pk = *x++;
}
#endif

template < class PF_TYPE >
bool
Nupack<PF_TYPE>::
load_parameters(const char* file)
{
  int p[] = { PAIR_AU, PAIR_CG, PAIR_GC, PAIR_UA, PAIR_GU, PAIR_UG };
  int b[] = { BASE_A-1, BASE_C-1, BASE_G-1, BASE_U-1 };
  
  std::ifstream is(file);
  if (!is) return false;

  std::string line;

  // stack
  std::getline(is, line);
  while (line[0]=='>')
    std::getline(is, line);
  for (int i=0; i!=6; ++i)
  {
    std::istringstream ss(line);
    for (int j=0; j!=6; ++j)
    {
      int v; ss >> v;
      stack37[p[i]][p[j]] = v/100.0;
    }
    std::getline(is, line);
  }
  
  // hairpin
  std::getline(is, line);
  while (line[0]=='>')
    std::getline(is, line);
  {
    int i, j, v;
    std::istringstream ss(line);
    for (i=0; i<30 && ss; ++i)
    {
      ss >> v;
      hairpin37[i] = v/100.0;
    }
    for (j=i; j<30; ++j)
      hairpin37[j] = hairpin37[i-1]+loop_greater30*LOG((j+1)/(1.0*i));
  }

  // bulge
  std::getline(is, line);
  while (line[0]=='>')
    std::getline(is, line);
  {
    int i, j, v;
    std::istringstream ss(line);
    for (i=0; i<30 && ss; ++i)
    {
      ss >> v;
      bulge37[i] = v/100.0;
    }
    for (j=i; j<30; ++j)
      bulge37[j] = bulge37[i-1]+loop_greater30*LOG((j+1)/(1.0*i));
  }


  // interior
  std::getline(is, line);
  while (line[0]=='>')
    std::getline(is, line);
  {
    int i, j, v;
    std::istringstream ss(line);
    for (i=0; i<30 && ss; ++i)
    {
      ss >> v;
      interior37[i] = v/100.0;
    }
    for (j=i; j<30; ++j)
      interior37[j] = interior37[i-1]+loop_greater30*LOG((j+1)/(1.0*i));
  }

  // asymmetry panelties
  std::getline(is, line);
  while (line[0]=='>')
    std::getline(is, line);
  {
    int v;
    std::istringstream ss(line);
    for (int i=0; i<4; ++i)
    {
      ss >> v;
      asymmetry_penalty[i] = v/100.0;
    }
    ss >> v; max_asymmetry = v/100.0;
  }

  // triloops
  //std::fill(triloop37.data(), triloop37.data()+triloop37.num_elements(), 0.0);
  std::fill(&triloop37[0][0][0][0][0], &triloop37[0][0][0][0][0]+4*4*4*4*4, 0.0);
  std::getline(is, line);
  while (line[0]=='>')
    std::getline(is, line);
  while (line[0]!='>')
  {
    int v;
    char loop[256];
    sscanf(line.c_str(), "%s %d", loop, &v);
    std::vector<int> idx(5);
    for (int i=0; i!=5; ++i) idx[i]=base(loop[i])-1;
    //triloop37(idx) = v/100.0;
    triloop37[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]] = v/100.0;
    std::getline(is, line);
  }

  // tloops
  //std::fill(tloop37.data(), tloop37.data()+tloop37.num_elements(), 0.0);
  std::fill(&tloop37[0][0][0][0][0][0], &tloop37[0][0][0][0][0][0]+4*4*4*4*4*4, 0.0);
  std::getline(is, line);
  while (line[0]=='>')
    std::getline(is, line);
  while (line[0]!='>')
  {
    int v;
    char loop[256];
    sscanf(line.c_str(), "%s %d", loop, &v);
    std::vector<int> idx(6);
    for (int i=0; i!=6; ++i) idx[i]=base(loop[i])-1;
    //tloop37(idx) = v/100.0;
    tloop37[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]] = v/100.0;
    std::getline(is, line);    
  }

  // mismatch hairpin
  while (line[0]=='>') std::getline(is, line);
  for (int i=0; i!=4; ++i)
  {
    for (int j=0; j!=4; ++j)
    {
      std::istringstream ss(line);
      for (int k=0; k!=6; ++k)
      {
        int v; ss >> v;
        mismatch_hairpin37[b[i]][b[j]][p[k]] = v/100.0;
      }
      std::getline(is, line);
    }
  }      

  // mismatch interior
  while (line[0]=='>') std::getline(is, line);
  for (int i=0; i!=4; ++i)
  {
    for (int j=0; j!=4; ++j)
    {
      std::istringstream ss(line);
      for (int k=0; k!=6; ++k)
      {
        int v; ss >> v;
        mismatch_interior37[b[i]][b[j]][p[k]] = v/100.0;
      }
      std::getline(is, line);
    }
  }

  // dangle5
  while (line[0]=='>') std::getline(is, line);
  for (int i=0; i!=6; ++i)
  {
    std::istringstream ss(line);
    for (int j=0; j!=4; ++j)
    {
      int v; ss >> v;
      dangle5_37[p[i]][b[j]] = v/100.0;
    }
    std::getline(is, line);
  }

  // dangle3
  while (line[0]=='>') std::getline(is, line);
  for (int i=0; i!=6; ++i)
  {
    std::istringstream ss(line);
    for (int j=0; j!=4; ++j)
    {
      int v; ss >> v;
      dangle3_37[p[i]][b[j]] = v/100.0;
    }
    std::getline(is, line);
  }

  // multiloop penalties
  while (line[0]=='>')
    std::getline(is, line);
  {
    int v;
    std::istringstream ss(line);
    ss >> v; multiloop_penalty = v/100.0;
    ss >> v; multiloop_paired_penalty = v/100.0;
    ss >> v; multiloop_unpaired_penalty = v/100.0;
    std::getline(is, line);
  }
  
  // AT terminate penalties
  while (line[0]=='>')
    std::getline(is, line);
  {
    int v;
    std::istringstream ss(line);
    ss >> v; at_penalty = v/100.0;
    std::getline(is, line);
  }

  // interior loops 1x1
  while (line[0]=='>')
    std::getline(is, line);
  for (int i=0; i!=6; ++i)
  {
    for (int j=0; j!=6; ++j)
    {
      std::getline(is, line);   // header
      for (int k=0; k!=4; ++k)
      {
        std::istringstream ss(line);
        for (int l=0; l!=4; ++l)
        {
          int v; ss >> v;
          int11_37[p[i]][p[j]][b[k]][b[l]] = v/100.0;
        }
        std::getline(is, line);
      }
    }
  }

  // interior loops 2x2
  while (line[0]=='>')
    std::getline(is, line);
  for (int i=0; i!=6; ++i)
  {
    for (int j=0; j!=6; ++j)
    {
      for (int m=0; m!=4; ++m)
      {
        for (int n=0; n!=4; ++n)
        {
          std::getline(is, line);   // header
          for (int k=0; k!=4; ++k)
          {
            std::istringstream ss(line);
            for (int l=0; l!=4; ++l)
            {
              int v; ss >> v;
              int22_37[p[i]][p[j]][b[m]][b[l]][b[n]][b[k]] = v/100.0;
            }
            std::getline(is, line);
          }
        }
      }
    }
  }

  // interior loops 1x2
  while (line[0]=='>')
    std::getline(is, line);
  for (int i=0; i!=6; ++i)
  {
    for (int j=0; j!=6; ++j)
    {
      for (int m=0; m!=4; ++m)
      {
        std::getline(is, line);   // header
        for (int k=0; k!=4; ++k)
        {
          std::istringstream ss(line);
          for (int l=0; l!=4; ++l)
          {
            int v; ss >> v;
            int21_37[p[i]][b[k]][b[m]][p[j]][b[l]] = v/100.0;
          }
          std::getline(is, line);
        }
      }
    }
  }

  // polyC hairpin parameters
  while (line[0]=='>')
    std::getline(is, line);
  {
    int v;
    std::istringstream ss(line);
    ss >> v; polyC_penalty = v/100.0;
    ss >> v; polyC_slope = v/100.0;
    ss >> v; polyC_int = v/100.0;
    std::getline(is, line);
  }

  // pseudoknot energy parameters
  while (line[0]=='>')
    std::getline(is, line);
  {
    int v;
    std::istringstream ss(line);
    ss >> v; pk_penalty = v/100.0;
    ss >> v; pk_paired_penalty = v/100.0;
    ss >> v; pk_unpaired_penalty = v/100.0;
    ss >> v; pk_multiloop_penalty = v/100.0;
    ss >> v; pk_pk_penalty = v/100.0;
    std::getline(is, line);
  }
  pk_band_penalty = 0.0;
  pk_stack_span = 1.0;
  pk_interior_span = 1.0;
  multiloop_penalty_pk = multiloop_penalty;
  multiloop_paired_penalty_pk = multiloop_paired_penalty;
  multiloop_unpaired_penalty_pk = multiloop_unpaired_penalty;

  // BIMOLECULAR TERM
  while (line[0]=='>')
    std::getline(is, line);
  {
    int v;
    std::istringstream ss(line);
    ss >> v; intermolecular_initiation = v/100.0;
    std::getline(is, line);
  }

  return true;
}

template < class PF_TYPE >
void
Nupack<PF_TYPE>::
load_default_parameters()
{
  int p[] = { PAIR_AU, PAIR_CG, PAIR_GC, PAIR_UA, PAIR_GU, PAIR_UG };
  int b[] = { BASE_A-1, BASE_C-1, BASE_G-1, BASE_U-1 };

  const int* v = &default_parameters[0];
  
  // stack
  for (int i=0; i!=6; ++i)
    for (int j=0; j!=6; ++j)
      stack37[p[i]][p[j]] = *(v++)/100.0;
  
  // hairpin
  for (int i=0; i<30; ++i)
    hairpin37[i] = *(v++)/100.0;

  // bulge
  for (int i=0; i<30; ++i)
    bulge37[i] = *(v++)/100.0;

  // interior
  for (int i=0; i<30; ++i)
    interior37[i] = *(v++)/100.0;

  // asymmetry panelties
  for (int i=0; i<4; ++i)
    asymmetry_penalty[i] = *(v++)/100.0;

  // mismatch hairpin
  for (int i=0; i!=4; ++i)
    for (int j=0; j!=4; ++j)
      for (int k=0; k!=6; ++k)
        mismatch_hairpin37[b[i]][b[j]][p[k]] = *(v++)/100.0;

  // mismatch interior
  for (int i=0; i!=4; ++i)
    for (int j=0; j!=4; ++j)
      for (int k=0; k!=6; ++k)
        mismatch_interior37[b[i]][b[j]][p[k]] = *(v++)/100.0;

  // dangle5
  for (int i=0; i!=6; ++i)
    for (int j=0; j!=4; ++j)
      dangle5_37[p[i]][b[j]] = *(v++)/100.0;

  // dangle3
  for (int i=0; i!=6; ++i)
    for (int j=0; j!=4; ++j)
      dangle3_37[p[i]][b[j]] = *(v++)/100.0;

  // multiloop penalties
  multiloop_penalty = *(v++)/100.0;
  multiloop_paired_penalty = *(v++)/100.0;
  multiloop_unpaired_penalty = *(v++)/100.0;
  
  // AT terminate penalties
  at_penalty = *(v++)/100.0;

  // interior loops 1x1
  for (int i=0; i!=6; ++i)
    for (int j=0; j!=6; ++j)
      for (int k=0; k!=4; ++k)
        for (int l=0; l!=4; ++l)
          int11_37[p[i]][p[j]][b[k]][b[l]] = *(v++)/100.0;

  // interior loops 2x2
  for (int i=0; i!=6; ++i)
    for (int j=0; j!=6; ++j)
      for (int m=0; m!=4; ++m)
        for (int n=0; n!=4; ++n)
          for (int k=0; k!=4; ++k)
            for (int l=0; l!=4; ++l)
              int22_37[p[i]][p[j]][b[m]][b[l]][b[n]][b[k]] = *(v++)/100.0;

  // interior loops 1x2
  for (int i=0; i!=6; ++i)
    for (int j=0; j!=6; ++j)
      for (int m=0; m!=4; ++m)
        for (int k=0; k!=4; ++k)
          for (int l=0; l!=4; ++l)
            int21_37[p[i]][b[k]][b[m]][p[j]][b[l]] = *(v++)/100.0;

  // polyC hairpin parameters
  polyC_penalty = *(v++)/100.0;
  polyC_slope = *(v++)/100.0;
  polyC_int = *(v++)/100.0;

  // pseudoknot energy parameters
  pk_penalty = *(v++)/100.0;
  pk_paired_penalty = *(v++)/100.0;
  pk_unpaired_penalty = *(v++)/100.0;
  pk_multiloop_penalty = *(v++)/100.0;
  pk_pk_penalty = *(v++)/100.0;
  pk_band_penalty = 0.0;
  pk_stack_span = 1.0;
  pk_interior_span = 1.0;
  multiloop_penalty_pk = multiloop_penalty;
  multiloop_paired_penalty_pk = multiloop_paired_penalty;
  multiloop_unpaired_penalty_pk = multiloop_unpaired_penalty;

  // BIMOLECULAR TERM
  intermolecular_initiation = *(v++)/100.0;

  // triloops
  //std::fill(triloop37.data(), triloop37.data()+triloop37.num_elements(), 0.0);
  std::fill(&triloop37[0][0][0][0][0], &triloop37[0][0][0][0][0]+4*4*4*4*4, 0.0);
  for (int i=0; default_triloops[i].s!=NULL; ++i)
  {
    int v=default_triloops[i].e;
    const char* loop=default_triloops[i].s;
    std::vector<int> idx(5);
    for (int i=0; i!=5; ++i) idx[i]=base(loop[i])-1;
    //triloop37(idx) = v/100.0;
    triloop37[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]] = v/100.0;
  }

  // tloops
  //std::fill(tloop37.data(), tloop37.data()+tloop37.num_elements(), 0.0);
  std::fill(&tloop37[0][0][0][0][0][0], &tloop37[0][0][0][0][0][0]+4*4*4*4*4*4, 0.0);
  for (int i=0; default_tloops[i].s!=NULL; ++i)
  {
    int v=default_tloops[i].e;
    const char* loop=default_tloops[i].s;
    std::vector<int> idx(6);
    for (int i=0; i!=6; ++i) idx[i]=base(loop[i])-1;
    //tloop37(idx) = v/100.0;
    tloop37[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]] = v/100.0;
  }
}

template < class PF_TYPE >
void
Nupack<PF_TYPE>::
dump_parameters(std::ostream& os) const
{
  // stack
  for (int i=0; i!=6; ++i)
    for (int j=0; j!=6; ++j)
      os << "stack37[" << i << "][" << j << "]=" << stack37[i][j] << std::endl;

  // hairpin
  for (int i=0; i!=30; ++i)
    os << "hairpin37[" << i << "]=" << hairpin37[i] << std::endl;

  // bulge
  for (int i=0; i!=30; ++i)
    os << "bulge37[" << i << "]=" << bulge37[i] << std::endl;

  // interior
  for (int i=0; i!=30; ++i)
    os << "interior37[" << i << "]=" << interior37[i] << std::endl;

  // asymmetry
  for (int i=0; i!=4; ++i)
    os << "asymmetry_penalty[" << i << "]=" << asymmetry_penalty[i] << std::endl;
  os << "max_asymmetry=" << max_asymmetry << std::endl;

  // triloops
  for (int i0=0; i0!=4; ++i0)
    for (int i1=0; i1!=4; ++i1)
      for (int i2=0; i2!=4; ++i2)
        for (int i3=0; i3!=4; ++i3)
          for (int i4=0; i4!=4; ++i4)
            if (triloop37[i0][i1][i2][i3][i4]!=0.0)
              os << "triloop37["<<i0<<"]["<<i1<<"]["<<i2<<"]["<<i3<<"]["<<i4<<"]="
                 << triloop37[i0][i1][i2][i3][i4] <<std::endl;

  // tloops
  for (int i0=0; i0!=4; ++i0)
    for (int i1=0; i1!=4; ++i1)
      for (int i2=0; i2!=4; ++i2)
        for (int i3=0; i3!=4; ++i3)
          for (int i4=0; i4!=4; ++i4)
            for (int i5=0; i5!=4; ++i5)
            if (tloop37[i0][i1][i2][i3][i4][i5]!=0.0)
              os << "tloop37["<<i0<<"]["<<i1<<"]["<<i2<<"]["<<i3<<"]["<<i4<<"]["<<i5<<"]="
                 << tloop37[i0][i1][i2][i3][i4][i5] <<std::endl;

  // mismatch hairpin
  for (int i=0; i!=4; ++i)
    for (int j=0; j!=4; ++j)
      for (int k=0; k!=6; ++k)
        os << "mismatch_hairpin37["<<i<<"]["<<j<<"]["<<k<<"]="
           << mismatch_hairpin37[i][j][k] << std::endl;

  // mismatch interior37
  for (int i=0; i!=4; ++i)
    for (int j=0; j!=4; ++j)
      for (int k=0; k!=6; ++k)
        os << "mismatch_interior37["<<i<<"]["<<j<<"]["<<k<<"]="
           << mismatch_interior37[i][j][k] << std::endl;

  // dangle5
  for (int i=0; i!=6; ++i)
    for (int j=0; j!=4; ++j)
      os << "dangle5_37["<<i<<"]["<<j<<"]=" << dangle5_37[i][j] << std::endl;

  // dangle3
  for (int i=0; i!=6; ++i)
    for (int j=0; j!=4; ++j)
      os << "dangle3_37["<<i<<"]["<<j<<"]=" << dangle3_37[i][j] << std::endl;

  // multiloop penalties
  os << "multiloop_penalty=" << multiloop_penalty << std::endl
     << "multiloop_paired_penalty=" << multiloop_paired_penalty << std::endl
     << "multiloop_unpaired_penalty=" << multiloop_unpaired_penalty << std::endl;

  // AT terminate penalties
  os << "at_penalty=" << at_penalty << std::endl;

  // interior loops 1x1
  for (int i=0; i!=6; ++i)
    for (int j=0; j!=6; ++j)
      for (int k=0; k!=4; ++k)
        for (int l=0; l!=4; ++l)
          os << "int11_37["<<i<<"]["<<j<<"]["<<k<<"]["<<l<<"]="
             << int11_37[i][j][k][l] << std::endl;

  // interior loops 2x2
  for (int i=0; i!=6; ++i)
    for (int j=0; j!=6; ++j)
      for (int m=0; m!=4; ++m)
        for (int n=0; n!=4; ++n)
          for (int k=0; k!=4; ++k)
            for (int l=0; l!=4; ++l)
              os << "int22_37["<<i<<"]["<<j<<"]["<<m<<"]["<<l<<"]["<<n<<"]["<<k<<"]="
                 << int22_37[i][j][m][l][n][k] << std::endl;

  // interior loops 1x2
  for (int i=0; i!=6; ++i)
    for (int j=0; j!=6; ++j)
      for (int m=0; m!=4; ++m)
        for (int k=0; k!=4; ++k)
          for (int l=0; l!=4; ++l)
            os << "int21_37["<<i<<"]["<<k<<"]["<<m<<"]["<<j<<"]["<<l<<"]="
               << int21_37[i][k][m][j][l] << std::endl;

  // polyC hairpin parameters
  os << "polyC_penalty=" << polyC_penalty << std::endl
     << "polyC_slope=" << polyC_slope << std::endl
     << "polyC_int=" << polyC_int << std::endl;

  // pseudoknot energy parameters
  os << "pk_penalty=" << pk_penalty << std::endl
     << "pk_paired_penalty=" << pk_paired_penalty << std::endl
     << "pk_unpaired_penalty=" << pk_unpaired_penalty << std::endl
     << "pk_multiloop_penalty=" << pk_multiloop_penalty << std::endl
     << "pk_pk_penalty=" << pk_pk_penalty << std::endl
     << "pk_band_penalty=" << pk_band_penalty << std::endl
     << "pk_stack_span=" << pk_stack_span << std::endl
     << "pk_interior_span=" << pk_interior_span << std::endl
     << "multiloop_penalty_pk=" << multiloop_penalty_pk << std::endl
     << "multiloop_paired_penalty_pk=" << multiloop_paired_penalty_pk << std::endl
     << "multiloop_unpaired_penalty_pk" << multiloop_unpaired_penalty_pk << std::endl;

  // BIMOLECULAR TERM
  os << "intermolecular_initiation=" << intermolecular_initiation << std::endl;

  // misc
  os << "loop_greater30=" << loop_greater30 << std::endl
     << "hairpin_GGG=" << hairpin_GGG << std::endl;

}


template < class PF_TYPE >
typename Nupack<PF_TYPE>::pf_type
Nupack<PF_TYPE>::
calculate_partition_function()
{
  Q.resize(N);    Q.fill(0.0);
  Qb.resize(N);   Qb.fill(0.0);
  Qm.resize(N);   Qm.fill(0.0);
  Qp.resize(N);   Qp.fill(0.0);
  Qz.resize(N);   Qz.fill(0.0);
  Qg.resize(N);   Qg.fill(0.0);
  Qgl.resize(N);  Qgl.fill(0.0);
  Qgr.resize(N);  Qgr.fill(0.0);
  Qgls.resize(N); Qgls.fill(0.0);
  Qgrs.resize(N); Qgrs.fill(0.0);
  for (int i=0; i!=N; ++i) Q(i,i-1) = Qz(i,i-1) = 1.0;
  DPTableX<PF_TYPE> Qx, Qx1, Qx2;

  for (int l=1; l<=N; ++l)
  {
    Qx.swap(Qx1);
    Qx1.swap(Qx2);
    Qx2.resize(l+1, N);
    Qx2.fill(0.0);

    for (int i=0; i+l<=N; ++i)
    {
      int j=i+l-1;

      // Qb recursion
      if (allow_paired(i,j))
      {
        Qb(i,j) = EXP( -score_hairpin(i,j)/RT );

        for (int d=i+1; d<=j-5; ++d) // all possible rightmost pairs d-e
        {
          for (int e=d+4; e<=j-1; ++e)
          {
            if (allow_paired(d,e))
            {
              Qb(i,j) += Qb(d,e) *
                EXP( -score_interior(i,d,e,j,false)/RT );

              if (d>=i+6 && wc_pair(d,e) && wc_pair(i,j))
              {
                Qb(i,j) += Qm(i+1,d-1) * Qb(d,e) *
                  EXP( -( score_multiloop(false) +
                          score_multiloop_paired(2,false) +
                          score_multiloop_unpaired(j-e-1,false) +
                          score_at_penalty(i,j) +
                          score_at_penalty(d,e) +
                          score_dangle(e+1,j-1) )/RT ) ;
              }
            }
          }
        }

        if (wc_pair(i,j))
        {
          for (int d=i+1; d<=j-6; ++d) // all possible rightmost pseudoknots filling [d,e]
          {
            for (int e=d+5; e<=j-1; ++e)
            {
              Qb(i,j) += Qp(d,e) *
                EXP( -( score_multiloop(false) +
                        score_pk_multiloop() +
                        score_multiloop_paired(3,false) +
                        score_multiloop_unpaired(j-e-1 + d-i-1,false) +
                        score_at_penalty(i,j) +
                        score_dangle(e+1,j-1) +
                        score_dangle(i+1,d-1) )/RT );
              
              Qb(i,j) += Qm(i+1,d-1) * Qp(d,e) *
                EXP( -( score_multiloop(false) +
                        score_pk_multiloop() +
                        score_multiloop_paired(3,false) +
                        score_multiloop_unpaired(j-e-1,false) +
                        score_at_penalty(i,j) +
                        score_dangle(e+1,j-1) )/RT );
            }
          }
        }
      }

      // Qg recursion
      if (allow_paired(i,j))
      {
        // case 0: only 1 pair
        Qg(i,i,j,j) = 1;

        // case 1: terminal inner pair
        for (int d=i+1; d<=j-5; ++d)
        {
          for (int e=d+4; e<=j-1; ++e)
          {
            if (allow_paired(d,e))
              Qg(i,d,e,j) += EXP( -score_interior(i,d,e,j,true)/RT );
          }
        }
      }

      fastiloops(i, j, Qg, Qx, Qx2);

      if (allow_paired(i,j) && wc_pair(i,j))
      {

        // case 2: multiloop left
        for (int d=i+6; d<=j-5; ++d)
        {
          for (int e=d+4; e<=j-1; ++e)
          {
            if (allow_paired(d,e) && wc_pair(d,e))
            {
              Qg(i,d,e,j) += Qm(i+1,d-1) *
                EXP( -( score_multiloop(true) +
                        score_multiloop_paired(2,true) +
                        score_multiloop_unpaired(j-e-1,true) + 
                        score_at_penalty(i,j) +
                        score_at_penalty(d,e) +
                        score_dangle(e+1, j-1) )/RT );
            }
          }
        }

        // case 3: multiloop right
        for (int d=i+1; d<=j-10; ++d)
        {
          for (int e=d+4; e<=j-6; ++e)
          {
            if (allow_paired(d,e) && wc_pair(d,e))
            {
              Qg(i,d,e,j) += Qm(e+1,j-1) *
                EXP( -( score_multiloop(true) +
                        score_multiloop_paired(2,true) +
                        score_multiloop_unpaired(d-i-1,true) +
                        score_at_penalty(i,j) +
                        score_at_penalty(d,e) +
                        score_dangle(i+1, d-1) )/RT );
            }
          }
        }

        // case 4: multiloop both sides
        for (int d=i+6; d<=j-10; ++d)
        {
          for (int e=d+4; e<=j-6; ++e)
          {
            if (allow_paired(d,e) && wc_pair(d,e))
            {
              Qg(i,d,e,j) += Qm(i+1,d-1) * Qm(e+1,j-1) *
                EXP( -( score_multiloop(true) +
                        score_multiloop_paired(2,true) +
                        score_at_penalty(i,j) +
                        score_at_penalty(d,e) )/RT );
            }
          }
        }

        // case 5: interior loop + multi left
        for (int d=i+7; d<=j-6; ++d)
        {
          for (int e=d+4; e<=j-2; ++e)
          {
            if (allow_paired(d,e))
            {
              for (int f=e+1; f<=j-1; ++f)
              {
                Qg(i,d,e,j) += Qgls(i+1,d,e,f) *
                  EXP( -( score_multiloop(true) +
                          score_multiloop_paired(1,true) +
                          score_multiloop_unpaired(j-f-1,true) +
                          score_at_penalty(i,j) +
                          score_dangle(f+1,j-1) )/RT );
              }
            }
          }
        }

        // case 6: interior loop + multi right
        for (int d=i+2; d<=j-11; ++d)
        {
          for (int e=d+4; e<=j-7; ++e)
          {
            if (allow_paired(d,e))
            {
              for (int c=i+1; c<=d-1; ++c)
              {
                Qg(i,d,e,j) += Qgrs(c,d,e,j-1) *
                  EXP( -( score_multiloop(true) +
                          score_multiloop_paired(1,true) +
                          score_multiloop_unpaired(c-i-1,true) +
                          score_at_penalty(i,j) +
                          score_dangle(i+1,c-1) )/RT );
              }
            }
          }
        }

        // case 7: interior loop + multi both sides
        for (int d=i+7; d<=j-11; ++d)
        {
          for (int e=d+4; e<=j-7; ++e)
          {
            if (allow_paired(d,e))
            {
              for (int c=i+6; c<=d-1; ++c)
              {
                Qg(i,d,e,j) += Qm(i+1,c-1) * Qgrs(c,d,e,j-1) *
                  EXP( -( score_multiloop(true) +
                          score_multiloop_paired(1,true) +
                          score_at_penalty(i,j) )/RT );
              }
            }
          }
        }
      }

      // Qgls recursion
      for (int c=i+5; c<=j-6; ++c)
      {
        if (allow_paired(c,j) && wc_pair(c,j))
        {
          for (int d=c+1; d<=j-5; ++d)
          {
            for (int e=d+4; e<=j-1; ++e)
            {
              if (allow_paired(d,e))
              {
                Qgls(i,d,e,j) += Qm(i,c-1) * Qg(c,d,e,j) *
                  EXP( -( score_multiloop_paired(1,true) +
                          score_at_penalty(c,j) )/RT );
              }
            }
          }
        }
      }

      // Qgrs recursion
      for (int d=i+1; d<=j-10; ++d)
      {
        for (int e=d+4; e<=j-6; ++e)
        {
          if (allow_paired(d,e))
          {
            for (int f=e+1; f<=j-5; ++f)
            {
              if (allow_paired(i,f) && wc_pair(i,f))
              {
                Qgrs(i,d,e,j) += Qg(i,d,e,f) * Qm(f+1,j) *
                  EXP( -( score_multiloop_paired(1,true) +
                          score_at_penalty(i,f) )/RT );
              }
            }
          }
        }
      }

      // Qgl recursions
      for (int d=i+1; d<=j-5; ++d)
      {
        for (int f=d+4; f<=j-1; ++f)
        {
          if (allow_paired(d,f) && wc_pair(d,f))
          {
            for (int e=d; e<=f-2; ++e) // f-3???
            {
              Qgl(i,e,f,j) += Qg(i,d,f,j) * Qz(d+1,e) *
                EXP( -( score_pk_paired(1) +
                        score_at_penalty(d,f) )/RT );
            }
          }
        }
      }

      // Qgr recursion
      for (int d=i+1; d<=j-3; ++d)
      {
        for (int e=d+2; e<=j-1; ++e)
        {
          for (int f=e; f<=j-1; ++f)
          {
            Qgr(i,d,e,j) += Qgl(i,d,f,j) * Qz(e,f-1);
          }
        }
      }

      // Qp recursion
      // case 1: both Qg are exactly 1 pair
      // first case is exactly 1 pair per Og
      if (j-i>4)
      {
        int a=i;
        int f=j;
        for (int b=a+1; b<=j-4; ++b)
        {
          if (allow_paired(b,j) && wc_pair(b,j))
          {
            int c=b;
            for (int d=std::max(c+1,a+4); d<=j-1; ++d)
            {
              if (allow_paired(a,d) && wc_pair(a,d))
              {
                int e=d;
                Qp(i,j) += Qg(i,a,d,e) * Qg(b,c,f,j) * Qz(e+1,f-1) * Qz(c+1,d-1) * Qz(a+1,b-1) *
                  EXP( -( score_pk_paired(2) +
                          score_pk_band(2) +
                          score_at_penalty(a,d) +
                          score_at_penalty(c,f) +
                          score_at_penalty(i,e) +
                          score_at_penalty(b,j) )/RT );
              }
            }
          }
        }
      }

      if (j-i>6)
      {
        // case 2 left Og is exactly 1 pair, right is 2+
        for (int d=i+1; d<=j-6; ++d)
        {
          if (allow_paired(d,j) && wc_pair(d,j))
          {
            for (int e=std::max(d+2,i+4); e<=j-2; ++e)
            {
              int f=e;
              if (allow_paired(i,f) && wc_pair(i,f))
              {
                Qp(i,j) += Qg(i,i,e,f) * Qz(i+1,d-1) * Qgr(d,e-1,f+1,j) *
                  EXP( -( score_pk_paired(1) +
                          score_pk_band(2) +
                          score_at_penalty(d,j) +
                          score_at_penalty(i,f)*2 )/RT );
              }
            }
          }
        }

        // case 2 left Qg is 2+ pairs, right is 1
        for (int d=i+2; d<=j-4; ++d)
        {
          if (allow_paired(d,j) && wc_pair(d,j))
          {
            for (int e=std::max(d+1,i+4); e<=j-2; ++e)
            {
              for (int f=e+1; f<=j-1; ++f)
              {
                if (allow_paired(i,f) && wc_pair(i,f))
                {
                  Qp(i,j) += Qgl(i,d-1,e,f) * Qg(d,d,j,j) * Qz(d+1,e-1) * Qz(f+1,j-1) *
                    EXP( -( score_pk_paired(1) +
                            score_pk_band(2) +
                            score_at_penalty(d,j)*2 +
                            score_at_penalty(i,f) )/RT );
                }
              }
            }
          }
        }
      }

      // otherwise
      if (j-i>7)
      {
        for (int d=i+2; d<=j-4; ++d)
        {
          if (allow_paired(d,j) && wc_pair(d,j))
          {
            for (int e=std::max(d+2,i+5); e<=j-3; ++e)
            {
              for (int f=e+1; f<=j-2; ++f)
              {
                if (allow_paired(i,f) && wc_pair(i,f))
                {
                  Qp(i,j) += Qgl(i,d-1,e,f) * Qgr(d,e-1,f+1,j) *
                    EXP( -( score_pk_band(2) +
                            score_at_penalty(d,j) +
                            score_at_penalty(i,j) )/RT );
                }
              }
            }
          }
        }
      }

      // Q, Qm, Qz recursions
      Q(i,j) = EXP( -score_dangle(i,j)/RT ); // empty recursion

      if (i!=0 && j!=N-1)
      {
        Qz(i,j) =
          EXP( -( score_dangle(i,j) +
                  score_pk_unpaired(j-i+1) )/RT );
      }

      for (int d=i; d<=j-4; ++d)
      {
        for (int e=d+4; e<=j; ++e)
        {
          if (allow_paired(d,e) && wc_pair(d,e))
          {
            Q(i,j) += Q(i,d-1) * Qb(d,e) *
              EXP( -( score_at_penalty(d,e) +
                      score_dangle(e+1,j) )/RT );

            if (i!=0 && j!=N-1)
            {
              Qm(i,j) += Qb(d,e) *
                EXP( -( score_multiloop_paired(1,false) +
                        score_multiloop_unpaired(d-i + j-e,false) +
                        score_at_penalty(d,e) +
                        score_dangle(e+1,j) +
                        score_dangle(i,d-1) )/RT );

              if (d>=i+5)
              {
                Qm(i,j) += Qm(i,d-1) * Qb(d,e) *
                  EXP( -( score_multiloop_paired(1,false) +
                          score_multiloop_unpaired(j-e,false) +
                          score_at_penalty(d,e) +
                          score_dangle(e+1,j) )/RT );
              }
                
              Qz(i,j) += Qz(i,d-1) * Qb(d,e) *
                EXP( -( score_pk_paired(1) +
                        score_pk_unpaired(j-e) +
                        score_at_penalty(d,e) +
                        score_dangle(e+1,j) )/RT );
            }
          }
        }
      }

      for (int d=i; d<=j-5; ++d)
      {
        for (int e=d+5; e<=j; ++e)
        {
          Q(i,j) += Q(i,d-1) * Qp(d,e) *
            EXP( -( score_pk() +
                    score_dangle(e+1,j) )/RT );

          if (i!=0 && j!=N-1)
          {
            Qm(i,j) += Qp(d,e) *
              EXP( -( score_pk_multiloop() +
                      score_multiloop_paired(2,false) +
                      score_multiloop_unpaired(d-i + j-e,false) +
                      score_dangle(e+1,j) +
                      score_dangle(i,d-1) )/RT );

            if (d>=i+5)
            {
              Qm(i,j) += Qm(i,d-1) * Qp(d,e) *
                EXP( -( score_pk_multiloop() +
                        score_multiloop_paired(2,false) +
                        score_multiloop_unpaired(j-e,false) +
                        score_dangle(e+1,j) )/RT );
            }

            Qz(i,j) += Qz(i,d-1) * Qp(d,e) *
              EXP( -( score_pk_pk() +
                      score_pk_paired(2) +
                      score_pk_unpaired(j-e) +
                      score_dangle(e+1,j) )/RT );
          }
        }
      }
      //printf("%d,%d: %Lf\n", i, j, Q(i,j));
    }
  }

  return Q(0,N-1);
}

template < class PF_TYPE >
void
Nupack<PF_TYPE>::
fastiloops(int i, int j, DPTable4<PF_TYPE>& Qg, DPTableX<PF_TYPE>& Qx, DPTableX<PF_TYPE>& Qx2)
{
  int l=j-i+1;
  if (l>=17)   // smallest subsequence not added to Qg as special case
  {
    for (int d=i+6; d<=j-10; ++d)
    {
      for (int e=d+4; e<=j-6; ++e)
      {
        if (allow_paired(d,e))
        {
          int l1=4;               // explicitly add in terms for l1=4, l2>=4
          int c=i+l1+1;
          for (int l2=4; l2<=j-e-2; ++l2)
          {
            int s=l1+l2;
            int f=j-l2-1;
            if (allow_paired(c,f))
            {
              Qx(i,d,e,s) += Qg(c,d,e,f) *
                EXP( -( score_interior_asymmetry(l1, l2) +
                        score_interior_mismatch(f,c,f+1,c-1) )/RT );
            }
          }

          if (d>=i+7)
          {
            int l2=4;             // explicitly add in terms of l1>=5, l2=4
            int f=j-l2-1;
            for (int l1=5; l1<=d-i-2; ++l1)
            {
              int s=l1+l2;
              int c=i+l1+1;
              if (allow_paired(c,f))
              {
                Qx(i,d,e,s) += Qg(c,d,e,f) *
                  EXP( -( score_interior_asymmetry(l1, l2) +
                          score_interior_mismatch(f,c,f+1,c-1) )/RT );
              }
            }
          }
        }
      }
    }
  }

  for (int d=i+1; d<=j-5; ++d)
  {
    for (int e=d+4; e<=j-1; ++e) 
    {
      if (allow_paired(d,e))
      {
        // convert Qx into interior loop energies
        if (l>=17 && allow_paired(i,j))
        {
          for (int s=8; s<=l-9; ++s)
          {
            Qg(i,d,e,j) += Qx(i,d,e,s) *
              EXP( -score_interior_mismatch(i,j,i+1,j-1)/RT );
          }
        }

        // extend loops for future use
        if (i!=0 && j!=N-1)
        {
          for (int s=8; s<=l-9; ++s)
          {
            Qx2(i-1,d,e,s+2) = Qx(i,d,e,s) *
              EXP( -(score_loop(s+2)-score_loop(s))/RT );
          }
        }
      
        if (allow_paired(i,j))
        {
          // Add small inextensible interior loops to Qg as special cases
          for (int l1=0; l1<=std::min(3,d-i-2); ++l1)
          {
            int c=i+l1+1;
            for (int l2=0; l2<=std::min(3,j-e-2); ++l2)
            {
              int f=j-l2-1;
              if (allow_paired(c,f))
              {
                Qg(i,d,e,j) += Qg(c,d,e,f) *
                  EXP( -score_interior(i,c,f,j,true)/RT );
              }
            }
          }
          // Add bulge loops and large asymmetric loops as special cases
          for (int l1=0; l1<=std::min(3,d-i-2); ++l1) // cases l1=0,1,2,3, l2>=4
          {
            int c=i+l1+1;
            for (int l2=4; l2<=j-e-2; ++l2)
            {
              int f=j-l2-1;
              if (allow_paired(c,f))
              {
                Qg(i,d,e,j) += Qg(c,d,e,f) *
                  EXP( -score_interior(i,c,f,j,true)/RT );
              }
            }
          }
          for (int l2=0; l2<=std::min(3,j-e-2); ++l2)
          {
            int f=j-l2-1;
            for (int l1=4; l1<=d-i-2; ++l1)
            {
              int c=i+l1+1;
              if (allow_paired(c,f))
              {
                Qg(i,d,e,j) += Qg(c,d,e,f) *
                  EXP( -score_interior(i,c,f,j,true)/RT );
              }
            }
          }
        }
      }
    }
  }
}

template < class PF_TYPE >
void
Nupack<PF_TYPE>::
calculate_posterior()
{
  P.resize(N);    P.fill(0.0);
  Pb.resize(N);   Pb.fill(0.0);
  Pm.resize(N);   Pm.fill(0.0);
  Pp.resize(N);   Pp.fill(0.0);
  Pz.resize(N);   Pz.fill(0.0);
  Pbg.resize(N);  Pbg.fill(0.0);
  Pg.resize(N);   Pg.fill(0.0);
  Pgl.resize(N);  Pgl.fill(0.0);
  Pgr.resize(N);  Pgr.fill(0.0);
  Pgls.resize(N); Pgls.fill(0.0);
  Pgrs.resize(N); Pgrs.fill(0.0);
  P(0,N-1)=1.0;
  DPTableX<PF_TYPE> Qx, Qx1, Qx2;
  Qx.resize(N,N); Qx1.resize(N-1,N); Qx2.resize(N-2,N);
  DPTableX<DBL_TYPE> Px, Px1, Px2;
  Px.resize(N,N); Px1.resize(N-1,N); Px2.resize(N-2,N);
#ifdef ENABLE_RECALCULATE
  DPTableX<float> Precx, Precx1, Precx2;
  Precx.resize(N,N); Precx1.resize(N-1,N); Precx2.resize(N-2,N);
#endif

  for (int l=N; l>=1; --l)
  {
    Qx.swap(Qx1);
    Qx1.swap(Qx2);
    Qx2.resize(l-3, N);
    Qx2.fill(0.0);
    Px.swap(Px1);
    Px1.swap(Px2);
    Px2.resize(l-3, N);
    Px2.fill(0.0);
#ifdef ENABLE_RECALCULATE
    Precx.swap(Precx1);
    Precx1.swap(Precx2);
    Precx2.resize(l-3, N);
    Precx2.fill(0.0);
#endif

    for (int i=0; i+l<=N; ++i)
    {
      int j=i+l-1;

      // P, Pm, Pz recursions
      for (int d=i; d<=j-4; ++d)
      {
        for (int e=d+4; e<=j; ++e)
        {
          if (allow_paired(d,e) && wc_pair(d,e))
          {
            if (P(i,j)>0.0)
            {
              DBL_TYPE p = Q(i,d-1) * Qb(d,e) *
                EXP( -( score_at_penalty(d,e) +
                        score_dangle(e+1,j) )/RT )
                / Q(i,j) * P(i,j); 
              P(i,d-1) += p;
              Pb(d,e) += p;
              assert(!std::isnan(p));
            }
            
            if (i!=0 && j!=N-1)
            {
              if (Pm(i,j)>0.0)
              {
                DBL_TYPE p = Qb(d,e) *
                  EXP( -( score_multiloop_paired(1,false) +
                          score_multiloop_unpaired(d-i + j-e,false) +
                          score_at_penalty(d,e) +
                          score_dangle(e+1,j) +
                          score_dangle(i,d-1) )/RT )
                  / Qm(i,j) * Pm(i,j);
                Pb(d,e) += p;

                if (d>=i+5)
                {
                  p = Qm(i,d-1) * Qb(d,e) *
                    EXP( -( score_multiloop_paired(1,false) +
                            score_multiloop_unpaired(j-e,false) +
                            score_at_penalty(d,e) +
                            score_dangle(e+1,j) )/RT )
                    / Qm(i,j) * Pm(i,j);
                  Pm(i,d-1) += p;
                  Pb(d,e) += p;
                  assert(!std::isnan(p));
                }
              }

              if (Pz(i,j)>0.0)
              {
                DBL_TYPE p = Qz(i,d-1) * Qb(d,e) *
                  EXP( -( score_pk_paired(1) +
                          score_pk_unpaired(j-e) +
                          score_at_penalty(d,e) +
                          score_dangle(e+1,j) )/RT )
                  / Qz(i,j) * Pz(i,j);
                Pz(i,d-1) += p;
                Pb(d,e) += p;
                assert(!std::isnan(p));
              }
            }
          }
        }
      }

      for (int d=i; d<=j-5; ++d)
      {
        for (int e=d+5; e<=j; ++e)
        {
          if (P(i,j)>0.0)
          {
            DBL_TYPE p = Q(i,d-1) * Qp(d,e) *
              EXP( -( score_pk() +
                      score_dangle(e+1,j) )/RT )
              / Q(i,j) * P(i,j);
            P(i,d-1) += p;
            Pp(d,e) += p;
            assert(!std::isnan(p));
          }

          if (i!=0 && j!=N-1)
          {
            if (Pm(i,j)>0.0)
            {
              DBL_TYPE p = Qp(d,e) *
                EXP( -( score_pk_multiloop() +
                        score_multiloop_paired(2,false) +
                        score_multiloop_unpaired(d-i + j-e,false) +
                        score_dangle(e+1,j) +
                        score_dangle(i,d-1) )/RT )
                / Qm(i,j) * Pm(i,j);
              Pp(d,e) += p;
              assert(!std::isnan(p));

              if (d>=i+5)
              {
                p = Qm(i,d-1) * Qp(d,e) *
                  EXP( -( score_pk_multiloop() +
                          score_multiloop_paired(2,false) +
                          score_multiloop_unpaired(j-e,false) +
                          score_dangle(e+1,j) )/RT )
                  / Qm(i,j) * Pm(i,j);
                Pm(i,d-1) += p;
                Pp(d,e) += p;
                assert(!std::isnan(p));
              }
            }

            if (Pz(i,j)>0.0)
            {
              DBL_TYPE p = Qz(i,d-1) * Qp(d,e) *
                EXP( -( score_pk_pk() +
                        score_pk_paired(2) +
                        score_pk_unpaired(j-e) +
                        score_dangle(e+1,j) )/RT )
                / Qz(i,j) * Pz(i,j);
              Pz(i,d-1) += p;
              Pp(d,e) += p;
              assert(!std::isnan(p));
            }
          }
        }
      }

      // Pp recursion
      // case 1: both Qg are exactly 1 pair
      // first case is exactly 1 pair per Og
      if (j-i>4)
      {
        int a=i;
        int f=j;
        for (int b=a+1; b<=j-4; ++b)
        {
          if (allow_paired(b,j) && wc_pair(b,j))
          {
            int c=b;
            for (int d=std::max(c+1,a+4); d<=j-1; ++d)
            {
              if (allow_paired(a,d) && wc_pair(a,d) && Pp(i,j)>0.0)
              {
                int e=d;
                DBL_TYPE p = Qg(i,a,d,e) * Qg(b,c,f,j) * Qz(e+1,f-1) * Qz(c+1,d-1) * Qz(a+1,b-1) *
                  EXP( -( score_pk_paired(2) +
                          score_at_penalty(a,d) +
                          score_at_penalty(c,f) +
                          score_at_penalty(i,e) +
                          score_at_penalty(b,j) )/RT )
                  / Qp(i,j) * Pp(i,j);
                Pg(i,a,d,e) += p;
                Pg(b,c,f,j) += p;
                Pz(e+1,f-1) += p;
                Pz(c+1,d-1) += p;
                Pz(a+1,b-1) += p;
                assert(!std::isnan(p));
              }
            }
          }
        }
      }

      if (j-i>6)
      {
        // case 2 left Og is exactly 1 pair, right is 2+
        for (int d=i+1; d<=j-6; ++d)
        {
          if (allow_paired(d,j) && wc_pair(d,j))
          {
            for (int e=std::max(d+2,i+4); e<=j-2; ++e)
            {
              int f=e;
              if (allow_paired(i,f) && wc_pair(i,f) && Pp(i,j)>0.0)
              {
                DBL_TYPE p = Qg(i,i,e,f) * Qz(i+1,d-1) * Qgr(d,e-1,f+1,j) *
                  EXP( -( score_pk_paired(1) +
                          score_at_penalty(d,j) +
                          score_at_penalty(i,f)*2 )/RT )
                  / Qp(i,j) * Pp(i,j);
                Pg(i,i,e,f) += p;
                Pz(i+1,d-1) += p;
                Pgr(d,e-1,f+1,j) += p;
                assert(!std::isnan(p));
              }
            }
          }
        }

        // case 2 left Qg is 2+ pairs, right is 1
        for (int d=i+2; d<=j-4; ++d)
        {
          if (allow_paired(d,j) && wc_pair(d,j))
          {
            for (int e=std::max(d+1,i+4); e<=j-2; ++e)
            {
              for (int f=e+1; f<=j-1; ++f)
              {
                if (allow_paired(i,f) && wc_pair(i,f) && Pp(i,j)>0.0)
                {
                  DBL_TYPE p = Qgl(i,d-1,e,f) * Qg(d,d,j,j) * Qz(d+1,e-1) * Qz(f+1,j-1) *
                    EXP( -( score_pk_paired(1) +
                            score_at_penalty(d,j)*2 +
                            score_at_penalty(i,f) )/RT )
                    / Qp(i,j) * Pp(i,j);
                  Pgl(i,d-1,e,f) += p;
                  Pg(d,d,j,j) += p;
                  Pz(d+1,e-1) += p;
                  Pz(f+1,j-1) += p;
                  assert(!std::isnan(p));
                }
              }
            }
          }
        }
      }

      // otherwise
      if (j-i>7)
      {
        for (int d=i+2; d<=j-4; ++d)
        {
          if (allow_paired(d,j) && wc_pair(d,j))
          {
            for (int e=std::max(d+2,i+5); e<=j-3; ++e)
            {
              for (int f=e+1; f<=j-2; ++f)
              {
                if (allow_paired(i,f) && wc_pair(i,f) && Pp(i,j)>0.0)
                {
                  DBL_TYPE p = Qgl(i,d-1,e,f) * Qgr(d,e-1,f+1,j) *
                    EXP( -( score_at_penalty(d,j) +
                            score_at_penalty(i,j) )/RT )
                    / Qp(i,j) * Pp(i,j);
                  Pgl(i,d-1,e,f) += p;
                  Pgr(d,e-1,f+1,j) += p;
                  assert(!std::isnan(p));
                }
              }
            }
          }
        }
      }

      // Pgr recursion
      for (int d=i+1; d<=j-3; ++d)
      {
        for (int e=d+2; e<=j-1; ++e)
        {
          for (int f=e; f<=j-1; ++f)
          {
            if (Pgr(i,d,e,j)>0.0)
            {
              DBL_TYPE p = Qgl(i,d,f,j) * Qz(e,f-1) / Qgr(i,d,e,j) * Pgr(i,d,e,j);
              Pgl(i,d,f,j) += p;
              Pz(e,f-1) += p;
              assert(!std::isnan(p));
            }
          }
        }
      }
      
      // Pgl recursions
      for (int d=i+1; d<=j-5; ++d)
      {
        for (int f=d+4; f<=j-1; ++f)
        {
          if (allow_paired(d,f) && wc_pair(d,f))
          {
            for (int e=d; e<=f-2; ++e) // f-3???
            {
              if (Pgl(i,e,f,j)>0.0)
              {
                DBL_TYPE p = Qg(i,d,f,j) * Qz(d+1,e) *
                  EXP( -( score_pk_paired(1) +
                          score_at_penalty(d,f) )/RT )
                  / Qgl(i,e,f,j) * Pgl(i,e,f,j);
                Pg(i,d,f,j) += p;
                Pz(d+1,e) += p;
                Pbg(d,f) += p; // Pbg inner gap-spanning base-pairing prob
                assert(!std::isnan(p));
              }
            }
          }
        }
      }

      // Pgrs recursion
      for (int d=i+1; d<=j-10; ++d)
      {
        for (int e=d+4; e<=j-6; ++e)
        {
          if (allow_paired(d,e))
          {
            for (int f=e+1; f<=j-5; ++f)
            {
              if (allow_paired(i,f) && wc_pair(i,f) && Pgrs(i,d,e,j)>0.0)
              {
                DBL_TYPE p = Qg(i,d,e,f) * Qm(f+1,j) *
                  EXP( -( score_multiloop_paired(1,true) +
                          score_at_penalty(i,f) )/RT )
                  / Qgrs(i,d,e,j) * Pgrs(i,d,e,j);
                Pg(i,d,e,f) += p;
                Pm(f+1,j) += p;
                assert(!std::isnan(p));
              }
            }
          }
        }
      }

      // Pgls recursion
      for (int c=i+5; c<=j-6; ++c)
      {
        if (allow_paired(c,j) && wc_pair(c,j))
        {
          for (int d=c+1; d<=j-5; ++d)
          {
            for (int e=d+4; e<=j-1; ++e)
            {
              if (allow_paired(d,e) && Pgls(i,d,e,j)>0.0)
              {
                DBL_TYPE p = Qm(i,c-1) * Qg(c,d,e,j) *
                  EXP( -( score_multiloop_paired(1,true) +
                          score_at_penalty(c,j) )/RT )
                  / Qgls(i,d,e,j) * Pgls(i,d,e,j);
                Pm(i,c-1) += p;
                Pg(c,d,e,j) += p;
                assert(!std::isnan(p));
              }
            }
          }
        }
      }

      // Pg recursion
      fastiloops_pr(i, j, Qg, Qx, Qx2, Pg, Px, Px2);

      if (allow_paired(i,j) && wc_pair(i,j))
      {
        // case 2: multiloop left
        for (int d=i+6; d<=j-5; ++d)
        {
          for (int e=d+4; e<=j-1; ++e)
          {
            if (allow_paired(d,e) && wc_pair(d,e) && Pg(i,d,e,j)>0.0)
            {
              DBL_TYPE p = Qm(i+1,d-1) *
                EXP( -( score_multiloop(true) +
                        score_multiloop_paired(2,true) +
                        score_multiloop_unpaired(j-e-1,true) + 
                        score_at_penalty(i,j) +
                        score_at_penalty(d,e) +
                        score_dangle(e+1, j-1) )/RT )
                / Qg(i,d,e,j) * Pg(i,d,e,j);
              Pm(i+1,d-1) += p;
              assert(!std::isnan(p));
            }
          }
        }

        // case 3: multiloop right
        for (int d=i+1; d<=j-10; ++d)
        {
          for (int e=d+4; e<=j-6; ++e)
          {
            if (allow_paired(d,e) && wc_pair(d,e) && Pg(i,d,e,j)>0.0)
            {
              DBL_TYPE p = Qm(e+1,j-1) *
                EXP( -( score_multiloop(true) +
                        score_multiloop_paired(2,true) +
                        score_multiloop_unpaired(d-i-1,true) +
                        score_at_penalty(i,j) +
                        score_at_penalty(d,e) +
                        score_dangle(i+1, d-1) )/RT )
                / Qg(i,d,e,j) * Pg(i,d,e,j);
              Pm(e+1,j-1) += p;
              assert(!std::isnan(p));
            }
          }
        }

        // case 4: multiloop both sides
        for (int d=i+6; d<=j-10; ++d)
        {
          for (int e=d+4; e<=j-6; ++e)
          {
            if (allow_paired(d,e) && wc_pair(d,e) && Pg(i,d,e,j)>0.0)
            {
              DBL_TYPE p = Qm(i+1,d-1) * Qm(e+1,j-1) *
                EXP( -( score_multiloop(true) +
                        score_multiloop_paired(2,true) +
                        score_at_penalty(i,j) +
                        score_at_penalty(d,e) )/RT )
                / Qg(i,d,e,j) * Pg(i,d,e,j);
              Pm(i+1,d-1) += p;
              Pm(e+1,j-1) += p;
              assert(!std::isnan(p));
            }
          }
        }

        // case 5: interior loop + multi left
        for (int d=i+7; d<=j-6; ++d)
        {
          for (int e=d+4; e<=j-2; ++e)
          {
            if (allow_paired(d,e))
            {
              for (int f=e+1; f<=j-1; ++f)
              {
                if (Pg(i,d,e,j)>0.0)
                {
                  DBL_TYPE p = Qgls(i+1,d,e,f) *
                    EXP( -( score_multiloop(true) +
                            score_multiloop_paired(1,true) +
                            score_multiloop_unpaired(j-f-1,true) +
                            score_at_penalty(i,j) +
                            score_dangle(f+1,j-1) )/RT )
                    / Qg(i,d,e,j) * Pg(i,d,e,j);
                  Pgls(i+1,d,e,f) += p;
                  assert(!std::isnan(p));
                }
              }
            }
          }
        }

        // case 6: interior loop + multi right
        for (int d=i+2; d<=j-11; ++d)
        {
          for (int e=d+4; e<=j-7; ++e)
          {
            if (allow_paired(d,e))
            {
              for (int c=i+1; c<=d-1; ++c)
              {
                if (Pg(i,d,e,j)>0.0)
                {
                  DBL_TYPE p = Qgrs(c,d,e,j-1) *
                    EXP( -( score_multiloop(true) +
                            score_multiloop_paired(1,true) +
                            score_multiloop_unpaired(c-i-1,true) +
                            score_at_penalty(i,j) +
                            score_dangle(i+1,c-1) )/RT )
                    / Qg(i,d,e,j) * Pg(i,d,e,j);
                  Pgrs(c,d,e,j-1) += p;
                  assert(!std::isnan(p));
                }
              }
            }
          }
        }

        // case 7: interior loop + multi both sides
        for (int d=i+7; d<=j-11; ++d)
        {
          for (int e=d+4; e<=j-7; ++e)
          {
            if (allow_paired(d,e))
            {
              for (int c=i+6; c<=d-1; ++c)
              {
                if (Pg(i,d,e,j)>0.0)
                {
                  DBL_TYPE p = Qm(i+1,c-1) * Qgrs(c,d,e,j-1) *
                    EXP( -( score_multiloop(true) +
                            score_multiloop_paired(1,true) +
                            score_at_penalty(i,j) )/RT )
                    / Qg(i,d,e,j) * Pg(i,d,e,j);
                  Pm(i+1,c-1) += p;
                  Pgrs(c,d,e,j-1) += p;
                  assert(!std::isnan(p));
                }
              }
            }
          }
        }
      }

      // Pbg outer gap-spanning base-pairing prob
      for (int d=i+1; d<=j-5; ++d)
        for (int e=d+4; e<=j-1; ++e)
          Pbg(i,j) += Pg(i,d,e,j);

      // Pb recursion
      if (allow_paired(i,j))
      {
        for (int d=i+1; d<=j-5; ++d) // all possible rightmost pairs d-e
        {
          for (int e=d+4; e<=j-1; ++e)
          {
            if (allow_paired(d,e))
            {
              if (Pb(i,j)>0.0)
              {
                DBL_TYPE p = Qb(d,e) *
                  EXP( -score_interior(i,d,e,j,false)/RT )
                  / Qb(i,j) * Pb(i,j);
                Pb(d,e) += p;
                assert(!std::isnan(p));

                if (d>=i+6 && wc_pair(d,e) && wc_pair(i,j))
                {
                  DBL_TYPE p = Qm(i+1,d-1) * Qb(d,e) *
                    EXP( -( score_multiloop(false) +
                            score_multiloop_paired(2,false) +
                            score_multiloop_unpaired(j-e-1,false) +
                            score_at_penalty(i,j) +
                            score_at_penalty(d,e) +
                            score_dangle(e+1,j-1) )/RT )
                    / Qb(i,j) * Pb(i,j);
                  Pm(i+1,d-1) += p;
                  Pb(d,e) += p;
                  assert(!std::isnan(p));
                }
              }
            }
          }
        }

        if (wc_pair(i,j))
        {
          for (int d=i+1; d<=j-6; ++d) // all possible rightmost pseudoknots filling [d,e]
          {
            for (int e=d+5; e<=j-1; ++e)
            {
              if (Pb(i,j)>0.0)
              {
                DBL_TYPE p;
                p = Qp(d,e) *
                  EXP( -( score_multiloop(false) +
                          score_pk_multiloop() +
                          score_multiloop_paired(3,false) +
                          score_multiloop_unpaired(j-e-1 + d-i-1,false) +
                          score_at_penalty(i,j) +
                          score_dangle(e+1,j-1) +
                          score_dangle(i+1,d-1) )/RT )
                  / Qb(i,j) * Pb(i,j);
                Pp(d,e) += p;
                assert(!std::isnan(p));
              
                p = Qm(i+1,d-1) * Qp(d,e) *
                  EXP( -( score_multiloop(false) +
                          score_pk_multiloop() +
                          score_multiloop_paired(3,false) +
                          score_multiloop_unpaired(j-e-1,false) +
                          score_at_penalty(i,j) +
                          score_dangle(e+1,j-1) )/RT )
                  / Qb(i,j) * Pb(i,j);
                Pm(i+1,d-1) += p;
                Pp(d,e) += p;
                assert(!std::isnan(p));
              }
            }
          }
        }
      }
    }
  }
}

template < class PF_TYPE >
void
Nupack<PF_TYPE>::
fastiloops_pr(int i, int j,
#ifdef ENABLE_RECALCULATE
              DPTableX<float>& precX, DPTableX<float>& precX2,
#endif
              DPTable4<PF_TYPE>& Qg, DPTableX<PF_TYPE>& Qx, DPTableX<PF_TYPE>& Qx2,
              DPTable4<DBL_TYPE>& Pg, DPTableX<DBL_TYPE>& Px, DPTableX<DBL_TYPE>& Px2)
{
  int l=j-i+1;

  if (allow_paired(i,j))
  {
    for (int d=i+1; d<=j-5; ++d)
    {
      for (int e=d+4; e<=j-1; ++e) 
      {
        if (allow_paired(d,e))
        {
          // Add small inextensible interior loops to Qg as special cases
          for (int l1=0; l1<=std::min(3,d-i-2); ++l1)
          {
            int c=i+l1+1;
            for (int l2=0; l2<=std::min(3,j-e-2); ++l2)
            {
              int f=j-l2-1;
              if (allow_paired(c,f) && Pg(i,d,e,j)>0.0)
              {
                DBL_TYPE p = Qg(c,d,e,f) *
                  EXP( -score_interior(i,c,f,j,true)/RT )
                  / Qg(i,d,e,j) * Pg(i,d,e,j);
                Pg(c,d,e,f) += p;
                assert(!std::isnan(p));
              }
            }
          }
          // Add bulge loops and large asymmetric loops as special cases
          for (int l1=0; l1<=std::min(3,d-i-2); ++l1) // cases l1=0,1,2,3, l2>=4
          {
            int c=i+l1+1;
            for (int l2=4; l2<=j-e-2; ++l2)
            {
              int f=j-l2-1;
              if (allow_paired(c,f) && Pg(i,d,e,j)>0.0)
              {
                DBL_TYPE p = Qg(c,d,e,f) *
                  EXP( -score_interior(i,c,f,j,true)/RT )
                  / Qg(i,d,e,j) * Pg(i,d,e,j);
                Pg(c,d,e,f) += p;
                assert(!std::isnan(p));
              }
            }
          }
          for (int l2=0; l2<=std::min(3,j-e-2); ++l2)
          {
            int f=j-l2-1;
            for (int l1=4; l1<=d-i-2; ++l1)
            {
              int c=i+l1+1;
              if (allow_paired(c,f) && Pg(i,d,e,j)>0.0)
              {
                DBL_TYPE p = Qg(c,d,e,f) *
                  EXP( -score_interior(i,c,f,j,true)/RT )
                  / Qg(i,d,e,j) * Pg(i,d,e,j);
                Pg(c,d,e,f) += p;
                assert(!std::isnan(p));
              }
            }
          }
        }
      }
    }
  }

  // Add cases that are at an end with l1>=4, l2>=4
  if ((i==0 || j==N-1) && l>=17)
  {
    for (int d=i+6; d<=j-10; ++d)
    {
      for (int e=d+4; e<=j-6; ++e)
      {
        if (allow_paired(d,e))
        {
          for (int c=i+5; c<=d-1; ++c)
          {
            for (int f=e+1; f<=j-5; ++f)
            {
              if (allow_paired(c,f))
              {
                int l1=c-i-1;
                int l2=j-f-1;
                int s=l1+l2;
                Qx(i,d,e,s) += Qg(c,d,e,f) *
                  EXP( -( score_interior_asymmetry(l1, l2) +
                          score_interior_mismatch(f,c,f+1,c-1) )/RT );
              }
            }
          }
        }
      }
    }
  }

  // Use Qx to finish calculation of Px
  if (allow_paired(i,j))
  {
    for (int d=i+1; d<=j-5; ++d)
    {
      for (int e=d+4; e<=j-1; ++e)
      {
        if (allow_paired(d,e))
        {
          for (int s=8; s<=l-9; ++s)
          {
            DBL_TYPE p = Qx(i,d,e,s) *
              EXP( -score_interior_mismatch(i,j,i+1,j-1)/RT )
              / Qg(i,d,e,j) * Pg(i,d,e,j);
            Px(i,d,e,s) += p;
            assert(!std::isnan(p));
          }
        }
      }
    }
  }

  // Calculate Pg contribution using Qx and Px
  if (l>=17)
  {
    for (int d=i+6; d<=j-10; ++d)
    {
      for (int e=d+4; e<=j-6; ++e)
      {
        if (allow_paired(d,e))
        {
          int l1=4;               // explicitly add in terms for l1=4, l2>=4
          int c=i+l1+1;
          for (int l2=4; l2<=j-e-2; ++l2)
          {
            int s=l1+l2;
            int f=j-l2-1;
            if (allow_paired(c,f) && Qx(i,d,e,s)>0.0)
            {
              PF_TYPE temp = Qg(c,d,e,f) *
                EXP( -( score_interior_asymmetry(l1, l2) +
                        score_interior_mismatch(f,c,f+1,c-1) )/RT );
              DBL_TYPE p = temp / Qx(i,d,e,s) * Px(i,d,e,s);
              Pg(c,d,e,f) += p;
              Px(i,d,e,s) -= p;
              assert(!std::isnan(p));
              if (temp>Qx(i,d,e,s))
              {
                temp = 0.0;
                int l1min = 5;
                int l2min = 4;
                for (int c=i+l1min+1; c<=d-1; ++c)
                {
                  int f=c-i+j-3;
                  if (j-f-1>l2min && f>=e+1 && allow_paired(c,f))
                  {
                    temp += Qg(c,d,e,f) *
                      EXP ( -score_interior(i,j,c,f,true)/RT );
                  }
                }
                Qx(i,d,e,s) = temp;
              }
              else
              {
                Qx(i,d,e,s) -= temp;
              }
            }
          }

          if (d>=i+7)
          {
            int l2=4;             // explicitly add in terms of l1>=5, l2=4
            int f=j-l2-1;
            for (int l1=5; l1<=d-i-2; ++l1)
            {
              int s=l1+l2;
              int c=i+l1+1;
              if (allow_paired(c,f) && Qx(i,d,e,s)>0.0)
              {
                PF_TYPE temp = Qg(c,d,e,f) *
                  EXP( -( score_interior_asymmetry(l1, l2) +
                          score_interior_mismatch(f,c,f+1,c-1) )/RT );
                DBL_TYPE p = temp / Qx(i,d,e,s) * Px(i,d,e,s);
                Pg(c,d,e,f) += p;
                Px(i,d,e,s) -= p;
                assert(!std::isnan(p));
                if (temp>Qx(i,d,e,s))
                {
                  temp = 0.0;
                  int l1min = 5;
                  int l2min = 5;
                  for (int c=i+l1min+1; c<=d-1; ++c)
                  {
                    int f=c-i+j-4;
                    if (j-f-1>l2min && f>=e+1 && allow_paired(c,f))
                    {
                      temp += Qg(c,d,e,f) *
                        EXP ( -score_interior(i,j,c,f,true)/RT );
                    }
                  }
                  Qx(i,d,e,s) = temp;
                }
                else
                {
                  Qx(i,d,e,s) -= temp;
                }
              }
            }
          }
        }

        // Store partial values for Qx2 and Px2
        for (int s=10; s<=l-9; ++s)
        {
          Qx2(i+1,d,e,s-2) = Qx(i,d,e,s) *
            EXP( -(score_loop(s+2)-score_loop(s))/RT );
          Px2(i+1,d,e,s-2) = Px(i,d,e,s);
        }
      }
    }
  }
}

template < class PF_TYPE >
void
Nupack<PF_TYPE>::
get_posterior(std::vector<float>& bp1, std::vector<float>& bp2, std::vector<int>& offset) const
{
  bp1.resize((N+1)*(N+2)/2);
  bp2.resize((N+1)*(N+2)/2);
  offset.resize(N+1);
  for (int i=0; i<=N; ++i)
    offset[i] = i*((N+1)+(N+1)-i-1)/2;
  for (int i=0; i!=N-1; ++i)
    for (int j=i+1; j!=N; ++j)
    {
      bp1[offset[i+1]+(j+1)] = Pb(i,j);
      bp2[offset[i+1]+(j+1)] = Pbg(i,j);
    }
}

template < class PF_TYPE >
void
Nupack<PF_TYPE>::
get_posterior(std::vector<float>& bp, std::vector<int>& offset) const
{
  bp.resize((N+1)*(N+2)/2);
  offset.resize(N+1);
  for (int i=0; i<=N; ++i)
    offset[i] = i*((N+1)+(N+1)-i-1)/2;
  for (int i=0; i!=N-1; ++i)
    for (int j=i+1; j!=N; ++j)
      bp[offset[i+1]+(j+1)] = Pb(i,j) + Pbg(i,j);
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_hairpin(int i, int j) const
{
  energy_t e=0.0;
  bool polyC = true;
  for (int k=i+1; k<j; ++k)
  {
    if (seq[k]!=BASE_C)
    {
      polyC = false;
      break;
    }
  }

  int size=j-i-1;
#if 0
  if (size<3) return INF;
  if (!allow_paired(i,j)) return INF;
#else
  assert(size>=3);
  assert(allow_paired(i,j));
#endif

  e += size<=30 ?
    hairpin37[size-1] : 
    hairpin37[30 - 1] + loop_greater30*LOG(size/30.0);

  if (size==3)
  {
    e += score_at_penalty(i,j);
    e += triloop37[seq[i]-1][seq[i+1]-1][seq[i+2]-1][seq[j-1]-1][seq[j]-1];
    if (polyC) e += polyC_penalty;
    if (seq[i+1]==BASE_G && seq[i+2]==BASE_G && seq[j-1]==BASE_G)
      e += hairpin_GGG;
  }
  else if (size==4)
  {
    e += tloop37[seq[i]-1][seq[i+1]-1][seq[i+2]-1][seq[j-2]-1][seq[j-1]-1][seq[j]-1];
    e += mismatch_hairpin37[seq[i+1]-1][seq[j-1]-1][pair_type(i,j)];
    if (polyC) e += polyC_slope*size + polyC_int;
  }
  else /*if (size>4)*/
  {
    e += mismatch_hairpin37[seq[i+1]-1][seq[j-1]-1][pair_type(i,j)];
    if (polyC) e += polyC_slope*size + polyC_int;
  }
  return e;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_loop(int l) const
{
  return l<=30 ?
    interior37[l-1] :
    interior37[30-1]+loop_greater30*LOG(l/30.0);
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_interior(int i, int h, int m, int j, bool pk) const
{
  int l1 = h - i - 1;
  int l2 = j - m - 1;
  int size = l1 + l2;
  energy_t e = 0;

  // helix
  if (size==0)
  {
    return stack37[pair_type(i,j)][pair_type(h,m)]
      * (pk ? pk_stack_span : 1.0);
  }
  
  // bulge
  else if (l1==0 || l2==0)
  {
    e += size<=30 ?
      bulge37[size-1] :
      bulge37[30-1] + loop_greater30*LOG(size/30.0);

    if (l1+l2==1)           //single bulge...treat as a stacked region
    {
      e += stack37[pair_type(i,j)][pair_type(h,m)];
      e -= SALT_CORRECTION;
    }
    else
    {
      e += score_at_penalty(i,j);
      e += score_at_penalty(h,m);
    }
  }

  // interior loop
  else if (l1>0 && l2>0)
  {
    int asymmetry = std::abs(l1-l2);
    if (asymmetry>1 || size>4)
    {
      e += score_interior_asymmetry(l1, l2);
      if (l1>1 && l2>1)
      {
        e += score_interior_mismatch(m, h, m+1, h-1);
        e += score_interior_mismatch(i, j, i+1, j-1);
      }
      else if (l1==1 || l2==1)
      {
#if 1                           // assume AA terminal mismatch?
        e += score_interior_mismatch(m, h);
        e += score_interior_mismatch(i, j);
#else
        e += score_interior_mismatch(m, h, m+1, h-1);
        e += score_interior_mismatch(i, j, i+1, j-1);
#endif
      }
      else
      {
        assert(!"unclassified interior loop");
        exit(1);
      }
    }
    else if (l1==1 && l2==1)
      e += int11_37[pair_type(i,j)][pair_type(h,m)][seq[i+1]-1][seq[j-1]-1];
    else if (l1==2 && l2==2)
      e += int22_37[pair_type(i,j)][pair_type(h,m)][seq[i+1]-1][seq[j-1]-1][seq[i+2]-1][seq[j-2]-1];
    else if (l1==1 && l2==2)
      e += int21_37[pair_type(i,j)][seq[j-2]-1][seq[i+1]-1][pair_type(h,m)][seq[j-1]-1];
    else if (l1==2 && l2==1)
      e += int21_37[pair_type(m,h)][seq[i+1]-1][seq[j-1]-1][pair_type(j,i)][seq[i+2]-1];
    else
    {
      assert(!"error in tabulated interior loop");
      exit(1);
    }
  }
  else
  {
    assert(!"improperly classifed interior loop");
    exit(1);
  }
  return e * (pk ? pk_interior_span : 1.0);
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_interior_mismatch(int i, int j, int k, int l) const
{
  return mismatch_interior37[seq[k]-1][seq[l]-1][pair_type(i,j)];
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_interior_mismatch(int i, int j) const
{
  return mismatch_interior37[BASE_N][BASE_N][pair_type(i,j)];
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_interior_asymmetry(int l1, int l2) const
{
  energy_t e=0.0;
  int size = l1+l2;
  int asymmetry = std::abs(l1-l2);
  e += size<=30 ?
    interior37[size-1] :
    interior37[30-1] + loop_greater30*LOG(size/30.0);

  // asymmetry penalty
  e += std::min(max_asymmetry,
                asymmetry * asymmetry_penalty[std::min(4,std::min(l1,l2))-1]);

  return e;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_multiloop(bool pk) const
{
  return pk ? multiloop_penalty_pk : multiloop_penalty ;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_multiloop_paired(int n, bool pk) const
{
  return (pk ? multiloop_paired_penalty_pk : multiloop_paired_penalty) * n;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_multiloop_unpaired(int n, bool pk) const
{
  return (pk ? multiloop_unpaired_penalty_pk : multiloop_unpaired_penalty) * n;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_at_penalty(int i, int j) const
{
  return pair_type(i,j)==PAIR_AU || pair_type(i,j)==PAIR_UA ? at_penalty : 0;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_dangle(int i, int j) const
{
  energy_t d5=0.0, d3=0.0;

#if 0
  if (j!=N-1)
    d3 = dangle3_37[pair_type(i-1,j+1)][seq[j]-1];
  if (i!=0)
    d5 = dangle5_37[pair_type(i-1,j+1)][seq[i]-1];

#else  // Is this correct? This implementation is the same as the original nupack.

  //if( DANGLETYPE != 2) {
    if( (j == i - 1) || (j==-1 && i>0)) {
      return 0.0;
    }
  //}
  //else if( (j==-1 && i>0) || (j == i - 1 && (i == 0 || j == seqlength - 1)) ) {
    //return 1.0;
  //}

  if( (j==-1 && i>0) || (j==i-1 && (i==0 || j==N-1)) )
    return 0.0;

  if (j!=N-1)
    d3 = dangle3_37[3-pair_type(j+1)][seq[j]-1];
  if (i!=0)
    d5 = dangle5_37[pair_type(i-1)][seq[i]-1];
#endif

  if (i==j && i!=0 && j!=N-1 /* && DANGLETYPE!=2 */)
    return std::min(d3, d5);
  else
    return d3+d5;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_pk() const
{
  return pk_penalty;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_pk_multiloop() const
{
  return pk_multiloop_penalty;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_pk_pk() const
{
  return pk_pk_penalty;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_pk_paired(int n) const
{
  return pk_paired_penalty*n;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_pk_unpaired(int n) const
{
  return pk_unpaired_penalty*n;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_pk_band(int n) const
{
  return pk_band_penalty*n;
}

// instantiation
template class Nupack<long double>;

#if 0
int
main(int argc, char* argv[])
{
  std::string seq("GGGCUGUUUUUCUCGCUGACUUUCAGCCCCAAACAAAAAAUGUCAGCA");
  Nupack<long double> nu;
  nu.load_parameters(argv[1]);
  nu.load_sequence(seq);
  std::cout << nu.calculate_partition_function() << std::endl;
  nu.calculate_posterior();
  std::vector<float> bp1, bp2;
  std::vector<int> offset;
  nu.get_posterior(bp1, /*bp2,*/ offset);
  for (int i=1; i<=(int)seq.size(); ++i)
  {
    std::cout << i << " " << seq[i-1] << " ";
    for (int j=i+1; j<=(int)seq.size(); ++j)
      if (bp1[offset[i]+j]/*+bp2[offset[i]+j]*/>1e-5)
        std::cout << j << ":" << bp1[offset[i]+j] /*<< "," << bp2[offset[i]+j] */ << " " ;
    std::cout << std::endl;
  }
  return 0;
}
#endif
