// $Id:$

#ifndef __INC_NUPACK_H__
#define __INC_NUPACK_H__

#include <vector>
#include <string>
#include <boost/multi_array.hpp>

template <class T> class DPTable2;
template <class T> class DPTable4;

template < class PF_TYPE >
class Nupack
{
public:
  typedef PF_TYPE pf_type;

public:
  Nupack();
  void load_sequence(const std::string& s);
  void load_parameters();
  pf_type calculate_partition_function();
  void calcualt_posterior();
  
private:
  void fastiloop();
  energy_t score_hairpin(int i, int j) const;
  energy_t score_loop(int l) const;
  energy_t score_interior(int i, int d, int e, int j) const;
  energy_t score_interior_mismatch(int i, int j, int k, int l) const;
  energy_t score_interior_asymmetry(int d) const;
  energy_t score_multiloop() const;
  energy_t score_multiloop_paired(int n) const;
  energy_t score_multiloop_unpaired(int n) const;
  energy_t score_at_penalty(int i, int j) const;
  energy_t score_dangle(int i, int j) const;
  energy_t score_pk() const;  
  energy_t score_pk_multiloop() const;
  energy_t score_pk_pk() const;  
  energy_t score_pk_paired(int n) const;
  energy_t score_pk_unpaired(int n) const;

  bool allow_paired(int i, int j) const;
  bool wc_pair(int i, int j) const;
  int pair_type(int i, int j) const;

private:
  std::vector<int> seq;
  int N;
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

  // energy parameters
  std::vector<energy_t> hairpin37;
  std::vector<energy_t> bulge37;
  std::vector<energy_t> interior3;
  boost::multi_array<energy_t,2> stack37;
  boost::multi_array<energy_t,4> int11;
  boost::multi_array<energy_t,5> int12;
  boost::multi_array<energy_t,6> int22;
  boost::multi_array<energy_t,3> dangle3_37;
  boost::multi_array<energy_t,3> dangle5_37;
  boost::multi_array<energy_t,5> triloop37;
  boost::multi_array<energy_t,6> tetraloop37;
  boost::multi_array<energy_t,3> mismatch_hairpin37;
  boost::multi_array<energy_t,3> mismatch_interior37;
  std::vector<energy_t> asymmetry_penalty;
  energy_t polyC_penalty, polyC_slope, polyC_int;
  energy_t multiloop_penalty; // alpha1
  energy_t multiloop_paired_penalty; // alpha2
  energy_t multiloop_unpaired_penalty; // alpha3
  energy_t pk_penalty; // beta1
  energy_t pk_multiloop_penalty; // beta1m
  energy_t pk_pk_penalty; // beta1p
  energy_t pk_paired_penalty; // beta2
  energy_t pk_unpaired_penalty; // beta3
  energy_t at_penalty;

  energy_t max_asymmetry;
  energy_t SALT_CORRECTION;
};

#endif // __INC_NUPACK_H__

// Local Variables:
// mode: C++
// End:
