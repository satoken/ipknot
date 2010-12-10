// $Id:$

#ifndef __INC_NUPACK_H__
#define __INC_NUPACK_H__

#include <vector>
#include <string>

template < class PF_TYPE >
class Nupack
{
public:
  typedef PF_TYPE pf_type;

private:  
  class DPTable2;
  class DPTable4;
  class DPTableX;

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

private:
  std::vector<int> seq;
  int N;
  DPTable2 Q;
  DPTable2 Qb;
  DPTable2 Qm;
  DPTable2 Qp;
  DPTable2 Qz;
  DPTable4 Qg;
  DPTable4 Qgl;
  DPTable4 Qgr;
  DPTable4 Qgls;
  DPTable4 Qgrs;
};

#endif // __INC_NUPACK_H__

// Local Variables:
// mode: C++
// End:
