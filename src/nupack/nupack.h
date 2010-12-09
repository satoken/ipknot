// $Id:$

#ifndef __INC_NUPACK_H__
#define __INC_NUPACK_H__

#include <vector>

template < class PF_TYPE >
class Nupack
{
public:
  typedef PF_TYPE pf_type;
public:
  Nupack();
  void load_sequence();
  void load_parameters();
  pf_type calculate_partition_function();
  void calcualt_posterior();
  
private:
  void fastiloop();

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
