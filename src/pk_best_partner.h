#ifndef IPKNOT_PK_BEST_PARTNER_H
#define IPKNOT_PK_BEST_PARTNER_H
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

struct PKPartnerTerm { int group; double score, multiplier; };
struct PKPartnerState { double value = 0; int upper = 0; std::vector<int> lower; };
struct PKPartnerWorkspace {
  PKPartnerState state;
  std::vector<double> original,modified,force;
  std::vector<int> force_index,positive;
};

// Exact max over a,y with a<=sum(y), and an actually active partner group.
// Group indices are made contiguous once for each bounded support row.
// Scan all members once; optimizing each group's forced-nonempty choice
// takes O(d+g), without enumerating the 2^d subsets.
inline const PKPartnerState& pk_best_partner_factor(double upper_multiplier,
    const std::vector<PKPartnerTerm>& terms, int groups,PKPartnerWorkspace& work) {
  auto& result=work.state;result.value=0;result.upper=0;result.lower.resize(terms.size());
  auto& original=work.original;auto& modified=work.modified;auto& force=work.force;
  auto& force_index=work.force_index;auto& positive=work.positive;
  original.assign(groups,0);modified.assign(groups,0);force.assign(groups,-std::numeric_limits<double>::infinity());
  force_index.assign(groups,-1);positive.assign(groups,0);
  double total = 0;
  for (std::size_t k=0; k<terms.size(); ++k) {
    const auto& t=terms[k];
    if (t.group<0 || t.group>=groups || !std::isfinite(t.score) || !std::isfinite(t.multiplier))
      throw std::invalid_argument("Invalid best-partner factor");
    const double before=std::max(0.,t.multiplier), after=t.multiplier+t.score;
    total+=before; original[t.group]+=before; modified[t.group]+=std::max(0.,after);
    positive[t.group]+=after>0;
    if (after>force[t.group]) { force[t.group]=after; force_index[t.group]=k; }
    result.lower[k]=t.multiplier>0;
  }
  result.value=total;
  int winner=-1;
  for (int g=0;g<groups;++g) if (force_index[g]>=0) {
    const double value=upper_multiplier+total-original[g]+modified[g]+(positive[g]?0:force[g]);
    if (value>result.value) { result.value=value; winner=g; }
  }
  if (winner>=0) {
    result.upper=1;
    for (std::size_t k=0;k<terms.size();++k) if(terms[k].group==winner)
      result.lower[k]=terms[k].multiplier+terms[k].score>0;
    if (!positive[winner]) result.lower[force_index[winner]]=1;
  }
  return result;
}
inline PKPartnerState pk_best_partner_factor(double upper_multiplier,
    const std::vector<PKPartnerTerm>& terms,int groups) {
  PKPartnerWorkspace work;return pk_best_partner_factor(upper_multiplier,terms,groups,work);
}

inline double pk_best_partner_value(const std::vector<PKPartnerTerm>& terms,
    const std::vector<int>& lower, int groups) {
  std::vector<double> sums(groups); std::vector<int> counts(groups);
  for(std::size_t k=0;k<terms.size();++k) if(lower[k]) {
    sums[terms[k].group]+=terms[k].score; ++counts[terms[k].group];
  }
  double result=-std::numeric_limits<double>::infinity();
  for(int g=0;g<groups;++g) if(counts[g]) result=std::max(result,sums[g]);
  return result;
}
#endif
