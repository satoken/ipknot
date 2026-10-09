#ifndef IPKNOT_PK_RANKER_H
#define IPKNOT_PK_RANKER_H
#include "pk_learned.h"
#include "pk_ensemble.h"
#include <array>
#include <set>
#include <sstream>

constexpr std::size_t PK_RANK_FEATURES=8;
using PKRankFeatures=std::array<double,PK_RANK_FEATURES>;
inline PKRankFeatures pk_rank_features(const std::vector<int>& pairs,
    const PKPosteriorContext& context,double pf,double pfpk,double pk_score) {
  std::vector<std::pair<int,int>> physical;
  for(int i=0;i<static_cast<int>(pairs.size());++i) if(pairs[i]>i) physical.push_back({i,pairs[i]});
  std::vector<bool> crossed(physical.size());
  for(std::size_t a=0;a<physical.size();++a) for(std::size_t b=a+1;b<physical.size();++b) {
    const auto [i,j]=physical[a]; const auto [k,l]=physical[b];
    if(i<k && k<j && j<l) crossed[a]=crossed[b]=true;
  }
  int count=0; double support=0,margin=0;
  for(std::size_t k=0;k<physical.size();++k) if(crossed[k]) {
    const auto [i,j]=physical[k]; ++count;
    if(context.contains(i,j)) { support+=std::min(1.,context.evidence(i,j)); margin+=context.bounded_margin(i,j); }
  }
  const double n=std::max<std::size_t>(1,pairs.size());
  return {pf,pfpk,physical.size()/n,count/n,count?support/count:0,count?margin/count:0,
      pk_component_energy(pairs).kcal/(10*n),pk_score/n};
}

class PKRankModel {
  bool loaded_=false;
public:
  PKRankFeatures weights{};
  bool loaded() const {return loaded_;}
  void load(const std::string& filename) {
    std::ifstream input(filename); if(!input) throw std::runtime_error("Cannot read PK ranking model");
    std::string header; if(!(input>>header) || header!="IPKNOT_PK_RANK_V1")
      throw std::invalid_argument("Invalid PK ranking header");
    PKRankFeatures candidate{}; std::array<bool,PK_RANK_FEATURES> seen{};
    std::string key; double value;
    while(input>>key) {
      if(!(input>>value) || !std::isfinite(value) || key.size()!=2 || key[0]!='f' || key[1]<'0' || key[1]>'7')
        throw std::invalid_argument("Invalid PK ranking coefficient");
      const int k=key[1]-'0'; if(seen[k]) throw std::invalid_argument("Duplicate PK ranking coefficient");
      candidate[k]=value; seen[k]=true;
    }
    if(!input.eof() || !std::all_of(seen.begin(),seen.end(),[](bool x){return x;}))
      throw std::invalid_argument("Incomplete PK ranking model");
    weights=candidate; loaded_=true;
  }
  double score(const PKRankFeatures& features) const {
    double result=0; for(std::size_t k=0;k<weights.size();++k) result+=weights[k]*features[k];
    if(!std::isfinite(result)) throw std::invalid_argument("PK ranking overflow");
    return result;
  }
};

struct PKRankCandidate { std::vector<int> pairs; double base; PKRankFeatures features; };
inline void write_pk_rank_pool(const std::string& filename,const std::vector<PKRankCandidate>& pool,
    std::size_t selected) {
  std::ofstream out(filename,std::ios::app); if(!out) throw std::runtime_error("Cannot write PK rank pool");
  out<<std::setprecision(17)<<"{\"selected\":"<<selected<<",\"candidates\":[";
  for(std::size_t k=0;k<pool.size();++k) {
    if(k) out<<','; const auto& c=pool[k];
    out<<"{\"base\":"<<c.base<<",\"features\":[";
    for(std::size_t j=0;j<c.features.size();++j) {if(j)out<<',';out<<c.features[j];}
    out<<"],\"pairs\":["; bool first=true;
    for(std::size_t i=0;i<c.pairs.size();++i) if(c.pairs[i]>static_cast<int>(i)) {
      if(!first)out<<',';first=false;out<<'['<<i<<','<<c.pairs[i]<<']';
    }
    out<<"]}";
  }
  out<<"]}\n"; if(!out) throw std::runtime_error("Failed writing PK rank pool");
}
#endif
