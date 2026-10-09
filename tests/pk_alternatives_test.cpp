#include "pk_score.h"
#include "pk_best_partner.h"
#include <random>
#include <cstdio>
#include <iostream>

void require(bool condition) {if(!condition)throw std::runtime_error("PK alternative test failed");}
void close(double a,double b) {require(std::abs(a-b)<1e-9);}
int main(int argc,char** argv) {
  require(argc==2); const std::string path=argv[1];
  {std::ofstream out(path);out<<"IPKNOT_PK_EXCLUSION_V1\nintercept 0\nenergy_scale 0\ntemperature 37\n";}
  PKLearnedModel model;model.load(path);require(model.active());
  std::mt19937 rng(8831);std::uniform_real_distribution<double> u(0,1);
  for(int k=0;k<5000;++k) {
    const double a=u(rng)*.9,b=u(rng)*.9,z=1-a*b;
    const double p=a*(1-b)/z,q=b*(1-a)/z,empty=1-p-q;
    PKLearnedFeatureArray f{};f[0]=p;f[1]=q;f[2]=2;f[3]=3;f[4]=20;
    if(empty>=.01) close(model.score(f),2*(a-p)+3*(b-q));
    else close(model.score(f),0);
  }
  PKLearnedFeatureArray f{};f[0]=.7;f[1]=.4;f[2]=2;f[3]=2;close(model.score(f),0);
  f[0]=0;f[1]=.2;close(model.score(f),0);
  {std::ofstream out(path);out<<"IPKNOT_PK_EXCLUSION_V1\nintercept 1000000\nenergy_scale 0\ntemperature 37\n";}
  model.load(path);f[0]=.6;close(model.score(f),2*(1-.6)+2*(1-.2));
  {std::ofstream out(path);out<<"IPKNOT_PK_EXCLUSION_V1\nintercept 0\nenergy_scale -1\ntemperature 37\n";}
  bool invalid=false;try {model.load(path);}catch(const std::invalid_argument&){invalid=true;}require(invalid);
  close(model.score(f),2*(1-.6)+2*(1-.2));

  for(int width:{0,2,3})for(int n=1;n<100;++n) {
    int offset=0;for(auto [begin,size]:pk_core_segments(n,width)) {
      require(begin==offset && size>0);if(width)require(size<=width+1 && (n==1||size>=2));offset+=size;
    }require(offset==n);
  }
  for(int k=0;k<5000;++k) {
    const int d=1+rng()%9,g=1+rng()%d;const double upper=4*u(rng)-2;
    std::vector<PKPartnerTerm> terms;
    for(int j=0;j<d;++j)terms.push_back({int(rng()%g),4*u(rng)-2,4*u(rng)-2});
    const auto actual=pk_best_partner_factor(upper,terms,g);double best=-1e100;
    for(int a=0;a<2;++a)for(int mask=0;mask<(1<<d);++mask) {
      if(a&&!mask)continue;std::vector<int> selected(d);double value=a*upper;
      for(int j=0;j<d;++j) {selected[j]=(mask>>j)&1;value+=selected[j]*terms[j].multiplier;}
      if(a)value+=pk_best_partner_value(terms,selected,g);best=std::max(best,value);
    }
    close(actual.value,best);require(!actual.upper||std::any_of(actual.lower.begin(),actual.lower.end(),[](int x){return x;}));
    double value=actual.upper*upper;for(int j=0;j<d;++j)value+=actual.lower[j]*terms[j].multiplier;
    if(actual.upper)value+=pk_best_partner_value(terms,actual.lower,g);close(value,best);
  }
  {std::ofstream out(path);out<<"IPKNOT_PK_RANK_V1\n";for(int k=0;k<8;++k)out<<'f'<<k<<' '<<k<<'\n';}
  PKRankModel ranker;ranker.load(path);close(ranker.score({1,1,1,1,1,1,1,1}),28);
  {std::ofstream out(path);out<<"IPKNOT_PK_RANK_V1\nf0 1\nf0 2\n";}
  invalid=false;try{ranker.load(path);}catch(const std::invalid_argument&){invalid=true;}require(invalid);
  close(ranker.score({1,1,1,1,1,1,1,1}),28);
  std::remove(path.c_str());
  std::cout<<"5000 latent-prior inversions, 5000 exhaustive signed partner factors, disjoint cores, atomic rank loader passed\n";
}
