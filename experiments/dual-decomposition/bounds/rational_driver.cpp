#include "dd_bounds.h"
#include "dual_decomposition.h"
#include <iostream>
#include <iomanip>
int main(){int n,m,c,width,offset,lonely; while(std::cin>>n>>m>>c>>width>>offset>>lonely){std::vector<DDPair> p(m);std::vector<double>w(m);std::vector<unsigned char>a(m);for(int k=0;k<m;++k){int mask;std::cin>>p[k].left>>p[k].right>>p[k].level>>p[k].weight>>mask;w[k]=p[k].weight;a[k]=mask;}std::vector<DDScoredContact> cs(c);for(auto&x:cs)std::cin>>x.upper>>x.lower>>x.score;DDBlockBound b(n,p,0,width,offset,lonely);std::cout<<std::setprecision(17)<<b.evaluate(w,a)<<' '<<dd_global_bound(n,p,cs,a)<<'\n';}}
