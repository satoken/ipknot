#include "pk_learned.h"
#include "pk_ensemble.h"
#include <cstdio>
#include <random>
#include <iostream>

void require(bool condition) { if (!condition) throw std::runtime_error("DP conversion test failed"); }
void close(double a, double b) { if (std::abs(a-b) >= 1e-10) { std::cerr << std::setprecision(17) << a << " != " << b << '\n'; require(false); } }
int main(int argc, char** argv) {
  require(argc == 2);
  const std::string filename = argv[1];
  auto model = [&](bool local, double eta, double tau) {
    { std::ofstream f(filename); f << (local ? "IPKNOT_PK_DP_LOCAL_V1\n" : "IPKNOT_PK_DP_LINEAR_V1\n")
        << std::setprecision(17) << "intercept " << eta << "\nenergy_scale " << tau << "\ntemperature 37\n"; }
    PKLearnedModel result; result.load(filename); return result;
  };
  std::mt19937 random(48217); std::uniform_real_distribution<double> uniform(0,1);
  for (int k = 0; k < 3000; ++k) {
    const double a=uniform(random), b=uniform(random), q=30*uniform(random);
    const double eta=10*uniform(random)-5, tau=uniform(random);
    const int A=2+random()%12, B=2+random()%12;
    PKLearnedFeatureArray features{}; features[0]=a; features[1]=b; features[2]=A; features[3]=B; features[4]=q;
    const double factor=std::exp(eta-tau*q);
    const double z=1-a*b+a*b*factor;
    const double expected=A*(a*(1-b)/z+a*b*factor/z-a)+B*(b*(1-a)/z+a*b*factor/z-b);
    close(model(true,eta,tau).score(features), expected);
    close(model(false,eta,tau).score(features), (A+B)*(eta-tau*q));
    close(model(true,0,0).score(features), 0);
  }
  auto m=model(true,1000000,0);
  PKLearnedFeatureArray f{}; f[0]=.99; f[1]=.99; f[2]=5; f[3]=5; f[4]=20;
  close(m.score(f),.1);
  f[0]=1; f[1]=1; close(m.score(f),0);
  f[0]=0; close(m.score(f),0);
  { std::ofstream s(filename); s << "IPKNOT_PK_DP_LOCAL_V1\nintercept 0\nenergy_scale -1\ntemperature 37\n"; }
  bool failed=false; try { m.load(filename); } catch (const std::invalid_argument&) { failed=true; }
  require(failed); f[0]=.99; f[1]=.99; close(m.score(f),.1); // load is atomic
  PKPosteriorPairs posterior(12);
  posterior[1]={{8,.8f}}; posterior[2]={{7,.4f}};
  PKPosteriorContext context(posterior);
  close(context.mean_probability({0,7,2}),(double(.8f)+double(.4f))/2.0);
  auto c=m.features(context,{0,7,2},{0,7,2},{2,2,1,1,1});
  close(c.values[4],10.1/(.00198720425864083*310.15));

  auto candidate=[](std::vector<std::pair<int,int>> contacts, double utility, int n=14) {
    PKEnsembleCandidate c; c.pairs.assign(n,-1); c.levels.assign(n,-1); c.utility=utility;
    for (auto [i,j]:contacts) { c.pairs[i]=j; c.pairs[j]=i; c.levels[i]=c.levels[j]=0; }
    return c;
  };
  auto h=candidate({{0,7},{1,6},{3,10},{4,9}},0);
  auto energy=pk_component_energy(h.pairs); close(energy.kcal,10.1); require(energy.components==1);
  auto nested=candidate({{0,13},{1,8},{2,7},{4,11},{5,10}},0);
  close(pk_component_energy(nested.pairs).kcal,15.5);
  auto empty=candidate({},0);
  auto result=pk_finite_ensemble({h,empty,h},1,.1,0,.6163314008174533,.5);
  require(result.visits==3 && result.candidates.size()==2);
  for (int base=0;base<14;++base) { double sum=0; for(auto [ij,p]:result.marginals) if(ij.first==base||ij.second==base) sum+=p; require(sum<=1+1e-12); }
  auto duplicate=pk_finite_ensemble({empty,h,empty,h},1,.1,0,.6163314008174533,.5);
  require(result.probabilities==duplicate.probabilities && result.marginals==duplicate.marginals);
  auto one=pk_finite_ensemble({h,h},1,1,0,.6163314008174533,.5);
  require(one.candidates.size()==1 && one.candidates[one.selected].pairs==h.pairs);
  auto mode=candidate({{0,9}},std::log(.4));
  auto left=candidate({{1,8},{2,7}},std::log(.3));
  auto right=candidate({{1,8},{3,6}},std::log(.3));
  auto mea=pk_finite_ensemble({mode,left,right},1,0,0,.6163314008174533,.25);
  require(mea.candidates[mea.selected].pairs!=mode.pairs);
  for(std::size_t k=0;k<mea.candidates.size();++k) {
    double oracle=0;
    for(std::size_t j=0;j<mea.candidates.size();++j) {
      double utility=0; const auto& chosen=mea.candidates[k].pairs; const auto& truth=mea.candidates[j].pairs;
      for(std::size_t i=0;i<chosen.size();++i) if(chosen[i]>static_cast<int>(i)) utility+=(truth[i]==chosen[i]?1:0)-.25;
      oracle+=mea.probabilities[j]*utility;
    }
    close(mea.expected_gain[k],oracle);
  }
  std::remove(filename.c_str());
  std::cout << "3000 independent four-state oracles, finite ensemble MEA/capacity/dedup, original DP loop energies passed\n";
}
