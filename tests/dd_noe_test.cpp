#include "dd_noe.h"
#include "dd_constrained.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <random>

void check(bool x, const char* message) { if (!x) throw std::runtime_error(message); }
bool feasible(const IPModel& m, const std::vector<double>& x) {
  for (const auto& r : m.rows) {
    double v=0; for (const auto& [c,a] : r.terms) v+=a*x[c];
    if ((r.bound==IP::LO || r.bound==IP::DB || r.bound==IP::FX) && v<r.lower-1e-8) return false;
    if ((r.bound==IP::UP || r.bound==IP::DB) && v>r.upper+1e-8) return false;
    if (r.bound==IP::FX && v>r.lower+1e-8) return false;
  }
  return true;
}
int main() {
  std::mt19937 rng(616);
  for (int trial=0; trial<40; ++trial) {
    const bool proxy = trial%2;
    IPModel model; IP ip(model);
    std::vector<DDPair> pairs; std::vector<int> columns;
    for (int p=0;p<3;++p) for(int lv=0;lv<2;++lv) {
      double w=(int(rng()%21)-5)/13.0;
      columns.push_back(ip.make_variable(w)); pairs.push_back({p,p+5,lv,w});
    }
    for(int obs=0;obs<2;++obs) {
      ip.mark_noe_variable(ip.make_variable(0));
      ip.mark_noe_variable(ip.make_variable(0));
      ip.mark_noe_variable(ip.make_variable(-double(1+rng()%8)/5));
      int r=ip.make_constraint(IP::FX,1,1);
      for(int k=0;k<3;++k)ip.add_constraint(r,6+3*obs+k,1);
    }
    std::vector<int> proxies;
    for(int p=0;p<3;++p) {
      int capacity=ip.make_constraint(IP::UP,0,1);
      for(int lv=0;lv<2;++lv)ip.add_constraint(capacity,2*p+lv,1);
      if(proxy) {
        int col=ip.make_variable(0);proxies.push_back(col);
        int r=ip.make_constraint(IP::FX,0,0);ip.add_constraint(r,col,1);
        for(int lv=0;lv<2;++lv)ip.add_constraint(r,2*p+lv,-1);
      }
    }
    auto add_pair=[&](int r,int p,double a) {
      if(proxy)ip.add_constraint(r,proxies[p],a);
      else for(int lv=0;lv<2;++lv)ip.add_constraint(r,2*p+lv,a);
    };
    for(int p=0;p<3;++p) {
      int r=ip.make_constraint(IP::UP,0,0);ip.add_constraint(r,2*p+1,1);
      for(int q=0;q<3;++q)if(p!=q)ip.add_constraint(r,2*q,-1);
    }
    const std::vector<std::vector<int>> req{{6,9},{7},{10}};
    for(int p=0;p<3;++p) {
      int r=ip.make_constraint(IP::UP,0,0);
      for(int y:req[p])ip.add_constraint(r,y,1);
      add_pair(r,p,-1);
    }
    // Blocker uses real DP columns, even when the requirement uses a proxy.
    int blocked=ip.make_constraint(IP::UP,0,1);ip.add_constraint(blocked,7,1);
    for(int lv=0;lv<2;++lv)ip.add_constraint(blocked,4+lv,1);
    int face=ip.make_constraint(IP::UP,0,1);ip.add_constraint(face,6,1);ip.add_constraint(face,10,1);
    auto factor=dd_noe_factor(model,pairs,columns);
    check(factor.columns.size()==6,"RNA column entered NOE factor");
    check(factor.sharing_rows>=1 && factor.blocker_rows>=1,"Missing valid NOE projections");
    double optimum=-std::numeric_limits<double>::infinity();
    for(int bits=0;bits<(1<<12);++bits) {
      std::vector<double>x(model.variables.size());
      for(int c=0;c<12;++c)x[c]=(bits>>c)&1;
      if(proxy)for(int p=0;p<3;++p)x[proxies[p]]=x[2*p]+x[2*p+1];
      bool planar=true;
      for(int lv=0;lv<2;++lv)if(x[lv]+x[2+lv]+x[4+lv]>1)planar=false;
      if(!planar || !feasible(model,x))continue;
      std::vector<double>y;for(int c:factor.columns)y.push_back(x[c]);
      check(feasible(factor.model,y),"NOE projection excludes a full feasible assignment");
      double score=0;for(std::size_t c=0;c<x.size();++c)score+=model.variables[c].coefficient*x[c];
      optimum=std::max(optimum,score);
    }
    check(std::isfinite(optimum),"Invalid fixture");
    std::vector<double>lo(model.variables.size(),0),hi(model.variables.size(),1);
    DDNoEOracle oracle(factor,lo,hi);
    std::vector<double> empty(model.variables.size(),0);
    check(oracle.recover(model,empty) && feasible(model,empty),"Fixed-RNA NOE repair failed");
    check(empty[8]==1 && empty[11]==1,"Fixed empty RNA did not pay NOE violations");
    for(int c:columns)check(empty[c]==0,"Conditional NOE repair changed a RNA pair");
    std::vector<double>cost(model.variables.size()),selected(model.variables.size(),7.5);
    for(int sample=0;sample<3;++sample) {
      for(int c:factor.columns)cost[c]=(int(rng()%101)-50)/19.0;
      double best=-std::numeric_limits<double>::infinity();
      for(int bits=0;bits<64;++bits) {
        std::vector<double>y(6);double score=0;
        for(int k=0;k<6;++k){y[k]=(bits>>k)&1;score+=cost[factor.columns[k]]*y[k];}
        if(feasible(factor.model,y))best=std::max(best,score);
      }
      auto result=oracle.solve(cost,selected);
      check(std::abs(result.value-best)<1e-7,"NOE oracle differs from exhaustive enumeration");
      check(result.upper_bound+1e-7>=best,"Invalid NOE oracle bound");
      for(int c:columns)check(selected[c]==7.5,"NOE ILP changed a DP pair column");
      const auto calls=oracle.calls();oracle.solve(cost,selected);
      check(oracle.calls()==calls && oracle.cache_hits()>0,"Identical objective was not cached");
    }
    model.solution.assign(model.variables.size(),0);model.solution[8]=model.solution[11]=1;
    DDOptions options;options.noe_ilp=true;options.nussinov_dp=true;
    options.linear_constraints=false;options.constraint_states=0;options.max_iterations=25;
    auto r=solve_constrained_dd(8,pairs,columns,2,model,options);
    check(feasible(model,model.solution),"Hybrid DD returned an infeasible primal");
    check(std::abs(r.objective-optimum)<1e-7,"Hybrid DD repair differs from global enumeration");
    check(r.upper_bound+1e-7>=optimum,"Hybrid DD bound excludes global optimum");
    check(r.noe_ilp_variables==6 && r.noe_ilp_calls>0,"Joint factor was not executed");
    if (trial < 8) {
      double previous=-std::numeric_limits<double>::infinity();
      for (int iterations:{50,100}) {
        model.solution.assign(model.variables.size(),0);model.solution[8]=model.solution[11]=1;
        options.linear_constraints=true;options.constraint_states=20;
        options.constraint_recovery_every=10;options.max_iterations=iterations;
        auto finite=solve_constrained_dd(8,pairs,columns,2,model,options);
        check(finite.repair_states<=20,"NOE hybrid exceeded the shared DFS budget");
        check(finite.objective+1e-8>=previous,"Periodic NOE hybrid lost its feasible incumbent");
        check(finite.upper_bound+1e-8>=optimum,"Invalid finite-budget NOE hybrid bound");
        previous=finite.objective;
      }
    }
  }
  // Integer triangle gap is entirely inside NOE: the joint factor closes it.
  IPModel triangle;IP ip(triangle);
  for(int k=0;k<3;++k)ip.mark_noe_variable(ip.make_variable(1));
  for(int k=0;k<3;++k){int r=ip.make_constraint(IP::UP,0,1);ip.add_constraint(r,k,1);ip.add_constraint(r,(k+1)%3,1);}
  triangle.solution.assign(3,0);
  DDOptions o;o.linear_constraints=false;o.constraint_states=0;o.max_iterations=1;
  auto relaxed=solve_constrained_dd(0,{}, {},1,triangle,o);
  o.noe_ilp=true;auto joint=solve_constrained_dd(0,{}, {},1,triangle,o);
  check(std::abs(joint.upper_bound-1)<1e-7 && relaxed.upper_bound>joint.upper_bound+1,"NOE integer hull did not strengthen the bound");
  // Capacity projections require an actual global capacity certificate.
  IPModel uncapacitated;IP u(uncapacitated);
  u.make_variable(0);u.make_variable(0);
  u.mark_noe_variable(u.make_variable(1));u.mark_noe_variable(u.make_variable(1));
  int link=u.make_constraint(IP::UP,0,0);
  for(int c:{0,1})u.add_constraint(link,c,-1);
  for(int c:{2,3})u.add_constraint(link,c,1);
  std::vector<DDPair> two{{0,5,0,0},{0,5,1,0}};
  auto f=dd_noe_factor(uncapacitated,two,{0,1});
  check(f.sharing_rows==0,"Invented a global physical-pair capacity");
  DDNoEOracle free(f,std::vector<double>(4,0),std::vector<double>(4,1));
  std::vector<double> costs{0,0,1,1},selected(4);
  check(free.solve(costs,selected).value==2,"NOE projection excluded a feasible uncapacitated assignment");
  // Different level subsets are not the same RNA expression.
  IPModel subset;IP s(subset);s.make_variable(0);s.make_variable(0);
  s.mark_noe_variable(s.make_variable(1));s.mark_noe_variable(s.make_variable(1));
  int cap=s.make_constraint(IP::UP,0,1);s.add_constraint(cap,0,1);s.add_constraint(cap,1,1);
  int req=s.make_constraint(IP::UP,0,0);s.add_constraint(req,2,1);s.add_constraint(req,0,-1);
  int block=s.make_constraint(IP::UP,0,1);s.add_constraint(block,3,1);s.add_constraint(block,1,1);
  auto sf=dd_noe_factor(subset,two,{0,1});
  check(sf.blocker_rows==0,"Confused different level expressions in a blocker projection");
  check(feasible(subset,{1,0,1,1}),"Invalid subset fixture");
  // An infeasible NOE assignment for one RNA proposal is not global infeasibility.
  IPModel hard;IP h(hard);h.make_variable(0);h.mark_noe_variable(h.make_variable(0));
  int mandatory=h.make_constraint(IP::FX,1,1);h.add_constraint(mandatory,1,1);
  int support=h.make_constraint(IP::UP,0,0);h.add_constraint(support,1,1);h.add_constraint(support,0,-1);
  DDNoEOracle ho(dd_noe_factor(hard,{{0,5,0,0}},{0}),{0,0},{1,1});
  std::vector<double> missing{0,0};check(!ho.recover(hard,missing),"Accepted unsupported hard NOE");
  missing[0]=1;check(ho.recover(hard,missing) && missing[1]==1,"Rejected supported hard NOE");
  IPModel continuous;IP c(continuous);c.mark_noe_variable(c.make_variable(1));c.make_continuous_variable(1,0,1);
  int cr=c.make_constraint(IP::UP,0,1);c.add_constraint(cr,0,1);c.add_constraint(cr,1,1);
  continuous.solution={0,0};o.max_iterations=5;
  auto mixed=solve_constrained_dd(0,{}, {},1,continuous,o);
  check(std::abs(mixed.objective-1)<1e-7 && mixed.upper_bound>=1-1e-7,"Mixed continuous/NOE repair was not exact");
  std::cout<<"NOE factor: exhaustive oracles, valid projections, pair isolation, repeated objectives, global certificates and integer gap passed\n";
}
