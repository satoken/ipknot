#include "dd_constrained.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <random>

void check(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
bool feasible(const IPModel& model, const std::vector<double>& value) {
  for (const auto& row : model.rows) {
    double sum = 0;
    for (const auto& [col, a] : row.terms) sum += a * value[col];
    if ((row.bound == IP::LO || row.bound == IP::DB || row.bound == IP::FX) && sum < row.lower - 1e-8) return false;
    if ((row.bound == IP::UP || row.bound == IP::DB) && sum > row.upper + 1e-8) return false;
    if (row.bound == IP::FX && sum > row.lower + 1e-8) return false;
  }
  return true;
}
bool planar(const std::vector<DDPair>& pairs, const std::vector<int>& columns,
            const std::vector<double>& value) {
  for(std::size_t i=0;i<pairs.size();++i) if(value[columns[i]]>.5)
    for(std::size_t j=0;j<i;++j) if(value[columns[j]]>.5 && pairs[i].level==pairs[j].level) {
      const auto& a=pairs[i];const auto& b=pairs[j];
      if(a.left==b.left || a.left==b.right || a.right==b.left || a.right==b.right ||
          (a.left<b.left && b.left<a.right && a.right<b.right) ||
          (b.left<a.left && a.left<b.right && b.right<a.right)) return false;
    }
  return true;
}
void audit(IPModel& model, const std::vector<DDPair>& pairs,
           const std::vector<int>& columns, double expected, bool linear = false) {
  DDOptions options; options.linear_constraints = linear;
  options.max_iterations = 40; options.patience = 0; options.constraint_states = linear ? 4096 : 0;
  options.constraint_passes = 1;
  for (int beam : {0, 1, 100}) {
    options.beam = beam;
    try {
      const auto result = solve_constrained_dd(16, pairs, columns, 1, model, options);
      check(std::isfinite(expected), "Returned an infeasible model");
      check(feasible(model, model.solution), "Primal violates a recorded row");
      check(planar(pairs,columns,model.solution), "Primal violates implicit planarity");
      double actual = 0;
      for (std::size_t col = 0; col < model.variables.size(); ++col)
        actual += model.variables[col].coefficient * model.solution[col];
      check(std::abs(actual - result.objective) < 1e-8, "Primal score mismatch");
      check(std::abs(actual - expected) < 1e-8, "Unlimited repair differs from independent enumeration");
      check(std::isfinite(result.upper_bound) && result.upper_bound + 1e-8 >= expected, "Invalid constrained DD certificate");
    } catch (const DDInfeasible&) {
      check(!std::isfinite(expected), "Rejected an independently feasible model");
    }
  }
}
int main() {
  std::mt19937 random(1701);
  for (int trial = 0; trial < 100; ++trial) {
    IPModel model; IP ip(model);
    std::vector<DDPair> pairs; std::vector<int> columns;
    for (int col = 0; col < 6; ++col) {
      const auto weight = (int(random() % 17) - 8) / 7.0;
      const auto id = ip.make_variable(weight);
      if (col < 4) { pairs.push_back({col, 15-col, 0, weight}); columns.push_back(id); }
    }
    // Include bounded integer slacks and signed/repeated row terms.
    ip.make_variable(-.3, 0, 2);
    for (int r = 0; r < 5; ++r) {
      const auto bnd = static_cast<IP::BoundType>(1 + random() % 4);
      const double lo = int(random() % 5) - 2;
      const auto row = ip.make_constraint(bnd, lo, lo + random() % 3);
      for (int col = 0; col < 7; ++col) if (random() % 2)
        ip.add_constraint(row, col, int(random() % 5) - 2);
      if (trial % 3 == 0) { ip.add_constraint(row, 0, 1); ip.add_constraint(row, 0, -1); }
    }
    double expected = -std::numeric_limits<double>::infinity();
    for (int bits = 0; bits < 64; ++bits) for (int slack = 0; slack < 3; ++slack) {
      std::vector<double> value(7); double score = 0;
      for (int col = 0; col < 6; ++col) value[col] = (bits >> col) & 1;
      value[6] = slack;
      if (feasible(model, value)) {
        for (int col = 0; col < 7; ++col) score += model.variables[col].coefficient * value[col];
        expected = std::max(expected, score);
      }
    }
    audit(model, pairs, columns, expected);
  }
  // Arbitrary long/crossing pairs with no pair-pair exclusion rows. Compare
  // the implicit linear sweep and bounded propagation against an independent
  // O(m^2) matching/planarity check over every complete assignment.
  for(int trial=0;trial<120;++trial) {
    IPModel model;IP ip(model);std::vector<DDPair> pairs;std::vector<int> columns;
    for(int id=0;id<7;++id) {
      int left=random()%11,right=left+1+random()%(15-left),level=random()%2;
      bool duplicate=false;for(const auto& p:pairs) duplicate|=p.left==left && p.right==right && p.level==level;
      if(duplicate) {--id;continue;}
      const double weight=(int(random()%13)-6)/5.;
      columns.push_back(ip.make_variable(weight));pairs.push_back({left,right,level,weight});
    }
    for(int i=0;i<16;++i) {
      int row=ip.make_constraint(IP::UP,0,1);
      for(int id=0;id<7;++id) if(pairs[id].left==i || pairs[id].right==i) ip.add_constraint(row,columns[id],1);
    }
    const int row=ip.make_constraint(IP::FX,trial%4,trial%4);
    for(int col:columns) ip.add_constraint(row,col,1);
    double expected=-std::numeric_limits<double>::infinity();
    for(int bits=0;bits<128;++bits) {
      std::vector<double> value(7);double score=0;
      for(int id=0;id<7;++id) {value[id]=(bits>>id)&1;score+=value[id]*model.variables[id].coefficient;}
      if(feasible(model,value) && planar(pairs,columns,value)) expected=std::max(expected,score);
    }
    DDOptions options;options.constraint_states=4096;options.constraint_passes=1;
    options.max_iterations=10;options.patience=0;options.beam=1;
    try {
      const auto result=solve_constrained_dd(16,pairs,columns,2,model,options);
      check(feasible(model,model.solution) && planar(pairs,columns,model.solution),"Invalid linear repair structure");
      check(std::abs(result.objective-expected)<1e-8,"Linear repair disagrees with independent planar enumeration");
      check(result.upper_bound+1e-8>=expected,"Linear DD certificate below exact planar optimum");
      check(result.propagation_work<=2*(result.repair_states+1)*options.constraint_passes*result.nonzeros,"Unbounded linear propagation work");
    } catch(const DDInfeasible&) {check(!std::isfinite(expected),"Linear repair pruned a feasible planar structure");}
  }
  // A wide count row must update all unique terms in a fixed number of scans.
  IPModel wide;IP count(wide);
  int empty=count.make_constraint(IP::FX,0,0);
  for(int i=0;i<1024;++i) count.add_constraint(empty,count.make_variable(.1),1);
  DDOptions linear;linear.constraint_states=8;linear.constraint_passes=2;linear.max_iterations=2;
  const auto wide_result=solve_constrained_dd(0,{}, {},1,wide,linear);
  check(wide_result.objective==0 && wide_result.propagation_work<=8192,"Wide count row did not use bounded sweeps");
  // Negative forced pairs: zero is not a feasible lower bound.
  IPModel forced; IP ip(forced); const int pair = ip.make_variable(-2.5);
  const int row = ip.make_constraint(IP::FX, 1, 1); ip.add_constraint(row, pair, 1);
  audit(forced, {{0,15,0,-2.5}}, {pair}, -2.5);
  // A generic continuous signed product z = -.3*x*y, independent of integer
  // witnesses. Exhaustively compare its exact vertex value.
  IPModel pk; IP product(pk);
  const int x = product.make_variable(.2), y = product.make_variable(.4);
  const int z = product.make_continuous_variable(1, -.3, 0);
  auto add = [&](IP::BoundType b, double lower, double upper,
                 std::initializer_list<std::pair<int,double>> terms) {
    int id = product.make_constraint(b, lower, upper);
    for (const auto& [col,a] : terms) product.add_constraint(id,col,a);
  };
  add(IP::UP,0,0,{{z,1}});
  add(IP::LO,0,0,{{z,1},{x,.3}});
  add(IP::LO,0,0,{{z,1},{y,.3}});
  add(IP::UP,.3,.3,{{z,1},{x,.3},{y,.3}});
  audit(pk,{{0,15,0,.2},{1,14,0,.4}},{x,y},.4);
  // A tiny budget is explicitly reported rather than called infeasibility.
  IPModel limited; IP limit(limited);
  for (int i=0;i<8;++i) limit.make_variable(.1);
  DDOptions options; options.constraint_states=1; options.max_iterations=1;
  int exactly=limit.make_constraint(IP::FX,4,4);
  for(int i=0;i<8;++i) limit.add_constraint(exactly,i,1);
  bool reported=false;
  try { solve_constrained_dd(0,{}, {},1,limited,options); }
  catch(const DDInfeasible&) { throw std::runtime_error("Budget exhaustion claimed infeasibility"); }
  catch(const std::runtime_error& error) { reported=std::string(error.what()).find("repair budget")!=std::string::npos; }
  check(reported,"Missing explicit repair budget diagnostic");
  // A triangle of mutually exclusive nested pairs has integer optimum 1
  // and dual optimum 1.5. Symmetric dual oracles propose all or no pairs,
  // so improvement must come from the interrupted integer search. One-node
  // slices must resume the frontier instead of visiting the root repeatedly.
  const auto triangle = [] {
    IPModel m; IP ip(m);
    for (int i=0;i<3;++i) ip.make_variable(1);
    for (int i=0;i<3;++i) for (int j=0;j<i;++j) {
      int row=ip.make_constraint(IP::UP,0,1);
      ip.add_constraint(row,i,1); ip.add_constraint(row,j,1);
    }
    m.solution.assign(3,0);
    return m;
  };
  std::vector<DDPair> nested{{0,15,0,1},{1,14,0,1},{2,13,0,1}};
  DDOptions resumed; resumed.beam=1; resumed.nussinov_dp=true; resumed.constraint_states=5;
  resumed.constraint_recovery_every=10; resumed.patience=0;
  for (bool target_best : {false,true}) {
  resumed.recovery_target_best=target_best;
  double previous=-std::numeric_limits<double>::infinity();
  for (int iterations : {50,100,500}) {
    auto m=triangle(); resumed.max_iterations=iterations;
    const auto r=solve_constrained_dd(16,nested,{0,1,2},1,m,resumed);
    check(r.objective==1 && feasible(m,m.solution),"Resumed search lost the feasible one-pair completion");
    check(r.repair_states<=5 && r.periodic_repair_calls>=2 && r.repair_calls==r.repair_states,
          "Periodic recovery restarted the root or exceeded the shared state budget");
    check(r.objective>=previous,"More DD iterations discarded a checkpoint incumbent");
    check(r.dp_pruned_states==0 && r.upper_bound>=r.objective,"Exact constrained DP bound is invalid");
    previous=r.objective;
  }
  }
  // Exhausting the search tree is separate from consuming the budget. The
  // same resumed search must still return the independently known optimum.
  for (bool full : {false,true}) {
    auto m=triangle(); resumed.max_iterations=50;
    resumed.linear_constraints=!full; resumed.constraint_states=full ? 0 : 100;
    const auto r=solve_constrained_dd(16,nested,{0,1,2},1,m,resumed);
    check(r.objective==1 && !r.repair_budget_exhausted,"Completed periodic repair was called budget exhaustion");
  }
  // Dense input longer than the default beam. Independent cardinality
  // bound: n/2 pairs is attainable using adjacent pairs. The DP selector
  // must disable pruning even when the configured beam is one.
  const int dense_length=180;
  IPModel dense; IP dense_ip(dense);
  std::vector<DDPair> dense_pairs; std::vector<int> dense_columns;
  for(int i=0;i<dense_length;++i) for(int j=i+1;j<dense_length;++j) {
    dense_columns.push_back(dense_ip.make_variable(1));
    dense_pairs.push_back({i,j,0,1});
  }
  dense.solution.assign(dense.variables.size(),0);
  DDOptions exact_dense; exact_dense.nussinov_dp=true; exact_dense.beam=1;
  const auto exact=solve_constrained_dd(dense_length,dense_pairs,dense_columns,1,dense,exact_dense);
  check(exact.objective==dense_length/2 && exact.dp_pruned_states==0,
        "Explicit Nussinov mode pruned dense interval states or missed the cardinality optimum");
  check(planar(dense_pairs,dense_columns,dense.solution),"Dense exact DP traceback is not a matching");
  exact_dense.nussinov_dp=false;
  const auto beam=solve_constrained_dd(dense_length,dense_pairs,dense_columns,1,dense,exact_dense);
  check(beam.dp_pruned_states>0,"Dense regression input did not distinguish beam from exact DP");
  std::cout << "Constrained DD row, implicit planarity, linear work, primal, auxiliary products, bounds and resumed recovery audits passed\n";
}
