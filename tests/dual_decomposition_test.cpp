#include "dual_decomposition.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <unordered_map>

void check(bool v, const char* message) { if (!v) throw std::runtime_error(message); }
bool crossing(const DDPair& a, const DDPair& b) {
  return (a.left < b.left && b.left < a.right && a.right < b.right) ||
         (b.left < a.left && a.left < b.right && b.right < a.right);
}
bool feasible(int n, const std::vector<DDPair>& pairs, const std::vector<int>& ids, bool lonely) {
  std::vector<int> mate(n, -1);
  for (int id : ids) {
    const auto& p = pairs[id];
    if (mate[p.left] >= 0 || mate[p.right] >= 0) return false;
    mate[p.left] = p.right; mate[p.right] = p.left;
    for (int other : ids) if (id != other && p.level == pairs[other].level && crossing(p, pairs[other])) return false;
    for (int lower = 0; lower < p.level; ++lower) {
      bool found = false;
      for (int other : ids) if (pairs[other].level == lower && crossing(p, pairs[other])) found = true;
      if (!found) return false;
    }
    if (lonely) {
      bool found = false;
      for (int other : ids) if (pairs[other].level == p.level &&
          ((pairs[other].left == p.left + 1 && pairs[other].right == p.right - 1) ||
           (pairs[other].left == p.left - 1 && pairs[other].right == p.right + 1))) found = true;
      if (!found) return false;
    }
  }
  return true;
}
std::vector<int> ids_from_result(const std::vector<DDPair>& p, const DDResult& result) {
  std::vector<int> ids;
  for (int id = 0; id < static_cast<int>(p.size()); ++id)
    if (result.bpseq[p[id].left] == p[id].right && result.levels[p[id].left] == p[id].level) ids.push_back(id);
  return ids;
}
// Enumerate physical matchings independently of any interval DP recurrence.
double oracle(int n, const std::vector<DDPair>& pairs, bool lonely,
              const std::function<double(const std::vector<int>&)>& score) {
  std::vector<int> ids, used(n);
  std::vector<std::vector<int>> by_left(n);
  for (int id = 0; id < static_cast<int>(pairs.size()); ++id) by_left[pairs[id].left].push_back(id);
  double best = 0;
  std::function<void(int)> visit = [&](int i) {
    if (i == n) { if (feasible(n, pairs, ids, lonely)) best = std::max(best, score(ids)); return; }
    visit(i + 1);
    if (used[i]) return;
    for (int id : by_left[i]) if (!used[pairs[id].right]) {
      used[i] = used[pairs[id].right] = 1; ids.push_back(id);
      visit(i + 1);
      ids.pop_back(); used[i] = used[pairs[id].right] = 0;
    }
  };
  visit(0); return best;
}

void check_large_improved_column() {
  // More than 1024 starts exercises the fixed radix sorting path, rather
  // than only the small-column implementation checked by exhaustive tests.
  constexpr int n = 1301;
  std::vector<DDPair> pairs;
  std::vector<double> weights;
  for (int left = 0; left < n - 1; ++left) {
    const double weight = 2.0 - left / 1000.0;
    pairs.push_back({left, n - 1, 0, weight});
    weights.push_back(weight);
  }
  const std::vector<unsigned char> allowed(pairs.size(), 1);
  for (int beam : {0, 100}) {
    DDNussinov decoder(n, pairs, 0, beam, false, true);
    std::vector<int> selected;
    check(decoder.decode(weights, allowed, selected) == 2.0 &&
          selected == std::vector<int>{0}, "Large improved column lost its best pair");
    check((decoder.pruned_states() > 0) == (beam > 0),
          "Exact dominance was counted as approximate beam pruning");
    for (auto& weight : weights) weight = -weight;
    check(decoder.decode(weights, allowed, selected) == 0 && selected.empty(),
          "Large reused improved DP selected a negative component");
    for (auto& weight : weights) weight = -weight;
  }
}


int main() {
  check_large_improved_column();
  std::mt19937 random(314159);
  for (int trial = 0; trial < 150; ++trial) {
    const int n = 2 + random() % 7;
    std::vector<DDPair> pairs;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j)
      if (random() % 3) pairs.push_back({i, j, 0, (static_cast<int>(random() % 13) - 6) / 3.0});
    std::vector<double> weights;
    std::vector<unsigned char> mask;
    for (const auto& p : pairs) { weights.push_back(p.weight); mask.push_back(random() % 5 != 0); }
    for (bool lonely : {false, true}) {
      const auto objective = [&](const std::vector<int>& ids) {
        double score = 0;
        for (int id : ids) { if (!mask[id]) return -1e9; score += weights[id]; }
        return score;
      };
      const double expected = oracle(n, pairs, lonely, objective);
      for (bool improved : {false, true}) for (int beam : {0, 100, 2}) {
        DDNussinov decoder(n, pairs, 0, beam, lonely, improved);
        std::vector<int> selected;
        const double actual = decoder.decode(weights, mask, selected);
        if (!decoder.pruned_states()) check(std::abs(actual-expected)<1e-9, "Unpruned oracle was not exact");
        check(feasible(n, pairs, selected, lonely), "DP returned an invalid structure");
        check(std::abs(actual - objective(selected)) < 1e-9, "DP traceback score mismatch");
        if (beam != 2) check(std::abs(actual - expected) < 1e-9, "Nussinov disagrees with exhaustive matching oracle");
        else check(actual <= expected + 1e-9, "Beam score exceeds exact optimum");
        const std::vector<unsigned char> empty_mask(mask.size(), 0);
        check(decoder.decode(weights, empty_mask, selected) == 0 && selected.empty(), "Reused DP leaked a masked pair");
        DDNussinov fresh(n,pairs,0,beam,lonely,improved);
        fresh.decode(weights,empty_mask,selected);
        check(decoder.pruned_states()==fresh.pruned_states(), "Pruning count leaked across oracle calls");
      }
    }
  }
  // The unbounded metric sample equals a direct all-crossing oracle.
  std::vector<std::vector<std::pair<unsigned int, float>>> evidence(10);
  for (int i=1;i<10;++i) for (int j=i+1;j<10;++j)
    if (random()%3==0) evidence[i].push_back({static_cast<unsigned>(j),.1f*(random()%9+1)});
  auto contacts=dd_crossing_evidence(evidence,0);
  double metric_sum=0, metric_oracle=0;
  for(const auto& c:contacts) metric_sum+=c.product;
  for(int i=1;i<10;++i) for(const auto& [j,p]:evidence[i])
    for(int k=i+1;k<static_cast<int>(j);++k) for(const auto& [l,q]:evidence[k])
      if(j<l) metric_oracle+=static_cast<double>(p)*q;
  check(std::abs(metric_sum-metric_oracle)<1e-8,"Crossing evidence differs from complete enumeration");
  DDOptions options;
  options.beam = options.crossing_beam = options.witnesses = options.patience = 0;
  options.max_iterations = 150;
  for (int trial = 0; trial < 60; ++trial) {
    const int n = 7, levels = trial % 3 + 1;
    std::vector<DDPair> pairs;
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j)
      if (random() % 6 == 0) for (int level = 0; level < levels; ++level)
        pairs.push_back({i, j, level, (static_cast<int>(random() % 13) - 4) / 5.0});
    for (bool lonely : {false, true}) {
      const auto objective = [&](const std::vector<int>& ids) {
        double score = 0; for (int id : ids) score += pairs[id].weight; return score;
      };
      const double optimum = oracle(n, pairs, lonely, objective);
      for (int beam : {0, 3}) {
        options.beam = beam;
        auto result = solve_dual_decomposition(n, pairs, levels, lonely, options);
        auto selected = ids_from_result(pairs, result);
        check(feasible(n, pairs, selected, lonely), "DD recovery violates structural constraints");
        check(std::abs(result.objective - objective(selected)) < 1e-8, "DD returned objective differs from its structure");
        check(result.objective <= optimum + 1e-8, "DD exceeds exhaustive optimum");
        check(result.upper_bound + 1e-8 >= optimum, "DD bound is below exhaustive optimum");
      }
    }
  }
  // A narrow beam remains a lower-valued oracle even with the new bound.
  options.unpruned_bound = true;
  for (int trial=0;trial<40;++trial) {
    const int n=8; std::vector<DDPair> pairs;
    for(int i=0;i<n;++i) for(int j=i+1;j<n;++j) if(random()%4==0)
      for(int l=0;l<2;++l) pairs.push_back({i,j,l,(int(random()%13)-4)/5.0});
    for (int beam : {1,3,100}) for (bool lonely : {false,true}) {
      options.beam=beam;
      auto result=solve_dual_decomposition(n,pairs,2,lonely,options);
      double opt=oracle(n,pairs,lonely,[&](const std::vector<int>& ids){double v=0;for(int id:ids)v+=pairs[id].weight;return v;});
      check(result.upper_bound>=opt-1e-8,"Unpruned certificate below integer optimum");
      check(result.objective<=opt+1e-8,"Unpruned stopping exceeds optimum");
    }
  }
  options.unpruned_bound = false;
  std::vector<DDPair> h{{0,7,0,.4},{1,6,0,.4},{3,10,1,.4},{4,9,1,.4}};
  options.beam = 0;
  const auto ordinary = solve_dual_decomposition(11,h,2,true,options);
  check(std::abs(ordinary.objective - 1.6) < 1e-8,
        "Ordinary crossing structure changed its pair-only objective");
  // Alternative layer assignments: compare the returned pair-only primal
  // and bound with an independently enumerated structural objective.
  for (int trial=0;trial<30;++trial) {
    std::vector<DDPair> alternatives;
    for (const auto& pair:h) for (int level=0;level<2;++level)
      alternatives.push_back({pair.left,pair.right,level,(static_cast<int>(random()%9)-3)/5.0});
    const auto score=[&](const std::vector<int>& ids) {
      double value=0;
      for(int id:ids) value+=alternatives[id].weight;
      return value;
    };
    const double optimum=oracle(11,alternatives,true,score);
    for(bool projected:{false,true}) {
      options.projected_norm=projected;
      auto value=solve_dual_decomposition(11,alternatives,2,true,options);
      auto ids=ids_from_result(alternatives,value);
      check(feasible(11,alternatives,ids,true),"Product-factor recovery invalid");
      check(std::abs(value.objective-score(ids))<1e-8,"Product-factor primal mismatch");
      check(value.upper_bound+1e-8>=optimum,"Product-factor bound below pair-only optimum");
    }
    // A whole physical window couples both level structures. Its
    // certificate and feasible exchange meet at
    // the independently enumerated integer optimum, even for a narrow beam.
    for (bool optimized : {false, true}) {
      DDOptions coupled;
      coupled.beam=2; coupled.crossing_beam=coupled.witnesses=0;
      coupled.max_iterations=30; coupled.patience=0;
      coupled.joint_bound_width=12; coupled.joint_bound_states=0;
      coupled.exchange_width=12; coupled.exchange_passes=2; coupled.exchange_states=0;
      coupled.recovery_every=7; coupled.recovery_cache=coupled.recovery_share=optimized;
      const auto result=solve_dual_decomposition(11,alternatives,2,true,coupled);
      const auto selected=ids_from_result(alternatives,result);
      check(feasible(11,alternatives,selected,true),"Joint certificate/exchange returned infeasible H structure");
      check(std::abs(result.objective-score(selected))<1e-8,"Joint exchange used proposal instead of original pair weights");
      check(std::abs(result.objective-optimum)<1e-8,"Whole-window H exchange failed independent integer optimum");
      check(std::abs(result.upper_bound-optimum)<1e-8,"Whole-window H certificate differs from integer optimum");
      check(result.stop_reason=="bound_gap","Whole-window integer bounds did not certify optimal stopping");
    }
  }
  // Additional primal proposals leave the baseline Polyak trajectory in
  // place. They must dominate the baseline under both a fixed budget and
  // patience, including signed pair weights and genuinely pruned oracles.
  for (int trial=0;trial<40;++trial) {
    std::vector<DDPair> alternatives;
    for (const auto& pair:h) for (int level=0;level<2;++level)
      alternatives.push_back({pair.left,pair.right,level,(int(random()%13)-5)/4.0});
    const auto score=[&](const std::vector<int>& ids) {
      double value=0;for(int id:ids)value+=alternatives[id].weight;
      return value;
    };
    const double optimum=oracle(11,alternatives,true,score);
    for(int patience:{0,3}) {
      DDOptions base; base.beam=2;base.crossing_beam=base.witnesses=0;
      base.max_iterations=25;base.patience=patience;base.unpruned_bound=true;
      const auto original=solve_dual_decomposition(11,alternatives,2,true,base);
      for(int mode=0;mode<5;++mode) {
        auto improved=base;
        improved.global_bound=mode==0 || mode>=3;
        improved.bound_block=mode==1 || mode>=3 ? 4 : 0;
        improved.bound_every=7;
        improved.recovery_every=mode>=2 ? 7 : 0;
        improved.recovery_target_best=mode==4;
        const auto result=solve_dual_decomposition(11,alternatives,2,true,improved);
        const auto ids=ids_from_result(alternatives,result);
        check(feasible(11,alternatives,ids,true),"Bound/recovery integration is infeasible");
        check(std::abs(result.objective-score(ids))<1e-8,"Bound/recovery integration score mismatch");
        check(result.upper_bound>=optimum-1e-8 && result.objective<=optimum+1e-8,"Integrated certificate failed exact oracle");
        if(mode<4)check(result.objective>=original.objective-1e-8,"Baseline-target recovery worsened output");
      }
      if (!patience) {
        DDResult previous;
        for (bool optimized : {false, true}) {
          auto enhanced=base;
          enhanced.joint_bound_width=4; enhanced.joint_bound_shift=true;
          enhanced.exchange_width=8; enhanced.exchange_passes=2;
          enhanced.recovery_every=7;
          enhanced.recovery_cache=enhanced.recovery_share=optimized;
          const auto result=solve_dual_decomposition(11,alternatives,2,true,enhanced);
          const auto ids=ids_from_result(alternatives,result);
          check(feasible(11,alternatives,ids,true),"Pruned integer-window integration violated feasibility");
          check(result.objective>=original.objective-1e-8,"Baseline-target joint exchange worsened fixed-budget output");
          check(result.upper_bound>=optimum-1e-8 && result.objective<=optimum+1e-8,"Integer-window integration failed pair-only oracle");
          if (optimized)
            check(result.objective==previous.objective && result.bpseq==previous.bpseq && result.levels==previous.levels,
                  "Recovery caching/sharing changed integer-window prediction");
          previous=result;
        }
      }
    }
  }
  options.projected_norm=true;
  // A low-score inner pair must be retained to support a profitable outer pair.
  std::vector<DDPair> negative_stack{{0,5,0,3},{1,4,0,-1}};
  auto stack = solve_dual_decomposition(6,negative_stack,1,true,options);
  check(std::abs(stack.objective - 2) < 1e-8, "Stack DP discarded necessary negative support");
  // Finite coefficients can overflow the baseline accumulation before a
  // later cancelling term. An infinity must never become a proof of optimum.
  std::vector<DDPair> overflow{{0,5,0,1e308},{2,3,0,1e308},{1,4,0,-1e308}};
  bool overflow_rejected=false;
  try { auto unused=solve_dual_decomposition(6,overflow,1,true,options); (void)unused; }
  catch(const std::overflow_error&) { overflow_rejected=true; }
  check(overflow_rejected,"Nonfinite primal was returned as a certified optimum");
  // Fixed crossing budgets may restrict feasibility, but recovery must remain
  // feasible in the original model; their bound is only for their own graph.
  options.beam=4; options.crossing_beam=2; options.witnesses=1;
  std::vector<DDPair> many;
  for(int i=0;i<20;++i) for(int level=0;level<3;++level) many.push_back({i,i+25,level,.3});
  auto limited=solve_dual_decomposition(45,many,3,false,options);
  check(feasible(45,many,ids_from_result(many,limited),false), "Bounded witness recovery invalid");
  check(limited.contacts <= limited.support_rows && limited.crossing_beam_drops > 0 && limited.witness_drops > 0, "Crossing budgets not enforced");
  std::vector<DDPair> long_pairs;
  for(int i=0;i<20000;i+=6) long_pairs.push_back({i,i+4,0,1});
  DDNussinov long_decoder(20005,long_pairs,0,4,false);
  std::vector<double> weights(long_pairs.size(),1);
  std::vector<unsigned char> mask(long_pairs.size(),1);
  std::vector<int> selected;
  check(long_decoder.decode(weights,mask,selected)==long_pairs.size(), "Long iterative traceback failed");
  std::cout << "Exhaustive Nussinov/DD feasibility, bounds, pair-only objectives and long traceback passed\n";
}
