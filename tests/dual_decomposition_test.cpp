#include "dual_decomposition.h"
#include "ip.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <chrono>
#include <filesystem>
#include <fstream>
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
void check_pk_contact_parity() {
  const auto path = std::filesystem::temp_directory_path() /
      ("ipknot-dd-pk-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  {
    std::ofstream output(path);
    output << "IPKNOT_PK_LINEAR_V1\n";
    int index = 0;
    for (const auto name : pk_learned_feature_names())
      output << name << ' ' << (index++ < 8 ? .1 * index - .45 : 0) << '\n';
  }
  PKLearnedModel learned;
  learned.load(path.string());
  {
    std::ofstream output(path);
    output << "IPKNOT_PK_BOUNDED_V1\nbias .2\nsupport 1\ncompetition .5\nloop_cost .002\ncap .05\n";
  }
  PKLearnedModel bounded;
  bounded.load(path.string());
  {
    std::ofstream output(path);
    output << "IPKNOT_PK_LINEAR_AB_V1\n";
    for (std::size_t i=0;i<learned.weights.size();++i)
      output << pk_learned_feature_names()[i] << ' ' << learned.weights[i] << '\n';
  }
  PKLearnedModel block_linear;
  block_linear.load(path.string());
  {std::ofstream output(path);output<<"IPKNOT_PK_EXCLUSION_V1\nintercept .2\nenergy_scale .03\ntemperature 37\n";}
  PKLearnedModel exclusion;exclusion.load(path.string());
  {std::ofstream output(path);output<<"IPKNOT_PK_DP_LOCAL_V1\nintercept .2\nenergy_scale .03\ntemperature 37\n";}
  PKLearnedModel local;local.load(path.string());
  std::filesystem::remove(path);
  for (int levels : {2, 3}) {
    IPModel recording;
    IP ip(recording);
    PKLevelPairs native(levels, std::vector<std::vector<std::pair<unsigned int, int>>>(15));
    PKPosteriorPairs evidence(16);
    std::vector<DDPair> pairs;
    int index = 0;
    for (const auto coordinate : std::vector<std::pair<int, int>>{
        {0,9},{1,8},{2,7},{4,14},{5,13},{6,12},{3,11},{4,10}}) {
      evidence[coordinate.first+1].push_back({static_cast<unsigned>(coordinate.second+1),
                                             coordinate.first < 3 ? .6f : .2f});
      for (int level = 0; level < levels; ++level) if ((index + level) % 3 != 1) {
        const int variable = ip.make_variable(0);
        native[level][coordinate.first].push_back({static_cast<unsigned>(coordinate.second),variable});
        check(variable == static_cast<int>(pairs.size()), "Native and DD pair ids differ");
        pairs.push_back({coordinate.first,coordinate.second,level,.1});
      }
      ++index;
    }
    PKPosteriorContext posterior(evidence);
    for (int core : {0,2,3}) for (int mode = 0; mode < 12; ++mode) {
      PKScoreOptions pk;
      pk.core_width=core;
      pk.crossing = pk.fixed_blocks = true;
      if (mode < 2) pk.intercept = mode ? -1.2 : 1.2;
      else if (mode < 4) {
        pk.energy.model = mode == 2 ? PKLoopEnergyModel::DP : PKLoopEnergyModel::CC;
        pk.energy_scale = .01;
        pk.energy_intercept = .2;
      } else {
        pk.learned = mode == 8 ? bounded : mode == 9 ? block_linear : mode==10 ? exclusion : mode==11 ? local : learned;
        pk.learned_scale = mode == 7 ? 0 : .05;
        pk.hybrid_shape = mode >= 5 && mode < 8;
        if (pk.hybrid_shape) {
          pk.intercept = -.05;
          pk.stem_reward = .025;
          pk.loop_penalty = .0075;
          if (mode == 6) pk.coax_bonus = .2;
        }
      }
      const auto expected = add_pk_h_score(ip,native,pk,&posterior).crossing_scores;
      for (int width : {0, 2, 100}) {
        DDOptions options;
        options.crossing_beam = width;
        options.witnesses = width ? 2 : 0;
        const auto actual = dd_bounded_graph(15,pairs,levels,options,pk,&posterior);
        for (const auto& row : actual.rows) for (const auto& contact : row.contacts) {
          double value = 0;
          const auto upper = expected.find(row.upper);
          if (upper != expected.end()) {
            const auto lower = upper->second.find(contact.first);
            if (lower != upper->second.end()) value = lower->second;
          }
          check(std::abs(contact.second-value) < 1e-12,
                "DD retained contact differs from native shape, energy or learned/hybrid allocation");
        }
        if (!width) for (const auto& [upper,contacts] : expected)
          for (const auto& [lower,value] : contacts) {
            bool found = false;
            for (const auto& row : actual.rows) if (row.upper == upper)
              for (const auto& contact : row.contacts) found |= contact.first == lower;
            check(found, "Unbounded DD omitted a scored native contact");
          }
      }
      pk.crossing = false;
      pk.projected = true;
      std::vector<double> expected_projection(pairs.size());
      for (const auto& [variable, coefficient] : add_pk_h_score(ip,native,pk,&posterior).terms)
        expected_projection[variable] += coefficient;
      for (int width : {0, 1, 100}) {
        DDOptions options;
        options.crossing_beam = width;
        options.witnesses = 1;
        const auto actual = dd_bounded_graph(15,pairs,levels,options,pk,&posterior);
        const auto base = dd_bounded_graph(15,pairs,levels,options,PKScoreOptions(),nullptr);
        check(actual.rows.size() == base.rows.size(), "Projection changed support rows");
        for (std::size_t r = 0; r < actual.rows.size(); ++r)
          check(actual.rows[r].upper == base.rows[r].upper &&
                actual.rows[r].lower_level == base.rows[r].lower_level &&
                actual.rows[r].contacts == base.rows[r].contacts,
                "Projection changed support witnesses or added product factors");
        if (!width) for (std::size_t id = 0; id < pairs.size(); ++id) {
          const double coefficient = actual.projected_coefficients.empty() ? 0 : actual.projected_coefficients[id];
          check(std::abs(coefficient-expected_projection[id]) < 1e-12,
                "Unbounded DD projection differs from native signed full-partner coefficients");
        }
        options.witnesses = 0;
        const auto all_witnesses = dd_bounded_graph(15,pairs,levels,options,pk,&posterior);
        check(actual.projected_coefficients == all_witnesses.projected_coefficients,
              "Witness trimming reduced full-partner projection");
      }
    }
  }
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

void check_best_partner_model() {
  std::mt19937 rng(57213);std::uniform_real_distribution<double> uniform(-.5,.5);
  for(int trial=0;trial<40;++trial) {
    std::vector<DDPair> pairs;
    for(auto [i,j]:std::vector<std::pair<int,int>>{{0,9},{1,8},{3,12},{4,11},{6,15},{7,14}})
      for(int level=0;level<2;++level)pairs.push_back({i,j,level,uniform(rng)});
    DDOptions options;options.beam=0;options.crossing_beam=options.witnesses=0;options.max_iterations=50;
    PKScoreOptions pk;pk.best_partner=pk.crossing=pk.fixed_blocks=true;pk.intercept=trial%2?-.8:.8;
    const auto graph=dd_bounded_graph(16,pairs,2,options,pk,nullptr);
    auto score=[&](const std::vector<int>& ids) {
      std::vector<bool> selected(pairs.size());double value=0;
      for(int id:ids){selected[id]=true;value+=pairs[id].weight;}
      for(const auto& row:graph.rows)if(selected[row.upper]) {
        std::map<int,double> groups;
        for(auto [lower,w]:row.contacts)if(selected[lower])groups[lower/4]+=w;
        check(!groups.empty(),"Best-partner structure lacks a witness");
        double best=-1e100;for(auto [g,w]:groups)best=std::max(best,w);value+=best;
      }
      return value;
    };
    const double exact=oracle(16,pairs,true,score);
    const auto result=solve_dual_decomposition(16,pairs,2,true,options,pk);
    const auto ids=ids_from_result(pairs,result);
    check(feasible(16,pairs,ids,true),"Best-partner DD recovery infeasible");
    check(std::abs(result.objective-score(ids))<1e-10,"Best-partner primal objective differs");
    check(result.objective<=exact+1e-9 && result.upper_bound>=exact-1e-9,"Best-partner DD certificate invalid");
  }
}

int main() {
  check_best_partner_model();
  check_pk_contact_parity();
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
  PKPosteriorPairs evidence(10);
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
  PKScoreOptions pk; pk.crossing = pk.fixed_blocks = true; pk.intercept = 8;
  options.beam = 0;
  auto rewarded = solve_dual_decomposition(11,h,2,true,options,pk);
  check(std::abs(rewarded.objective - 9.6) < 1e-8 && std::abs(rewarded.pk_score - 8) < 1e-8, "PK block contact allocation differs from independent product score");
  pk.intercept = -8;
  auto penalized = solve_dual_decomposition(11,h,2,true,options,pk);
  check(std::abs(penalized.objective - .8) < 1e-8 && penalized.pk_score == 0, "Negative PK correction was dropped or charged to absent contacts");
  pk.intercept = 0; pk.energy.model = PKLoopEnergyModel::DP; pk.energy_intercept = 10; pk.energy_scale = .01;
  auto energy = solve_dual_decomposition(11,h,2,true,options,pk);
  check(std::abs(energy.pk_score - pk.score({2,2,1,1,1})) < 1e-8, "DD PK energy differs from physical block geometry");
  PKScoreOptions projected;
  projected.projected = projected.fixed_blocks = true;
  projected.intercept = 1;
  const std::vector<DDPair> partial{{0,7,0,.4},{1,6,0,.4},{2,5,0,-.1},
                                  {3,10,1,.4},{4,9,1,.4}};
  options.crossing_beam = 0;
  options.witnesses = 1;
  const auto projection = solve_dual_decomposition(11,partial,2,true,options,projected);
  check(std::abs(projection.objective-2.6) < 1e-8 &&
        std::abs(projection.pk_score-1) < 1e-8 && projection.bpseq[2] < 0,
        "Projection did not charge full potential partner when only two of three lower pairs are selected");
  check(projection.scored_contacts == 0 && projection.projected_pairs == 2,
        "Projection unexpectedly created product factors");
  projected.intercept = -1;
  const auto negative_projection = solve_dual_decomposition(11,partial,2,true,options,projected);
  check(std::abs(negative_projection.objective-.8) < 1e-8 && negative_projection.pk_score == 0,
        "Negative projected score was clipped to zero or charged to absent upper pairs");
  options.witnesses = 0;
  // Signed product factors with alternative layer assignments: compare the
  // returned primal and bound with a separately enumerated contact objective.
  for (int trial=0;trial<30;++trial) {
    std::vector<DDPair> alternatives;
    for (const auto& pair:h) for (int level=0;level<2;++level)
      alternatives.push_back({pair.left,pair.right,level,(static_cast<int>(random()%9)-3)/5.0});
    PKScoreOptions contact_pk;
    contact_pk.crossing=contact_pk.fixed_blocks=true;
    contact_pk.intercept=(static_cast<int>(random()%9)-4)/2.0;
    const auto score=[&](const std::vector<int>& ids) {
      double value=0;
      for(int id:ids) value+=alternatives[id].weight;
      for(int a:ids) for(int b:ids)
        if(alternatives[a].level>alternatives[b].level && crossing(alternatives[a],alternatives[b]))
          value+=contact_pk.intercept/4;
      return value;
    };
    const double optimum=oracle(11,alternatives,true,score);
    for(bool projected:{false,true}) {
      options.projected_norm=projected;
      auto value=solve_dual_decomposition(11,alternatives,2,true,options,contact_pk);
      auto ids=ids_from_result(alternatives,value);
      check(feasible(11,alternatives,ids,true),"Product-factor recovery invalid");
      check(std::abs(value.objective-score(ids))<1e-8,"Product-factor primal mismatch");
      check(value.upper_bound+1e-8>=optimum,"Product-factor bound below signed contact optimum");
    }
    // A whole physical window couples the original signed H contacts and
    // both level structures. Its certificate and feasible exchange meet at
    // the independently enumerated integer optimum, even for a narrow beam.
    for (bool optimized : {false, true}) {
      DDOptions coupled;
      coupled.beam=2; coupled.crossing_beam=coupled.witnesses=0;
      coupled.max_iterations=30; coupled.patience=0;
      coupled.joint_bound_width=12; coupled.joint_bound_states=0;
      coupled.exchange_width=12; coupled.exchange_passes=2; coupled.exchange_states=0;
      coupled.recovery_every=7; coupled.recovery_cache=coupled.recovery_share=optimized;
      const auto result=solve_dual_decomposition(11,alternatives,2,true,coupled,contact_pk);
      const auto selected=ids_from_result(alternatives,result);
      check(feasible(11,alternatives,selected,true),"Joint certificate/exchange returned infeasible H structure");
      check(std::abs(result.objective-score(selected))<1e-8,"Joint exchange used proposal instead of original H score");
      check(std::abs(result.objective-optimum)<1e-8,"Whole-window H exchange failed independent integer optimum");
      check(std::abs(result.upper_bound-optimum)<1e-8,"Whole-window H certificate differs from integer optimum");
      check(result.stop_reason=="bound_gap","Whole-window integer bounds did not certify optimal stopping");
    }
  }
  // Additional primal proposals leave the baseline Polyak trajectory in
  // place. They must dominate the baseline under both a fixed budget and
  // patience, including signed contacts and genuinely pruned oracles.
  for (int trial=0;trial<40;++trial) {
    std::vector<DDPair> alternatives;
    for (const auto& pair:h) for (int level=0;level<2;++level)
      alternatives.push_back({pair.left,pair.right,level,(int(random()%13)-5)/4.0});
    PKScoreOptions contact_pk; contact_pk.crossing=contact_pk.fixed_blocks=true;
    contact_pk.intercept=(int(random()%11)-5)/2.0;
    const auto score=[&](const std::vector<int>& ids) {
      double value=0;for(int id:ids)value+=alternatives[id].weight;
      for(int a:ids)for(int b:ids)if(alternatives[a].level>alternatives[b].level && crossing(alternatives[a],alternatives[b]))
        value+=contact_pk.intercept/4;
      return value;
    };
    const double optimum=oracle(11,alternatives,true,score);
    for(int patience:{0,3}) {
      DDOptions base; base.beam=2;base.crossing_beam=base.witnesses=0;
      base.max_iterations=25;base.patience=patience;base.unpruned_bound=true;
      const auto original=solve_dual_decomposition(11,alternatives,2,true,base,contact_pk);
      for(int mode=0;mode<5;++mode) {
        auto improved=base;
        improved.global_bound=mode==0 || mode>=3;
        improved.bound_block=mode==1 || mode>=3 ? 4 : 0;
        improved.bound_every=7;
        improved.recovery_every=mode>=2 ? 7 : 0;
        improved.recovery_target_best=mode==4;
        const auto result=solve_dual_decomposition(11,alternatives,2,true,improved,contact_pk);
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
          const auto result=solve_dual_decomposition(11,alternatives,2,true,enhanced,contact_pk);
          const auto ids=ids_from_result(alternatives,result);
          check(feasible(11,alternatives,ids,true),"Pruned integer-window integration violated feasibility");
          check(result.objective>=original.objective-1e-8,"Baseline-target joint exchange worsened fixed-budget output");
          check(result.upper_bound>=optimum-1e-8 && result.objective<=optimum+1e-8,"Integer-window integration failed signed H oracle");
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
  std::cout << "Exhaustive Nussinov/DD feasibility, bounds, signed PK scores and long traceback passed\n";
}
