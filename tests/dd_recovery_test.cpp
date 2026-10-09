#include "dd_recovery.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>

namespace {
void check(bool value, const char* message) {
  if (!value) throw std::runtime_error(message);
}
bool crosses(const DDPair& a, const DDPair& b) {
  return (a.left < b.left && b.left < a.right && a.right < b.right) ||
         (b.left < a.left && a.left < b.right && b.right < a.right);
}
bool adjacent(const DDPair& a, const DDPair& b) {
  return a.level == b.level && ((a.left == b.left + 1 && a.right == b.right - 1) ||
                              (b.left == a.left + 1 && b.right == a.right - 1));
}
bool feasible(int n, const std::vector<DDPair>& pairs,
              const std::vector<DDRecoveryRow>& rows,
              const std::vector<unsigned char>& allowed,
              const std::vector<unsigned char>& selected, bool no_lonely) {
  if (selected.size() != pairs.size()) return false;
  std::vector<unsigned char> used(n);
  for (int id = 0; id < static_cast<int>(pairs.size()); ++id) if (selected[id]) {
    const auto& p = pairs[id];
    if (!allowed[id] || used[p.left] || used[p.right]) return false;
    used[p.left] = used[p.right] = 1;
    bool stacked = false;
    for (int other = 0; other < static_cast<int>(pairs.size()); ++other)
      if (selected[other] && id != other) {
        if (p.level == pairs[other].level && crosses(p, pairs[other])) return false;
        stacked |= adjacent(p, pairs[other]);
      }
    if (no_lonely && !stacked) return false;
    for (int lower = 0; lower < p.level; ++lower) {
      bool supported = false;
      for (const auto& row : rows) if (row.upper == id)
        for (const auto& contact : row.contacts)
          if (pairs[contact.first].level == lower && selected[contact.first]) supported = true;
      if (!supported) return false;
    }
  }
  return true;
}
double score(const std::vector<DDPair>& pairs, const std::vector<DDRecoveryRow>& rows,
             const std::vector<unsigned char>& selected) {
  double value = 0;
  for (std::size_t id = 0; id < pairs.size(); ++id) if (selected[id]) value += pairs[id].weight;
  for (const auto& row : rows) if (selected[row.upper])
    for (const auto& c : row.contacts) if (selected[c.first]) value += c.second;
  return value;
}
double enumerate(int n, const std::vector<DDPair>& pairs,
                 const std::vector<DDRecoveryRow>& rows,
                 const std::vector<unsigned char>& allowed, bool no_lonely) {
  std::vector<std::vector<int>> by_left(n);
  for (int id = 0; id < static_cast<int>(pairs.size()); ++id) if (allowed[id])
    by_left[pairs[id].left].push_back(id);
  std::vector<unsigned char> selected(pairs.size()), used(n);
  double best = 0;
  std::function<void(int)> visit = [&](int base) {
    if (base == n) {
      if (feasible(n, pairs, rows, allowed, selected, no_lonely))
        best = std::max(best, score(pairs, rows, selected));
      return;
    }
    visit(base + 1);
    if (used[base]) return;
    for (int id : by_left[base]) if (!used[pairs[id].right]) {
      selected[id] = used[base] = used[pairs[id].right] = 1;
      visit(base + 1);
      selected[id] = used[base] = used[pairs[id].right] = 0;
    }
  };
  visit(0);
  return best;
}
// Build random signed contacts and independently propagate missing support.
std::vector<DDRecoveryRow> graph(const std::vector<DDPair>& pairs,
                               std::vector<unsigned char>& allowed,
                               bool no_lonely, std::mt19937& random) {
  std::vector<DDRecoveryRow> rows;
  for (int upper = 0; upper < static_cast<int>(pairs.size()); ++upper)
    for (int level = 0; level < pairs[upper].level; ++level) {
      DDRecoveryRow row{upper, {}};
      for (int lower = 0; lower < static_cast<int>(pairs.size()); ++lower)
        if (pairs[lower].level == level && crosses(pairs[upper], pairs[lower]))
          row.contacts.push_back({lower, (static_cast<int>(random() % 11) - 5) / 2.0});
      rows.push_back(std::move(row));
    }
  bool changed = true;
  while (changed) {
    changed = false;
    for (int id = 0; id < static_cast<int>(pairs.size()); ++id) if (allowed[id]) {
      bool possible = true;
      for (const auto& row : rows) if (row.upper == id) {
        bool support = false;
        for (const auto& c : row.contacts) support |= bool(allowed[c.first]);
        possible &= support;
      }
      if (no_lonely) {
        bool support = false;
        for (int other = 0; other < static_cast<int>(pairs.size()); ++other)
          if (allowed[other] && adjacent(pairs[id], pairs[other])) support = true;
        possible &= support;
      }
      if (!possible) { allowed[id] = 0; changed = true; }
    }
  }
  for (auto& row : rows) {
    if (!allowed[row.upper]) row.contacts.clear();
    else row.contacts.erase(std::remove_if(row.contacts.begin(), row.contacts.end(),
         [&](const auto& c) { return !allowed[c.first]; }), row.contacts.end());
  }
  return rows;
}
}

int main() {
  // Level zero prefers an endpoint-conflicting alternative and cannot see a
  // valuable upper helix. A lookahead seed escapes that feasible local choice.
  std::vector<DDPair> pairs{{0,7,0,.4},{1,6,0,.4},{0,4,0,2},{1,3,0,2},
                          {3,10,1,.4},{4,9,1,.4}};
  std::vector<DDRecoveryRow> rows{{4,{{0,2},{1,2}}},{5,{{0,2},{1,2}}}};
  std::vector<unsigned char> allowed(pairs.size(),1), dual(pairs.size());
  dual[2]=dual[3]=dual[4]=dual[5]=1;
  DDPrimalRecovery original(11,pairs,2,0,true,allowed,rows,8,1);
  auto local = original.propose(dual);
  check(std::abs(local.objective-4)<1e-9,"Original seed adversarial baseline mismatch");
  DDPrimalRecovery recovery(11,pairs,2,0,true,allowed,rows);
  std::vector<double> adjusted(pairs.size(),0);
  adjusted[0]=adjusted[1]=3;
  recovery.observe_adjusted(adjusted);
  auto improved = recovery.propose(dual);
  check(std::abs(improved.objective-9.6)<1e-9,"Upper-aware/mean seed missed rewarding motif");
  check(improved.objective > local.objective && feasible(11,pairs,rows,allowed,improved.selected,true),
        "Recovery improvement was infeasible");
  DDPrimalRecovery lookahead(11,pairs,2,0,true,allowed,rows,8,4);
  check(std::abs(lookahead.propose(dual).objective-9.6)<1e-9,"Upper lookahead alone failed");
  // Three alternative witnesses dilute an upper reward. Ordinary support
  // lookahead still prefers a root using the upper pair's left endpoint;
  // the opportunity seed makes room for that upper pair and a true witness.
  std::vector<DDPair> endpoint_pairs{{0,2,0,3},{0,3,0,1},{1,3,0,-.8},{1,4,0,-.8},{2,5,1,5}};
  std::vector<DDRecoveryRow> endpoint_rows{{4,{{1,0},{2,0},{3,0}}}};
  DDPrimalRecovery endpoint(6,endpoint_pairs,2,0,false,{1,1,1,1,1},endpoint_rows,8,4);
  auto freed=endpoint.propose({1,0,0,0,1});
  check(std::abs(freed.objective-6)<1e-9 && freed.method=="upper_opportunity",
        "Endpoint opportunity seed did not free the valuable upper pair");
  check(feasible(6,endpoint_pairs,endpoint_rows,{1,1,1,1,1},freed.selected,false),
        "Endpoint opportunity seed violated graph feasibility");
  // An endpoint penalty remains optimistic about signed interactions; exact
  // final scoring must reject an upper pair with a large negative product.
  for(auto& c:endpoint_rows[0].contacts) c.second=-10;
  DDPrimalRecovery signed_endpoint(6,endpoint_pairs,2,0,false,{1,1,1,1,1},endpoint_rows);
  auto signed_freed=signed_endpoint.propose({1,0,0,0,1});
  check(std::abs(signed_freed.objective-3)<1e-9,
        "Opportunity costs hid a negative signed product during evaluation");
  // Negative corrections must remain in the actual objective. No positive
  // contact filtering is permitted during completion or final evaluation.
  for (auto& row : rows) for (auto& c : row.contacts) c.second=-2;
  DDPrimalRecovery signed_recovery(11,pairs,2,0,true,allowed,rows);
  signed_recovery.observe_adjusted(adjusted);
  auto signed_result=signed_recovery.propose(dual);
  check(std::abs(signed_result.objective-4)<1e-9,"Signed PK contacts changed feasible objective");
  check(std::abs(signed_result.objective-score(pairs,rows,signed_result.selected))<1e-9,
        "Negative contact accounting mismatch");

  // The ring must contain the latest fixed-width window, with no older
  // coefficient retained after wrap-around. Mean-only provides a clean check.
  std::vector<DDPair> alternatives{{0,3,0,1},{0,4,0,2}};
  std::vector<unsigned char> a(2,1), s(2,0);
  DDPrimalRecovery mean(5,alternatives,1,0,false,a,{},2,2);
  mean.observe_adjusted({100,0});
  mean.observe_adjusted({0,2});
  mean.observe_adjusted({0,2});
  check(mean.propose(s).selected == std::vector<unsigned char>({0,1}),
        "Recovery window retained expired history");
  // A negative original seed can be favored by proposal coefficients; the
  // returned lower bound must still retain the feasible empty structure.
  std::vector<DDPair> negative{{0,3,0,-2}};
  DDPrimalRecovery negative_mean(4,negative,1,0,false,{1},{},8,2);
  negative_mean.observe_adjusted({3});
  auto empty=negative_mean.propose({0});
  check(empty.objective==0 && empty.selected==std::vector<unsigned char>({0}),
        "Negative original proposal was used as the lower bound");

  std::mt19937 random(390112);
  for (int trial=0;trial<80;++trial) {
    int n=8, levels=trial%3+1;
    std::vector<DDPair> p;
    for (int i=0;i<n;++i) for (int j=i+1;j<n;++j) if(random()%5==0)
      for (int level=0;level<levels;++level) if(random()%3)
        p.push_back({i,j,level,(static_cast<int>(random()%13)-6)/3.0});
    for(bool no_lonely:{false,true}) {
      std::vector<unsigned char> permitted(p.size(),1);
      auto r=graph(p,permitted,no_lonely,random);
      const double optimum=enumerate(n,p,r,permitted,no_lonely);
      for(int beam:{0,2}) for(unsigned methods:{1u,2u,4u,7u}) {
        DDPrimalRecovery engine(n,p,levels,beam,no_lonely,permitted,r,3,methods);
        DDPrimalRecovery uncached(n,p,levels,beam,no_lonely,permitted,r,3,methods,false);
        std::vector<std::unique_ptr<DDNussinov>> decoders;
        std::vector<DDNussinov*> shared;
        for(int level=0;level<levels;++level) {
          decoders.push_back(std::make_unique<DDNussinov>(n,p,level,beam,no_lonely));
          shared.push_back(decoders.back().get());
        }
        DDPrimalRecovery sharing(n,p,levels,beam,no_lonely,permitted,r,3,methods,true,shared);
        for (int iteration=0;iteration<4;++iteration) {
          std::vector<double> q(p.size());
          std::vector<unsigned char> decoded(p.size());
          for(std::size_t id=0;id<p.size();++id) {
            q[id]=(static_cast<int>(random()%21)-10)/4.0;
            decoded[id]=random()%2;
          }
          engine.observe_adjusted(q);
          uncached.observe_adjusted(q); sharing.observe_adjusted(q);
          auto result=engine.propose(decoded);
          const auto control=uncached.propose(decoded), reused=sharing.propose(decoded);
          check(result.selected==control.selected && result.objective==control.objective,
                "Static proposal caching changed the incumbent");
          check(result.selected==reused.selected && result.objective==reused.objective,
                "Shared decoder buffers changed the incumbent");
          check(result.proposals<=control.proposals,
                "Cached recovery increased proposal completions");
          check(feasible(n,p,r,permitted,result.selected,no_lonely),
                "Recovery violates structural/support/allowed constraints");
          check(std::abs(result.objective-score(p,r,result.selected))<1e-9,
                "Recovery objective differs from independent signed score");
          check(result.objective<=optimum+1e-9,"Recovery exceeds exhaustive optimum");
          check(result.proposals<=5,"Recovery exceeded fixed proposal budget");
        }
      }
    }
  }
  // Cleared impossible rows are valid, but silently missing required rows are
  // an integration error; reject them instead of emitting an infeasible seed.
  bool rejected=false;
  try { DDPrimalRecovery invalid(6,{{0,3,0,1},{2,5,1,1}},2,0,false,{1,1},{}); }
  catch(const std::invalid_argument&) { rejected=true; }
  check(rejected,"Missing witness rows were accepted");
  std::vector<DDPair> impossible{{0,3,0,1},{2,5,1,1}};
  DDPrimalRecovery disabled(6,impossible,2,0,false,{1,0},{{1,{}}});
  check(disabled.propose({0,0}).objective==1,"Cleared impossible row rejected");
  std::cout << "DD recovery original/mean/lookahead seeds, signed objectives and exhaustive feasibility passed\n";
}
