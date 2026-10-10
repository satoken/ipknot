#include "linear_nupack.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <unordered_map>

namespace {
constexpr double negative_infinity=-std::numeric_limits<double>::infinity();
double log_add(double a, double b)
{
  if (a==negative_infinity) return b;
  if (b==negative_infinity) return a;
  if (a<b) std::swap(a,b);
  return a+std::log1p(std::exp(b-a));
}

using Key=std::array<int,4>;
using Chart=std::map<Key,std::size_t>;
struct Edge
{
  double weight;
  std::array<std::size_t,2> children;
  int count;
};
struct Node
{
  double alpha=negative_infinity, beta=negative_infinity;
  std::vector<Edge> edges;
  int pair_left=-1, pair_right=-1;
};

// Only finalized nodes can be children. Recycle rejected candidates before
// processing the next beam, rather than retaining an unpruned forest.
struct Forest
{
  std::vector<Node> nodes;
  std::vector<std::size_t> free, order;
  std::size_t make(int i=-1, int j=-1)
  {
    std::size_t id;
    if (free.empty()) { id=nodes.size(); nodes.emplace_back(); }
    else { id=free.back(); free.pop_back(); nodes[id]=Node{}; }
    nodes[id].pair_left=i; nodes[id].pair_right=j;
    return id;
  }
  void add(std::size_t id, double weight,
           std::initializer_list<std::size_t> children={})
  {
    double alpha=weight;
    Edge e{weight,{},static_cast<int>(children.size())};
    int k=0;
    for (auto c:children) { e.children[k++]=c; alpha+=nodes[c].alpha; }
    if (alpha==negative_infinity) return;
    if (!std::isfinite(alpha)) throw std::runtime_error("non-finite LinearNUPACK energy");
    nodes[id].alpha=log_add(nodes[id].alpha,alpha);
    nodes[id].edges.push_back(e);
  }
  std::size_t get(Chart& chart, Key key, int i=-1, int j=-1)
  {
    auto p=chart.find(key);
    if (p!=chart.end()) return p->second;
    auto id=make(i,j); chart.emplace(key,id); return id;
  }
  std::size_t leaf(double weight)
  {
    auto id=make(); add(id,weight); order.push_back(id); return id;
  }
};
}

LinearNupack::LinearNupack(unsigned beam, int max_loop, bool pseudoknots)
  : beam_(beam), max_loop_(max_loop), pseudoknots_(pseudoknots)
{
  if (max_loop<0)
    throw std::invalid_argument("LinearNUPACK bounds must be nonnegative");
  energy_.load_default_parameters();
}

double LinearNupack::calculate(const std::string& sequence, const std::string& constraints)
{
  energy_.load_sequence(sequence);
  energy_.prepare_scoring();
  const int n=energy_.N;
  statistics_=Statistics{};
  posterior_.clear();
  if (!constraints.empty() && constraints.size()!=sequence.size())
    throw std::invalid_argument("LinearNUPACK constraint length differs from sequence");

  // -1 is unconstrained, -2 unpaired, -3 left, -4 right, -5 paired.
  std::vector<int> cons(n,-1), stack, blocked(n+1,0);
  for (int i=0; i<static_cast<int>(constraints.size()); ++i)
  {
    switch (constraints[i])
    {
      case '?': break;
      case '.': cons[i]=-2; break;
      case '<': cons[i]=-3; break;
      case '>': cons[i]=-4; break;
      case '|': cons[i]=-5; break;
      case '(': stack.push_back(i); break;
      case ')':
        if (stack.empty()) throw std::invalid_argument("unbalanced LinearNUPACK constraint");
        cons[i]=stack.back(); cons[stack.back()]=i; stack.pop_back(); break;
      default: throw std::invalid_argument("invalid LinearNUPACK constraint character");
    }
  }
  if (!stack.empty()) throw std::invalid_argument("unbalanced LinearNUPACK constraint");
  for (int i=0; i<n; ++i) blocked[i+1]=blocked[i]+(cons[i]!=-1 && cons[i]!=-2);
  auto unpaired=[&](int i,int j) { return j<i || blocked[j+1]==blocked[i]; };
  auto paired=[&](int i,int j) {
    return j-i>=4 && energy_.pair_type(i,j)>=0 &&
      (cons[i]==-1 || cons[i]==-3 || cons[i]==-5 || cons[i]==j) &&
      (cons[j]==-1 || cons[j]==-4 || cons[j]==-5 || cons[j]==i);
  };
  auto wc=[&](int i,int j) { return energy_.wc_pair(i,j); };
  auto weight=[&](double e) { return -e/energy_.RT; };
  auto dangle=[&](int i,int j) { return energy_.score_dangle(i,j); };
  auto tail_delta=[&](int tail,int j) {
    if (tail==0) return dangle(j,j);
    if (tail==1) return dangle(j-1,j)-dangle(j-1,j-1);
    return dangle(0,j)-dangle(0,j-1);
  };

  // Next complementary position strictly after each position; O(n) storage.
  std::vector<std::array<int,5>> next(n+1);
  next[n].fill(-1);
  for (int j=n-1; j>=0; --j)
  {
    next[j]=next[j+1];
    for (int b=1; b<=4; ++b)
      if (energy_.pair_map[b][energy_.seq[j]]>=0) next[j][b]=j;
  }

  Forest f;
  const auto one=f.leaf(0);
  std::vector<Chart> p(n), g(n), k(n+1), r(n), pk(n), m(n), z(n), c(n);
  std::vector<std::unordered_map<int,double>> h(n);
  std::vector<double> prefix(n+1,0);

  auto finish=[&](Chart& chart) {
    std::vector<std::pair<double,Key>> ranked;
    for (const auto& item:chart)
      ranked.emplace_back(prefix[item.first[0]]+f.nodes[item.second].alpha,item.first);
    std::sort(ranked.begin(),ranked.end(),[](const auto& a,const auto& b) {
      return a.first!=b.first ? a.first>b.first : a.second<b.second;
    });
    const auto keep=beam_ ? std::min<std::size_t>(beam_,ranked.size()) : ranked.size();
    for (std::size_t q=keep; q<ranked.size(); ++q)
    {
      auto it=chart.find(ranked[q].second);
      f.nodes[it->second]=Node{}; f.free.push_back(it->second); chart.erase(it);
    }
    statistics_.peak_beam=std::max(statistics_.peak_beam,chart.size());
    for (const auto& item:chart)
    {
      f.order.push_back(item.second);
      ++statistics_.retained_states;
      statistics_.retained_edges+=f.nodes[item.second].edges.size();
    }
  };

  // A gap always has its all-unpaired derivation available, independent of
  // pruning of the nonempty Z chart. No dense interval table is allocated.
  auto gap=[&](int i,int j) {
    std::vector<std::pair<std::size_t,double>> choices;
    if (j<i) { choices.emplace_back(one,0); return choices; }
    if (unpaired(i,j))
      choices.emplace_back(one,weight(dangle(i,j)+energy_.score_pk_unpaired(j-i+1)));
    for (int tail=0; tail<3; ++tail)
    {
      auto it=z[j].find({i,tail,0,0});
      if (it!=z[j].end()) choices.emplace_back(it->second,0);
    }
    return choices;
  };

  auto branch=[&](int i,int j,std::size_t node,bool pseudoknot) {
    const double at=pseudoknot ? 0 : energy_.score_at_penalty(i,j);
    const double external=weight(pseudoknot ? energy_.score_pk() : at);
    auto dest=f.get(c[j],{0,0,0,0});
    if (i==0) f.add(dest,external,{node});
    else for (const auto& left:c[i-1]) f.add(dest,external,{left.second,node});
    if (i==0 || j==n-1) return;
    const double multi=weight(pseudoknot ? energy_.score_pk_multiloop()+
      energy_.score_multiloop_paired(2,false) : energy_.score_multiloop_paired(1,false)+at);
    f.add(f.get(m[j],{i,0,pseudoknot ? 2 : 1,0}),multi,{node});
    if (i>0) for (const auto& left:m[i-1])
      f.add(f.get(m[j],{left.first[0],0,2,0}),multi,{left.second,node});
    if (!pseudoknots_) return;
    const double knot=weight(pseudoknot ? energy_.score_pk_pk()+energy_.score_pk_paired(2) :
      energy_.score_pk_paired(1)+at);
    // A completed component can also fill the third pseudoknot gap. The
    // remembered right endpoint becomes this component's end, so following
    // unpaired extensions apply its dangling-end energy.
    if (j+1<n) for (const auto& left:k[i])
    {
      auto key=left.first; key[1]=j;
      f.add(f.get(k[j+1],key),knot,{left.second,node});
    }
    // Like the left multiloop padding in LinearPartition, leading gap padding
    // is bounded. Subsequent rightward extension has no span limit.
    for (int a=std::max(1,i-max_loop_); a<=i; ++a)
      if (unpaired(a,i-1))
        f.add(f.get(z[j],{a,0,0,0}),
          knot+weight(dangle(a,i-1)+energy_.score_pk_unpaired(i-a)),{node});
    if (i>0) for (const auto& left:z[i-1])
      f.add(f.get(z[j],{left.first[0],0,0,0}),knot,{left.second,node});
  };

  for (int j=0; j<n; ++j)
  {
    // Three tail categories suffice: zero, one, and at least two unpaired
    // bases. Only the one-base case uses min(dangle3,dangle5).
    if (unpaired(j,j))
    {
      if (j==0) f.add(f.get(c[j],{0,1,0,0}),weight(tail_delta(0,j)),{one});
      else
      {
        for (const auto& item:c[j-1])
          f.add(f.get(c[j],{0,std::min(2,item.first[1]+1),0,0}),
            weight(tail_delta(item.first[1],j)),{item.second});
        for (const auto& item:m[j-1])
          f.add(f.get(m[j],{item.first[0],std::min(2,item.first[1]+1),item.first[2],0}),
            weight(tail_delta(item.first[1],j)+energy_.score_multiloop_unpaired(1,false)),{item.second});
        for (const auto& item:z[j-1])
          f.add(f.get(z[j],{item.first[0],std::min(2,item.first[1]+1),0,0}),
            weight(tail_delta(item.first[1],j)+energy_.score_pk_unpaired(1)),{item.second});
      }
    }

    if (j+4<n)
    {
      const int q=next[j+4][energy_.seq[j]];
      if (q>=0) h[q][j]=weight(energy_.score_hairpin(j,q));
    }
    std::vector<std::pair<double,int>> hairpins;
    for (const auto& item:h[j]) hairpins.emplace_back(prefix[item.first]+item.second,item.first);
    const auto better=[](const auto& a,const auto& b) {
      return a.first!=b.first ? a.first>b.first : a.second<b.second;
    };
    // A distant complementary base can collect O(n) hairpin candidates.
    // Select in expected linear time before sorting the bounded survivors.
    if (beam_ && hairpins.size()>beam_)
    {
      std::nth_element(hairpins.begin(),hairpins.begin()+beam_,hairpins.end(),better);
      hairpins.resize(beam_);
    }
    std::sort(hairpins.begin(),hairpins.end(),better);
    for (const auto& item:hairpins)
    {
      const int i=item.second;
      if (paired(i,j))
      {
        if (unpaired(i+1,j-1))
          f.add(f.get(p[j],{i,0,0,0},i,j),weight(energy_.score_hairpin(i,j)));
        if (pseudoknots_) f.add(f.get(g[j],{i,i,j,0},i,j),0);
      }
      const int q=next[j+1][energy_.seq[i]];
      if (q>=0) h[q][i]=weight(energy_.score_hairpin(i,q));
    }
    h[j].clear();

    // K remembers the first band's inner left endpoint and the second band's
    // inner left endpoint. Its third gap can grow arbitrarily far right.
    finish(k[j]);
    for (const auto& item:k[j])
    {
      const int a=item.first[0], end=item.first[1], d=item.first[2], inner=item.first[3];
      if (paired(inner,j) && wc(inner,j))
        f.add(f.get(r[j],{a,d,inner,0},inner,j),
          weight(energy_.score_at_penalty(inner,j)),{item.second});
      if (j+1<n && unpaired(j,j))
        f.add(f.get(k[j+1],item.first),weight(energy_.score_pk_unpaired(1)+
          dangle(end+1,j)-dangle(end+1,j-1)),{item.second});
    }

    // Grow the second band outward exactly as G grows the first one. Both
    // bands allow bounded bulges/interior loops and arbitrary stem length.
    finish(r[j]);
    for (const auto& item:r[j])
    {
      const int a=item.first[0], d=item.first[1], b=item.first[2];
      if (wc(b,j))
        for (const auto& left:gap(d+1,b-1))
          f.add(f.get(pk[j],{a,0,0,0}),
            weight(energy_.score_at_penalty(b,j))+left.second,{item.second,left.first});
      for (int outer=std::max(d+1,b-max_loop_-1); outer<b; ++outer)
      {
        if (!unpaired(outer+1,b-1)) continue;
        for (int q=j+1; q<n && (b-outer-1)+(q-j-1)<=max_loop_; ++q)
          if (paired(outer,q) && unpaired(j+1,q-1))
            f.add(f.get(r[q],{a,d,outer,0},outer,q),
              weight(energy_.score_interior(outer,b,j,q,true)),{item.second});
      }
    }

    finish(pk[j]);
    for (const auto& item:pk[j]) branch(item.first[0],j,item.second,true);
    finish(p[j]);
    for (const auto& item:p[j])
    {
      const int i=item.first[0];
      if (wc(i,j)) branch(i,j,item.second,false);
      for (int a=std::max(0,i-max_loop_-1); a<i; ++a)
      {
        if (!unpaired(a+1,i-1)) continue;
        for (int q=j+1; q<n && (i-a-1)+(q-j-1)<=max_loop_; ++q)
          if (paired(a,q) && unpaired(j+1,q-1))
            f.add(f.get(p[q],{a,0,0,0},a,q),
              weight(energy_.score_interior(a,i,j,q,false)),{item.second});
      }
    }
    finish(m[j]);
    if (j+1<n)
      for (const auto& item:m[j])
        if (item.first[2]==2)
        {
          const int i=item.first[0], q=j+1;
          for (int a=std::max(0,i-max_loop_-1); a<i; ++a)
            if (paired(a,q) && wc(a,q) && unpaired(a+1,i-1))
              f.add(f.get(p[q],{a,0,0,0},a,q),weight(energy_.score_multiloop(false)+
                energy_.score_multiloop_paired(1,false)+energy_.score_multiloop_unpaired(i-a-1,false)+
                energy_.score_at_penalty(a,q)+dangle(a+1,i-1)),{item.second});
        }
    finish(z[j]);

    finish(g[j]);
    for (const auto& item:g[j])
    {
      const int a=item.first[0], d=item.first[1], e=item.first[2];
      // G contains a band with bounded interior loops; multiloops inside
      // bands are excluded. Both terminal pairs must be Watson-Crick.
      if (j+1<n && wc(a,j) && wc(d,e))
      {
        std::vector<int> starts;
        // Bound all-unpaired middle padding, while also considering the
        // starts of retained structured gaps of arbitrary length.
        for (int inner=std::max(d+1,e-1-max_loop_); inner<e; ++inner)
          starts.push_back(inner);
        for (const auto& middle:z[e-1])
          if (middle.first[0]-1>d && middle.first[0]-1<e)
            starts.push_back(middle.first[0]-1);
        std::sort(starts.begin(),starts.end());
        starts.erase(std::unique(starts.begin(),starts.end()),starts.end());
        for (int inner:starts)
          for (const auto& middle:gap(inner+1,e-1))
          {
              f.add(f.get(k[j+1],{a,j,d,inner}),weight(energy_.score_pk_paired(2)+
                energy_.score_pk_band(2)+energy_.score_at_penalty(a,j)+energy_.score_at_penalty(d,e))+
                middle.second,{item.second,middle.first});
          }
      }
      for (int outer=std::max(0,a-max_loop_-1); outer<a; ++outer)
      {
        if (!unpaired(outer+1,a-1)) continue;
        for (int q=j+1; q<n && (a-outer-1)+(q-j-1)<=max_loop_; ++q)
          if (paired(outer,q) && unpaired(j+1,q-1))
            f.add(f.get(g[q],{outer,d,e,0},outer,q),
              weight(energy_.score_interior(outer,a,j,q,true)),{item.second});
      }
    }
    // The external prefix has only three states and must not be pruned.
    double total=negative_infinity;
    for (const auto& item:c[j])
    {
      f.order.push_back(item.second);
      total=log_add(total,f.nodes[item.second].alpha);
    }
    prefix[j+1]=total;
  }

  auto root=f.make();
  if (n==0) f.add(root,0,{one});
  else for (const auto& item:c[n-1]) f.add(root,0,{item.second});
  f.order.push_back(root);
  log_z_=f.nodes[root].alpha;
  if (!std::isfinite(log_z_))
    throw std::runtime_error("LinearNUPACK beam contains no structure satisfying the constraints");
  f.nodes[root].beta=0;
  std::unordered_map<std::uint64_t,double> probabilities;
  for (auto it=f.order.rbegin(); it!=f.order.rend(); ++it)
  {
    const auto id=*it;
    const auto& node=f.nodes[id];
    if (node.beta==negative_infinity) continue;
    if (node.pair_left>=0)
      probabilities[(std::uint64_t(node.pair_left)<<32)|std::uint32_t(node.pair_right)]+=
        std::exp(node.alpha+node.beta-log_z_);
    for (const auto& edge:node.edges)
    {
      for (int q=0; q<edge.count; ++q)
      {
        // Avoid subtraction of -infinity for an unreachable derivation.
        double outside=node.beta+edge.weight;
        for (int s=0; s<edge.count; ++s)
          if (s!=q) outside+=f.nodes[edge.children[s]].alpha;
        auto& child=f.nodes[edge.children[q]];
        child.beta=log_add(child.beta,outside);
      }
    }
  }
  for (const auto& item:probabilities)
    posterior_.emplace_back(static_cast<int>(item.first>>32),
      static_cast<int>(item.first&0xffffffffu),std::min(1.0,item.second));
  // Stable counting sorts preserve coordinate order without an O(n log n)
  // sort of all reported pairs (one nucleotide may have O(n) partners).
  std::vector<std::tuple<int,int,double>> sorted(posterior_.size());
  for (int coordinate:{1,0})
  {
    std::vector<std::size_t> counts(n+1,0);
    auto key=[&](const auto& value) { return coordinate ? std::get<1>(value) : std::get<0>(value); };
    for (const auto& item:posterior_) ++counts[key(item)+1];
    for (int i=1; i<=n; ++i) counts[i]+=counts[i-1];
    for (const auto& item:posterior_) sorted[counts[key(item)]++]=item;
    posterior_.swap(sorted);
  }
  return log_z_;
}
