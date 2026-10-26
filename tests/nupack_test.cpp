#include "nupack/nupack.h"
#include "nupack/linear_nupack.h"
#include <cassert>
#include <cmath>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <iostream>

static void close(double actual,double expected,double tolerance=3e-6)
{
  assert(std::isfinite(actual));
  if (std::abs(actual-expected)>tolerance*(1+std::abs(expected)))
    std::cerr << "actual=" << actual << " expected=" << expected << '\n';
  assert(std::abs(actual-expected)<=tolerance*(1+std::abs(expected)));
}

static void check_probabilities(const LinearNupack& parser,int n)
{
  std::vector<double> paired(n,0);
  for (auto [i,j,p]:parser.posterior())
  {
    assert(i>=0 && i<j && j<n && j-i>=4);
    assert(std::isfinite(p) && p>=0 && p<=1);
    paired[i]+=p; paired[j]+=p;
  }
  for (double p:paired) assert(p<=1+1e-8);
}

// The old polynomial index is replaced by a lexicographic combinatorial
// rank. Enumerating all tuples catches aliasing independently of the formula.
static void tables()
{
  for (int n=0; n<=32; ++n)
  {
    DPTable4<int> table;
    table.resize(n); table.fill(-1);
    int value=0;
    for (int a=0; a<n; ++a) for (int b=a+1; b<n; ++b)
      for (int c=b+1; c<n; ++c) for (int d=c+1; d<n; ++d)
        table(a,b,c,d)=value++;
    if (n) table(0,0,n-1,n-1)=123456;
    value=0;
    for (int a=0; a<n; ++a) for (int b=a+1; b<n; ++b)
      for (int c=b+1; c<n; ++c) for (int d=c+1; d<n; ++d)
        assert(table(a,b,c,d)==value++);
    if (n) assert(table(n-1,n-1,n-1,n-1)==123456);
  }
}

static void defaults()
{
  Nupack<long double> parser;
  parser.load_default_parameters();
  std::ostringstream dump;
  parser.dump_parameters(dump);
  for (const auto* parameter:{"max_asymmetry=3\n","multiloop_penalty=4.6\n",
      "multiloop_paired_penalty=0.1\n","multiloop_unpaired_penalty=0.4\n",
      "pk_penalty=9.6\n","pk_multiloop_penalty=15\n","at_penalty=0\n",
      "intermolecular_initiation=4.09\n"})
    assert(dump.str().find(parameter)!=std::string::npos);
  parser.load_sequence("");
  close(parser.calculate_partition_function(),1);
  parser.calculate_posterior();
  std::vector<float> bp; std::vector<int> off;
  parser.get_posterior(bp,off);
  assert(bp.size()==1 && bp[0]==0 && off==std::vector<int>{0});
  parser.load_sequence("GAAAC");
  parser.load_constraints({4,-1,-1,-1,0});
  assert(parser.calculate_partition_function()>1);
  // Sequence reuse clears old constraints and scoring caches.
  parser.load_sequence("GAAAAC");
  assert(parser.calculate_partition_function()>1);
  bool failed=false;
  try { parser.load_constraints({1}); } catch (const std::invalid_argument&) { failed=true; }
  assert(failed);
}

// Up to twelve bases all crossing bands fit the linear grammar. Compare both
// inside and outside against the full pseudoknot DP, with pruning disabled.
static void short_oracle()
{
  std::mt19937 rng(1947);
  for (int length=0; length<=12; ++length)
    for (int sample=0; sample<128; ++sample)
    {
      std::string seq(length,'A');
      for (auto& b:seq) b=sample<64 ? "ACGU"[rng()%4] : "ACGUN"[rng()%5];
      Nupack<long double> exact;
      exact.load_default_parameters(); exact.load_sequence(seq);
      auto z=exact.calculate_partition_function();
      exact.calculate_posterior();
      std::vector<float> bp; std::vector<int> off;
      exact.get_posterior(bp,off);
      LinearNupack linear(0);
      const auto logz=linear.calculate(seq);
      if (std::abs(logz-std::log(z))>3e-6*(1+std::abs(std::log(z))))
        std::cerr << "sequence=" << seq << '\n';
      close(logz,std::log(z));
      std::map<std::pair<int,int>,double> probabilities;
      for (auto [i,j,p]:linear.posterior()) probabilities[{i,j}]=p;
      for (int i=0; i<length; ++i) for (int j=i+1; j<length; ++j)
      {
        if (std::abs(probabilities[{i,j}]-bp[off[i+1]+j+1])>3e-6)
          std::cerr << "posterior sequence=" << seq << " pair=" << i << ',' << j << '\n';
        close(probabilities[{i,j}],bp[off[i+1]+j+1]);
      }
      check_probabilities(linear,length);
    }
}

// Möbius inversion extracts the weight of one fully paired motif from the
// exact partition functions of all allowed-pair subsets. This checks the
// multiloop and dangle recurrences without sharing their implementation.
static void motif(const std::string& seq,const std::vector<std::pair<int,int>>& pairs)
{
  long double reference=0;
  for (unsigned mask=0; mask<(1u<<pairs.size()); ++mask)
  {
    std::vector<int> allowed(seq.size(),-1);
    unsigned count=0;
    for (unsigned q=0; q<pairs.size(); ++q) if (mask&(1u<<q))
    {
      auto [i,j]=pairs[q]; allowed[i]=j; allowed[j]=i; ++count;
    }
    Nupack<long double> exact;
    exact.load_default_parameters(); exact.load_sequence(seq); exact.load_constraints(allowed);
    reference+=((pairs.size()-count)%2 ? -1 : 1)*exact.calculate_partition_function();
  }
  assert(reference>0);
  std::string constraints(seq.size(),'.');
  for (auto [i,j]:pairs) { constraints[i]='('; constraints[j]=')'; }
  LinearNupack linear(0);
  close(linear.calculate(seq,constraints),std::log(reference),1e-5);
  check_probabilities(linear,seq.size());
  assert(linear.posterior().size()==pairs.size());
  for (auto [i,j,p]:linear.posterior()) close(p,1,1e-8);
}

static void beams_and_constraints()
{
  std::mt19937 rng(42);
  for (int length:{80,160,320})
  {
    std::string seq(length,'A');
    for (auto& b:seq) b="ACGU"[rng()%4];
    LinearNupack linear(8,6);
    linear.calculate(seq);
    check_probabilities(linear,length);
    assert(linear.statistics().peak_beam<=8);
    // Seven pruned state families (P,G,K,R,PK,M,Z) per input position.
    assert(linear.statistics().retained_states<=7u*8u*length);
    const auto first=linear.posterior(); const auto logz=linear.log_partition_function();
    close(linear.calculate(seq),logz,0);
    assert(linear.posterior()==first);
  }
  LinearNupack linear(0);
  close(linear.calculate("GAAAAC","......"),0,1e-10);
  const auto dna=linear.calculate("aAAAAt","(....)");
  close(dna,linear.calculate("AAAAAU","(....)"));
  for (auto [i,j,p]:linear.posterior()) close(p,1,1e-8);
  for (const auto& constraint:{"??", "(?????", ")?????", "x?????", "(....("})
  {
    bool failed=false;
    try { linear.calculate("GAAAAC",constraint); } catch (const std::invalid_argument&) { failed=true; }
    assert(failed);
  }
  bool failed=false;
  try { linear.calculate("AAAAAA","(....)"); } catch (const std::runtime_error&) { failed=true; }
  assert(failed);
  linear.calculate("GNAAC");
  for (auto [i,j,p]:linear.posterior()) assert(i!=1 && j!=1);
  // Explicit crossing pairs must carry nonzero posterior probability.
  linear.calculate("GGAACC");
  std::map<std::pair<int,int>,double> pp;
  for (auto [i,j,p]:linear.posterior()) pp[{i,j}]=p;
  assert(pp[std::make_pair(0,4)]>0 && pp[std::make_pair(1,5)]>0);
  LinearNupack nested(0,30,0);
  assert(linear.log_partition_function()>nested.calculate("GGAACC"));
  // Long forced hairpins remain possible: no sliding sequence window.
  std::string long_seq="G"+std::string(1000,'A')+"C";
  std::string long_cons="("+std::string(1000,'.')+")";
  linear.calculate(long_seq,long_cons);
  assert(linear.posterior().size()==1);
  close(std::get<2>(linear.posterior()[0]),1,1e-8);
  // Both crossing stems can grow beyond eight pairs, without a span window.
  const auto crossing=std::string(9,'G')+std::string(9,'A')+
    std::string(9,'C')+std::string(9,'U');
  LinearNupack two_bands(0,0);
  two_bands.calculate(crossing,std::string(36,'|'));
  assert(two_bands.posterior().size()==18);
  for (auto [i,j,p]:two_bands.posterior()) close(p,1,1e-8);
}

static void ambiguous_bases()
{
  // Include every ambiguous/nonstandard alphabet symbol in the official set.
  // They all behave as N, preserving coordinates and never forming pairs.
  const std::string symbols="NYMRSWPOKVXIDHBZ";
  for (const auto& prototype:{"G?AAC", "GA?AAC", "G?AAAAAC", "GG?AAACC",
      "G?GAAAC?C", "G??GAAAC?C", "G?GAAAC??C", "G??GAAAC??C"})
  {
    std::vector<float> expected;
    double expected_z=0;
    for (char symbol:symbols)
    {
      std::string seq=prototype;
      for (auto& base:seq) if (base=='?') base=symbol;
      Nupack<long double> exact;
      exact.load_default_parameters(); exact.load_sequence(seq);
      const auto z=exact.calculate_partition_function();
      exact.calculate_posterior();
      std::vector<float> bp; std::vector<int> off;
      exact.get_posterior(bp,off);
      if (symbol=='N') { expected=bp; expected_z=std::log(z); }
      else { assert(bp==expected); close(std::log(z),expected_z,0); }
      LinearNupack linear(0);
      close(linear.calculate(seq),std::log(z));
      std::map<std::pair<int,int>,double> probabilities;
      for (auto [i,j,p]:linear.posterior())
      {
        assert(seq[i]!=symbol && seq[j]!=symbol);
        probabilities[{i,j}]=p;
      }
      for (int i=0; i<int(seq.size()); ++i) for (int j=i+1; j<int(seq.size()); ++j)
      {
        close(probabilities[{i,j}],bp[off[i+1]+j+1]);
        if (seq[i]==symbol || seq[j]==symbol) assert(bp[off[i+1]+j+1]==0);
      }
      check_probabilities(linear,seq.size());
    }
  }
  LinearNupack linear(100);
  close(linear.calculate(symbols),0,0);
  assert(linear.posterior().empty());
}

int main()
{
  tables(); defaults(); short_oracle(); ambiguous_bases();
  motif("GAAAC",{{0,4}});
  motif("GGAAACC",{{0,6},{1,5}});
  motif("GAGAAACC",{{0,7},{2,6}});
  motif("GGAAACGAAACC",{{0,11},{1,5},{6,10}});
  motif("GAGAAACAGAAACAC",{{0,14},{2,6},{8,12}});
  beams_and_constraints();
  return 0;
}
