// ScoringModel

#include "ScoringModel.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>


template <class ProbT, class CountT>
//static
const std::string 
ScoringModel<ProbT,CountT>::
alphabet = "ACGU";    // allowed symbols -- all other letters ignored

template <class ProbT, class CountT>
ScoringModel<ProbT,CountT>::
ScoringModel()
{
  MakeMapping();
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
LoadParameters(const std::string& fname, bool as_log_scale /*=false*/)
{
  std::ifstream is(fname.c_str());
  std::string k;
  float v;
  while (is >> k >> v)
  {
    std::vector<std::pair<ProbT,CountT>*>& e=rev_[index_[k]];
    typename std::vector<std::pair<ProbT,CountT>*>::iterator x;
    for (x=e.begin(); x!=e.end(); ++x)
      **x = std::pair<ProbT,CountT>(as_log_scale ? exp(v) : v, 0);
  }

  InitializeCache();
}


template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
SaveParameters(const std::string& fname, bool as_log_scale /*=false*/) const
{
  std::ofstream os(fname.c_str());
  for (unsigned int i=0; i!=name_.size(); ++i)
  {
    if (as_log_scale)
      os << name_[i] << " " << log(rev_[i].front()->first) << std::endl;
    else
      os << name_[i] << " " << log(rev_[i].front()->first) << std::endl;
  }
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
SetParameters(const std::vector<ProbT>& param)
{
  assert(param.size()==rev_.size());
  for (unsigned int i=0; i!=rev_.size(); ++i)
  {
    typename std::vector<std::pair<ProbT,CountT>*>::iterator x;
    for (x=rev_[i].begin(); x!=rev_[i].end(); ++x)
      **x = std::pair<ProbT,CountT>(param[i], 0);
  }

  InitializeCache();
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
GetParameters(std::vector<ProbT>& param) const
{
  param.resize(rev_.size());
  for (unsigned int i=0; i!=rev_.size(); ++i)
    param[i] = rev_[i].front()->first;
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
GetCounts(std::vector<CountT>& cnt)
{
  FinalizeCounts();

  cnt.resize(rev_.size());
  for (unsigned int i=0; i!=rev_.size(); ++i)
  {
    cnt[i] = 0;
    typename std::vector<std::pair<ProbT,CountT>*>::const_iterator x;
    for (x=rev_[i].begin(); x!=rev_[i].end(); ++x)
      cnt[i] += (*x)->second;
  }
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
ClearCounts()
{
  for (unsigned int i=0; i!=rev_.size(); ++i)
  {
    typename std::vector<std::pair<ProbT,CountT>*>::const_iterator x;
    for (x=rev_[i].begin(); x!=rev_[i].end(); ++x)
      (*x)->second = 0;
  }
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreFCHairpin(const std::vector<int>& seq, int i, int j) const
{
  return ScoreHairpin(seq, i, j);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountFCHairpin(const std::vector<int>& seq, int i, int j, CountT c)
{
  CountHairpin(seq, i, j, c);
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreFCSingle(const std::vector<int>& seq, int i, int j, int p, int q) const
{
  if (p == i && q == j) // helix stacking
    return ScoreBasePair(seq, i+1, j) * ScoreHelixStacking(seq, i, j+1);
  else                // interior loops
    return ScoreSingle(seq, i, j, p, q);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountFCSingle(const std::vector<int>& seq, int i, int j, int p, int q, CountT c)
{
  if (p == i && q == j) // helix stacking
  {
    CountBasePair(seq, i+1, j, c);
    CountHelixStacking(seq, i, j+1, c);
  }
  else                // interior loops
    CountSingle(seq, i, j, p, q, c);
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreFCBifurcation(const std::vector<int>& seq, int i, int j) const
{
  return ScoreJunctionA(seq, i, j) * ScoreMultiPaired() * ScoreMultiBase();
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountCBifurcation(const std::vector<int>& seq, int i, int j, CountT c)
{
  CountJunctionA(seq, i, j, c);
  CountMultiPaired(c);
  CountMultiBase(c);
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreFM1Paired(const std::vector<int>& seq, int i, int j) const
{
  return ScoreJunctionA(seq, j, i) * ScoreMultiPaired() * ScoreBasePair(seq, i+1, j);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountFM1Paired(const std::vector<int>& seq, int i, int j, CountT c)
{
  CountJunctionA(seq, j, i, c);
  CountMultiPaired(c);
  CountBasePair(seq, i+1, j, c);
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreFM1Unpaired(const std::vector<int>& seq, int i, int j) const
{
  return ScoreMultiUnpaired(seq, i+1);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountFM1Unpaired(const std::vector<int>& seq, int i, int j, CountT c)
{
  CountMultiUnpaired(seq, i+1, c);
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreFMBifurcation(const std::vector<int>& seq, int i, int j) const
{
  return ProbT(1);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountFMBifurcation(const std::vector<int>& seq, int i, int j, CountT c)
{
}
  
template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreFMUnpaired(const std::vector<int>& seq, int i, int j) const
{
  return ScoreMultiUnpaired(seq, j);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountFMUnpaired(const std::vector<int>& seq, int i, int j, CountT c)
{
  CountMultiUnpaired(seq, j, c);
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreFMtoFM1(const std::vector<int>& seq, int i, int j) const
{
  return ProbT(1);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountFMtoFM1(const std::vector<int>& seq, int i, int j, CountT c)
{
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreF5Unpaired(const std::vector<int>& seq, int j) const
{
  return ScoreExternalUnpaired(seq, j);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountF5Unpaired(const std::vector<int>& seq, int j, CountT c)
{
  CountExternalUnpaired(seq, j, c);
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreF5Bifurcation(const std::vector<int>& seq, int k, int j) const
{
  return ScoreExternalPaired() * ScoreBasePair(seq, k+1, j) * ScoreJunctionA(seq, j, k);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountF5Bifurcation(const std::vector<int>& seq, int k, int j, CountT c)
{
  CountExternalPaired(c);
  CountBasePair(seq, k+1, j, c);
  CountJunctionA(seq, j, k, c);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
MakeMapping()
{
  int n=0;

  char buf1[1000], buf2[1000];
  // base_pair
  for (int i=0; i<=M; ++i)
    for (int j=0; j<=M; ++j)
      if (i==M || j==M)
        score_base_pair[i][j] = std::pair<ProbT,CountT>(1, 0);
      else
      {
        sprintf(buf1, "base_pair_%c%c", alphabet[i], alphabet[j]);
        sprintf(buf2, "base_pair_%c%c", alphabet[j], alphabet[i]);
        if (index_.count(buf1)==0)
        {
          index_[buf1] = index_[buf2] = n;
          rev_.resize(n+1);
          rev_[n].push_back(&score_base_pair[i][j]);
          rev_[n].push_back(&score_base_pair[j][i]);
          name_.push_back(strcmp(buf1, buf2)<0 ? buf1 : buf2);
          n++;
        } 
      }

  // terminal_mismatch
  for (int i1=0; i1<=M; i1++)
    for (int j1=0; j1<=M; j1++)
      for (int i2=0; i2<=M; i2++)
        for (int j2=0; j2<=M; j2++)
          if (i1==M || j1==M || i2==M || j2==M)
            score_terminal_mismatch[i1][j1][i2][j2] = std::pair<ProbT,CountT>(1, 0);
          else
          {
            sprintf(buf1, "terminal_mismatch_%c%c%c%c", alphabet[i1], alphabet[j1], alphabet[i2], alphabet[j2]);
            assert(index_.count(buf1)==0);
            index_[buf1] = n;
            rev_.resize(n+1);
            rev_[n].push_back(&score_terminal_mismatch[i1][j1][i2][j2]);
            name_.push_back(buf1);
            n++;
          }

  // hairpin_length_at_least
  for (int i=0; i<=D_MAX_HAIRPIN_LENGTH; i++)
  {
    sprintf(buf1, "hairpin_length_at_least_%d", i);
    assert(index_.count(buf1)==0);
    index_[buf1] = n;
    rev_.resize(n+1);
    rev_[n].push_back(&score_hairpin_length_at_least[i]);
    name_.push_back(buf1);
    n++;
  }

  // internal_explicit
  for (int i=0; i<=D_MAX_INTERNAL_EXPLICIT_LENGTH; i++)
    for (int j=0; j<=D_MAX_INTERNAL_EXPLICIT_LENGTH; j++)
      if (i==0 || j==0)
        score_internal_explicit[i][j] = std::pair<ProbT,CountT>(1, 0);
      else
      {
        sprintf(buf1, "internal_explicit_%d_%d", std::min(i, j), std::max(i, j));
        if (index_.count(buf1)==0)
        {
          index_[buf1] = n;
          rev_.resize(n+1);
          rev_[n].push_back(&score_internal_explicit[i][j]);
          name_.push_back(buf1);
          n++;
        }
      }

  // bulge_length_at_least
  for (int i=0; i<=D_MAX_BULGE_LENGTH; i++)
    if (i==0)
      score_bulge_length_at_least[i] = std::pair<ProbT,CountT>(1, 0);
    else
    {
      sprintf(buf1, "bulge_length_at_least_%d", i);
      assert(index_.count(buf1)==0);
      index_[buf1] = n;
      rev_.resize(n+1);
      rev_[n].push_back(&score_bulge_length_at_least[i]);
      name_.push_back(buf1);
      n++;
    }

  // internal_length_at_least
  for (int i=0; i<=D_MAX_INTERNAL_LENGTH; i++)
    if (i<2)
      score_internal_length_at_least[i] = std::pair<ProbT,CountT>(1, 0);
    else
    {
      sprintf(buf1, "internal_length_at_least_%d", i);
      assert(index_.count(buf1)==0);
      index_[buf1] = n;
      rev_.resize(n+1);
      rev_[n].push_back(&score_internal_length_at_least[i]);
      name_.push_back(buf1);
      n++;
    }

  // internal_symmetric_length_at_least
  for (int i=0; i<=D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
    if (i==0)
      score_internal_symmetric_length_at_least[i] = std::pair<ProbT,CountT>(1, 0);
    else
    {
      sprintf(buf1, "internal_symmetric_length_at_least_%d", i);
      assert(index_.count(buf1)==0);
      index_[buf1] = n;
      rev_.resize(n+1);
      rev_[n].push_back(&score_internal_symmetric_length_at_least[i]);
      name_.push_back(buf1);
      n++;
    }

  // internal_asymmetry_at_least
  for (int i=0; i<=D_MAX_INTERNAL_ASYMMETRY; i++)
    if (i == 0)
      score_internal_asymmetry_at_least[i] = std::pair<ProbT,CountT>(1, 0);
    else
    {
      sprintf(buf1, "internal_asymmetry_at_least_%d", i);
      assert(index_.count(buf1)==0);
      index_[buf1] = n;
      rev_.resize(n+1);
      rev_[n].push_back(&score_internal_asymmetry_at_least[i]);
      name_.push_back(buf1);
      n++;
    }

  // bulge_0x1_nucleotides
  for (int i1=0; i1<=M; i1++)
    if (i1==M)
    {
      score_bulge_0x1_nucleotides[i1] = std::pair<ProbT,CountT>(1, 0);
      score_bulge_1x0_nucleotides[i1] = std::pair<ProbT,CountT>(1, 0);
    }
    else
    {
      sprintf(buf1, "bulge_0x1_nucleotides_%c", alphabet[i1]);
      assert(index_.count(buf1)==0);
      index_[buf1] = n;
      rev_.resize(n+1);
      rev_[n].push_back(&score_bulge_0x1_nucleotides[i1]);
      rev_[n].push_back(&score_bulge_1x0_nucleotides[i1]);
      name_.push_back(buf1);
      n++;
    }

  // internal_1x1_nucleotides
  for (int i1=0; i1<=M; i1++)
    for (int i2=0; i2<=M; i2++)
      if (i1==M || i2==M)
        score_internal_1x1_nucleotides[i1][i2] = std::pair<ProbT,CountT>(1, 0);
      else 
      {          
        sprintf(buf1, "internal_1x1_nucleotides_%c%c", alphabet[i1], alphabet[i2]);
        sprintf(buf2, "internal_1x1_nucleotides_%c%c", alphabet[i2], alphabet[i1]);
        if (index_.count(buf1)==0)
        {
          index_[buf1] = index_[buf2] = n;
          rev_.resize(n+1);
          rev_[n].push_back(&score_internal_1x1_nucleotides[i1][i2]);
          rev_[n].push_back(&score_internal_1x1_nucleotides[i2][i1]);
          name_.push_back(strcmp(buf1, buf2)<0 ? buf1 : buf2);
          n++;
        }
      }

  // helix_stacking
  for (int i1=0; i1<=M; i1++)
    for (int j1=0; j1<=M; j1++)
      for (int i2=0; i2<=M; i2++)
        for (int j2=0; j2<=M; j2++)
          if (i1==M || j1==M || i2==M || j2==M)
            score_helix_stacking[i1][j1][i2][j2] = std::pair<ProbT,CountT>(1, 0);
          else
          {
            sprintf(buf1, "helix_stacking_%c%c%c%c", alphabet[i1], alphabet[j1], alphabet[i2], alphabet[j2]);
            sprintf(buf2, "helix_stacking_%c%c%c%c", alphabet[j2], alphabet[i2], alphabet[j1], alphabet[i1]);
            if (index_.count(buf1)==0)
            {
              index_[buf1] = index_[buf2] = n;
              rev_.resize(n+1);
              rev_[n].push_back(&score_helix_stacking[i1][j1][i2][j2]);
              rev_[n].push_back(&score_helix_stacking[j2][i2][j1][i1]);
              name_.push_back(strcmp(buf1, buf2)<0 ? buf1 : buf2);
              n++;
            }
          }

  // helix_closing
  for (int i=0; i<=M; i++)
    for (int j=0; j<=M; j++)
      if (i==M || j==M)
        score_helix_closing[i][j] = std::pair<ProbT,CountT>(1, 0);
      else
      {
        sprintf(buf1, "helix_closing_%c%c", alphabet[i], alphabet[j]);
        assert(index_.count(buf1)==0);
        index_[buf1] = n;
        rev_.resize(n+1);
        rev_[n].push_back(&score_helix_closing[i][j]);
        name_.push_back(buf1);
        n++;
      }

  // multi_base
  assert(index_.count("multi_base")==0);
  index_["multi_base"] = n;
  rev_.resize(n+1);
  rev_[n].push_back(&score_multi_base);
  name_.push_back("multi_base");
  n++;

  // multi_unpaired
  assert(index_.count("multi_unpaired")==0);
  index_["multi_unpaired"] = n;
  rev_.resize(n+1);
  rev_[n].push_back(&score_multi_unpaired);
  name_.push_back("multi_unpaired");
  n++;

  // multi_paired
  assert(index_.count("multi_paired")==0);
  index_["multi_paired"] = n;
  rev_.resize(n+1);
  rev_[n].push_back(&score_multi_paired);
  name_.push_back("multi_paired");
  n++;

  // dangle_left
  for (int i1=0; i1<=M; i1++)
    for (int j1=0; j1<=M; j1++)
      for (int i2=0; i2<=M; i2++)
        if (i1==M || j1==M || i2==M)
          score_dangle_left[i1][j1][i2] = std::pair<ProbT,CountT>(1, 0);
        else
        {
          sprintf(buf1, "dangle_left_%c%c%c", alphabet[i1], alphabet[j1], alphabet[i2]);
          assert(index_.count(buf1)==0);
          index_[buf1] = n;
          rev_.resize(n+1);
          rev_[n].push_back(&score_dangle_left[i1][j1][i2]);
          name_.push_back(buf1);
          n++;
        }

  // dangle_right
  for (int i1=0; i1<=M; i1++)
    for (int j1=0; j1<=M; j1++)
      for (int j2=0; j2<=M; j2++)
        if (i1==M || j1==M || j2==M)
          score_dangle_right[i1][j1][j2] = std::pair<ProbT,CountT>(1, 0);
        else
        {
          sprintf(buf1, "dangle_right_%c%c%c", alphabet[i1], alphabet[j1], alphabet[j2]);
          assert(index_.count(buf1)==0);
          index_[buf1] = n;
          rev_.resize(n+1);
          rev_[n].push_back(&score_dangle_right[i1][j1][j2]);
          name_.push_back(buf1);
          n++;
        }

  // external_unpaired
  assert(index_.count("external_unpaired")==0);
  index_["external_unpaired"] = n;
  rev_.resize(n+1);
  rev_[n].push_back(&score_external_unpaired);
  name_.push_back("external_unpaired");
  n++;
    
  // external_paired
  assert(index_.count("external_paired")==0);
  index_["external_paired"] = n;
  rev_.resize(n+1);
  rev_[n].push_back(&score_external_paired);
  name_.push_back("external_paired");
  n++;

  assert(rev_.size()==name_.size());
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreHairpin(const std::vector<int>& seq, int i, int j) const
{
  return ScoreUnpaired(seq, i, j) * ScoreJunctionB(seq, i, j)
    * cache_score_hairpin_length[std::min(j-i, int(D_MAX_HAIRPIN_LENGTH))].first;
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountHairpin(const std::vector<int>& seq, int i, int j, CountT c)
{
  CountUnpaired(seq, i, j, c);
  CountJunctionB(seq, i, j, c);
  cache_score_hairpin_length[std::min(j-i, int(D_MAX_HAIRPIN_LENGTH))].second += c;
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreBasePair(const std::vector<int>& seq, int i, int j) const
{
  return score_base_pair[seq[i]][seq[j]].first;
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountBasePair(const std::vector<int>& seq, int i, int j, CountT c)
{
  score_base_pair[seq[i]][seq[j]].second += c;
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreHelixStacking(const std::vector<int>& seq, int i, int j) const
{
  return score_helix_stacking[seq[i]][seq[j]][seq[i+1]][seq[j-1]].first;
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountHelixStacking(const std::vector<int>& seq, int i, int j, CountT c)
{
  score_helix_stacking[seq[i]][seq[j]][seq[i+1]][seq[j-1]].second += c;
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreSingle(const std::vector<int>& seq, int i, int j, int p, int q) const
{
  return 
    cache_score_single[p-i][j-q].first
    * ScoreBasePair(seq, p+1, q)
    * ScoreJunctionB(seq, i, j) 
    * ScoreJunctionB(seq, q, p)
    * ScoreSingleNucleotides(seq, i, j, p, q);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountSingle(const std::vector<int>& seq, int i, int j, int p, int q, CountT c)
{
  cache_score_single[p-i][j-q].second += c;
  CountBasePair(seq, p+1, q, c);
  CountJunctionB(seq, i, j, c); 
  CountJunctionB(seq, q, p, c);
  CountSingleNucleotides(seq, i, j, p, q, c);
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreSingleNucleotides(const std::vector<int>& seq, int i, int j, int p, int q) const
{
  const int l1 = p-i;
  const int l2 = j-q;
  return 
    ScoreUnpaired(seq, i,p)
    * ScoreUnpaired(seq, q,j)
    * (l1==0 && l2==1 ? score_bulge_0x1_nucleotides[seq[j]].first : ProbT(1))
    * (l1==1 && l2==0 ? score_bulge_1x0_nucleotides[seq[i+1]].first : ProbT(1))
    * (l1==1 && l2==1 ? score_internal_1x1_nucleotides[seq[i+1]][seq[j]].first : ProbT(1));
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountSingleNucleotides(const std::vector<int>& seq, int i, int j, int p, int q, CountT c)
{
  const int l1 = p-i;
  const int l2 = j-q;
  CountUnpaired(seq, i, p, c);
  CountUnpaired(seq, q, j, c);
  if (l1==0 && l2==1) score_bulge_0x1_nucleotides[seq[j]].second +=c;
  if (l1==1 && l2==0) score_bulge_1x0_nucleotides[seq[i+1]].second += c;
  if (l1==1 && l2==1) score_internal_1x1_nucleotides[seq[i+1]][seq[j]].second += c;
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreJunctionA(const std::vector<int>& seq, int i, int j) const
{
  const int L=seq.size()-1;
  return score_helix_closing[seq[i]][seq[j+1]].first
    * (i<L ? score_dangle_left[seq[i]][seq[j+1]][seq[i+1]].first : ProbT(1))
    * (j>0 ? score_dangle_right[seq[i]][seq[j+1]][seq[j]].first : ProbT(1));
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountJunctionA(const std::vector<int>& seq, int i, int j, CountT c)
{
  const int L=seq.size()-1;
  score_helix_closing[seq[i]][seq[j+1]].second += c;
  if (i<L) score_dangle_left[seq[i]][seq[j+1]][seq[i+1]].second += c;
  if (j>0) score_dangle_right[seq[i]][seq[j+1]][seq[j]].second += c;
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreJunctionB(const std::vector<int>& seq, int i, int j) const
{
  return score_helix_closing[seq[i]][seq[j+1]].first
    * score_terminal_mismatch[seq[i]][seq[j+1]][seq[i+1]][seq[j]].first;
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountJunctionB(const std::vector<int>& seq, int i, int j, CountT c)
{
  score_helix_closing[seq[i]][seq[j+1]].second += c;
  score_terminal_mismatch[seq[i]][seq[j+1]][seq[i+1]][seq[j]].second += c;
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreMultiPaired() const
{
  return score_multi_paired.first;
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountMultiPaired(CountT c)
{
  score_multi_paired.second += c;
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreMultiBase() const
{
  return score_multi_base.first;
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountMultiBase(CountT c)
{
  score_multi_base.second += c;
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreMultiUnpaired(const std::vector<int>& seq, int i) const
{
  return score_multi_unpaired.first * ScoreUnpairedPosition(seq, i);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountMultiUnpaired(const std::vector<int>& seq, int i, CountT c )
{
  score_multi_unpaired.second += c;
  CountUnpairedPosition(seq, i, c);
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreExternalPaired() const
{
  return score_external_paired.first;
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountExternalPaired(CountT c)
{
  score_external_paired.second += c;
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreExternalUnpaired(const std::vector<int>& seq, int i) const
{
  return score_external_unpaired.first * ScoreUnpairedPosition(seq, i);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountExternalUnpaired(const std::vector<int>& seq, int i, CountT c)
{
  score_external_unpaired.second += c;
  CountUnpairedPosition(seq, i, c);
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreUnpaired(const std::vector<int>& seq, int i, int j) const
{
  return ProbT(1);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountUnpaired(const std::vector<int>& seq, int i, int j, CountT c)
{
}

template <class ProbT, class CountT>
ProbT
ScoringModel<ProbT,CountT>::
ScoreUnpairedPosition(const std::vector<int>& seq, int i) const
{
  return ProbT(1);
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
CountUnpairedPosition(const std::vector<int>& seq, int i, CountT c) const
{
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
InitializeCache()
{
  cache_score_hairpin_length[0].first = score_hairpin_length_at_least[0].first;
  for (int i=1; i <= D_MAX_HAIRPIN_LENGTH; i++)
    cache_score_hairpin_length[i].first = cache_score_hairpin_length[i-1].first * score_hairpin_length_at_least[i].first;

  ProbT temp_cache_score_bulge_length[D_MAX_BULGE_LENGTH+1];
  temp_cache_score_bulge_length[0] = score_bulge_length_at_least[0].first;
  for (int i=1; i<=D_MAX_BULGE_LENGTH; i++)
    temp_cache_score_bulge_length[i] = temp_cache_score_bulge_length[i-1] * score_bulge_length_at_least[i].first;
    
  ProbT temp_cache_score_internal_length[D_MAX_INTERNAL_LENGTH+1];
  temp_cache_score_internal_length[0] = score_internal_length_at_least[0].first;
  for (int i=1; i<=D_MAX_INTERNAL_LENGTH; i++)
    temp_cache_score_internal_length[i] = temp_cache_score_internal_length[i-1] * score_internal_length_at_least[i].first;
    
  ProbT temp_cache_score_internal_symmetric_length[D_MAX_INTERNAL_SYMMETRIC_LENGTH+1];
  temp_cache_score_internal_symmetric_length[0] = score_internal_symmetric_length_at_least[0].first;
  for (int i=1; i<=D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
    temp_cache_score_internal_symmetric_length[i] = temp_cache_score_internal_symmetric_length[i-1] * score_internal_symmetric_length_at_least[i].first;
    
  ProbT temp_cache_score_internal_asymmetry[D_MAX_INTERNAL_ASYMMETRY+1];
  temp_cache_score_internal_asymmetry[0] = score_internal_asymmetry_at_least[0].first;
  for (int i=1; i<=D_MAX_INTERNAL_ASYMMETRY; i++)
    temp_cache_score_internal_asymmetry[i] = temp_cache_score_internal_asymmetry[i-1] * score_internal_asymmetry_at_least[i].first;
    
  // precompute score for single-branch loops of length l1 and l2
  for (int l1=0; l1<=C_MAX_SINGLE_LENGTH; l1++)
  {
    for (int l2=0; l1+l2<=C_MAX_SINGLE_LENGTH; l2++)
    {
      cache_score_single[l1][l2].first = ProbT(1);

      // skip over stacking pairs
      if (l1==0 && l2==0) continue;
        
      // consider bulge loops
      if (l1==0 ||l2 ==0)
      {
        cache_score_single[l1][l2].first *= temp_cache_score_bulge_length[std::min(int(D_MAX_BULGE_LENGTH), l1+l2)];
      }
      // consider internal loops
      else
      {
        if (l1<=D_MAX_INTERNAL_EXPLICIT_LENGTH && l2<=D_MAX_INTERNAL_EXPLICIT_LENGTH)
          cache_score_single[l1][l2].first *= score_internal_explicit[l1][l2].first;
        cache_score_single[l1][l2].first *= temp_cache_score_internal_length[std::min(int(D_MAX_INTERNAL_LENGTH), l1+l2)];
        if (l1==l2)
          cache_score_single[l1][l2].first *= temp_cache_score_internal_symmetric_length[std::min(int(D_MAX_INTERNAL_SYMMETRIC_LENGTH), l1)];
        cache_score_single[l1][l2].first *= temp_cache_score_internal_asymmetry[std::min(int(D_MAX_INTERNAL_ASYMMETRY), abs(l1-l2))];
      }
    }
  }
}

template <class ProbT, class CountT>
void
ScoringModel<ProbT,CountT>::
FinalizeCounts()
{
  for (int i = 0; i <= D_MAX_HAIRPIN_LENGTH; i++)
    for (int j = i; j <= D_MAX_HAIRPIN_LENGTH; j++)
      score_hairpin_length_at_least[i].second += cache_score_hairpin_length[j].second;


  CountT temp_cache_counts_bulge_length[D_MAX_BULGE_LENGTH+1];
  std::fill(temp_cache_counts_bulge_length, temp_cache_counts_bulge_length + D_MAX_BULGE_LENGTH+1, CountT(0));
  CountT temp_cache_counts_internal_length[D_MAX_INTERNAL_LENGTH+1];
  std::fill(temp_cache_counts_internal_length, temp_cache_counts_internal_length + D_MAX_INTERNAL_LENGTH+1, CountT(0));

  CountT temp_cache_counts_internal_symmetric_length[D_MAX_INTERNAL_SYMMETRIC_LENGTH+1];
  std::fill(temp_cache_counts_internal_symmetric_length, temp_cache_counts_internal_symmetric_length + D_MAX_INTERNAL_SYMMETRIC_LENGTH+1, CountT(0));

  CountT temp_cache_counts_internal_asymmetry[D_MAX_INTERNAL_ASYMMETRY+1];
  std::fill(temp_cache_counts_internal_asymmetry, temp_cache_counts_internal_asymmetry + D_MAX_INTERNAL_ASYMMETRY+1, CountT(0));

  // compute contributions
  for (int l1=0; l1<=C_MAX_SINGLE_LENGTH; l1++)
  {
    for (int l2=0; l1+l2<=C_MAX_SINGLE_LENGTH; l2++)
    {
      // skip over stacking pairs
      if (l1==0 && l2==0) continue;
        
      // consider bulge loops
      if (l1==0 || l2==0)
      {
        temp_cache_counts_bulge_length[std::min(int(D_MAX_BULGE_LENGTH), l1+l2)] += cache_score_single[l1][l2].second;
      }
      // consider internal loops
      else
      {
        if (l1<=D_MAX_INTERNAL_EXPLICIT_LENGTH && l2<=D_MAX_INTERNAL_EXPLICIT_LENGTH)
          score_internal_explicit[l1][l2].second += cache_score_single[l1][l2].second;
        temp_cache_counts_internal_length[std::min(int(D_MAX_INTERNAL_LENGTH), l1+l2)] += cache_score_single[l1][l2].second;
        if (l1==l2)
          temp_cache_counts_internal_symmetric_length[std::min(int(D_MAX_INTERNAL_SYMMETRIC_LENGTH), l1)] += cache_score_single[l1][l2].second;
        temp_cache_counts_internal_asymmetry[std::min(int(D_MAX_INTERNAL_ASYMMETRY), abs(l1-l2))] += cache_score_single[l1][l2].second;
      }
    }
  }

  for (int i=0; i<=D_MAX_BULGE_LENGTH; i++)
    for (int j = i; j <= D_MAX_BULGE_LENGTH; j++)
      score_bulge_length_at_least[i].second += temp_cache_counts_bulge_length[j];
    
  for (int i=0; i<=D_MAX_INTERNAL_LENGTH; i++)
    for (int j=i; j<=D_MAX_INTERNAL_LENGTH; j++)
      score_internal_length_at_least[i].second += temp_cache_counts_internal_length[j];
    
  for (int i=0; i<=D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
    for (int j=i; j<=D_MAX_INTERNAL_SYMMETRIC_LENGTH; j++)
      score_internal_symmetric_length_at_least[i].second += temp_cache_counts_internal_symmetric_length[j];
    
  for (int i=0; i<=D_MAX_INTERNAL_ASYMMETRY; i++)
    for (int j=i; j<=D_MAX_INTERNAL_ASYMMETRY; j++)
      score_internal_asymmetry_at_least[i].second += temp_cache_counts_internal_asymmetry[j];
}

// instantiation
#include "log_value.h"
typedef LogValue<float> ProbT;
typedef float CountT;
template class ScoringModel<ProbT,CountT>;
