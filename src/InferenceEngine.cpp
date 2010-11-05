// InferenceEngine.cpp

#include <iostream>

#include <string>
#include <vector>
#include <queue>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "InferenceEngine.h"

typedef unsigned char BYTE;

// struct triple
template<typename T1, typename T2, typename T3>
struct triple {
  T1 first;
  T2 second;
  T3 third;

  // constructors
  triple()
    : first(), second(), third()
  {}
  triple(const T1 &first, const T2 &second, const T3 &third)
    : first(first), second(second), third(third)
  {}
  triple(const triple &rhs)
    : first(rhs.first), second(rhs.second), third(rhs.third)
  {}    
};

// utility function for making triples
template<typename T1, typename T2, typename T3> 
inline triple<T1,T2,T3> 
make_triple(T1 first, T2 second, T3 third)
{
  return triple<T1,T2,T3>(first, second, third);
}

// constructor
template < class ProbT, class CountT, class ScoringModel >
InferenceEngine<ProbT,CountT,ScoringModel>::
InferenceEngine(bool allow_noncomplementary)
  : ScoringModel(),
    allow_noncomplementary_(allow_noncomplementary),
    char_mapping_(256),
    is_complementary_(M+1)
{
  // precompute mapping from characters to index representation
  std::fill(char_mapping_.begin(), char_mapping_.end(), BYTE(alphabet.size()));
  for (size_t i=0; i<alphabet.size(); i++)
  {
    char_mapping_[BYTE(tolower(alphabet[i]))] = 
      char_mapping_[BYTE(toupper(alphabet[i]))] = i;
  }

  // precompute complementary pairings
  for (int i=0; i<=M; i++)
  {
    is_complementary_[i].resize(M+1);
    for (int j=0; j<=M; j++)
      is_complementary_[i][j] = 0;
  }
  
  is_complementary_[char_mapping_[BYTE('A')]][char_mapping_[BYTE('U')]] = 
    is_complementary_[char_mapping_[BYTE('U')]][char_mapping_[BYTE('A')]] = 
    is_complementary_[char_mapping_[BYTE('G')]][char_mapping_[BYTE('U')]] = 
    is_complementary_[char_mapping_[BYTE('U')]][char_mapping_[BYTE('G')]] = 
    is_complementary_[char_mapping_[BYTE('C')]][char_mapping_[BYTE('G')]] = 
    is_complementary_[char_mapping_[BYTE('G')]][char_mapping_[BYTE('C')]] = 1;
}

template < class ProbT, class CountT, class ScoringModel >
void
InferenceEngine<ProbT, CountT, ScoringModel>::
LoadSequence(const std::string& rna)
{
  const int L = rna.size();
  const int SIZE = (L+1)*(L+2)/2;

  seq_.resize(L+1);
  offset_.resize(L+1);
  allow_unpaired_position_.resize(L+1);
  allow_unpaired_.resize(SIZE);
  allow_paired_.resize(SIZE);

  // convert sequences to index representation
  seq_[0] = BYTE(alphabet.size());
  for (int i=1; i<=L; i++)
    seq_[i] = char_mapping_[BYTE(rna[i-1])];
  
  // compute indexing scheme for upper triangular arrays;
  // also allow each position to be unpaired by default, and
  // set the loss for each unpaired position to zero
  for (int i=0; i<=L; i++)
  {
    offset_[i] = ComputeRowOffset(i,L+1);
    allow_unpaired_position_[i] = 1;
  }

  // allow all ranges to be unpaired, and all pairs of letters
  // to be paired; set the respective losses to zero    
  for (int i=0; i<SIZE; i++)
  {
    allow_unpaired_[i] = 1;
    allow_paired_[i] = 1;
  }

  // prevent the non-letter before each sequence from pairing with anything;
  // also prevent each letter from pairing with itself
  for (int i=0; i<=L; i++)
  {
    allow_paired_[offset_[0]+i] = 0;
    allow_paired_[offset_[i]+i] = 0;
  }

  // enforce complementarity of base-pairings
  if (!allow_noncomplementary_)
  {
    // for each pair of non-complementary letters in the sequence, disallow the pairing
    for (int i=1; i<=L; i++)
    {
      for (int j=i+1; j<=L; j++)
      {
        if (!IsComplementary(i,j))
          allow_paired_[offset_[i]+j] = 0;
      }
    }
  }
}

template<class ProbT, class CountT, class ScoringModel>
bool
InferenceEngine<ProbT,CountT,ScoringModel>::
IsComplementary(int i, int j) const
{
  return is_complementary_[seq_[i]][seq_[j]];
}

//////////////////////////////////////////////////////////////////////
// ComputeRowOffset()
//
// Consider an N x N upper triangular matrix whose elements are
// stored in a one-dimensional flat array using the following
// row-major indexing scheme:
//
//     0  1  2  3     <-- row 0
//        4  5  6     <-- row 1
//           7 [8]    <-- row 2
//              9     <-- row 3
//
// Assuming 0-based indexing, this function computes offset[i]
// for the ith row such that offset[i]+j is the index of the
// (i,j)th element of the upper triangular matrix in the flat
// array.
//
// For example, offset[2] = 5, so the (2,3)th element of the
// upper triangular matrix (marked in the picture above) can be 
// found at position offset[2]+3 = 5+3 = 8 in the flat array.
//////////////////////////////////////////////////////////////////////

template <class ProbT, class CountT, class ScoringModel>
int
InferenceEngine<ProbT,CountT,ScoringModel>::
ComputeRowOffset(int i, int N) const
{
    // equivalent to:
    //   return N*(N+1)/2 - (N-i)*(N-i+1)/2 - i;
    return i*(N+N-i-1)/2;
}

template <class ProbT, class CountT, class ScoringModel>
void
InferenceEngine<ProbT,CountT,ScoringModel>::
UseConstraints(const std::vector<int> &true_mapping)
{
  const int L = seq_.size()-1;
  assert(L+1==(int)true_mapping.size());
  
  // determine whether we allow each position to be unpaired
  for (int i=1; i<=L; i++)
  {
    allow_unpaired_position_[i] =
      (true_mapping[i]==-1 ||   // unknown
       true_mapping[i]==0);   // unpaired
  }

  // determine whether we allow ranges of positions to be unpaired;
  // also determine which base-pairings we allow
  for (int i=0; i<=L; i++)
  {
    allow_unpaired_[offset_[i]+i] = 1;
    allow_paired_[offset_[i]+i] = 0;
    for (int j=i+1; j<=L; j++)
    {
      allow_unpaired_[offset_[i]+j] = 
        allow_unpaired_[offset_[i]+j-1] && 
        allow_unpaired_position_[j];
      allow_paired_[offset_[i]+j] =
        (i>0 &&
         (true_mapping[i]==-1 || true_mapping[i]==j) &&
         (true_mapping[j]==-1 || true_mapping[j]==i) &&
         (allow_noncomplementary_ || IsComplementary(i,j)));
    }
  }
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::EncodeTraceback()
// InferenceEngine::DecodeTraceback()
//
// Encode and decode traceback as an integer.  Here, i encodes
// a traceback type, and j encodes a length.
//////////////////////////////////////////////////////////////////////

template <class ProbT, class CountT, class ScoringModel>
inline int
InferenceEngine<ProbT,CountT,ScoringModel>::
EncodeTraceback(int i, int j) const
{
  return (j * NUM_TRACEBACK_TYPES) + i;
}

template <class ProbT, class CountT, class ScoringModel>
inline std::pair<int,int>
InferenceEngine<ProbT,CountT,ScoringModel>::
DecodeTraceback(int s) const
{
  return std::make_pair (s % NUM_TRACEBACK_TYPES, s / NUM_TRACEBACK_TYPES);
}

template <class ProbT, class CountT, class ScoringModel>
void
InferenceEngine<ProbT,CountT,ScoringModel>::
UpdateViterbi(ProbT& bs, int& bt, ProbT s, int t)
{
  if (bs<s)
  {
    bs=s;
    bt=t;
  }
}

template < class ProbT, class CountT, class ScoringModel >
ProbT
InferenceEngine<ProbT, CountT,ScoringModel>::
ComputeViterbi()
{
  const int L = seq_.size()-1;
  const int SIZE = (L+1)*(L+2)/2;
  
  F5v_.clear(); F5v_.resize(L+1, ProbT(0));
  FCv_.clear(); FCv_.resize(SIZE, ProbT(0));
  FMv_.clear(); FMv_.resize(SIZE, ProbT(0));
  FM1v_.clear(); FM1v_.resize(SIZE, ProbT(0));

  F5t_.clear(); F5t_.resize(L+1, -1);
  FCt_.clear(); FCt_.resize(SIZE, -1);
  FMt_.clear(); FMt_.resize(SIZE, -1);
  FM1t_.clear(); FM1t_.resize(SIZE, -1);

  for (int i=L; i>=0; --i)
  {
    for (int j=i; j<=L; ++j)
    {
      // FM2[i,j] = MAX(i<k<j : FM1[i,k] * FM[k,j])
      ProbT FM2v = ProbT(0);
      int FM2t = -1;
      for (int k=i+1; k<j; ++k)
        UpdateViterbi(FM2v, FM2t, FM1v_[offset_[i]+k] * FMv_[offset_[k]+j], k);

      // FC[i,j] = MAX [ScoreHairpin(i,j),
      //                MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1]),
      //                ScoreJunctionA(i,j) * a * c * MAX (i<k<j : FM1[i,k] * FM[k,j])]
      if (0<i && j<L && allow_paired_[offset_[i]+j+1])
      {
        // compute ScoreHairpin(i,j)
        if (allow_unpaired_[offset_[i]+j] && j-i>=C_MIN_HAIRPIN_LENGTH)
          UpdateViterbi(FCv_[offset_[i]+j], FCt_[offset_[i]+j],
                        ScoringModel::ScoreFCHairpin(seq_,i,j),
                        EncodeTraceback(TB_FC_HAIRPIN,0));
        
        // compute MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1])
        for (int p=i; p<=std::min(i+C_MAX_SINGLE_LENGTH,j); ++p)
        {
          if (p>i && !allow_unpaired_position_[p]) break;
          int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
          for (int q=j; q>=q_min; --q)
          {
            if (q<j && !allow_unpaired_position_[q+1]) break;
            if (!allow_paired_[offset_[p+1]+q]) continue;
            UpdateViterbi(FCv_[offset_[i]+j], FCt_[offset_[i]+j],
                          FCv_[offset_[p+1]+q-1] * ScoringModel::ScoreFCSingle(seq_,i,j,p,q),
                          EncodeTraceback(TB_FC_SINGLE,(p-i)*(C_MAX_SINGLE_LENGTH+1)+j-q));
          }
        }

        // compute MAX (i<k<j : FM1[i,k] * FM[k,j] * ScoreJunctionA(i,j) * a * c)
        UpdateViterbi(FCv_[offset_[i]+j], FCt_[offset_[i]+j],
                      FM2v * ScoringModel::ScoreFCBifurcation(seq_,i,j),
                      EncodeTraceback(TB_FC_BIFURCATION, FM2t));
      }

      // FM1[i,j] = MAX [FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)  if i+2<=j,
      //                 FM1[i+1,j] * b                                          if i+2<=j]
      if (0<i && i+2<=j && j<L)
      {
        // compute FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)
        if (allow_paired_[offset_[i+1]+j])
        {
          UpdateViterbi(FM1v_[offset_[i]+j], FM1t_[offset_[i]+j], 
                        FCv_[offset_[i+1]+j-1] * ScoringModel::ScoreFM1Paired(seq_,i,j),
                        EncodeTraceback(TB_FM1_PAIRED, 0));
        }

        // compute FM1[i+1,j] * b
        if (allow_unpaired_position_[i+1])
        {
          UpdateViterbi(FM1v_[offset_[i]+j], FM1t_[offset_[i]+j], 
                        FM1v_[offset_[i+1]+j] * ScoringModel::ScoreFM1Unpaired(seq_,i,j),
                        EncodeTraceback(TB_FM1_UNPAIRED,0));
        }
      }

      // FM[i,j] = MAX [MAX (i<k<j : FM1[i,k] * FM[k,j]),
      //                FM[i,j-1] * b,
      //                FM1[i,j]]
      if (0<i && i+2<=j && j<L)
      {
        // compute MAX (i<k<j : FM1[i,k] * FM[k,j])
        UpdateViterbi(FMv_[offset_[i]+j], FMt_[offset_[i]+j],
                      FM2v * ScoringModel::ScoreFMBifurcation(seq_,i,j),
                      EncodeTraceback(TB_FM_BIFURCATION,FM2t));

        // compute FM[i,j-1] * b
        if (allow_unpaired_position_[j])
        {
          UpdateViterbi(FMv_[offset_[i]+j], FMt_[offset_[i]+j],
                        FMv_[offset_[i]+j-1] * ScoringModel::ScoreFMUnpaired(seq_,i,j),
                        EncodeTraceback(TB_FM_UNPAIRED,0));
        }

        // compute FM1[i,j]
        UpdateViterbi(FMv_[offset_[i]+j], FMt_[offset_[i]+j],
                      FM1v_[offset_[i]+j] * ScoringModel::ScoreFMtoFM1(seq_,i,j),
                      EncodeTraceback(TB_FM_FM1,0));
      }
    }
  }

  F5v_[0] = ProbT(1);
  F5t_[0] = EncodeTraceback(TB_F5_ZERO,0);
  for (int j=1; j<=L; ++j)
  {
    // F5[j] = MAX [F5[j-1] * ScoreExternalUnpaired(),
    //              MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k)

    // compute F5[j-1] * ScoreExternalUnpaired()
    if (allow_unpaired_position_[j])
    {
      UpdateViterbi(F5v_[j], F5t_[j], 
                    F5v_[j-1] * ScoringModel::ScoreF5Unpaired(seq_,j),
                    EncodeTraceback(TB_F5_UNPAIRED,0));
    }
    
    // compute MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k))
    for (int k = 0; k < j; k++)
    {
      if (allow_paired_[offset_[k+1]+j])
      {
        UpdateViterbi(F5v_[j], F5t_[j], 
                      F5v_[k] * FCv_[offset_[k+1]+j-1] * ScoringModel::ScoreF5Bifurcation(seq_,k,j),
                      EncodeTraceback(TB_F5_BIFURCATION,k));
      }
    }
  }

  return F5v_[L];
}

template < class ProbT, class CountT, class ScoringModel >
std::vector<int>
InferenceEngine<ProbT, CountT, ScoringModel>::
DecodeViterbi()
{
  const int L = seq_.size()-1;

  std::vector<int> solution(L+1, 0);
  solution[0] = -1;

  std::queue<triple<int *,int,int> > traceback_queue;
  traceback_queue.push(make_triple(&F5t_[0], 0, L));
    
  while (!traceback_queue.empty())
  {
    triple<int *,int,int> t = traceback_queue.front();
    traceback_queue.pop();
    const int *V = t.first;
    const int i = t.second;
    const int j = t.third;
    
    std::pair<int,int> traceback = DecodeTraceback(V == &F5t_[0] ? V[j] : V[offset_[i]+j]);
    
    switch (traceback.first)
    {
      case TB_FC_HAIRPIN: 
        break;
      case TB_FC_SINGLE: 
      {
        const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
        const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);
        solution[p+1] = q;
        solution[q] = p+1;
        traceback_queue.push(make_triple(&FCt_[0], p+1, q-1));
      }
      break;
      case TB_FC_BIFURCATION:
      {
        const int k = traceback.second;
        traceback_queue.push(make_triple(&FM1t_[0], i, k));
        traceback_queue.push(make_triple(&FMt_[0], k, j));
      }
      break;
      case TB_FM1_PAIRED:
      {
        solution[i+1] = j;
        solution[j] = i+1;
        traceback_queue.push(make_triple(&FCt_[0], i+1, j-1));
      }
      break;
      case TB_FM1_UNPAIRED:
      {
        traceback_queue.push(make_triple(&FM1t_[0], i+1, j));
      }
      break;
      case TB_FM_BIFURCATION:
      {
        const int k = traceback.second;
        traceback_queue.push(make_triple(&FM1t_[0], i, k));
        traceback_queue.push(make_triple(&FMt_[0], k, j));
      }
      break;
      case TB_FM_UNPAIRED:
      {
        traceback_queue.push(make_triple(&FMt_[0], i, j-1));
      }
      break;
      case TB_FM_FM1: 
      {
        traceback_queue.push(make_triple(&FM1t_[0], i, j));
      }
      break;
      case TB_F5_ZERO:
        break;
      case TB_F5_UNPAIRED:
      {
        traceback_queue.push(make_triple(&F5t_[0], 0, j-1));
      }
      break;
      case TB_F5_BIFURCATION:
      {
        const int k = traceback.second;
        solution[k+1] = j;
        solution[j] = k+1;
        traceback_queue.push(make_triple(&F5t_[0], 0, k));
        traceback_queue.push(make_triple(&FCt_[0], k+1, j-1));
      }
      break;
      default:
        assert(!"Bad traceback.");
    }
  }
    
  return solution;
}

template < class ProbT, class CountT, class ScoringModel >
ProbT
InferenceEngine<ProbT, CountT, ScoringModel>::
ComputeInside()
{
  const int L = seq_.size()-1;
  const int SIZE = (L+1)*(L+2)/2;

  F5i_.clear(); F5i_.resize(L+1, ProbT(0));
  FCi_.clear(); FCi_.resize(SIZE, ProbT(0));
  FMi_.clear(); FMi_.resize(SIZE, ProbT(0));
  FM1i_.clear(); FM1i_.resize(SIZE, ProbT(0));

  for (int i=L; i>=0; --i)
  {
    for (int j=i; j<=L; ++j)
    {
      // FM2[i,j] = MAX(i<k<j : FM1[i,k] * FM[k,j])
      ProbT FM2i = ProbT(0);
      for (int k=i+1; k<j; ++k)
        FM2i += FM1i_[offset_[i]+k] * FMi_[offset_[k]+j];

      // FC[i,j] = MAX [ScoreHairpin(i,j),
      //                MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1]),
      //                ScoreJunctionA(i,j) * a * c * MAX (i<k<j : FM1[i,k] * FM[k,j])]
      if (0<i && j<L && allow_paired_[offset_[i]+j+1])
      {
        // compute ScoreHairpin(i,j)
        if (allow_unpaired_[offset_[i]+j] && j-i>=C_MIN_HAIRPIN_LENGTH)
          FCi_[offset_[i]+j] += ScoringModel::ScoreFCHairpin(seq_,i,j);
        
        // compute MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1])
        for (int p=i; p<=std::min(i+C_MAX_SINGLE_LENGTH,j); ++p)
        {
          if (p>i && !allow_unpaired_position_[p]) break;
          int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
          for (int q=j; q>=q_min; --q)
          {
            if (q<j && !allow_unpaired_position_[q+1]) break;
            if (!allow_paired_[offset_[p+1]+q]) continue;
            FCi_[offset_[i]+j] += FCi_[offset_[p+1]+q-1] * ScoringModel::ScoreFCSingle(seq_,i,j,p,q);
          }
        }

        // compute MAX (i<k<j : FM1[i,k] * FM[k,j] * ScoreJunctionA(i,j) * a * c)
        FCi_[offset_[i]+j] += FM2i * ScoringModel::ScoreFCBifurcation(seq_,i,j);
      }

      // FM1[i,j] = MAX [FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)  if i+2<=j,
      //                 FM1[i+1,j] * b                                          if i+2<=j]
      if (0<i && i+2<=j && j<L)
      {
        // compute FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)
        if (allow_paired_[offset_[i+1]+j])
          FM1i_[offset_[i]+j] += FCi_[offset_[i+1]+j-1] * ScoringModel::ScoreFM1Paired(seq_,i,j);

        // compute FM1[i+1,j] * b
        if (allow_unpaired_position_[i+1])
          FM1i_[offset_[i]+j] += FM1i_[offset_[i+1]+j] * ScoringModel::ScoreFM1Unpaired(seq_,i,j);
      }

      // FM[i,j] = MAX [MAX (i<k<j : FM1[i,k] * FM[k,j]),
      //                FM[i,j-1] * b,
      //                FM1[i,j]]
      if (0<i && i+2<=j && j<L)
      {
        // compute MAX (i<k<j : FM1[i,k] * FM[k,j])
        FMi_[offset_[i]+j] += FM2i * ScoringModel::ScoreFMBifurcation(seq_,i,j);

        // compute FM[i,j-1] * b
        if (allow_unpaired_position_[j])
          FMi_[offset_[i]+j] += FMi_[offset_[i]+j-1] * ScoringModel::ScoreFMUnpaired(seq_,i,j);

        // compute FM1[i,j]
        FMi_[offset_[i]+j] += FM1i_[offset_[i]+j] * ScoringModel::ScoreFMtoFM1(seq_,i,j);
      }
    }
  }

  F5i_[0] = ProbT(1);
  for (int j=1; j<=L; ++j)
  {
    // F5[j] = MAX [F5[j-1] * ScoreExternalUnpaired(),
    //              MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k)

    // compute F5[j-1] * ScoreExternalUnpaired()
    if (allow_unpaired_position_[j])
      F5i_[j] += F5i_[j-1] * ScoringModel::ScoreF5Unpaired(seq_,j);
    
    // compute MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k))
    for (int k = 0; k < j; k++)
    {
      if (allow_paired_[offset_[k+1]+j])
        F5i_[j] += F5i_[k] * FCi_[offset_[k+1]+j-1] * ScoringModel::ScoreF5Bifurcation(seq_,k,j);
    }
  }

  return F5i_[L];
}

template < class ProbT, class CountT, class ScoringModel >
ProbT
InferenceEngine<ProbT, CountT, ScoringModel>::
ComputeInside2(const std::vector<CountT>& weight)
{
  const int L = seq_.size()-1;
  const int SIZE = (L+1)*(L+2)/2;

  F5j_.clear(); F5j_.resize(L+1, ProbT(0));
  FCj_.clear(); FCj_.resize(SIZE, ProbT(0));
  FMj_.clear(); FMj_.resize(SIZE, ProbT(0));
  FM1j_.clear(); FM1j_.resize(SIZE, ProbT(0));

  for (int i=L; i>=0; --i)
  {
    for (int j=i; j<=L; ++j)
    {
      // FM2[i,j] = MAX(i<k<j : FM1[i,k] * FM[k,j])
      ProbT FM2j = ProbT(0);
      for (int k=i+1; k<j; ++k)
        FM2j += FM1j_[offset_[i]+k] * FMi_[offset_[k]+j] +
          FM1i_[offset_[i]+k] * FMj_[offset_[k]+j];
 
      // FC[i,j] = MAX [ScoreHairpin(i,j),
      //                MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1]),
      //                ScoreJunctionA(i,j) * a * c * MAX (i<k<j : FM1[i,k] * FM[k,j])]
      if (0<i && j<L && allow_paired_[offset_[i]+j+1])
      {
        // compute ScoreHairpin(i,j)
        
        // compute MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1])
        for (int p=i; p<=std::min(i+C_MAX_SINGLE_LENGTH,j); ++p)
        {
          if (p>i && !allow_unpaired_position_[p]) break;
          int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
          for (int q=j; q>=q_min; --q)
          {
            if (q<j && !allow_unpaired_position_[q+1]) break;
            if (!allow_paired_[offset_[p+1]+q]) continue;
            FCj_[offset_[i]+j] +=
              (FCj_[offset_[p+1]+q-1] + FCi_[offset_[p+1]+q-1] * weight[offset_[p+1]+q])
              * ScoringModel::ScoreFCSingle(seq_,i,j,p,q);
          }
        }

        // compute MAX (i<k<j : FM1[i,k] * FM[k,j] * ScoreJunctionA(i,j) * a * c)
        FCj_[offset_[i]+j] += FM2j * ScoringModel::ScoreFCBifurcation(seq_,i,j);
      }

      // FM1[i,j] = MAX [FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)  if i+2<=j,
      //                 FM1[i+1,j] * b                                          if i+2<=j]
      if (0<i && i+2<=j && j<L)
      {
        // compute FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)
        if (allow_paired_[offset_[i+1]+j])
          FM1j_[offset_[i]+j] +=
            (FCj_[offset_[i+1]+j-1] + FCi_[offset_[i+1]+j-1] * weight[offset_[i+1]+j])
            * ScoringModel::ScoreFM1Paired(seq_,i,j);

        // compute FM1[i+1,j] * b
        if (allow_unpaired_position_[i+1])
          FM1j_[offset_[i]+j] += FM1j_[offset_[i+1]+j] * ScoringModel::ScoreFM1Unpaired(seq_,i,j);
      }

      // FM[i,j] = MAX [MAX (i<k<j : FM1[i,k] * FM[k,j]),
      //                FM[i,j-1] * b,
      //                FM1[i,j]]
      if (0<i && i+2<=j && j<L)
      {
        // compute MAX (i<k<j : FM1[i,k] * FM[k,j])
        FMj_[offset_[i]+j] += FM2j * ScoringModel::ScoreFMBifurcation(seq_,i,j);

        // compute FM[i,j-1] * b
        if (allow_unpaired_position_[j])
          FMj_[offset_[i]+j] += FMj_[offset_[i]+j-1] * ScoringModel::ScoreFMUnpaired(seq_,i,j);

        // compute FM1[i,j]
        FMj_[offset_[i]+j] += FM1j_[offset_[i]+j] * ScoringModel::ScoreFMtoFM1(seq_,i,j);
      }
    }
  }

  F5j_[0] = ProbT(0);
  for (int j=1; j<=L; ++j)
  {
    // F5[j] = MAX [F5[j-1] * ScoreExternalUnpaired(),
    //              MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k)

    // compute F5[j-1] * ScoreExternalUnpaired()
    if (allow_unpaired_position_[j])
      F5j_[j] += F5j_[j-1] * ScoringModel::ScoreF5Unpaired(seq_,j);
    
    // compute MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k))
    for (int k = 0; k < j; k++)
    {
      if (allow_paired_[offset_[k+1]+j])
        F5j_[j] +=
          (F5j_[k] * FCi_[offset_[k+1]+j-1] + F5i_[k] * FCj_[offset_[k+1]+j-1] +
           F5i_[k] * FCi_[offset_[k+1]+j-1] * weight[offset_[k+1]+j])
          * ScoringModel::ScoreF5Bifurcation(seq_,k,j);
    }
  }

  return F5j_[L];
}

template < class ProbT, class CountT, class ScoringModel >
ProbT
InferenceEngine<ProbT, CountT, ScoringModel>::
ComputeOutside()
{
  const int L = seq_.size()-1;
  const int SIZE = (L+1)*(L+2)/2;

  F5o_.clear(); F5o_.resize(L+1, ProbT(0));
  FCo_.clear(); FCo_.resize(SIZE, ProbT(0));
  FMo_.clear(); FMo_.resize(SIZE, ProbT(0));
  FM1o_.clear(); FM1o_.resize(SIZE, ProbT(0));

  F5o_[L] = ProbT(1);
  for (int j=L; j>=1; --j)
  {
    // F5[j] = MAX [F5[j-1] * ScoreExternalUnpaired(),
    //              MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k)

    // compute F5[j-1] * ScoreExternalUnpaired()
    if (allow_unpaired_position_[j])
      F5o_[j-1] += F5o_[j] * ScoringModel::ScoreF5Unpaired(seq_,j);
    
    // compute MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k))
    for (int k = 0; k < j; k++)
    {
      if (allow_paired_[offset_[k+1]+j])
      {
        ProbT temp = F5o_[j] * ScoringModel::ScoreF5Bifurcation(seq_,k,j);
        F5o_[k] += temp * FCi_[offset_[k+1]+j-1];
        FCo_[offset_[k+1]+j-1] += temp * F5i_[k];
      }
    }
  }

  for (int i=0; i<=L; ++i)
  {
    for (int j=L; j>=i; --j)
    {
      ProbT FM2o = ProbT(0);

      // FM[i,j] = MAX [MAX (i<k<j : FM1[i,k] * FM[k,j]),
      //                FM[i,j-1] * b,
      //                FM1[i,j]]
      if (0<i && i+2<=j && j<L)
      {
        // compute MAX (i<k<j : FM1[i,k] * FM[k,j])
        FM2o += FMo_[offset_[i]+j] * ScoringModel::ScoreFMBifurcation(seq_,i,j);

        // compute FM[i,j-1] * b
        if (allow_unpaired_position_[j])
          FMo_[offset_[i]+j-1] += FMo_[offset_[i]+j] * ScoringModel::ScoreFMUnpaired(seq_,i,j);

        // compute FM1[i,j]
        FM1o_[offset_[i]+j] += FMo_[offset_[i]+j] * ScoringModel::ScoreFMtoFM1(seq_,i,j);
      }

      // FM1[i,j] = MAX [FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)  if i+2<=j,
      //                 FM1[i+1,j] * b                                          if i+2<=j]
      if (0<i && i+2<=j && j<L)
      {
        // compute FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)
        if (allow_paired_[offset_[i+1]+j])
          FCo_[offset_[i+1]+j-1] += FM1o_[offset_[i]+j] * ScoringModel::ScoreFM1Paired(seq_,i,j);

        // compute FM1[i+1,j] * b
        if (allow_unpaired_position_[i+1])
          FM1o_[offset_[i+1]+j] += FM1o_[offset_[i]+j] * ScoringModel::ScoreFM1Unpaired(seq_,i,j);
      }

      // FC[i,j] = MAX [ScoreHairpin(i,j),
      //                MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1]),
      //                ScoreJunctionA(i,j) * a * c * MAX (i<k<j : FM1[i,k] * FM[k,j])]
      if (0<i && j<L && allow_paired_[offset_[i]+j+1])
      {
        // compute ScoreHairpin(i,j) -- do nothing

        // compute MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1])
        for (int p=i; p<=std::min(i+C_MAX_SINGLE_LENGTH,j); ++p)
        {
          if (p>i && !allow_unpaired_position_[p]) break;
          int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
          for (int q=j; q>=q_min; --q)
          {
            if (q<j && !allow_unpaired_position_[q+1]) break;
            if (!allow_paired_[offset_[p+1]+q]) continue;
            FCo_[offset_[p+1]+q-1] += FCo_[offset_[i]+j] * ScoringModel::ScoreFCSingle(seq_,i,j,p,q);
          }
        }

        // compute MAX (i<k<j : FM1[i,k] * FM[k,j] * ScoreJunctionA(i,j) * a * c)
        FM2o += FCo_[offset_[i]+j] * ScoringModel::ScoreFCBifurcation(seq_,i,j);
      }

      // FM2[i,j] = SUM (i<k<j : FM1[i,k] * FM[k,j])
      for (int k = i+1; k < j; k++)
      {
        FM1o_[offset_[i]+k] += FM2o * FMi_[offset_[k]+j];
        FMo_[offset_[k]+j] += FM2o * FM1i_[offset_[i]+k];
      }
    }
  }

  return F5o_[0];
}

template < class ProbT, class CountT, class ScoringModel >
ProbT
InferenceEngine<ProbT, CountT, ScoringModel>::
ComputeOutside2(const std::vector<CountT>& weight)
{
  const int L = seq_.size()-1;
  const int SIZE = (L+1)*(L+2)/2;

  F5p_.clear(); F5p_.resize(L+1, ProbT(0));
  FCp_.clear(); FCp_.resize(SIZE, ProbT(0));
  FMp_.clear(); FMp_.resize(SIZE, ProbT(0));
  FM1p_.clear(); FM1p_.resize(SIZE, ProbT(0));

  F5p_[L] = ProbT(0);
  for (int j=L; j>=1; --j)
  {
    // F5[j] = MAX [F5[j-1] * ScoreExternalUnpaired(),
    //              MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k)

    // compute F5[j-1] * ScoreExternalUnpaired()
    if (allow_unpaired_position_[j])
      F5p_[j-1] += F5p_[j] * ScoringModel::ScoreF5Unpaired(seq_,j);
    
    // compute MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k))
    for (int k = 0; k < j; k++)
    {
      if (allow_paired_[offset_[k+1]+j])
      {
        F5p_[k] +=
          (F5p_[j] * FCi_[offset_[k+1]+j-1] + F5o_[j] * FCj_[offset_[k+1]+j-1] +
           F5o_[j] * FCi_[offset_[k+1]+j-1] * weight[offset_[k+1]+j])
          * ScoringModel::ScoreF5Bifurcation(seq_,k,j);
        FCp_[offset_[k+1]+j-1] +=
          (F5p_[j] * F5i_[k] + F5o_[j] * F5j_[k] +  F5o_[j] * F5i_[k] * weight[offset_[k+1]+j])
          * ScoringModel::ScoreF5Bifurcation(seq_,k,j);
      }
    }
  }

  for (int i=0; i<=L; ++i)
  {
    for (int j=L; j>=i; --j)
    {
      ProbT FM2o = ProbT(0);
      ProbT FM2p = ProbT(0);

      // FM[i,j] = MAX [MAX (i<k<j : FM1[i,k] * FM[k,j]),
      //                FM[i,j-1] * b,
      //                FM1[i,j]]
      if (0<i && i+2<=j && j<L)
      {
        // compute MAX (i<k<j : FM1[i,k] * FM[k,j])
        FM2o += FMo_[offset_[i]+j] * ScoringModel::ScoreFMBifurcation(seq_,i,j);
        FM2p += FMp_[offset_[i]+j] * ScoringModel::ScoreFMBifurcation(seq_,i,j);

        // compute FM[i,j-1] * b
        if (allow_unpaired_position_[j])
          FMp_[offset_[i]+j-1] += FMp_[offset_[i]+j] * ScoringModel::ScoreFMUnpaired(seq_,i,j);

        // compute FM1[i,j]
        FM1p_[offset_[i]+j] += FMp_[offset_[i]+j] * ScoringModel::ScoreFMtoFM1(seq_,i,j);
      }

      // FM1[i,j] = MAX [FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)  if i+2<=j,
      //                 FM1[i+1,j] * b                                          if i+2<=j]
      if (0<i && i+2<=j && j<L)
      {
        // compute FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)
        if (allow_paired_[offset_[i+1]+j])
          FCp_[offset_[i+1]+j-1] +=
            (FM1p_[offset_[i]+j] + FM1o_[offset_[i]+j] * weight[offset_[i+1]+j])
            * ScoringModel::ScoreFM1Paired(seq_,i,j);

        // compute FM1[i+1,j] * b
        if (allow_unpaired_position_[i+1])
          FM1p_[offset_[i+1]+j] += FM1p_[offset_[i]+j] * ScoringModel::ScoreFM1Unpaired(seq_,i,j);
      }

      // FC[i,j] = MAX [ScoreHairpin(i,j),
      //                MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1]),
      //                ScoreJunctionA(i,j) * a * c * MAX (i<k<j : FM1[i,k] * FM[k,j])]
      if (0<i && j<L && allow_paired_[offset_[i]+j+1])
      {
        // compute ScoreHairpin(i,j) -- do nothing

        // compute MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1])
        for (int p=i; p<=std::min(i+C_MAX_SINGLE_LENGTH,j); ++p)
        {
          if (p>i && !allow_unpaired_position_[p]) break;
          int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
          for (int q=j; q>=q_min; --q)
          {
            if (q<j && !allow_unpaired_position_[q+1]) break;
            if (!allow_paired_[offset_[p+1]+q]) continue;
            FCp_[offset_[p+1]+q-1] +=
              (FCp_[offset_[i]+j] + FCo_[offset_[i]+j] * weight[offset_[p+1]+q])
              * ScoringModel::ScoreFCSingle(seq_,i,j,p,q);
          }
        }

        // compute MAX (i<k<j : FM1[i,k] * FM[k,j] * ScoreJunctionA(i,j) * a * c)
        FM2o += FCo_[offset_[i]+j] * ScoringModel::ScoreFCBifurcation(seq_,i,j);
        FM2p += FCp_[offset_[i]+j] * ScoringModel::ScoreFCBifurcation(seq_,i,j);
      }

      // FM2[i,j] = SUM (i<k<j : FM1[i,k] * FM[k,j])
      for (int k = i+1; k < j; k++)
      {
        FM1p_[offset_[i]+k] += FM2o * FMj_[offset_[k]+j] + FM2p * FMi_[offset_[k]+j];
        FMp_[offset_[k]+j] += FM2o * FM1j_[offset_[i]+k] + FM2p * FM1i_[offset_[i]+k];
      }
    }
  }

  return F5p_[0];
}

template < class ProbT, class CountT, class ScoringModel >
void
InferenceEngine<ProbT, CountT, ScoringModel>::
ComputePosterior()
{
  const int L = seq_.size()-1;
  const int SIZE = (L+1)*(L+2)/2;

  posterior_.clear();
  posterior_.resize(SIZE, CountT(0));

  const ProbT Z = ComputePartitionCoefficient();

  for (int i=L; i>=0; --i)
  {
    for (int j=i; j<=L; ++j)
    {
      // FM2[i,j] = MAX(i<k<j : FM1[i,k] * FM[k,j])
#if 0 // FM2 is not required
      ProbT FM2i = ProbT(0);
      for (int k=i+1; k<j; ++k)
        FM2i += FM1i_[offset_[i]+k] * FMi_[offset_[k]+j];
#endif

      // FC[i,j] = MAX [ScoreHairpin(i,j),
      //                MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1]),
      //                ScoreJunctionA(i,j) * a * c * MAX (i<k<j : FM1[i,k] * FM[k,j])]
      if (0<i && j<L && allow_paired_[offset_[i]+j+1])
      {
        // compute ScoreHairpin(i,j)
        
        // compute MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1])
        for (int p=i; p<=std::min(i+C_MAX_SINGLE_LENGTH,j); ++p)
        {
          if (p>i && !allow_unpaired_position_[p]) break;
          int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
          for (int q=j; q>=q_min; --q)
          {
            if (q<j && !allow_unpaired_position_[q+1]) break;
            if (!allow_paired_[offset_[p+1]+q]) continue;
            posterior_[offset_[i+1]+j] += FCi_[offset_[p+1]+q-1]
              * ScoringModel::ScoreFCSingle(seq_,i,j,p,q) * FCo_[offset_[i]+j] / Z;
          }
        }

        // compute MAX (i<k<j : FM1[i,k] * FM[k,j] * ScoreJunctionA(i,j) * a * c)
      }

      // FM1[i,j] = MAX [FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)  if i+2<=j,
      //                 FM1[i+1,j] * b                                          if i+2<=j]
      if (0<i && i+2<=j && j<L)
      {
        // compute FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)
        if (allow_paired_[offset_[i+1]+j])
          posterior_[offset_[i+1]+j] += FCi_[offset_[i+1]+j-1]
            * ScoringModel::ScoreFM1Paired(seq_,i,j) * FM1o_[offset_[i]+j] / Z;

        // compute FM1[i+1,j] * b
      }

      // FM[i,j] = MAX [MAX (i<k<j : FM1[i,k] * FM[k,j]),
      //                FM[i,j-1] * b,
      //                FM1[i,j]]

      // compute MAX (i<k<j : FM1[i,k] * FM[k,j])

      // compute FM[i,j-1] * b

      // compute FM1[i,j]
    }
  }

  for (int j=1; j<=L; ++j)
  {
    // F5[j] = MAX [F5[j-1] * ScoreExternalUnpaired(),
    //              MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k)

    // compute F5[j-1] * ScoreExternalUnpaired()
    
    // compute MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k))
    for (int k = 0; k < j; k++)
    {
      if (allow_paired_[offset_[k+1]+j])
        posterior_[offset_[k+1]+j] += F5i_[k] * FCi_[offset_[k+1]+j-1]
          * ScoringModel::ScoreF5Bifurcation(seq_,k,j) * F5o_[j] / Z;
    }
  }

  for (int i = 1; i <= L; i++)
    for (int j = i+1; j <= L; j++)
      posterior_[offset_[i]+j] = std::min(std::max(posterior_[offset_[i]+j], CountT(0)), CountT(1));
}

template < class ProbT, class CountT, class ScoringModel >
void
InferenceEngine<ProbT, CountT, ScoringModel>::
ComputeExpectation()
{
  const int L = seq_.size()-1;

  const ProbT Z = ComputePartitionCoefficient();

  for (int i=L; i>=0; --i)
  {
    for (int j=i; j<=L; ++j)
    {
      // FM2[i,j] = MAX(i<k<j : FM1[i,k] * FM[k,j])
      ProbT FM2i = ProbT(0);
      for (int k=i+1; k<j; ++k)
        FM2i += FM1i_[offset_[i]+k] * FMi_[offset_[k]+j];

      // FC[i,j] = MAX [ScoreHairpin(i,j),
      //                MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1]),
      //                ScoreJunctionA(i,j) * a * c * MAX (i<k<j : FM1[i,k] * FM[k,j])]
      if (0<i && j<L && allow_paired_[offset_[i]+j+1])
      {
        // compute ScoreHairpin(i,j)
        if (allow_unpaired_[offset_[i]+j] && j-i>=C_MIN_HAIRPIN_LENGTH)
        {
          CountT value = ScoringModel::ScoreFCHairpin(seq_,i,j) * FCo_[offset_[i]+j] / Z;
          ScoringModel::CountFCHairpin(seq_,i,j,value);
        }
        
        // compute MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1])
        for (int p=i; p<=std::min(i+C_MAX_SINGLE_LENGTH,j); ++p)
        {
          if (p>i && !allow_unpaired_position_[p]) break;
          int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
          for (int q=j; q>=q_min; --q)
          {
            if (q<j && !allow_unpaired_position_[q+1]) break;
            if (!allow_paired_[offset_[p+1]+q]) continue;
            CountT value = FCi_[offset_[p+1]+q-1] * ScoringModel::ScoreFCSingle(seq_,i,j,p,q) * FCo_[offset_[i]+j] / Z;
            ScoringModel::CountFCSingle(seq_,i,j,p,q,value);
          }
        }

        // compute MAX (i<k<j : FM1[i,k] * FM[k,j] * ScoreJunctionA(i,j) * a * c)
        {
          CountT value = FM2i * ScoringModel::ScoreFM1Paired(seq_,i,j) * FCo_[offset_[i]+j] / Z;
          ScoringModel::CountFM1Paired(seq_,i,j,value);
        }
      }

      // FM1[i,j] = MAX [FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)  if i+2<=j,
      //                 FM1[i+1,j] * b                                          if i+2<=j]
      if (0<i && i+2<=j && j<L)
      {
        // compute FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)
        if (allow_paired_[offset_[i+1]+j])
        {
          CountT value = FCi_[offset_[i+1]+j-1] * ScoringModel::ScoreFM1Paired(seq_,i,j) * FM1o_[offset_[i]+j] / Z;
          ScoringModel::CountFM1Paired(seq_,i,j,value);
        }

        // compute FM1[i+1,j] * b
        if (allow_unpaired_position_[i+1])
        {
          CountT value = FM1i_[offset_[i+1]+j] * ScoringModel::ScoreFM1Unpaired(seq_,i,j) * FM1o_[offset_[i]+j] / Z;
          ScoringModel::CountFM1Unpaired(seq_,i,j,value);
        }
      }

      // FM[i,j] = MAX [MAX (i<k<j : FM1[i,k] * FM[k,j]),
      //                FM[i,j-1] * b,
      //                FM1[i,j]]
      if (0<i && i+2<=j && j<L)
      {
        // compute MAX (i<k<j : FM1[i,k] * FM[k,j])
        CountT value = FM2i * ScoringModel::ScoreFMBifurcation(seq_,i,j) * FMo_[offset_[i]+j] / Z;
        ScoringModel::CountFMBifurcation(seq_,i,j,value);

        // compute FM[i,j-1] * b
        if (allow_unpaired_position_[j])
        {
          CountT value = FMi_[offset_[i]+j-1] * ScoringModel::ScoreFMUnpaired(seq_,i,j) * FMo_[offset_[i]+j] / Z;
          ScoringModel::CountFMUnpaired(seq_,i,j,value);
        }

        // compute FM1[i,j]
        value = FM1i_[offset_[i]+j] * ScoringModel::ScoreFMtoFM1(seq_,i,j) * FMo_[offset_[i]+j] / Z;
        ScoringModel::CountFMtoFM1(seq_,i,j,value);
      }
    }
  }

  for (int j=1; j<=L; ++j)
  {
    // F5[j] = MAX [F5[j-1] * ScoreExternalUnpaired(),
    //              MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k)

    // compute F5[j-1] * ScoreExternalUnpaired()
    if (allow_unpaired_position_[j])
    {
      CountT value = F5i_[j-1] * ScoringModel::ScoreF5Unpaired(seq_,j) * F5o_[j] / Z;
      ScoringModel::CountF5Unpaired(seq_,j,value);
    }
    
    // compute MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k))
    for (int k = 0; k < j; k++)
    {
      if (allow_paired_[offset_[k+1]+j])
      {
        CountT value = F5i_[k] * FCi_[offset_[k+1]+j-1] * ScoringModel::ScoreF5Bifurcation(seq_,k,j) * F5o_[j] / Z;
        ScoringModel::CountF5Bifurcation(seq_,k,j,value);
      }
    }
  }
}

template < class ProbT, class CountT, class ScoringModel >
void
InferenceEngine<ProbT, CountT, ScoringModel>::
ComputeExpectation2(const std::vector<CountT>& weight)
{
  const int L = seq_.size()-1;

  const ProbT Z = ComputePartitionCoefficient();

  for (int i=L; i>=0; --i)
  {
    for (int j=i; j<=L; ++j)
    {
      // FM2[i,j] = MAX(i<k<j : FM1[i,k] * FM[k,j])
      ProbT FM2i = ProbT(0);
      ProbT FM2j = ProbT(0);
      for (int k=i+1; k<j; ++k)
      {
        FM2i += FM1i_[offset_[i]+k] * FMi_[offset_[k]+j];
        FM2j += FM1i_[offset_[i]+k] * FMj_[offset_[k]+j] +
          FM1j_[offset_[i]+k] * FMi_[offset_[k]+j];
      }

      // FC[i,j] = MAX [ScoreHairpin(i,j),
      //                MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1]),
      //                ScoreJunctionA(i,j) * a * c * MAX (i<k<j : FM1[i,k] * FM[k,j])]
      if (0<i && j<L && allow_paired_[offset_[i]+j+1])
      {
        // compute ScoreHairpin(i,j)
        if (allow_unpaired_[offset_[i]+j] && j-i>=C_MIN_HAIRPIN_LENGTH)
        {
          CountT value = ScoringModel::ScoreFCHairpin(seq_,i,j) * FCp_[offset_[i]+j] / Z;
          ScoringModel::CountFCHairpin(seq_,i,j,value);
        }
        
        // compute MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) * FC[p+1,q-1])
        for (int p=i; p<=std::min(i+C_MAX_SINGLE_LENGTH,j); ++p)
        {
          if (p>i && !allow_unpaired_position_[p]) break;
          int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
          for (int q=j; q>=q_min; --q)
          {
            if (q<j && !allow_unpaired_position_[q+1]) break;
            if (!allow_paired_[offset_[p+1]+q]) continue;
            CountT value =
              (FCj_[offset_[p+1]+q-1] * FCo_[offset_[i]+j]
               + FCi_[offset_[p+1]+q-1] * FCp_[offset_[i]+j]
               + FCi_[offset_[p+1]+q-1] * FCo_[offset_[i]+j] * weight[offset_[p+1]+q])
              * ScoringModel::ScoreFCSingle(seq_,i,j,p,q) / Z;
            ScoringModel::CountFCSingle(seq_,i,j,p,q,value);
          }
        }

        // compute MAX (i<k<j : FM1[i,k] * FM[k,j] * ScoreJunctionA(i,j) * a * c)
        {
          CountT value = (FM2j * FCo_[offset_[i]+j] + FM2i * FCp_[offset_[i]+j])
            * ScoringModel::ScoreFM1Paired(seq_,i,j) / Z;
          ScoringModel::CountFM1Paired(seq_,i,j,value);
        }
      }

      // FM1[i,j] = MAX [FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)  if i+2<=j,
      //                 FM1[i+1,j] * b                                          if i+2<=j]
      if (0<i && i+2<=j && j<L)
      {
        // compute FC[i+1,j-1] * ScoreJunctionA(j,i) * c * ScoreBP(i+1,j)
        if (allow_paired_[offset_[i+1]+j])
        {
          CountT value =
            (FCj_[offset_[i+1]+j-1] * FM1o_[offset_[i]+j]
             + FCi_[offset_[i+1]+j-1] * FM1p_[offset_[i]+j]
             + FCi_[offset_[i+1]+j-1] * FM1o_[offset_[i]+j] * weight[offset_[i+1]+j])
            * ScoringModel::ScoreFM1Paired(seq_,i,j) / Z;
          ScoringModel::CountFM1Paired(seq_,i,j,value);
        }

        // compute FM1[i+1,j] * b
        if (allow_unpaired_position_[i+1])
        {
          CountT value = (FM1j_[offset_[i+1]+j] * FM1o_[offset_[i]+j] + FM1i_[offset_[i+1]+j] * FM1p_[offset_[i]+j])
            * ScoringModel::ScoreFM1Unpaired(seq_,i,j) / Z;
          ScoringModel::CountFM1Unpaired(seq_,i,j,value);
        }
      }

      // FM[i,j] = MAX [MAX (i<k<j : FM1[i,k] * FM[k,j]),
      //                FM[i,j-1] * b,
      //                FM1[i,j]]
      if (0<i && i+2<=j && j<L)
      {
        // compute MAX (i<k<j : FM1[i,k] * FM[k,j])
        CountT value = (FM2j * FMo_[offset_[i]+j] + FM2i * FMp_[offset_[i]+j])
          * ScoringModel::ScoreFMBifurcation(seq_,i,j) / Z;
        ScoringModel::CountFMBifurcation(seq_,i,j,value);

        // compute FM[i,j-1] * b
        if (allow_unpaired_position_[j])
        {
          CountT value = (FMj_[offset_[i]+j-1] * FMo_[offset_[i]+j] + FMi_[offset_[i]+j-1] * FMp_[offset_[i]+j])
            * ScoringModel::ScoreFMUnpaired(seq_,i,j) / Z;
          ScoringModel::CountFMUnpaired(seq_,i,j,value);
        }

        // compute FM1[i,j]
        value = (FM1j_[offset_[i]+j] * FMo_[offset_[i]+j] + FM1i_[offset_[i]+j] * FMp_[offset_[i]+j])
          * ScoringModel::ScoreFMtoFM1(seq_,i,j) / Z;
        ScoringModel::CountFMtoFM1(seq_,i,j,value);
      }
    }
  }

  for (int j=1; j<=L; ++j)
  {
    // F5[j] = MAX [F5[j-1] * ScoreExternalUnpaired(),
    //              MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k)

    // compute F5[j-1] * ScoreExternalUnpaired()
    if (allow_unpaired_position_[j])
    {
      CountT value = (F5j_[j-1] * F5o_[j] + F5i_[j-1] * F5p_[j])
        * ScoringModel::ScoreF5Unpaired(seq_,j) / Z;
      ScoringModel::CountF5Unpaired(seq_,j,value);
    }
    
    // compute MAX (0<=k<j : F5[k] * FC[k+1,j-1] * ScoreExternalPaired() * ScoreBP(k+1,j) * ScoreJunctionA(j,k))
    for (int k = 0; k < j; k++)
    {
      if (allow_paired_[offset_[k+1]+j])
      {
        CountT value =
          (F5j_[k] * FCi_[offset_[k+1]+j-1] * F5o_[j]
           + F5i_[k] * FCj_[offset_[k+1]+j-1] * F5o_[j]
           + F5i_[k] * FCi_[offset_[k+1]+j-1] * F5p_[j]
           + F5i_[k] * FCi_[offset_[k+1]+j-1] * F5o_[j] * weight[offset_[k+1]+j])
          * ScoringModel::ScoreF5Bifurcation(seq_,k,j) / Z;
        ScoringModel::CountF5Bifurcation(seq_,k,j,value);
      }
    }
  }
}

// instantiation
#include "ScoringModel.h"
#include "log_value.h"

typedef LogValue<float> ProbT;
typedef float CountT;
template class InferenceEngine<ProbT,CountT,ScoringModel<ProbT,CountT> >;
