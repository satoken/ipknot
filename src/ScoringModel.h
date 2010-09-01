// ScoringModel

#ifndef __INC_SCORING_MODEL_H__
#define __INC_SCORING_MODEL_H__

#include <vector>
#include <string>
#include <map>

template <class ProbT, class CountT>
class ScoringModel
{
protected:
  enum {
    C_MIN_HAIRPIN_LENGTH = 0,
    C_MAX_SINGLE_LENGTH = 30,
    D_MAX_HAIRPIN_LENGTH = 30,
    D_MAX_BULGE_LENGTH = 30,
    D_MAX_INTERNAL_LENGTH = 30,
    D_MAX_INTERNAL_SYMMETRIC_LENGTH = 15,
    D_MAX_INTERNAL_ASYMMETRY = 28,
    D_MAX_INTERNAL_EXPLICIT_LENGTH = 4,
  };
  static const std::string alphabet/* = "ACGU"*/;    // allowed symbols -- all other letters ignored
  static const int M = 4;                        // number of alphabet symbols

public:
  ScoringModel();

  void LoadParameters(const std::string& fname, bool as_log_scale=false);
  void SaveParameters(const std::string& fname, bool as_log_scale=false) const;
  void SetParameters(const std::vector<ProbT>& param);
  void GetParameters(std::vector<ProbT>& param) const;
  void GetCounts(std::vector<CountT>& cnt);
  void ClearCounts();
  std::vector<std::string> GetParameterNames() const
  {
    return name_;
  }

protected:
  ProbT ScoreFCHairpin(const std::vector<int>& seq, int i, int j) const;
  void CountFCHairpin(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreFCSingle(const std::vector<int>& seq, int i, int j, int p, int q) const;
  void CountFCSingle(const std::vector<int>& seq, int i, int j, int p, int q, CountT c);

  ProbT ScoreFCBifurcation(const std::vector<int>& seq, int i, int j) const;
  void CountCBifurcation(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreFM1Paired(const std::vector<int>& seq, int i, int j) const;
  void CountFM1Paired(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreFM1Unpaired(const std::vector<int>& seq, int i, int j) const;
  void CountFM1Unpaired(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreFMBifurcation(const std::vector<int>& seq, int i, int j) const;
  void CountFMBifurcation(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreFMUnpaired(const std::vector<int>& seq, int i, int j) const;
  void CountFMUnpaired(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreFMtoFM1(const std::vector<int>& seq, int i, int j) const;
  void CountFMtoFM1(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreF5Unpaired(const std::vector<int>& seq, int j) const;
  void CountF5Unpaired(const std::vector<int>& seq, int j, CountT c);

  ProbT ScoreF5Bifurcation(const std::vector<int>& seq, int j, int k) const;
  void CountF5Bifurcation(const std::vector<int>& seq, int j, int k, CountT c);

//private:
  void MakeMapping();
  void InitializeCache();
  void FinalizeCounts();

  ProbT ScoreHairpin(const std::vector<int>& seq, int i, int j) const;
  void CountHairpin(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreBasePair(const std::vector<int>& seq, int i, int j) const;
  void CountBasePair(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreHelixStacking(const std::vector<int>& seq, int i, int j) const;
  void CountHelixStacking(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreSingle(const std::vector<int>& seq, int i, int j, int p, int q) const;
  void CountSingle(const std::vector<int>& seq, int i, int j, int p, int q, CountT c);

  ProbT ScoreSingleNucleotides(const std::vector<int>& seq, int i, int j, int p, int q) const;
  void CountSingleNucleotides(const std::vector<int>& seq, int i, int j, int p, int q, CountT c);

  ProbT ScoreJunctionA(const std::vector<int>& seq, int i, int j) const;
  void CountJunctionA(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreJunctionB(const std::vector<int>& seq, int i, int j) const;
  void CountJunctionB(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreMultiPaired() const;
  void CountMultiPaired(CountT c);

  ProbT ScoreMultiBase() const;
  void CountMultiBase(CountT c);

  ProbT ScoreMultiUnpaired(const std::vector<int>& seq, int i) const;
  void CountMultiUnpaired(const std::vector<int>& seq, int i, CountT c);

  ProbT ScoreExternalPaired() const;
  void CountExternalPaired(CountT c);

  ProbT ScoreExternalUnpaired(const std::vector<int>& seq, int i) const;
  void CountExternalUnpaired(const std::vector<int>& seq, int i, CountT c);

  ProbT ScoreUnpaired(const std::vector<int>& seq, int i, int j) const;
  void CountUnpaired(const std::vector<int>& seq, int i, int j, CountT c);

  ProbT ScoreUnpairedPosition(const std::vector<int>& seq, int i) const;
  void CountUnpairedPosition(const std::vector<int>& seq, int i, CountT c) const;
  
  // parameters
  std::pair<ProbT,CountT> score_base_pair[M+1][M+1];
  std::pair<ProbT,CountT> score_terminal_mismatch[M+1][M+1][M+1][M+1];
  std::pair<ProbT,CountT> score_hairpin_length_at_least[D_MAX_HAIRPIN_LENGTH+1];
  std::pair<ProbT,CountT> score_bulge_length_at_least[D_MAX_BULGE_LENGTH+1];
  std::pair<ProbT,CountT> score_bulge_0x1_nucleotides[M+1];
  std::pair<ProbT,CountT> score_bulge_1x0_nucleotides[M+1];
  std::pair<ProbT,CountT> score_internal_explicit[D_MAX_INTERNAL_EXPLICIT_LENGTH+1][D_MAX_INTERNAL_EXPLICIT_LENGTH+1];
  std::pair<ProbT,CountT> score_internal_length_at_least[D_MAX_INTERNAL_LENGTH+1];
  std::pair<ProbT,CountT> score_internal_symmetric_length_at_least[D_MAX_INTERNAL_SYMMETRIC_LENGTH+1];
  std::pair<ProbT,CountT> score_internal_asymmetry_at_least[D_MAX_INTERNAL_ASYMMETRY+1];
  std::pair<ProbT,CountT> score_internal_1x1_nucleotides[M+1][M+1];
  std::pair<ProbT,CountT> score_helix_stacking[M+1][M+1][M+1][M+1];
  std::pair<ProbT,CountT> score_helix_closing[M+1][M+1];
  std::pair<ProbT,CountT> score_multi_base;
  std::pair<ProbT,CountT> score_multi_unpaired;
  std::pair<ProbT,CountT> score_multi_paired;
  std::pair<ProbT,CountT> score_dangle_left[M+1][M+1][M+1];
  std::pair<ProbT,CountT> score_dangle_right[M+1][M+1][M+1];
  std::pair<ProbT,CountT> score_external_unpaired;
  std::pair<ProbT,CountT> score_external_paired;

  // cached values
  std::pair<ProbT,CountT> cache_score_hairpin_length[D_MAX_HAIRPIN_LENGTH+1];
  std::pair<ProbT,CountT> cache_score_single[C_MAX_SINGLE_LENGTH+1][C_MAX_SINGLE_LENGTH+1];

  // index mapping
  std::map<std::string,int> index_;
  std::vector< std::vector<std::pair<ProbT,CountT>*> > rev_;
  std::vector<std::string> name_;
};


#endif  //  __INC_SCORING_MODEL_H__

// Local Variables:
// mode: C++
// End:
