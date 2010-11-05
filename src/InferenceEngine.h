// InferenceEngine

#ifndef __INC_INFERENCE_ENGINE_H__
#define __INC_INFERENCE_ENGINE_H__

#include <vector>
#include <string>

template < class ProbT, class CountT, class ScoringModel >
class InferenceEngine : public ScoringModel
{
  using ScoringModel::C_MIN_HAIRPIN_LENGTH;
  using ScoringModel::C_MAX_SINGLE_LENGTH;
  using ScoringModel::M;
  using ScoringModel::alphabet;

public:
  InferenceEngine(bool allow_noncomplementary);
  ~InferenceEngine() { };

  void LoadSequence(const std::string& rna);
  void UseConstraints(const std::vector<int> &true_mapping);

  ProbT ComputeViterbi();
  std::vector<int> DecodeViterbi();
  ProbT ComputeInside();
  ProbT ComputeOutside();
  void ComputePosterior();
  std::vector<int> DecodePosterior();
  void ComputeExpectation();
  ProbT ComputeInside2(const std::vector<CountT>& weight);
  ProbT ComputeOutside2(const std::vector<CountT>& weight);
  void ComputeExpectation2(const std::vector<CountT>& weight);
  ProbT ComputePartitionCoefficient() const
  {
    return F5i_[seq_.size()-1];
  }
  void GetPosterior(std::vector<CountT>& p) const
  {
    p = posterior_;
  }
  void GetPosterior(std::vector<CountT>& p, std::vector<int>& offset) const
  {
    p = posterior_;
    offset = offset_;
  }

private:
  int ComputeRowOffset(int i, int N) const;
  bool IsComplementary(int i, int j) const;
  void UpdateViterbi(ProbT& bs, int& bt, ProbT s, int t);
  int EncodeTraceback(int i, int j) const;
  std::pair<int,int> DecodeTraceback(int s) const;

  bool allow_noncomplementary_;
  std::vector<int> char_mapping_;
  std::vector<std::vector<int> > is_complementary_;

  // sequence data
  std::vector<int> seq_;
  std::vector<int> offset_;
  std::vector<int> allow_unpaired_position_;
  std::vector<int> allow_unpaired_;
  std::vector<int> allow_paired_;

  // dynamic programming matrices
  std::vector<int> FCt_, F5t_, FMt_, FM1t_;   // traceback
  std::vector<ProbT> FCv_, F5v_, FMv_, FM1v_; // Viterbi  
  std::vector<ProbT> FCi_, F5i_, FMi_, FM1i_; // inside
  std::vector<ProbT> FCj_, F5j_, FMj_, FM1j_; // inside2
  std::vector<ProbT> FCo_, F5o_, FMo_, FM1o_; // outside
  std::vector<ProbT> FCp_, F5p_, FMp_, FM1p_; // outside2
  std::vector<CountT> posterior_;             // posterior

  enum TRACEBACK_TYPE {
    TB_FC_HAIRPIN,
    TB_FC_SINGLE,
    TB_FC_BIFURCATION,
    TB_FM1_PAIRED,
    TB_FM1_UNPAIRED,
    TB_FM_BIFURCATION,
    TB_FM_UNPAIRED,
    TB_FM_FM1,
    TB_F5_ZERO,
    TB_F5_UNPAIRED,
    TB_F5_BIFURCATION,
    NUM_TRACEBACK_TYPES
  };
};

#endif  //  __INC_INFERENCE_ENGINE_H__

// Local Variables:
// mode: C++
// End:
